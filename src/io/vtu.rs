use std::str::FromStr;

use roxmltree::{Document, Node};

use crate::mesh::{BoundaryNodes, Mesh};

/// Prefix of the name of a point data array which gives one boundary:
/// `boundary__<name>`
pub const BOUNDARY_ARRAY_PREFIX: &str = "boundary__";

/// VTK cell type of a linear triangle
const VTK_TRIANGLE: u8 = 5;

/// Reads the fracture mesh from the text of a VTK XML UnstructuredGrid file
/// (`.vtu`, ASCII, one `Piece`).
///
/// All cells must be triangles. The cell data arrays `fracture_id`,
/// `aperture` and `permeability` are necessary. The point data arrays
/// `boundary__<name>` (optional) give the boundary nodes.
pub fn parse_mesh(text: &str) -> Result<Mesh, String> {
    let document = Document::parse(text).map_err(|e| e.to_string())?;
    let root = document.root_element();

    if !root.has_tag_name("VTKFile") || root.attribute("type") != Some("UnstructuredGrid") {
        return Err("the file is not a VTK UnstructuredGrid (.vtu) file".to_string());
    }

    let pieces: Vec<Node> = root
        .descendants()
        .filter(|node| node.has_tag_name("Piece"))
        .collect();

    match pieces.as_slice() {
        [piece] => read_piece(*piece),
        _ => Err(format!(
            "the file must have one Piece, but it has {}; the points of different pieces are not connected",
            pieces.len()
        )),
    }
}

fn read_piece(piece: Node) -> Result<Mesh, String> {
    let points_element = find_child(piece, "Points")?;
    let cells_element = find_child(piece, "Cells")?;

    let coordinates: Vec<f64> = read_values(find_child(points_element, "DataArray")?)?;
    let connectivity: Vec<usize> = read_values(find_data_array(cells_element, "connectivity")?)?;
    let offsets: Vec<usize> = read_values(find_data_array(cells_element, "offsets")?)?;
    let types: Vec<u8> = read_values(find_data_array(cells_element, "types")?)?;

    let (point_coordinates, remaining_coordinates) = coordinates.as_chunks::<3>();

    if !remaining_coordinates.is_empty() {
        return Err("the points must have 3 coordinates each".to_string());
    }
    if offsets.len() != types.len() {
        return Err("offsets and types must have one value for each cell".to_string());
    }

    let points = point_coordinates.to_vec();
    let mut triangles = Vec::with_capacity(types.len());
    let mut cell_start = 0;

    for (cell, (&cell_end, &cell_type)) in offsets.iter().zip(&types).enumerate() {
        if cell_type != VTK_TRIANGLE {
            return Err(format!(
                "cell {cell}: the cell type is {cell_type}, but all cells must be triangles (type {VTK_TRIANGLE})"
            ));
        }

        let point_ids = connectivity
            .get(cell_start..cell_end)
            .ok_or(format!("cell {cell}: the offset {cell_end} is not valid"))?;
        let triangle: [usize; 3] = point_ids
            .try_into()
            .map_err(|_| format!("cell {cell}: a triangle must have 3 points"))?;

        if let Some(point_id) = triangle.iter().find(|&&id| id >= points.len()) {
            return Err(format!("cell {cell}: there is no point {point_id}"));
        }
        if triangle[0] == triangle[1] || triangle[1] == triangle[2] || triangle[0] == triangle[2] {
            return Err(format!(
                "cell {cell}: the triangle uses one point two times"
            ));
        }

        triangles.push(triangle);
        cell_start = cell_end;
    }

    if triangles.is_empty() {
        return Err("the mesh has no cells".to_string());
    }

    let mut is_point_used = vec![false; points.len()];
    for &point_id in triangles.iter().flatten() {
        is_point_used[point_id] = true;
    }
    if let Some(point_id) = is_point_used.iter().position(|&is_used| !is_used) {
        return Err(format!("point {point_id} is not a point of a triangle"));
    }

    let cell_data = find_child(piece, "CellData")
        .map_err(|_| "the file has no CellData; the arrays fracture_id, aperture and permeability are necessary".to_string())?;
    let fracture_ids: Vec<i64> = read_cell_array(cell_data, "fracture_id", triangles.len())?;
    let apertures: Vec<f64> = read_cell_array(cell_data, "aperture", triangles.len())?;
    let permeabilities: Vec<f64> = read_cell_array(cell_data, "permeability", triangles.len())?;

    check_positive(&apertures, "aperture")?;
    check_positive(&permeabilities, "permeability")?;

    Ok(Mesh {
        boundaries: read_boundary_arrays(piece, points.len())?,
        points,
        triangles,
        fracture_ids,
        apertures,
        permeabilities,
    })
}

fn check_positive(values: &[f64], name: &str) -> Result<(), String> {
    match values
        .iter()
        .position(|value| !(*value > 0.0 && value.is_finite()))
    {
        Some(cell) => Err(format!(
            "cell {cell}: the {name} is {}, but it must be a finite value above 0",
            values[cell]
        )),
        None => Ok(()),
    }
}

/// Reads a cell data array with one value for each cell
fn read_cell_array<T: FromStr>(
    cell_data: Node,
    name: &str,
    cell_count: usize,
) -> Result<Vec<T>, String> {
    let data_array = find_data_array(cell_data, name)
        .map_err(|_| format!("the cell data array {name} is necessary"))?;
    let values = read_scalar_values(data_array, name)?;

    if values.len() != cell_count {
        return Err(format!(
            "the array {name} must have one value for each cell, but it has {} values for {cell_count} cells",
            values.len()
        ));
    }

    Ok(values)
}

/// Reads the point data arrays `boundary__<name>`, in the order of the file
fn read_boundary_arrays(piece: Node, point_count: usize) -> Result<Vec<BoundaryNodes>, String> {
    let Some(point_data) = piece.children().find(|node| node.has_tag_name("PointData")) else {
        return Ok(Vec::new());
    };

    let mut boundaries: Vec<BoundaryNodes> = Vec::new();

    for data_array in point_data
        .children()
        .filter(|node| node.has_tag_name("DataArray"))
    {
        let array_name = data_array.attribute("Name").unwrap_or_default();
        let Some(name) = array_name.strip_prefix(BOUNDARY_ARRAY_PREFIX) else {
            continue;
        };

        if name.is_empty() {
            return Err(format!(
                "the array {array_name} has no boundary name after the prefix"
            ));
        }
        if boundaries.iter().any(|boundary| boundary.name == name) {
            return Err(format!("the file has two arrays {array_name}"));
        }

        let values: Vec<f64> = read_scalar_values(data_array, array_name)?;

        if values.len() != point_count {
            return Err(format!(
                "the array {array_name} must have one value for each point, but it has {} values for {point_count} points",
                values.len()
            ));
        }

        boundaries.push(BoundaryNodes {
            name: name.to_string(),
            nodes: (0..point_count)
                .filter(|&point| values[point] > 0.0)
                .collect(),
        });
    }

    Ok(boundaries)
}

/// Reads a `DataArray` with one component
fn read_scalar_values<T: FromStr>(data_array: Node, name: &str) -> Result<Vec<T>, String> {
    let component_count = data_array.attribute("NumberOfComponents").unwrap_or("1");

    if component_count.trim() != "1" {
        return Err(format!(
            "the array {name} must have one component (NumberOfComponents=\"1\")"
        ));
    }

    read_values(data_array)
}

/// Finds the first child element with the tag name `name`
fn find_child<'a, 'input>(node: Node<'a, 'input>, name: &str) -> Result<Node<'a, 'input>, String> {
    node.children()
        .find(|child| child.has_tag_name(name))
        .ok_or(format!(
            "the element {} has no {name} element",
            node.tag_name().name()
        ))
}

/// Finds the child `DataArray` element with the attribute `Name="name"`
fn find_data_array<'a, 'input>(
    node: Node<'a, 'input>,
    name: &str,
) -> Result<Node<'a, 'input>, String> {
    node.children()
        .find(|child| child.has_tag_name("DataArray") && child.attribute("Name") == Some(name))
        .ok_or(format!(
            "the element {} has no DataArray {name}",
            node.tag_name().name()
        ))
}

/// Reads the values of an ASCII `DataArray` element
fn read_values<T: FromStr>(data_array: Node) -> Result<Vec<T>, String> {
    let name = data_array.attribute("Name").unwrap_or_default();

    if data_array.attribute("format") != Some("ascii") {
        return Err(format!(
            "the DataArray {name} is not ASCII (format=\"ascii\"); write the file in ASCII mode"
        ));
    }

    data_array
        .text()
        .unwrap_or_default()
        .split_whitespace()
        .map(|value| {
            value
                .parse()
                .map_err(|_| format!("the DataArray {name} has the invalid value \"{value}\""))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two triangles of the unit square in the plane z = 0
    fn unit_square(cell_types: &str, cell_data: &str, point_data: &str) -> String {
        format!(
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="2">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0 0 0  1 0 0  1 1 0  0 1 0
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 0 2 3</DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">3 6</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">{cell_types}</DataArray>
      </Cells>
      <CellData>{cell_data}</CellData>
      <PointData>{point_data}</PointData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>"#
        )
    }

    const CELL_DATA: &str = r#"
        <DataArray type="Int32" Name="fracture_id" format="ascii">1 1</DataArray>
        <DataArray type="Float64" Name="aperture" format="ascii">1e-3 1e-3</DataArray>
        <DataArray type="Float64" Name="permeability" format="ascii">1e-12 1e-12</DataArray>"#;

    const POINT_DATA: &str = r#"
        <DataArray type="Int32" Name="boundary__left" format="ascii">1 0 0 1</DataArray>
        <DataArray type="Int32" Name="other" format="ascii">1 1 1 1</DataArray>"#;

    #[test]
    fn reads_triangles_cell_arrays_and_boundaries() {
        let mesh = parse_mesh(&unit_square("5 5", CELL_DATA, POINT_DATA)).unwrap();

        assert_eq!(mesh.points.len(), 4);
        assert_eq!(mesh.triangles, vec![[0, 1, 2], [0, 2, 3]]);
        assert_eq!(mesh.fracture_ids, vec![1, 1]);
        assert_eq!(mesh.apertures, vec![1e-3, 1e-3]);
        assert_eq!(mesh.boundaries.len(), 1);
        assert_eq!(mesh.boundaries[0].name, "left");
        assert_eq!(mesh.boundaries[0].nodes, vec![0, 3]);
    }

    #[test]
    fn non_triangle_cell_is_an_error() {
        let error = parse_mesh(&unit_square("5 9", CELL_DATA, POINT_DATA)).unwrap_err();
        assert!(error.contains("cell 1: the cell type is 9"), "{error}");
    }

    #[test]
    fn missing_cell_array_is_an_error() {
        let cell_data = CELL_DATA.replace("permeability", "perm");
        let error = parse_mesh(&unit_square("5 5", &cell_data, POINT_DATA)).unwrap_err();
        assert!(error.contains("permeability is necessary"), "{error}");
    }

    #[test]
    fn negative_aperture_is_an_error() {
        let cell_data = CELL_DATA.replace("1e-3 1e-3", "1e-3 -1e-3");
        let error = parse_mesh(&unit_square("5 5", &cell_data, POINT_DATA)).unwrap_err();
        assert!(error.contains("cell 1: the aperture is -0.001"), "{error}");
    }

    #[test]
    fn binary_data_array_is_an_error() {
        let cell_data = CELL_DATA.replacen("format=\"ascii\"", "format=\"binary\"", 1);
        let error = parse_mesh(&unit_square("5 5", &cell_data, POINT_DATA)).unwrap_err();
        assert!(error.contains("is not ASCII"), "{error}");
    }
}
