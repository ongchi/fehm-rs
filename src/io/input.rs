use std::path::{Path, PathBuf};

use serde::Deserialize;

use tracing::warn;

use crate::FehmError;
use crate::eos;
use crate::io::vtu::{BOUNDARY_ARRAY_PREFIX, parse_mesh};
use crate::mesh::Mesh;

/// All data of one run
#[derive(Clone, Debug)]
pub struct Input {
    pub mesh: Mesh,
    /// Fluid temperature (°C)
    pub temperature: f64,
    /// True: gravity in the -z direction
    pub has_gravity: bool,
    /// In the order of the input file
    pub boundaries: Vec<PressureBoundary>,
    pub solver: SolverOptions,
}

/// Nodes with a fixed pressure
#[derive(Clone, Debug)]
pub struct PressureBoundary {
    /// `<name>` of the point data array `boundary__<name>`
    pub name: String,
    /// MPa
    pub pressure: f64,
    /// In increasing order
    pub nodes: Vec<usize>,
}

/// Options of the Newton iterations
#[derive(Clone, Debug)]
pub struct SolverOptions {
    /// Stop when the largest mass residual of a node is less than
    /// `tolerance` times the largest sum of the absolute mass fluxes of a node
    pub tolerance: f64,
    pub max_iterations: usize,
}

const DEFAULT_TOLERANCE: f64 = 1e-10;
const DEFAULT_MAX_ITERATIONS: usize = 20;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct InputFile {
    mesh: String,
    temperature: f64,
    gravity: Option<bool>,
    #[serde(default)]
    boundary: Vec<BoundarySection>,
    solver: Option<SolverSection>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BoundarySection {
    name: String,
    pressure: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SolverSection {
    tolerance: Option<f64>,
    max_iterations: Option<usize>,
}

/// Reads the input file and the mesh file which it names. A relative mesh
/// path is relative to the folder of the input file.
pub fn read_input(path: &Path) -> Result<Input, FehmError> {
    let input_text = read_file(path, "input file")?;
    let mesh_path = path
        .parent()
        .unwrap_or(Path::new(""))
        .join(read_mesh_file_name(&input_text)?);
    let mesh_text = read_file(&mesh_path, "mesh file")?;

    build_input_from_text(&input_text, &mesh_text).map_err(|error| match error {
        FehmError::InvalidMesh(message) => {
            FehmError::InvalidMesh(format!("{}: {message}", mesh_path.display()))
        }
        other => other,
    })
}

fn read_file(path: &Path, description: &str) -> Result<String, FehmError> {
    std::fs::read_to_string(path).map_err(|e| {
        FehmError::Io(std::io::Error::new(
            e.kind(),
            format!("cannot read the {description} {}: {e}", path.display()),
        ))
    })
}

/// The value of the `mesh` key of the input file
pub fn read_mesh_file_name(input_text: &str) -> Result<PathBuf, FehmError> {
    Ok(PathBuf::from(parse_input_file(input_text)?.mesh))
}

fn parse_input_file(input_text: &str) -> Result<InputFile, FehmError> {
    toml::from_str(input_text).map_err(|e| FehmError::InvalidInput(e.to_string()))
}

/// Builds the input from the text of the input file and of the mesh file
pub fn build_input_from_text(input_text: &str, mesh_text: &str) -> Result<Input, FehmError> {
    let input_file = parse_input_file(input_text)?;
    let mesh = parse_mesh(mesh_text).map_err(FehmError::InvalidMesh)?;
    let invalid = |message: String| Err(FehmError::InvalidInput(message));

    if !(eos::MIN_TEMPERATURE..=eos::MAX_TEMPERATURE).contains(&input_file.temperature) {
        return invalid(format!(
            "temperature = {} °C is outside the range of the water properties of FEHM ({} to {} °C)",
            input_file.temperature,
            eos::MIN_TEMPERATURE,
            eos::MAX_TEMPERATURE
        ));
    }
    if input_file.boundary.is_empty() {
        return invalid(
            "there is no [[boundary]]; at least one boundary with a fixed pressure is necessary"
                .to_string(),
        );
    }

    let mut boundaries: Vec<PressureBoundary> = Vec::new();

    for section in input_file.boundary {
        if boundaries
            .iter()
            .any(|boundary| boundary.name == section.name)
        {
            return invalid(format!(
                "two [[boundary]] tables have the name \"{}\"",
                section.name
            ));
        }
        if !(eos::MIN_PRESSURE..=eos::MAX_PRESSURE).contains(&section.pressure) {
            return invalid(format!(
                "boundary \"{}\": pressure = {} MPa is outside the range of the water properties of FEHM ({} to {} MPa)",
                section.name,
                section.pressure,
                eos::MIN_PRESSURE,
                eos::MAX_PRESSURE
            ));
        }

        let Some(boundary_nodes) = mesh
            .boundaries
            .iter()
            .find(|boundary_nodes| boundary_nodes.name == section.name)
        else {
            return invalid(format!(
                "boundary \"{}\": the mesh has no point data array {BOUNDARY_ARRAY_PREFIX}{}",
                section.name, section.name
            ));
        };
        if boundary_nodes.nodes.is_empty() {
            return invalid(format!(
                "boundary \"{}\": the array {BOUNDARY_ARRAY_PREFIX}{} has no point with a value above 0",
                section.name, section.name
            ));
        }

        boundaries.push(PressureBoundary {
            name: section.name,
            pressure: section.pressure,
            nodes: boundary_nodes.nodes.clone(),
        });
    }

    resolve_boundary_overlaps(&mut boundaries, mesh.points.len())?;

    let solver_section = input_file.solver.unwrap_or(SolverSection {
        tolerance: None,
        max_iterations: None,
    });
    let solver = SolverOptions {
        tolerance: solver_section.tolerance.unwrap_or(DEFAULT_TOLERANCE),
        max_iterations: solver_section
            .max_iterations
            .unwrap_or(DEFAULT_MAX_ITERATIONS),
    };

    if !(solver.tolerance > 0.0 && solver.tolerance < 1.0) {
        return invalid(format!(
            "solver.tolerance = {} must be above 0 and below 1",
            solver.tolerance
        ));
    }
    if solver.max_iterations == 0 {
        return invalid("solver.max_iterations must be 1 or more".to_string());
    }

    Ok(Input {
        mesh,
        temperature: input_file.temperature,
        has_gravity: input_file.gravity.unwrap_or(true),
        boundaries,
        solver,
    })
}

/// A node in two boundaries with different pressures keeps the pressure of
/// the first boundary (in the order of the input file): it is removed from
/// the later boundary, with a warning for each pair of boundaries. A node in
/// two boundaries with the same pressure stays in both. A boundary which
/// then has no node is an error.
fn resolve_boundary_overlaps(
    boundaries: &mut [PressureBoundary],
    point_count: usize,
) -> Result<(), FehmError> {
    let mut boundary_of_node: Vec<Option<usize>> = vec![None; point_count];

    for index in 0..boundaries.len() {
        let pressure = boundaries[index].pressure;
        let mut removed_by_first: Vec<(usize, usize)> = Vec::new();

        let mut nodes = std::mem::take(&mut boundaries[index].nodes);
        nodes.retain(|&node| match boundary_of_node[node] {
            Some(first) if boundaries[first].pressure != pressure => {
                match removed_by_first
                    .iter_mut()
                    .find(|(boundary, _)| *boundary == first)
                {
                    Some((_, count)) => *count += 1,
                    None => removed_by_first.push((first, 1)),
                }
                false
            }
            Some(_) => true,
            None => {
                boundary_of_node[node] = Some(index);
                true
            }
        });
        boundaries[index].nodes = nodes;

        for (first, count) in removed_by_first {
            warn!(
                "{count} points are in the boundaries \"{}\" and \"{}\", which have different pressures; they keep the pressure of \"{}\" (the first in the input file)",
                boundaries[first].name, boundaries[index].name, boundaries[first].name
            );
        }

        if boundaries[index].nodes.is_empty() {
            return Err(FehmError::InvalidInput(format!(
                "boundary \"{}\": each of its points is in an earlier boundary with another pressure, thus it has no effect",
                boundaries[index].name
            )));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two triangles of the unit square, with the boundaries left (points 0
    /// and 3), right (points 1 and 2) and corner (point 0)
    const MESH: &str = r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="2">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">0 0 0 1 0 0 1 1 0 0 1 0</DataArray>
      </Points>
      <Cells>
        <DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 0 2 3</DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">3 6</DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">5 5</DataArray>
      </Cells>
      <CellData>
        <DataArray type="Int32" Name="fracture_id" format="ascii">1 1</DataArray>
        <DataArray type="Float64" Name="aperture" format="ascii">1e-3 1e-3</DataArray>
        <DataArray type="Float64" Name="permeability" format="ascii">1e-12 1e-12</DataArray>
      </CellData>
      <PointData>
        <DataArray type="Int32" Name="boundary__left" format="ascii">1 0 0 1</DataArray>
        <DataArray type="Int32" Name="boundary__right" format="ascii">0 1 1 0</DataArray>
        <DataArray type="Int32" Name="boundary__corner" format="ascii">1 0 0 0</DataArray>
        <DataArray type="Int32" Name="boundary__none" format="ascii">0 0 0 0</DataArray>
      </PointData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>"#;

    const VALID_INPUT: &str = r#"
mesh = "mesh.vtu"
temperature = 50.0

[[boundary]]
name = "left"
pressure = 10.1

[[boundary]]
name = "right"
pressure = 10.0
"#;

    fn error_of(input_text: &str) -> String {
        build_input_from_text(input_text, MESH)
            .unwrap_err()
            .to_string()
    }

    #[test]
    fn valid_input_uses_the_defaults() {
        let input = build_input_from_text(VALID_INPUT, MESH).unwrap();

        assert!(input.has_gravity);
        assert_eq!(input.boundaries.len(), 2);
        assert_eq!(input.boundaries[0].nodes, vec![0, 3]);
        assert_eq!(input.boundaries[1].pressure, 10.0);
        assert_eq!(input.solver.tolerance, DEFAULT_TOLERANCE);
        assert_eq!(input.solver.max_iterations, DEFAULT_MAX_ITERATIONS);
    }

    #[test]
    fn unknown_key_is_an_error() {
        let error = error_of(&format!("{VALID_INPUT}\n[solver]\nmax_iteration = 5\n"));
        assert!(error.contains("unknown field `max_iteration`"), "{error}");
    }

    #[test]
    fn missing_temperature_is_an_error() {
        let error = error_of(&VALID_INPUT.replace("temperature = 50.0", ""));
        assert!(error.contains("missing field `temperature`"), "{error}");
    }

    #[test]
    fn unknown_boundary_name_is_an_error() {
        let error = error_of(&VALID_INPUT.replace("\"right\"", "\"top\""));
        assert!(
            error.contains("no point data array boundary__top"),
            "{error}"
        );
    }

    #[test]
    fn empty_boundary_is_an_error() {
        let error = error_of(&VALID_INPUT.replace("\"right\"", "\"none\""));
        assert!(
            error.contains("has no point with a value above 0"),
            "{error}"
        );
    }

    #[test]
    fn no_boundary_is_an_error() {
        let error = error_of("mesh = \"mesh.vtu\"\ntemperature = 50.0\n");
        assert!(error.contains("there is no [[boundary]]"), "{error}");
    }

    #[test]
    fn node_in_two_boundaries_with_different_pressures_keeps_the_first_boundary() {
        let input_text = VALID_INPUT
            .replace("\"left\"", "\"first\"")
            .replace("\"right\"", "\"left\"")
            .replace("\"first\"", "\"corner\"");
        let input = build_input_from_text(&input_text, MESH).unwrap();

        assert_eq!(input.boundaries[0].nodes, vec![0]);
        assert_eq!(input.boundaries[1].nodes, vec![3]);
    }

    #[test]
    fn boundary_with_only_points_of_an_earlier_boundary_is_an_error() {
        let error = error_of(&VALID_INPUT.replace("\"right\"", "\"corner\""));
        assert!(
            error.contains("boundary \"corner\": each of its points is in an earlier boundary"),
            "{error}"
        );
    }

    #[test]
    fn node_in_two_boundaries_with_the_same_pressure_is_valid() {
        let input_text = VALID_INPUT
            .replace("\"right\"", "\"corner\"")
            .replace("10.0", "10.1");
        assert!(build_input_from_text(&input_text, MESH).is_ok());
    }

    #[test]
    fn temperature_outside_the_eos_range_is_an_error() {
        let error = error_of(&VALID_INPUT.replace("50.0", "400.0"));
        assert!(error.contains("temperature = 400 °C is outside"), "{error}");
    }
}
