use std::fmt::Display;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::FehmError;
use crate::flow::Solution;
use crate::io::input::Input;
use crate::velocity::calculate_darcy_velocities;

/// Name of the file with the mesh, the pressures and the velocities
pub const PRESSURE_FILE: &str = "pressure.vtu";
/// Name of the file with the connections and their mass fluxes
pub const CONNECTIONS_FILE: &str = "connections.vtp";

/// Writes the output files into a folder. The folder is made if it does not exist.
pub fn write_output_to_folder(
    input: &Input,
    solution: &Solution,
    folder: &Path,
) -> Result<(), FehmError> {
    std::fs::create_dir_all(folder)?;
    write_output(input, solution, &mut |name| {
        Ok(Box::new(BufWriter::new(File::create(folder.join(name))?)))
    })
}

/// Writes the output files. `create_file` gives the writer of each file from its name.
pub fn write_output(
    input: &Input,
    solution: &Solution,
    create_file: &mut dyn FnMut(&str) -> std::io::Result<Box<dyn Write>>,
) -> Result<(), FehmError> {
    let mut writer = create_file(PRESSURE_FILE)?;
    write_pressure_file(&mut writer, input, solution)?;
    writer.flush()?;

    let mut writer = create_file(CONNECTIONS_FILE)?;
    write_connections_file(&mut writer, input, solution)?;
    writer.flush()?;

    Ok(())
}

/// ASCII VTK UnstructuredGrid: the triangles with the pressure and volume of
/// each point, and the Darcy velocity and the input arrays of each triangle
fn write_pressure_file(
    writer: &mut dyn Write,
    input: &Input,
    solution: &Solution,
) -> std::io::Result<()> {
    let mesh = &input.mesh;
    let velocities = calculate_darcy_velocities(input, &solution.pressures);

    writeln!(writer, r#"<?xml version="1.0"?>"#)?;
    writeln!(writer, r#"<VTKFile type="UnstructuredGrid" version="1.0">"#)?;
    writeln!(writer, "<UnstructuredGrid>")?;
    writeln!(
        writer,
        r#"<Piece NumberOfPoints="{}" NumberOfCells="{}">"#,
        mesh.points.len(),
        mesh.triangles.len()
    )?;
    write_points(writer, &mesh.points)?;

    writeln!(writer, "<Cells>")?;
    write_data_array(
        writer,
        "connectivity",
        "Int64",
        1,
        mesh.triangles.iter().flatten(),
    )?;
    write_data_array(
        writer,
        "offsets",
        "Int64",
        1,
        (1..=mesh.triangles.len()).map(|cell| 3 * cell),
    )?;
    write_data_array(
        writer,
        "types",
        "UInt8",
        1,
        mesh.triangles.iter().map(|_| 5),
    )?;
    writeln!(writer, "</Cells>")?;

    writeln!(writer, r#"<PointData Scalars="pressure">"#)?;
    write_data_array(
        writer,
        "pressure",
        "Float64",
        1,
        solution.pressures.iter().map(|&value| Exact(value)),
    )?;
    write_data_array(
        writer,
        "volume",
        "Float64",
        1,
        solution
            .control_volumes
            .volumes
            .iter()
            .map(|&value| Exact(value)),
    )?;
    writeln!(writer, "</PointData>")?;

    writeln!(writer, r#"<CellData Vectors="darcy_velocity">"#)?;
    write_data_array(
        writer,
        "darcy_velocity",
        "Float64",
        3,
        velocities.iter().flatten().map(|&value| Exact(value)),
    )?;
    write_data_array(writer, "fracture_id", "Int64", 1, mesh.fracture_ids.iter())?;
    write_data_array(
        writer,
        "aperture",
        "Float64",
        1,
        mesh.apertures.iter().map(|&value| Exact(value)),
    )?;
    write_data_array(
        writer,
        "permeability",
        "Float64",
        1,
        mesh.permeabilities.iter().map(|&value| Exact(value)),
    )?;
    writeln!(writer, "</CellData>")?;

    writeln!(writer, "</Piece>")?;
    writeln!(writer, "</UnstructuredGrid>")?;
    writeln!(writer, "</VTKFile>")
}

/// ASCII VTK PolyData: one line for each connection, with its mass flux and transmissibility
fn write_connections_file(
    writer: &mut dyn Write,
    input: &Input,
    solution: &Solution,
) -> std::io::Result<()> {
    let connections = &solution.control_volumes.connections;

    writeln!(writer, r#"<?xml version="1.0"?>"#)?;
    writeln!(writer, r#"<VTKFile type="PolyData" version="1.0">"#)?;
    writeln!(writer, "<PolyData>")?;
    writeln!(
        writer,
        r#"<Piece NumberOfPoints="{}" NumberOfVerts="0" NumberOfLines="{}" NumberOfStrips="0" NumberOfPolys="0">"#,
        input.mesh.points.len(),
        connections.len()
    )?;
    write_points(writer, &input.mesh.points)?;

    writeln!(writer, "<Lines>")?;
    write_data_array(
        writer,
        "connectivity",
        "Int64",
        1,
        connections.iter().flat_map(|connection| connection.nodes),
    )?;
    write_data_array(
        writer,
        "offsets",
        "Int64",
        1,
        (1..=connections.len()).map(|line| 2 * line),
    )?;
    writeln!(writer, "</Lines>")?;

    writeln!(writer, r#"<CellData Scalars="mass_flux">"#)?;
    write_data_array(
        writer,
        "mass_flux",
        "Float64",
        1,
        solution.mass_fluxes.iter().map(|&value| Exact(value)),
    )?;
    write_data_array(
        writer,
        "transmissibility",
        "Float64",
        1,
        connections
            .iter()
            .map(|connection| Exact(connection.transmissibility)),
    )?;
    writeln!(writer, "</CellData>")?;

    writeln!(writer, "</Piece>")?;
    writeln!(writer, "</PolyData>")?;
    writeln!(writer, "</VTKFile>")
}

fn write_points(writer: &mut dyn Write, points: &[[f64; 3]]) -> std::io::Result<()> {
    writeln!(writer, "<Points>")?;
    write_data_array(
        writer,
        "Points",
        "Float64",
        3,
        points.iter().flatten().map(|&value| Exact(value)),
    )?;
    writeln!(writer, "</Points>")
}

/// A float which is written with the shortest text that reads back to the same value
struct Exact(f64);

impl Display for Exact {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{:e}", self.0)
    }
}

/// Writes an ASCII `DataArray`, one tuple on each line
fn write_data_array<T: Display>(
    writer: &mut dyn Write,
    name: &str,
    data_type: &str,
    component_count: usize,
    values: impl Iterator<Item = T>,
) -> std::io::Result<()> {
    writeln!(
        writer,
        r#"<DataArray type="{data_type}" Name="{name}" NumberOfComponents="{component_count}" format="ascii">"#
    )?;
    for (index, value) in values.enumerate() {
        let separator = if (index + 1) % component_count == 0 {
            "\n"
        } else {
            " "
        };
        write!(writer, "{value}{separator}")?;
    }
    writeln!(writer, "</DataArray>")
}
