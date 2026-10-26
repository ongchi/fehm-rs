use std::path::{Path, PathBuf};

use serde::Deserialize;

use tracing::warn;

use crate::FehmError;
use crate::eos;
use crate::io::vtu::parse_mesh;
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
    toml::from_str(input_text).map_err(|e| FehmError::InvalidToml(e.to_string()))
}

/// Builds the input from the text of the input file and of the mesh file
pub fn build_input_from_text(input_text: &str, mesh_text: &str) -> Result<Input, FehmError> {
    let input_file = parse_input_file(input_text)?;
    let mesh = parse_mesh(mesh_text).map_err(FehmError::InvalidMesh)?;

    if !(eos::MIN_TEMPERATURE..=eos::MAX_TEMPERATURE).contains(&input_file.temperature) {
        return Err(FehmError::TemperatureRange {
            temperature: input_file.temperature,
        });
    }
    if input_file.boundary.is_empty() {
        return Err(FehmError::NoBoundaryDefined);
    }

    let mut boundaries: Vec<PressureBoundary> = Vec::new();

    for section in input_file.boundary {
        if boundaries
            .iter()
            .any(|boundary| boundary.name == section.name)
        {
            return Err(FehmError::DuplicateBoundary {
                boundary: section.name,
            });
        }
        if !(eos::MIN_PRESSURE..=eos::MAX_PRESSURE).contains(&section.pressure) {
            return Err(FehmError::PressureRange {
                boundary: section.name,
                pressure: section.pressure,
            });
        }

        let Some(boundary_nodes) = mesh
            .boundaries
            .iter()
            .find(|boundary_nodes| boundary_nodes.name == section.name)
        else {
            return Err(FehmError::UnknownBoundary {
                boundary: section.name,
            });
        };
        if boundary_nodes.nodes.is_empty() {
            return Err(FehmError::EmptyBoundary {
                boundary: section.name,
            });
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
        return Err(FehmError::SolverTolerance {
            tolerance: solver.tolerance,
        });
    }
    if solver.max_iterations == 0 {
        return Err(FehmError::SolverMaxIterations);
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
            return Err(FehmError::BoundaryWithoutEffect {
                boundary: boundaries[index].name.clone(),
            });
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
        input_error(input_text).to_string()
    }

    fn input_error(input_text: &str) -> FehmError {
        build_input_from_text(input_text, MESH).unwrap_err()
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
        let error = input_error(&VALID_INPUT.replace("\"right\"", "\"top\""));
        assert!(
            matches!(&error, FehmError::UnknownBoundary { boundary } if boundary == "top"),
            "{error}"
        );
    }

    #[test]
    fn empty_boundary_is_an_error() {
        let error = input_error(&VALID_INPUT.replace("\"right\"", "\"none\""));
        assert!(
            matches!(&error, FehmError::EmptyBoundary { boundary } if boundary == "none"),
            "{error}"
        );
    }

    #[test]
    fn no_boundary_is_an_error() {
        let error = input_error("mesh = \"mesh.vtu\"\ntemperature = 50.0\n");
        assert!(matches!(error, FehmError::NoBoundaryDefined), "{error}");
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
        let error = input_error(&VALID_INPUT.replace("\"right\"", "\"corner\""));
        assert!(
            matches!(&error, FehmError::BoundaryWithoutEffect { boundary } if boundary == "corner"),
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
        let error = input_error(&VALID_INPUT.replace("50.0", "400.0"));
        assert!(
            matches!(error, FehmError::TemperatureRange { temperature } if temperature == 400.0),
            "{error}"
        );
    }

    #[test]
    fn temperature_error_names_the_range() {
        let error = error_of(&VALID_INPUT.replace("50.0", "400.0"));
        assert!(
            error.contains("temperature = 400 °C is outside the range of the water properties of FEHM (15 to 360 °C)"),
            "{error}"
        );
    }

    #[test]
    fn two_boundaries_with_one_name_are_an_error() {
        let error = input_error(&VALID_INPUT.replace("\"right\"", "\"left\""));
        assert!(
            matches!(&error, FehmError::DuplicateBoundary { boundary } if boundary == "left"),
            "{error}"
        );
    }

    #[test]
    fn pressure_outside_the_eos_range_is_an_error() {
        let error = input_error(&VALID_INPUT.replace("10.1", "200.0"));
        assert!(
            matches!(&error, FehmError::PressureRange { boundary, pressure } if boundary == "left" && *pressure == 200.0),
            "{error}"
        );
    }

    #[test]
    fn solver_tolerance_of_one_is_an_error() {
        let error = input_error(&format!("{VALID_INPUT}\n[solver]\ntolerance = 1.0\n"));
        assert!(
            matches!(error, FehmError::SolverTolerance { tolerance } if tolerance == 1.0),
            "{error}"
        );
    }

    #[test]
    fn zero_solver_iterations_are_an_error() {
        let error = input_error(&format!("{VALID_INPUT}\n[solver]\nmax_iterations = 0\n"));
        assert!(matches!(error, FehmError::SolverMaxIterations), "{error}");
    }
}
