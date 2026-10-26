use crate::eos::{MAX_PRESSURE, MAX_TEMPERATURE, MIN_PRESSURE, MIN_TEMPERATURE};
use crate::io::vtu::BOUNDARY_ARRAY_PREFIX;

#[derive(thiserror::Error, Debug)]
pub enum FehmError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// The TOML text: syntax, an unknown or a missing key, a wrong type
    #[error("Invalid input file: {0}")]
    InvalidToml(String),
    #[error(
        "Invalid input file: temperature = {temperature} °C is outside the range of the water properties of FEHM ({MIN_TEMPERATURE} to {MAX_TEMPERATURE} °C)"
    )]
    TemperatureRange { temperature: f64 },
    #[error(
        "Invalid input file: there is no [[boundary]]; at least one boundary with a fixed pressure is necessary"
    )]
    NoBoundaryDefined,
    #[error("Invalid input file: two [[boundary]] tables have the name \"{boundary}\"")]
    DuplicateBoundary { boundary: String },
    #[error(
        "Invalid input file: boundary \"{boundary}\": pressure = {pressure} MPa is outside the range of the water properties of FEHM ({MIN_PRESSURE} to {MAX_PRESSURE} MPa)"
    )]
    PressureRange { boundary: String, pressure: f64 },
    #[error(
        "Invalid input file: boundary \"{boundary}\": the mesh has no point data array {BOUNDARY_ARRAY_PREFIX}{boundary}"
    )]
    UnknownBoundary { boundary: String },
    #[error(
        "Invalid input file: boundary \"{boundary}\": the array {BOUNDARY_ARRAY_PREFIX}{boundary} has no point with a value above 0"
    )]
    EmptyBoundary { boundary: String },
    #[error(
        "Invalid input file: boundary \"{boundary}\": each of its points is in an earlier boundary with another pressure, thus it has no effect"
    )]
    BoundaryWithoutEffect { boundary: String },
    #[error("Invalid input file: solver.tolerance = {tolerance} must be above 0 and below 1")]
    SolverTolerance { tolerance: f64 },
    #[error("Invalid input file: solver.max_iterations must be 1 or more")]
    SolverMaxIterations,
    /// The content of the mesh file
    #[error("Invalid mesh: {0}")]
    InvalidMesh(String),
    /// `cell` starts at 0, in the order of the mesh file
    #[error("Invalid mesh: cell {cell}: the triangle has no area")]
    DegenerateTriangle { cell: usize },
    #[error(
        "The fractures {} have no connection to a boundary, thus their pressure has no unique value. Add a boundary node to them, or remove them from the mesh.",
        comma_list(.fracture_ids)
    )]
    NoBoundary { fracture_ids: Vec<i64> },
    #[error(
        "The Newton iterations did not converge in {0} iterations. Increase solver.max_iterations or solver.tolerance."
    )]
    NoConvergence(usize),
    #[error("Linear solver: {0}")]
    LinearSolver(String),
}

fn comma_list(values: &[i64]) -> String {
    values
        .iter()
        .map(i64::to_string)
        .collect::<Vec<_>>()
        .join(", ")
}
