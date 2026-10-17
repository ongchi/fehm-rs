#[derive(thiserror::Error, Debug)]
pub enum FehmError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("Invalid input file: {0}")]
    InvalidInput(String),
    #[error("Invalid mesh: {0}")]
    InvalidMesh(String),
    #[error(
        "The fractures {0} have no connection to a boundary, thus their pressure has no unique value. Add a boundary node to them, or remove them from the mesh."
    )]
    NoBoundary(String),
    #[error(
        "The Newton iterations did not converge in {0} iterations. Increase solver.max_iterations or solver.tolerance."
    )]
    NoConvergence(usize),
    #[error("Linear solver: {0}")]
    LinearSolver(String),
}
