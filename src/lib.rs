pub mod coefficients;
pub mod eos;
pub mod error;
pub mod flow;
pub mod io;
pub mod mesh;
pub mod velocity;
pub mod voxel;

pub use error::FehmError;
pub use flow::{FlowProblem, NodeSolution, Solution, solve, solve_control_volumes};
pub use io::input::{
    Input, SolverOptions, build_input_from_mesh, build_input_from_text, read_input,
    read_mesh_file_name,
};
pub use io::output::{write_output, write_output_to_folder};
