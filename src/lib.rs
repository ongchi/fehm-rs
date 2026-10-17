pub mod coefficients;
pub mod eos;
pub mod error;
pub mod flow;
pub mod io;
pub mod mesh;
pub mod velocity;

pub use error::FehmError;
pub use flow::{Solution, solve};
pub use io::input::{Input, build_input_from_text, read_input, read_mesh_file_name};
pub use io::output::{write_output, write_output_to_folder};
