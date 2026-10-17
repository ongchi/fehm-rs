//! The output files read back with the same values.

use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::rc::Rc;

use fehm::io::vtu::parse_mesh;
use fehm::{build_input_from_text, solve, write_output};

const INPUT: &str = include_str!("data/crossing_fractures.toml");
const MESH: &str = include_str!("data/crossing_fractures.vtu");

type Buffer = Rc<RefCell<Vec<u8>>>;

/// A writer which keeps the bytes in a shared buffer
struct SharedBuffer(Buffer);

impl Write for SharedBuffer {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0.borrow_mut().write(bytes)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The values of the `DataArray` with the attribute `Name="name"`
fn read_data_array(text: &str, name: &str) -> Vec<f64> {
    let document = roxmltree::Document::parse(text).unwrap();
    let data_array = document
        .descendants()
        .find(|node| node.has_tag_name("DataArray") && node.attribute("Name") == Some(name))
        .unwrap();
    data_array
        .text()
        .unwrap()
        .split_whitespace()
        .map(|value| value.parse().unwrap())
        .collect()
}

#[test]
fn output_files_read_back_with_the_same_values() {
    let input = build_input_from_text(INPUT, MESH).unwrap();
    let solution = solve(&input).unwrap();

    let mut files: HashMap<String, Buffer> = HashMap::new();
    write_output(&input, &solution, &mut |name| {
        let buffer = Buffer::default();
        files.insert(name.to_string(), buffer.clone());
        Ok(Box::new(SharedBuffer(buffer)))
    })
    .unwrap();

    let pressure_text = String::from_utf8(files["pressure.vtu"].borrow().clone()).unwrap();
    let connections_text = String::from_utf8(files["connections.vtp"].borrow().clone()).unwrap();

    // pressure.vtu is a valid mesh for a new run
    let mesh = parse_mesh(&pressure_text).unwrap();
    assert_eq!(mesh.points, input.mesh.points);
    assert_eq!(mesh.triangles, input.mesh.triangles);
    assert_eq!(mesh.permeabilities, input.mesh.permeabilities);

    assert_eq!(
        read_data_array(&pressure_text, "pressure"),
        solution.pressures
    );
    assert_eq!(
        read_data_array(&pressure_text, "darcy_velocity").len(),
        3 * input.mesh.triangles.len()
    );
    assert_eq!(
        read_data_array(&connections_text, "mass_flux"),
        solution.mass_fluxes
    );
}
