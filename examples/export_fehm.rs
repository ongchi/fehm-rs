//! Writes a fehm-rs case as an input for the Fortran FEHM, for the
//! comparison of the two programs.
//!
//! Usage: `cargo run --example export_fehm -- <input.toml> <output_folder>`,
//! then `xfehm fehmn.files` in the output folder.
//!
//! The `.stor` file has the coefficients of fehm-rs: area coefficient =
//! transmissibility / REFERENCE_PERMEABILITY, and all nodes have the
//! permeability REFERENCE_PERMEABILITY. Thus FEHM uses the same
//! transmissibility for each connection, also where fractures with different
//! properties intersect. The run is heat and mass at one temperature, as in
//! dfnWorks, with time steps to a steady state. Each boundary is a `flow`
//! node with a large impedance.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use fehm::coefficients::{ControlVolumes, calculate_control_volumes};
use fehm::{Input, read_input};

const REFERENCE_PERMEABILITY: f64 = 1e-12;
/// Impedance of the boundary nodes (FEHM: kg/s/MPa · 10⁻⁶)
const BOUNDARY_IMPEDANCE: f64 = 1e3;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments: Vec<String> = std::env::args().collect();
    let [_, input_file, output_folder] = arguments.as_slice() else {
        return Err("usage: export_fehm <input.toml> <output_folder>".into());
    };

    let input = read_input(Path::new(input_file))?;
    let control_volumes = calculate_control_volumes(&input.mesh)?;
    let folder = PathBuf::from(output_folder);
    std::fs::create_dir_all(&folder)?;

    write_file(&folder, "fehmn.files", |writer| write_control_file(writer))?;
    write_file(&folder, "mesh.fehmn", |writer| write_grid(writer, &input))?;
    write_file(&folder, "mesh.stor", |writer| {
        write_stor(writer, &control_volumes)
    })?;
    write_file(&folder, "run.dat", |writer| write_macros(writer, &input))?;

    Ok(())
}

fn write_file(
    folder: &Path,
    name: &str,
    write: impl FnOnce(&mut dyn Write) -> std::io::Result<()>,
) -> std::io::Result<()> {
    let mut writer = BufWriter::new(File::create(folder.join(name))?);
    write(&mut writer)?;
    writer.flush()
}

fn write_control_file(writer: &mut dyn Write) -> std::io::Result<()> {
    writeln!(writer, "input: run.dat")?;
    writeln!(writer, "gridf: mesh.fehmn")?;
    writeln!(writer, "outp:  run.out")?;
    writeln!(writer, "rsto:  run.fin")?;
    writeln!(writer, "check: run.chk")?;
    writeln!(writer, "root:  run")?;
    writeln!(writer, "stor:  mesh.stor")?;
    writeln!(writer)?;
    writeln!(writer, "all")?;
    writeln!(writer, "0")
}

fn write_grid(writer: &mut dyn Write, input: &Input) -> std::io::Result<()> {
    let mesh = &input.mesh;
    writeln!(writer, "coor")?;
    writeln!(writer, "{:10}", mesh.points.len())?;
    for (index, [x, y, z]) in mesh.points.iter().enumerate() {
        writeln!(writer, "{:10} {x:e} {y:e} {z:e}", index + 1)?;
    }
    // A blank line ends the coordinates
    writeln!(writer)?;
    writeln!(writer, "elem")?;
    writeln!(writer, "{:10}{:10}", 3, mesh.triangles.len())?;
    for (index, [a, b, c]) in mesh.triangles.iter().enumerate() {
        writeln!(
            writer,
            "{:10}{:10}{:10}{:10}",
            index + 1,
            a + 1,
            b + 1,
            c + 1
        )?;
    }
    writeln!(writer)?;
    writeln!(writer, "stop")
}

/// ASCII `.stor` file in the LaGriT format (one area coefficient for each matrix entry)
fn write_stor(writer: &mut dyn Write, control_volumes: &ControlVolumes) -> std::io::Result<()> {
    let node_count = control_volumes.volumes.len();

    // The entries of each row: (column, area coefficient), with the diagonal
    let mut rows: Vec<Vec<(usize, f64)>> = (0..node_count).map(|node| vec![(node, 0.0)]).collect();
    for connection in &control_volumes.connections {
        let [first, second] = connection.nodes;
        // FEHM (and LaGriT) use −(area/length)
        let coefficient = -connection.transmissibility / REFERENCE_PERMEABILITY;
        rows[first].push((second, coefficient));
        rows[second].push((first, coefficient));
    }
    for row in &mut rows {
        row.sort_by_key(|&(column, _)| column);
    }

    let entry_count: usize = rows.iter().map(Vec::len).sum();
    let largest_row = rows.iter().map(Vec::len).max().unwrap_or(0);
    let pointer_count = node_count + 1 + entry_count;

    writeln!(
        writer,
        "fehmstor ascir8i4 fehm-rs Sparse Matrix Voronoi Coefficients"
    )?;
    writeln!(
        writer,
        " fehm-rs export_fehm 2-D Linear Diffusion Model (matbld2d_stor)"
    )?;
    writeln!(
        writer,
        "{entry_count:10}{node_count:10}{pointer_count:10}{:10}{largest_row:10}",
        1
    )?;

    write_values(
        writer,
        control_volumes.volumes.iter().map(|v| format!("{v:e}")),
    )?;

    // Row pointers: the entries of row i are at nelm(i)+1 ..= nelm(i+1)
    let mut row_pointers = vec![node_count + 1];
    for row in &rows {
        row_pointers.push(row_pointers.last().unwrap() + row.len());
    }
    write_values(writer, row_pointers.iter().map(usize::to_string))?;
    write_values(
        writer,
        rows.iter()
            .flatten()
            .map(|(column, _)| (column + 1).to_string()),
    )?;
    // Each entry has its own coefficient; then node_count + 1 zeros
    write_values(
        writer,
        (1..=entry_count)
            .chain(std::iter::repeat_n(0, node_count + 1))
            .map(|pointer| pointer.to_string()),
    )?;
    // Position of the diagonal entry of each row in the combined array
    write_values(
        writer,
        rows.iter().enumerate().map(|(node, row)| {
            let offset = row.iter().position(|&(column, _)| column == node).unwrap();
            (row_pointers[node] + offset + 1).to_string()
        }),
    )?;
    write_values(
        writer,
        rows.iter()
            .flatten()
            .map(|(_, coefficient)| format!("{coefficient:e}")),
    )
}

/// Five values on each line
fn write_values(
    writer: &mut dyn Write,
    values: impl Iterator<Item = String>,
) -> std::io::Result<()> {
    let mut count = 0;
    for value in values {
        write!(writer, " {value}")?;
        count += 1;
        if count % 5 == 0 {
            writeln!(writer)?;
        }
    }
    if count % 5 != 0 {
        writeln!(writer)?;
    }
    Ok(())
}

fn write_macros(writer: &mut dyn Write, input: &Input) -> std::io::Result<()> {
    let temperature = input.temperature;
    let initial_pressure = input
        .boundaries
        .iter()
        .map(|boundary| boundary.pressure)
        .sum::<f64>()
        / input.boundaries.len() as f64;
    let gravity_direction = if input.has_gravity { 3 } else { 0 };

    writeln!(writer, "fehm-rs export_fehm")?;
    writeln!(writer, "pres")?;
    writeln!(writer, "1 0 0 {initial_pressure:e} {temperature:e} 1")?;
    writeln!(writer)?;

    // One zone for each boundary; a node in two boundaries is in the first one
    writeln!(writer, "zone")?;
    let mut is_assigned = vec![false; input.mesh.points.len()];
    for (index, boundary) in input.boundaries.iter().enumerate() {
        let nodes: Vec<usize> = boundary
            .nodes
            .iter()
            .copied()
            .filter(|&node| !std::mem::replace(&mut is_assigned[node], true))
            .collect();
        writeln!(writer, "{}", index + 1)?;
        writeln!(writer, "nnum")?;
        writeln!(writer, "{}", nodes.len())?;
        write_values(writer, nodes.iter().map(|node| (node + 1).to_string()))?;
    }
    writeln!(writer)?;

    writeln!(writer, "perm")?;
    writeln!(
        writer,
        "1 0 0 {REFERENCE_PERMEABILITY:e} {REFERENCE_PERMEABILITY:e} {REFERENCE_PERMEABILITY:e}"
    )?;
    writeln!(writer)?;
    writeln!(writer, "rock")?;
    writeln!(writer, "1 0 0 2700. 1000. 0.1")?;
    writeln!(writer)?;
    writeln!(writer, "cond")?;
    writeln!(writer, "1 0 0 1. 1. 1.")?;
    writeln!(writer)?;

    // Hold each node at the temperature: heat flow = 1e6 MW/°C · (T − temperature)
    writeln!(writer, "hflx")?;
    writeln!(writer, "1 0 0 {temperature:e} 1.e6")?;
    writeln!(writer)?;

    // Fixed pressure and fixed temperature of the in-flow
    writeln!(writer, "flow")?;
    for (index, boundary) in input.boundaries.iter().enumerate() {
        writeln!(
            writer,
            "-{} 0 0 {:e} {:e} {BOUNDARY_IMPEDANCE:e}",
            index + 1,
            boundary.pressure,
            -temperature
        )?;
    }
    writeln!(writer)?;

    // Time steps (days) to a steady state
    writeln!(writer, "time")?;
    writeln!(writer, "1.e-4 1.e4 10000 1 2026 1")?;
    writeln!(writer)?;
    writeln!(writer, "ctrl")?;
    writeln!(writer, "-80 1.e-10 8 100 gmres")?;
    writeln!(writer, "1 0 0 4")?;
    writeln!(writer)?;
    writeln!(writer, "1.0 {gravity_direction} 1.")?;
    writeln!(writer, "80 2.0 1.e-8 1.e4")?;
    writeln!(writer, "0 +1")?;
    writeln!(writer, "iter")?;
    writeln!(writer, "1.e-10 1.e-10 1.e-10 -1.e-10 1.2")?;
    writeln!(writer, "0 0 0 0 1.e6")?;
    writeln!(writer, "stop")
}
