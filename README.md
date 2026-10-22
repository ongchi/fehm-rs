# fehm-rs

Steady single-phase water flow in a discrete fracture network (DFN), in Rust.
It replaces the FEHM step of the [dfnWorks](https://github.com/lanl/dfnWorks)
work flow for networks from [dfngen-rs](../dfngen-rs). The equations and the
water properties come from [FEHM](https://github.com/lanl/FEHM) (in
`../FEHM`); only the isothermal liquid flow of FEHM is ported.

## Usage

```bash
cargo run --release -- <input_file.toml> <output_folder>
```

Example: `cargo run --release -- tests/data/crossing_fractures.toml output`.

## Input file

One TOML file and one mesh file. An unknown key is an error.

```toml
mesh = "dfn_mesh.vtu"           # the mesh file; a relative path is relative to the folder of this file
temperature = 50.0              # °C, 15 to 360; the same temperature everywhere
gravity = true                  # optional, default true; gravity is in the -z direction

# One table for each boundary with a fixed pressure (at least one)
[[boundary]]
name = "inflow"                 # the point data array boundary__inflow of the mesh
pressure = 10.1                 # MPa, 0.001 to 110

[[boundary]]
name = "outflow"
pressure = 10.0

[solver]                        # optional section
tolerance = 1e-10               # optional, default 1e-10; see "Solution"
max_iterations = 20             # optional, default 20; Newton iterations
```

## Mesh file

An ASCII VTK UnstructuredGrid file (`.vtu`) with one `Piece`:

- All cells are triangles (VTK type 5). Each point must be a point of a triangle.
- Fractures which intersect must share the points on their intersection
  line. Points at the same position which are not the same point are not
  connected.
- Cell data arrays (necessary), one value for each triangle:
  - `fracture_id` (integer): the fracture of the triangle. It is used only in messages and in the output.
  - `aperture` (m, above 0)
  - `permeability` (m², above 0); for example the cubic law, b²/12.
- Point data arrays `boundary__<name>` (optional): a point with a value
  above 0 is on the boundary `<name>`. The input file gives the pressure of
  each boundary which it uses; the other arrays are ignored. A point on two
  boundaries with different pressures keeps the pressure of the first of them
  in the input file (a warning gives the number of such points); a boundary
  whose points are all in earlier boundaries with another pressure is an
  error.

Each group of connected fractures must have a point on a boundary; if not,
its pressure has no unique value, and the run stops with an error which
gives its `fracture_id` values.

## Solution

For each triangle, in its own plane, the connection coefficient of an edge
is ½·cot of the opposite angle, and each point gets its Voronoi area of the
triangle (the 2-D Voronoi coefficients of LaGriT, `matbld2d_stor`). The
transmissibility of a connection is the sum over its triangles of
coefficient · aperture · permeability (m³); the volume of a point is the sum
of Voronoi area · aperture.

The mass flux from point i to point j is

    F = T · 10⁶ · (ρ/μ)_upstream · ((P_i − P_j) + ½(ρ_i + ρ_j) · g · (z_i − z_j))

with P in MPa, g = 9.81 m/s², ρ(P, T) and μ(P, T) of FEHM (liquid water), and
ρ/μ of the point with the higher potential. The steady state (the sum of the
fluxes of each point without a fixed pressure is 0) is solved with Newton
iterations and a sparse LU decomposition (`faer`). The iterations stop when
the largest residual of a point is below `tolerance` times the largest sum of
the flux magnitudes of a point, or when the pressure change is at the
rounding level.

An obtuse triangle gives a negative coefficient and a negative area; the run
writes a warning with their number. A Delaunay mesh of each fracture has
none.

## Output files

All files are in `<output_folder>`.

| File | Content |
|------|---------|
| `pressure.vtu` | The mesh. Point data: `pressure` (MPa), `volume` (m³). Cell data: `darcy_velocity` (m/s, a vector in the plane of the triangle, from the pressure gradient in the triangle), `fracture_id`, `aperture`, `permeability`. It can be the mesh of a new run (the boundary arrays are not in it). |
| `connections.vtp` | One line for each connection. Cell data: `mass_flux` (kg/s, positive from the first point of the line to the second), `transmissibility` (m³). |

The log gives the mass flow into the network through each boundary and the
mass balance.

## Library

The crate is a library (`fehm`) and a binary. The library compiles for
`wasm32-unknown-unknown`: `build_input_from_text` takes the text of the files,
and `write_output` gets a writer for each output file. `read_input` and
`write_output_to_folder` are the file-system wrappers.

```rust
let input = fehm::build_input_from_text(&toml_text, &mesh_text)?;
let solution = fehm::solve(&input)?;
```

## Comparison with FEHM

`tests/fehm_reference.rs` compares the pressures with those of the Fortran
FEHM V3.6.3 for three cases (largest difference below 1e-8 of the pressure
range). To repeat the FEHM runs, build FEHM from `../FEHM/src` and run

```bash
scripts/compare_with_fehm.sh tests/data/crossing_fractures.toml <work_folder> <path to xfehm>
```

`examples/export_fehm.rs` writes the FEHM input of a case: the coefficients
of fehm-rs in a `.stor` file, each node held at the temperature (`hflx`), and
the boundaries as `flow` nodes with a large impedance.

## Test data

`scripts/make_test_meshes.py tests/data` makes the small meshes (one
fracture, two crossing fractures, and a fracture with no boundary).
`tests/data/lagrit_dfn` is the DFN mesh of the LaGriT test
`level01/recon1` with its `.stor` file, converted by `scripts/avs_to_vtu.py`.

## License

BSD 3-Clause, as FEHM. FEHM: LANL LA-CC-2012-083. LaGriT test data:
LANL LA-CC-15-069 (BSD 3-Clause).
