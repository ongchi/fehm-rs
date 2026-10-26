//! Steady single-phase liquid flow: the liquid, isothermal part of the mass
//! equation of FEHM (`geneq2.f`), with no storage term.
//!
//! The mass flux from node `i` to node `j` is
//! `F = T · 10⁶ · (ρ/μ)_upstream · ((P_i − P_j) + ½(ρ_i + ρ_j) · g · (z_i − z_j))`,
//! with `T` the transmissibility of the connection (m³), `P` in MPa and
//! `g = 9.81·10⁻⁶` MPa·m²/kg. The upstream node is the node with the higher
//! potential; for a potential difference of 0 the mean of the two nodes is
//! used (`fid = 0.5` in `geneq2.f`).

use faer::Mat;
use faer::linalg::solvers::Solve;
use faer::sparse::{SparseColMat, Triplet};
use tracing::{debug, info, warn};

use crate::FehmError;
use crate::coefficients::{ControlVolumes, calculate_control_volumes};
use crate::eos;
use crate::io::input::Input;

/// Gravity (m/s²), as `ctrl` in FEHM (`inctrl.f`)
pub const GRAVITY: f64 = 9.81;
/// MPa → Pa
const PASCALS_PER_MEGAPASCAL: f64 = 1e6;

/// Result of a run
#[derive(Clone, Debug)]
pub struct Solution {
    pub control_volumes: ControlVolumes,
    /// For each node (MPa)
    pub pressures: Vec<f64>,
    /// For each connection of `control_volumes`: mass flux from `nodes[0]`
    /// to `nodes[1]` (kg/s)
    pub mass_fluxes: Vec<f64>,
    /// For each boundary of the input: mass flow into the network (kg/s);
    /// negative for an out-flow. A node in two boundaries counts for the
    /// first one.
    pub boundary_flows: Vec<f64>,
    pub newton_iterations: usize,
}

/// Water properties at each node
struct NodeProperties {
    densities: Vec<eos::Property>,
    /// ρ/μ (kg/m³/(Pa·s)) and its pressure derivative (per MPa)
    mobilities: Vec<eos::Property>,
}

impl NodeProperties {
    fn new(pressures: &[f64], temperature: f64) -> Self {
        let densities: Vec<eos::Property> = pressures
            .iter()
            .map(|&pressure| eos::density(pressure, temperature))
            .collect();
        let mobilities = pressures
            .iter()
            .zip(&densities)
            .map(|(&pressure, density)| {
                let viscosity = eos::viscosity(pressure, temperature);
                eos::Property {
                    value: density.value / viscosity.value,
                    pressure_derivative: (density.pressure_derivative * viscosity.value
                        - density.value * viscosity.pressure_derivative)
                        / (viscosity.value * viscosity.value),
                }
            })
            .collect();

        Self {
            densities,
            mobilities,
        }
    }
}

/// Mass flux of one connection and its derivatives with respect to the two pressures
struct ConnectionFlux {
    flux: f64,
    /// The flux with the absolute values of its pressure term and gravity
    /// term: the scale of the rounding error of `flux`, also at hydrostatic
    /// equilibrium, where `flux` is 0
    flux_scale: f64,
    derivative_first: f64,
    derivative_second: f64,
}

/// Solves the steady flow equations with Newton iterations.
pub fn solve(input: &Input) -> Result<Solution, FehmError> {
    let control_volumes = calculate_control_volumes(&input.mesh)?;
    log_negative_coefficients(&control_volumes);

    let fixed_pressures = fixed_pressures(input);
    check_every_cluster_has_a_boundary(input, &control_volumes, &fixed_pressures)?;

    // Number of each node in the unknowns; `None` for a node with a fixed pressure
    let mut unknown_count = 0;
    let unknown_of_node: Vec<Option<usize>> = fixed_pressures
        .iter()
        .map(|fixed| {
            fixed.is_none().then(|| {
                unknown_count += 1;
                unknown_count - 1
            })
        })
        .collect();

    let mean_boundary_pressure = input
        .boundaries
        .iter()
        .map(|boundary| boundary.pressure)
        .sum::<f64>()
        / input.boundaries.len() as f64;
    let mut pressures: Vec<f64> = fixed_pressures
        .iter()
        .map(|fixed| fixed.unwrap_or(mean_boundary_pressure))
        .collect();
    let gravity = if input.has_gravity { GRAVITY } else { 0.0 };
    let largest_pressure = pressures.iter().fold(0.0, |max: f64, p| max.max(p.abs()));
    // Largest pressure change of the last Newton iteration (MPa)
    let mut largest_change = f64::INFINITY;

    for iteration in 0..=input.solver.max_iterations {
        let properties = NodeProperties::new(&pressures, input.temperature);
        let mut residuals = vec![0.0; unknown_count];
        // Sum of the flux scales of the connections of each unknown node: the scale of its residual
        let mut flux_scales = vec![0.0; unknown_count];
        let mut jacobian = Vec::new();

        for connection in &control_volumes.connections {
            let [first, second] = connection.nodes;
            let connection_flux = calculate_connection_flux(
                connection.transmissibility,
                [first, second],
                &pressures,
                &properties,
                input.mesh.points[first][2] - input.mesh.points[second][2],
                gravity,
            );

            for (node, sign) in [(first, 1.0), (second, -1.0)] {
                let Some(row) = unknown_of_node[node] else {
                    continue;
                };
                residuals[row] += sign * connection_flux.flux;
                flux_scales[row] += connection_flux.flux_scale;

                for (column_node, derivative) in [
                    (first, connection_flux.derivative_first),
                    (second, connection_flux.derivative_second),
                ] {
                    if let Some(column) = unknown_of_node[column_node] {
                        jacobian.push(Triplet::new(row, column, sign * derivative));
                    }
                }
            }
        }

        let largest_residual = residuals.iter().fold(0.0, |max: f64, r| max.max(r.abs()));
        let largest_scale = flux_scales.iter().fold(0.0, |max: f64, s| max.max(*s));
        debug!(
            "Newton iteration {iteration}: largest residual {largest_residual:e} kg/s, largest flux scale {largest_scale:e} kg/s, largest pressure change {largest_change:e} MPa"
        );

        // A pressure change at the rounding level of the pressures: no further iteration can decrease the residual
        let is_change_at_rounding_level = largest_change <= 4.0 * f64::EPSILON * largest_pressure;
        if largest_residual <= input.solver.tolerance * largest_scale || is_change_at_rounding_level
        {
            info!("Newton iterations: {iteration}");
            let mass_fluxes = calculate_mass_fluxes(input, &control_volumes, &pressures, gravity);
            let boundary_flows = calculate_boundary_flows(
                input,
                &control_volumes,
                &mass_fluxes,
                unknown_of_node.len(),
            );
            return Ok(Solution {
                control_volumes,
                pressures,
                mass_fluxes,
                boundary_flows,
                newton_iterations: iteration,
            });
        }
        if iteration == input.solver.max_iterations {
            break;
        }

        let pressure_changes = solve_linear_system(unknown_count, &jacobian, &residuals)?;
        for (node, unknown) in unknown_of_node.iter().enumerate() {
            if let Some(unknown) = *unknown {
                pressures[node] -= pressure_changes[unknown];
            }
        }
        largest_change = pressure_changes
            .iter()
            .fold(0.0, |max: f64, c| max.max(c.abs()));
    }

    Err(FehmError::NoConvergence(input.solver.max_iterations))
}

/// For each node: its fixed pressure, if it is in a boundary
fn fixed_pressures(input: &Input) -> Vec<Option<f64>> {
    let mut fixed_pressures = vec![None; input.mesh.points.len()];
    for boundary in &input.boundaries {
        for &node in &boundary.nodes {
            fixed_pressures[node] = Some(boundary.pressure);
        }
    }
    fixed_pressures
}

/// Mass flux from `nodes[0]` to `nodes[1]` (kg/s) and its pressure derivatives (per MPa)
fn calculate_connection_flux(
    transmissibility: f64,
    nodes: [usize; 2],
    pressures: &[f64],
    properties: &NodeProperties,
    elevation_difference: f64,
    gravity: f64,
) -> ConnectionFlux {
    let [first, second] = nodes;
    let density = [first, second].map(|node| properties.densities[node]);
    let mobility = [first, second].map(|node| properties.mobilities[node]);
    let gravity_term = 0.5 * gravity * elevation_difference / PASCALS_PER_MEGAPASCAL;
    let conductance = transmissibility * PASCALS_PER_MEGAPASCAL;

    let pressure_difference = pressures[first] - pressures[second];
    let gravity_difference = (density[0].value + density[1].value) * gravity_term;
    let potential_difference = pressure_difference + gravity_difference;
    let potential_derivative_first = 1.0 + density[0].pressure_derivative * gravity_term;
    let potential_derivative_second = -1.0 + density[1].pressure_derivative * gravity_term;

    // Upstream weights of the two nodes (`fid` in `geneq2.f`)
    let weights = if potential_difference > 0.0 {
        [1.0, 0.0]
    } else if potential_difference < 0.0 {
        [0.0, 1.0]
    } else {
        [0.5, 0.5]
    };
    let upstream_mobility = weights[0] * mobility[0].value + weights[1] * mobility[1].value;

    ConnectionFlux {
        flux: conductance * upstream_mobility * potential_difference,
        flux_scale: (conductance * upstream_mobility).abs()
            * (pressure_difference.abs() + gravity_difference.abs()),
        derivative_first: conductance
            * (upstream_mobility * potential_derivative_first
                + weights[0] * mobility[0].pressure_derivative * potential_difference),
        derivative_second: conductance
            * (upstream_mobility * potential_derivative_second
                + weights[1] * mobility[1].pressure_derivative * potential_difference),
    }
}

/// Mass flux of each connection (kg/s)
fn calculate_mass_fluxes(
    input: &Input,
    control_volumes: &ControlVolumes,
    pressures: &[f64],
    gravity: f64,
) -> Vec<f64> {
    let properties = NodeProperties::new(pressures, input.temperature);

    control_volumes
        .connections
        .iter()
        .map(|connection| {
            let [first, second] = connection.nodes;
            calculate_connection_flux(
                connection.transmissibility,
                connection.nodes,
                pressures,
                &properties,
                input.mesh.points[first][2] - input.mesh.points[second][2],
                gravity,
            )
            .flux
        })
        .collect()
}

/// Mass flow into the network through each boundary (kg/s)
fn calculate_boundary_flows(
    input: &Input,
    control_volumes: &ControlVolumes,
    mass_fluxes: &[f64],
    point_count: usize,
) -> Vec<f64> {
    // Net mass flux out of each node into the network
    let mut outflows = vec![0.0; point_count];
    for (connection, flux) in control_volumes.connections.iter().zip(mass_fluxes) {
        outflows[connection.nodes[0]] += flux;
        outflows[connection.nodes[1]] -= flux;
    }

    let mut is_counted = vec![false; point_count];
    input
        .boundaries
        .iter()
        .map(|boundary| {
            boundary
                .nodes
                .iter()
                .filter(|&&node| !std::mem::replace(&mut is_counted[node], true))
                .map(|&node| outflows[node])
                .sum()
        })
        .collect()
}

/// Solves `jacobian · x = residuals` with a sparse LU decomposition
fn solve_linear_system(
    size: usize,
    jacobian: &[Triplet<usize, usize, f64>],
    residuals: &[f64],
) -> Result<Vec<f64>, FehmError> {
    let matrix = SparseColMat::<usize, f64>::try_new_from_triplets(size, size, jacobian)
        .map_err(|e| FehmError::LinearSolver(format!("{e:?}")))?;
    let lu = matrix
        .sp_lu()
        .map_err(|e| FehmError::LinearSolver(e.to_string()))?;
    let right_hand_side = Mat::from_fn(size, 1, |row, _| residuals[row]);
    let solution = lu.solve(&right_hand_side);
    let solution: Vec<f64> = (0..size).map(|row| solution[(row, 0)]).collect();

    if solution.iter().any(|value| !value.is_finite()) {
        return Err(FehmError::LinearSolver(
            "the solution is not finite; the matrix is singular".to_string(),
        ));
    }

    Ok(solution)
}

fn log_negative_coefficients(control_volumes: &ControlVolumes) {
    let negative_connections = control_volumes.negative_connection_count();
    let negative_volumes = control_volumes.negative_volume_count();

    if negative_connections > 0 || negative_volumes > 0 {
        warn!(
            "Obtuse triangles give {negative_connections} connections with a negative coefficient and {negative_volumes} nodes with a negative volume. A Delaunay mesh in each fracture has none."
        );
    }
}

/// Each group of connected nodes must have a node with a fixed pressure:
/// else its pressure has no unique value.
fn check_every_cluster_has_a_boundary(
    input: &Input,
    control_volumes: &ControlVolumes,
    fixed_pressures: &[Option<f64>],
) -> Result<(), FehmError> {
    let point_count = input.mesh.points.len();
    let mut parents: Vec<usize> = (0..point_count).collect();

    fn find_root(parents: &mut [usize], mut node: usize) -> usize {
        while parents[node] != node {
            parents[node] = parents[parents[node]];
            node = parents[node];
        }
        node
    }

    for connection in &control_volumes.connections {
        if connection.transmissibility == 0.0 {
            continue;
        }
        let first = find_root(&mut parents, connection.nodes[0]);
        let second = find_root(&mut parents, connection.nodes[1]);
        parents[first.max(second)] = first.min(second);
    }

    let mut has_boundary = vec![false; point_count];
    for (node, fixed) in fixed_pressures.iter().enumerate() {
        if fixed.is_some() {
            let root = find_root(&mut parents, node);
            has_boundary[root] = true;
        }
    }

    let mut fracture_ids: Vec<i64> = Vec::new();
    for (cell, triangle) in input.mesh.triangles.iter().enumerate() {
        let root = find_root(&mut parents, triangle[0]);
        if !has_boundary[root] {
            fracture_ids.push(input.mesh.fracture_ids[cell]);
        }
    }

    if fracture_ids.is_empty() {
        return Ok(());
    }

    fracture_ids.sort_unstable();
    fracture_ids.dedup();
    Err(FehmError::NoBoundary { fracture_ids })
}
