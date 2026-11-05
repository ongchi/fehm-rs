//! Flow solutions on voxel grids, against analytic values.

use fehm::flow::GRAVITY;
use fehm::voxel::{VoxelGrid, calculate_voxel_control_volumes};
use fehm::{FehmError, FlowProblem, NodeSolution, SolverOptions, eos, solve_control_volumes};

const TEMPERATURE: f64 = 50.0;
const PERMEABILITY: f64 = 1e-12;

const SOLVER: SolverOptions = SolverOptions {
    tolerance: 1e-10,
    max_iterations: 20,
};

/// A row of voxels of 1 m along `axis` (0: x, 2: z)
fn row(axis: usize, permeabilities: Vec<Option<f64>>) -> VoxelGrid {
    let mut shape = [1, 1, 1];
    shape[axis] = permeabilities.len();
    VoxelGrid {
        shape,
        voxel_size: [1.0, 1.0, 1.0],
        origin: [0.0, 0.0, 0.0],
        permeabilities,
    }
}

/// Solves the flow on the grid with the fixed pressures `(node, pressure)`
fn solve_grid(
    grid: &VoxelGrid,
    fixed: &[(usize, f64)],
    has_gravity: bool,
) -> Result<NodeSolution, FehmError> {
    let voxels = calculate_voxel_control_volumes(grid)?;
    let mut fixed_pressures = vec![None; voxels.elevations.len()];
    for &(node, pressure) in fixed {
        fixed_pressures[node] = Some(pressure);
    }
    solve_control_volumes(&FlowProblem {
        control_volumes: &voxels.control_volumes,
        elevations: &voxels.elevations,
        fixed_pressures: &fixed_pressures,
        initial_pressure: 10.0,
        temperature: TEMPERATURE,
        has_gravity,
        solver: &SOLVER,
    })
}

/// ρ/μ (kg/m³/(Pa·s)) at a pressure (MPa)
fn mobility(pressure: f64) -> f64 {
    eos::density(pressure, TEMPERATURE).value / eos::viscosity(pressure, TEMPERATURE).value
}

#[test]
fn row_without_gravity_has_a_linear_pressure_profile() {
    let grid = row(0, vec![Some(PERMEABILITY); 10]);
    let solution = solve_grid(&grid, &[(0, 10.1), (9, 10.0)], false).unwrap();

    // ρ/μ changes a little over the pressure drop, thus the profile is not exactly linear
    for (node, pressure) in solution.pressures.iter().enumerate() {
        let linear = 10.1 - 0.1 * node as f64 / 9.0;
        assert!(
            (pressure - linear).abs() < 1e-5 * 0.1,
            "node {node}: pressure {pressure}, linear {linear}"
        );
    }
}

#[test]
fn row_flow_matches_darcy_law() {
    let grid = row(0, vec![Some(PERMEABILITY); 10]);
    let solution = solve_grid(&grid, &[(0, 10.1), (9, 10.0)], false).unwrap();

    // Q = ρ/μ · k · A · Δp / L, with A = 1 m² and L = 9 m between the end centers
    let expected = mobility(10.05) * PERMEABILITY * 0.1e6 / 9.0;
    for flux in &solution.mass_fluxes {
        assert!(
            (flux - expected).abs() < 1e-4 * expected,
            "flux {flux:e}, Darcy {expected:e}"
        );
    }
}

#[test]
fn two_layers_in_series_use_the_harmonic_mean_permeability() {
    let (low, high) = (PERMEABILITY, 4.0 * PERMEABILITY);
    let grid = row(0, vec![Some(low), Some(low), Some(high), Some(high)]);
    let solution = solve_grid(&grid, &[(0, 10.1), (3, 10.0)], false).unwrap();

    // Resistances of the three faces in series, each of 1 m: 1/k, 1/harmonic mean, 1/k
    let harmonic_mean = 2.0 * low * high / (low + high);
    let resistance = 1.0 / low + 1.0 / harmonic_mean + 1.0 / high;
    let expected = mobility(10.05) * 0.1e6 / resistance;
    for flux in &solution.mass_fluxes {
        assert!(
            (flux - expected).abs() < 1e-4 * expected,
            "flux {flux:e}, expected {expected:e}"
        );
    }
}

#[test]
fn hydrostatic_column_has_no_flow() {
    let grid = row(2, vec![Some(PERMEABILITY); 6]);
    let solution = solve_grid(&grid, &[(5, 10.0)], true).unwrap();

    let density = eos::density(10.0, TEMPERATURE).value;
    let flow_scale = mobility(10.0) * PERMEABILITY * density * GRAVITY;
    let largest_flux = solution
        .mass_fluxes
        .iter()
        .fold(0.0, |max: f64, flux| max.max(flux.abs()));
    assert!(
        largest_flux < 1e-12 * flow_scale,
        "largest flux {largest_flux:e}"
    );

    // Discrete hydrostatic pressure, 1 m for each step: P_lower = P_upper + ½(ρ_upper + ρ_lower) g · 1 m
    let mut expected = 10.0;
    for _ in 0..5 {
        let upper = expected;
        for _ in 0..50 {
            let mean_density = 0.5
                * (eos::density(upper, TEMPERATURE).value
                    + eos::density(expected, TEMPERATURE).value);
            expected = upper + mean_density * GRAVITY * 1e-6;
        }
    }
    let bottom_pressure = solution.pressures[0];
    assert!(
        (bottom_pressure - expected).abs() < 1e-9 * (expected - 10.0),
        "bottom pressure {bottom_pressure}, expected {expected}"
    );
}

#[test]
fn group_without_fixed_pressure_is_an_error() {
    // The voxel with no permeability separates the last voxel from the fixed pressure
    let grid = row(
        0,
        vec![
            Some(PERMEABILITY),
            Some(PERMEABILITY),
            None,
            Some(PERMEABILITY),
        ],
    );
    let error = solve_grid(&grid, &[(0, 10.0)], false).unwrap_err();
    assert!(
        matches!(&error, FehmError::NoBoundaryNodes { nodes } if nodes == &[2]),
        "{error}"
    );
}

#[test]
fn node_array_of_wrong_length_is_an_error() {
    let voxels = calculate_voxel_control_volumes(&row(0, vec![Some(PERMEABILITY); 3])).unwrap();
    let error = solve_control_volumes(&FlowProblem {
        control_volumes: &voxels.control_volumes,
        elevations: &voxels.elevations,
        fixed_pressures: &[Some(10.0)],
        initial_pressure: 10.0,
        temperature: TEMPERATURE,
        has_gravity: false,
        solver: &SOLVER,
    })
    .unwrap_err();
    assert!(
        matches!(
            error,
            FehmError::NodeArrayLength {
                array: "fixed_pressures",
                length: 1,
                node_count: 3
            }
        ),
        "{error}"
    );
}
