//! Flow solutions on the meshes of `scripts/make_test_meshes.py`.

use fehm::flow::GRAVITY;
use fehm::{Input, Solution, build_input_from_text, eos, solve};

const SINGLE_FRACTURE: &str = include_str!("data/single_fracture.vtu");
const CROSSING_FRACTURES: &str = include_str!("data/crossing_fractures.vtu");
const ISOLATED_FRACTURE: &str = include_str!("data/isolated_fracture.vtu");

const TEMPERATURE: f64 = 50.0;

/// Input with gravity on or off and boundaries `(name, pressure)`
fn input(mesh: &str, has_gravity: bool, boundaries: &[(&str, f64)]) -> Input {
    let mut text =
        format!("mesh = \"mesh.vtu\"\ntemperature = {TEMPERATURE}\ngravity = {has_gravity}\n");
    for (name, pressure) in boundaries {
        text += &format!("[[boundary]]\nname = \"{name}\"\npressure = {pressure}\n");
    }
    build_input_from_text(&text, mesh).unwrap()
}

fn largest_flux(solution: &Solution) -> f64 {
    solution
        .mass_fluxes
        .iter()
        .fold(0.0, |max: f64, flux| max.max(flux.abs()))
}

#[test]
fn single_fracture_without_gravity_has_a_linear_pressure_profile() {
    let input = input(SINGLE_FRACTURE, false, &[("west", 10.1), ("east", 10.0)]);
    let solution = solve(&input).unwrap();

    // ρ/μ changes by about 5e-5 over the pressure drop, thus the profile is
    // not exactly linear: the deviation is about 2e-6 of the drop
    for (point, pressure) in input.mesh.points.iter().zip(&solution.pressures) {
        let linear = 10.1 - 0.1 * point[0] / 10.0;
        assert!(
            (pressure - linear).abs() < 1e-5 * 0.1,
            "x = {}: pressure {pressure}, linear {linear}",
            point[0]
        );
    }
}

#[test]
fn single_fracture_flow_matches_darcy_law() {
    let input = input(SINGLE_FRACTURE, false, &[("west", 10.1), ("east", 10.0)]);
    let solution = solve(&input).unwrap();

    // Q = ρ/μ · k · b · H · Δp / L, with ρ and μ at the upstream side of each connection
    let aperture = 1e-4;
    let permeability = aperture * aperture / 12.0;
    let density = eos::density(10.05, TEMPERATURE).value;
    let viscosity = eos::viscosity(10.05, TEMPERATURE).value;
    let expected = density / viscosity * permeability * aperture * 5.0 * 0.1e6 / 10.0;

    let inflow = solution.boundary_flows[0];
    assert!(
        (inflow - expected).abs() < 1e-4 * expected,
        "in-flow {inflow:e}, Darcy {expected:e}"
    );
}

#[test]
fn hydrostatic_pressure_gives_no_flow() {
    let input = input(SINGLE_FRACTURE, true, &[("top", 10.0)]);
    let solution = solve(&input).unwrap();

    let density = eos::density(10.0, TEMPERATURE).value;
    let flow_scale = density / 5e-4 * 1e-4f64.powi(3) / 12.0 * density * GRAVITY;
    assert!(
        largest_flux(&solution) < 1e-12 * flow_scale,
        "largest flux {:e}",
        largest_flux(&solution)
    );

    // The bottom is 5 m below the top. Discrete hydrostatic pressure, 1 m
    // for each step: P_lower = P_upper + ½(ρ_upper + ρ_lower) g · 1 m
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
fn crossing_fractures_conserve_mass() {
    let input = input(
        CROSSING_FRACTURES,
        true,
        &[
            ("west", 10.2),
            ("east", 10.0),
            ("south", 10.1),
            ("north", 10.05),
        ],
    );
    let solution = solve(&input).unwrap();

    let net_flow: f64 = solution.boundary_flows.iter().sum();
    let inflow: f64 = solution
        .boundary_flows
        .iter()
        .filter(|flow| **flow > 0.0)
        .sum();
    assert!(inflow > 0.0);
    assert!(
        net_flow.abs() < 1e-10 * inflow,
        "net flow {net_flow:e}, in-flow {inflow:e}"
    );

    // Each node without a fixed pressure has no net flux
    let mut net_fluxes = vec![0.0; input.mesh.points.len()];
    for (connection, flux) in solution
        .control_volumes
        .connections
        .iter()
        .zip(&solution.mass_fluxes)
    {
        net_fluxes[connection.nodes[0]] += flux;
        net_fluxes[connection.nodes[1]] -= flux;
    }
    let fixed: Vec<usize> = input
        .boundaries
        .iter()
        .flat_map(|b| b.nodes.clone())
        .collect();
    for (node, net_flux) in net_fluxes.iter().enumerate() {
        if !fixed.contains(&node) {
            assert!(
                net_flux.abs() < 1e-10 * inflow,
                "node {node}: net flux {net_flux:e}"
            );
        }
    }
}

#[test]
fn crossing_fractures_exchange_flow_on_the_intersection() {
    // Flow enters only through fracture 1 (west) and leaves only through fracture 2 (north)
    let input = input(
        CROSSING_FRACTURES,
        false,
        &[("west", 10.1), ("north", 10.0)],
    );
    let solution = solve(&input).unwrap();

    assert!(solution.boundary_flows[0] > 0.0);
    assert!(
        (solution.boundary_flows[0] + solution.boundary_flows[1]).abs()
            < 1e-10 * solution.boundary_flows[0]
    );
}

#[test]
fn fracture_without_boundary_is_an_error() {
    let input = input(ISOLATED_FRACTURE, true, &[("west", 10.1), ("east", 10.0)]);
    let error = solve(&input).unwrap_err().to_string();
    assert!(
        error.starts_with("The fractures 3 have no connection to a boundary"),
        "{error}"
    );
}
