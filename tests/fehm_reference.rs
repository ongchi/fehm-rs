//! The pressures of fehm-rs against those of the Fortran FEHM V3.6.3.
//!
//! `tests/data/fehm_reference/<case>.fin` is the final state of FEHM for
//! the case, from `scripts/compare_with_fehm.sh` (the FEHM input of
//! `examples/export_fehm.rs`: the same transmissibilities, each node held at
//! the temperature, impedance boundaries, time steps to a steady state).

use fehm::{build_input_from_text, solve};

/// The values of the `pressure` block of a FEHM `.fin` file
fn read_fin_pressures(text: &str) -> Vec<f64> {
    let lines: Vec<&str> = text.lines().collect();
    let node_count: usize = lines[3].split_whitespace().next().unwrap().parse().unwrap();
    let start = lines
        .iter()
        .position(|line| line.trim() == "pressure")
        .unwrap();
    lines[start + 1..]
        .iter()
        .flat_map(|line| line.split_whitespace())
        .take(node_count)
        .map(|value| value.parse().unwrap())
        .collect()
}

/// The largest pressure difference must be below 1e-8 of the pressure range
fn assert_pressures_match_fehm(input_text: &str, mesh_text: &str, fin_text: &str) {
    let input = build_input_from_text(input_text, mesh_text).unwrap();
    let pressures = solve(&input).unwrap().pressures;
    let fehm_pressures = read_fin_pressures(fin_text);

    assert_eq!(pressures.len(), fehm_pressures.len());
    let range = pressures.iter().fold(f64::MIN, |max, &p| max.max(p))
        - pressures.iter().fold(f64::MAX, |min, &p| min.min(p));
    for (node, (pressure, fehm_pressure)) in pressures.iter().zip(&fehm_pressures).enumerate() {
        assert!(
            (pressure - fehm_pressure).abs() < 1e-8 * range,
            "node {node}: fehm-rs {pressure}, FEHM {fehm_pressure}"
        );
    }
}

#[test]
fn single_fracture_matches_fehm() {
    assert_pressures_match_fehm(
        include_str!("data/single_fracture.toml"),
        include_str!("data/single_fracture.vtu"),
        include_str!("data/fehm_reference/single_fracture.fin"),
    );
}

#[test]
fn crossing_fractures_match_fehm() {
    assert_pressures_match_fehm(
        include_str!("data/crossing_fractures.toml"),
        include_str!("data/crossing_fractures.vtu"),
        include_str!("data/fehm_reference/crossing_fractures.fin"),
    );
}

#[test]
fn lagrit_dfn_mesh_matches_fehm() {
    assert_pressures_match_fehm(
        include_str!("data/lagrit_dfn/output_dfn_recon1.toml"),
        include_str!("data/lagrit_dfn/output_dfn_recon1.vtu"),
        include_str!("data/fehm_reference/output_dfn_recon1.fin"),
    );
}
