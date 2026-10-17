use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use fehm::{FehmError, read_input, solve, write_output_to_folder};
use tracing::{error, info};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to input file (TOML)
    input_file: PathBuf,

    /// Path to output folder
    output_folder: PathBuf,
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    init_tracing();

    info!("Starting fehm-rs");

    match run(&cli) {
        Ok(()) => {
            info!("fehm-rs - Complete");
            ExitCode::SUCCESS
        }
        Err(e) => {
            error!("{}", e);
            ExitCode::FAILURE
        }
    }
}

/// Reads the input, solves the flow and writes the output files.
fn run(cli: &Cli) -> Result<(), FehmError> {
    let input = read_input(&cli.input_file)?;
    info!(
        "Mesh: {} points, {} triangles",
        input.mesh.points.len(),
        input.mesh.triangles.len()
    );

    let solution = solve(&input)?;

    let total_inflow: f64 = solution
        .boundary_flows
        .iter()
        .filter(|flow| **flow > 0.0)
        .sum();
    for (boundary, flow) in input.boundaries.iter().zip(&solution.boundary_flows) {
        info!(
            "Boundary {}: mass flow into the network {flow:e} kg/s",
            boundary.name
        );
    }
    let net_flow: f64 = solution.boundary_flows.iter().sum();
    info!(
        "Mass balance: net flow {net_flow:e} kg/s, {:e} of the in-flow",
        if total_inflow > 0.0 {
            net_flow.abs() / total_inflow
        } else {
            0.0
        }
    );

    write_output_to_folder(&input, &solution, &cli.output_folder)?;

    Ok(())
}

/// Logs to stdout. The level is `info`, unless the environment variable `RUST_LOG` gives another.
fn init_tracing() {
    let level_filter = tracing_subscriber::EnvFilter::try_from_default_env()
        .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info"));

    tracing_subscriber::fmt()
        .compact()
        .with_env_filter(level_filter)
        .init();
}
