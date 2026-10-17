use crate::coefficients::{cross, norm, subtract};
use crate::eos;
use crate::flow::GRAVITY;
use crate::io::input::Input;

/// Darcy velocity (m/s) of each triangle, in the plane of the triangle:
/// `q = −(k/μ) · (∇P + ρ g ∇z)`, with the gradient of the linear pressure
/// in the triangle, and ρ and μ at the mean pressure of its points.
pub fn calculate_darcy_velocities(input: &Input, pressures: &[f64]) -> Vec<[f64; 3]> {
    let gravity = if input.has_gravity { GRAVITY } else { 0.0 };

    input
        .mesh
        .triangles
        .iter()
        .enumerate()
        .map(|(cell, triangle)| {
            let vertices = triangle.map(|point| input.mesh.points[point]);
            // Pa
            let vertex_pressures = triangle.map(|point| pressures[point] * 1e6);
            let mean_pressure = triangle.iter().map(|&point| pressures[point]).sum::<f64>() / 3.0;
            let density = eos::density(mean_pressure, input.temperature).value;
            let viscosity = eos::viscosity(mean_pressure, input.temperature).value;

            let normal = cross(
                subtract(vertices[1], vertices[0]),
                subtract(vertices[2], vertices[0]),
            );
            let double_area = norm(normal);
            let unit_normal = normal.map(|component| component / double_area);

            // ∇P = Σ P_k · (n × (edge opposite k)) / (2A), edges counterclockwise about n
            let mut pressure_gradient = [0.0; 3];
            for (k, i, j) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
                let shape_gradient = cross(unit_normal, subtract(vertices[j], vertices[i]));
                for axis in 0..3 {
                    pressure_gradient[axis] +=
                        vertex_pressures[k] * shape_gradient[axis] / double_area;
                }
            }

            // ∇z in the plane: e_z minus its normal part
            let vertical_gradient = {
                let normal_part = unit_normal[2];
                [
                    -normal_part * unit_normal[0],
                    -normal_part * unit_normal[1],
                    1.0 - normal_part * unit_normal[2],
                ]
            };

            let factor = -input.mesh.permeabilities[cell] / viscosity;
            [0, 1, 2].map(|axis| {
                factor * (pressure_gradient[axis] + density * gravity * vertical_gradient[axis])
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::input::SolverOptions;
    use crate::mesh::Mesh;

    fn input(points: Vec<[f64; 3]>, has_gravity: bool) -> Input {
        Input {
            mesh: Mesh {
                points,
                triangles: vec![[0, 1, 2]],
                fracture_ids: vec![1],
                apertures: vec![1e-3],
                permeabilities: vec![1e-12],
                boundaries: Vec::new(),
            },
            temperature: 50.0,
            has_gravity,
            boundaries: Vec::new(),
            solver: SolverOptions {
                tolerance: 1e-10,
                max_iterations: 10,
            },
        }
    }

    #[test]
    fn velocity_follows_the_pressure_gradient() {
        let input = input(
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 5.0]],
            false,
        );
        // P = 10 MPa − 0.5 MPa/m · x
        let pressures = [10.0, 9.0, 10.0];
        let velocity = calculate_darcy_velocities(&input, &pressures)[0];
        let viscosity = eos::viscosity(29.0 / 3.0, 50.0).value;
        let expected = 1e-12 / viscosity * 0.5e6;

        assert!((velocity[0] - expected).abs() < 1e-12 * expected);
        assert!(velocity[1].abs() < 1e-12 * expected);
        assert!(velocity[2].abs() < 1e-12 * expected);
    }

    #[test]
    fn hydrostatic_pressure_gives_no_velocity() {
        let input = input(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 1.0]],
            true,
        );
        let density = eos::density(10.0, 50.0).value;
        // P = 10 MPa − ρ g z
        let pressures = [10.0, 10.0, 10.0 - density * GRAVITY * 1e-6];
        let velocity = calculate_darcy_velocities(&input, &pressures)[0];
        let scale = 1e-12 / 5e-4 * density * GRAVITY;

        for component in velocity {
            // The density at the mean pressure differs from ρ(10 MPa) by about 1e-6
            assert!(component.abs() < 1e-5 * scale, "{velocity:?}");
        }
    }
}
