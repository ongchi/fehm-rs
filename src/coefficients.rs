//! Control volumes and connection coefficients of the triangle mesh.
//!
//! This replaces the `.stor` file which LaGriT writes for FEHM
//! (`matbld2d_stor`: 2-D Voronoi volumes and area/length coefficients), with
//! the aperture correction of dfnWorks (`correct_stor_file` in
//! `pydfnworks/dfnFlow/fehm.py`) applied for each triangle.

use crate::FehmError;
use crate::mesh::Mesh;

/// One connection of two nodes: an edge of one or more triangles
#[derive(Clone, Debug, PartialEq)]
pub struct Connection {
    /// The two nodes, `nodes[0] < nodes[1]`
    pub nodes: [usize; 2],
    /// Sum over the triangles of the edge of ½·cot(angle opposite the edge)
    /// (the 2-D Voronoi face length / edge length; LaGriT writes it with a
    /// negative sign)
    pub area_coefficient: f64,
    /// Sum over the triangles of the edge of ½·cot(angle opposite the edge)
    /// · aperture (m): the aperture-corrected `.stor` coefficient of dfnWorks
    pub aperture_coefficient: f64,
    /// Sum over the triangles of the edge of ½·cot(angle opposite the edge)
    /// · aperture · permeability (m³). Mass flux = transmissibility · ρ/μ · Δpotential.
    pub transmissibility: f64,
}

/// Control volumes and connections of a mesh
#[derive(Clone, Debug)]
pub struct ControlVolumes {
    /// For each node: Voronoi area (m²) of the node in its triangles
    pub areas: Vec<f64>,
    /// For each node: Voronoi area · aperture of the triangles (m³)
    pub volumes: Vec<f64>,
    /// Sorted by `nodes`
    pub connections: Vec<Connection>,
}

impl ControlVolumes {
    /// Number of connections with a negative transmissibility (from obtuse triangles)
    pub fn negative_connection_count(&self) -> usize {
        self.connections
            .iter()
            .filter(|connection| connection.transmissibility < 0.0)
            .count()
    }

    /// Number of nodes with a negative volume (from obtuse triangles)
    pub fn negative_volume_count(&self) -> usize {
        self.volumes.iter().filter(|&&volume| volume < 0.0).count()
    }
}

/// Calculates the control volumes and the connections of the mesh.
///
/// For each triangle, in its own plane, the coefficient of the edge `ij` is
/// ½·cot of the angle at the third point `k`, and the Voronoi area of point
/// `i` is (|ij|²·cot k + |ik|²·cot j) / 8. These are the values of a
/// Voronoi (circumcentric) dual mesh; for an obtuse triangle some of them
/// are negative.
pub fn calculate_control_volumes(mesh: &Mesh) -> Result<ControlVolumes, FehmError> {
    let mut areas = vec![0.0; mesh.points.len()];
    let mut volumes = vec![0.0; mesh.points.len()];
    // One entry for each edge of each triangle:
    // (nodes, area coefficient, aperture coefficient, transmissibility)
    let mut edges = Vec::with_capacity(3 * mesh.triangles.len());

    for (cell, triangle) in mesh.triangles.iter().enumerate() {
        let vertices = triangle.map(|point| mesh.points[point]);
        let cotangents = cotangents(&vertices).ok_or(FehmError::DegenerateTriangle { cell })?;
        let thickness = mesh.apertures[cell];
        let conductance = thickness * mesh.permeabilities[cell];

        // Vertex k is opposite the edge (i, j)
        for (k, i, j) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
            let edge_coefficient = 0.5 * cotangents[k];
            let nodes = [triangle[i].min(triangle[j]), triangle[i].max(triangle[j])];
            edges.push((
                nodes,
                edge_coefficient,
                edge_coefficient * thickness,
                edge_coefficient * conductance,
            ));

            // Voronoi area of the edge (i, j): ¼·|ij|² · ½·cot k, half to i and half to j
            let half_area = 0.125 * squared_distance(vertices[i], vertices[j]) * cotangents[k];
            for node in [triangle[i], triangle[j]] {
                areas[node] += half_area;
                volumes[node] += half_area * thickness;
            }
        }
    }

    edges.sort_by_key(|&(nodes, _, _, _)| nodes);

    let mut connections: Vec<Connection> = Vec::new();
    for (nodes, area_coefficient, aperture_coefficient, transmissibility) in edges {
        match connections.last_mut() {
            Some(connection) if connection.nodes == nodes => {
                connection.area_coefficient += area_coefficient;
                connection.aperture_coefficient += aperture_coefficient;
                connection.transmissibility += transmissibility;
            }
            _ => connections.push(Connection {
                nodes,
                area_coefficient,
                aperture_coefficient,
                transmissibility,
            }),
        }
    }

    Ok(ControlVolumes {
        areas,
        volumes,
        connections,
    })
}

/// Cotangent of the angle at each vertex of a triangle; `None` for a triangle with no area
fn cotangents(vertices: &[[f64; 3]; 3]) -> Option<[f64; 3]> {
    let double_area = norm(cross(
        subtract(vertices[1], vertices[0]),
        subtract(vertices[2], vertices[0]),
    ));
    let longest_squared_edge = [(0, 1), (1, 2), (2, 0)]
        .map(|(a, b)| squared_distance(vertices[a], vertices[b]))
        .into_iter()
        .fold(0.0, f64::max);

    if double_area <= 1e-12 * longest_squared_edge {
        return None;
    }

    Some([(0, 1, 2), (1, 2, 0), (2, 0, 1)].map(|(k, i, j)| {
        let to_i = subtract(vertices[i], vertices[k]);
        let to_j = subtract(vertices[j], vertices[k]);
        dot(to_i, to_j) / double_area
    }))
}

pub(crate) fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn squared_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let difference = subtract(a, b);
    dot(difference, difference)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mesh(points: Vec<[f64; 3]>, triangles: Vec<[usize; 3]>) -> Mesh {
        let cell_count = triangles.len();
        Mesh {
            points,
            triangles,
            fracture_ids: vec![1; cell_count],
            apertures: vec![2.0; cell_count],
            permeabilities: vec![3.0; cell_count],
            boundaries: Vec::new(),
        }
    }

    #[test]
    fn unit_square_matches_lagrit() {
        // LaGriT test level01/pflotran_stor, out_tri_1x1m.stor
        let control_volumes = calculate_control_volumes(&mesh(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        ))
        .unwrap();

        assert_eq!(control_volumes.areas, vec![0.25; 4]);
        assert_eq!(control_volumes.volumes, vec![0.5; 4]);

        let coefficients: Vec<([usize; 2], f64)> = control_volumes
            .connections
            .iter()
            .map(|connection| (connection.nodes, connection.area_coefficient))
            .collect();
        let expected = [
            ([0, 1], 0.5),
            ([0, 2], 0.0),
            ([0, 3], 0.5),
            ([1, 2], 0.5),
            ([2, 3], 0.5),
        ];
        assert_eq!(coefficients.len(), expected.len());
        for ((nodes, coefficient), (expected_nodes, expected_coefficient)) in
            coefficients.iter().zip(expected)
        {
            assert_eq!(*nodes, expected_nodes);
            assert!((coefficient - expected_coefficient).abs() < 1e-15);
        }
        assert!((control_volumes.connections[0].aperture_coefficient - 1.0).abs() < 1e-15);
        assert!((control_volumes.connections[0].transmissibility - 3.0).abs() < 1e-15);
    }

    #[test]
    fn aperture_coefficient_sums_the_aperture_of_each_triangle() {
        // Two triangles on the edge (0, 1); the angle opposite it has cot = 3/4 in each
        let mut rhombus = mesh(
            vec![
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
                [1.0, -2.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 3, 1]],
        );
        rhombus.apertures = vec![1.0, 2.0];
        let control_volumes = calculate_control_volumes(&rhombus).unwrap();

        let shared_edge = &control_volumes.connections[0];
        assert_eq!(shared_edge.nodes, [0, 1]);
        // 0.375 · 1 + 0.375 · 2
        assert!((shared_edge.aperture_coefficient - 1.125).abs() < 1e-15);
    }

    #[test]
    fn coefficients_do_not_depend_on_the_plane_of_the_triangle() {
        let flat = calculate_control_volumes(&mesh(
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.5, 1.5, 0.0]],
            vec![[0, 1, 2]],
        ))
        .unwrap();
        // The same triangle, rotated about the x axis into the plane y = z
        let s = 1.5 / 2.0_f64.sqrt();
        let tilted = calculate_control_volumes(&mesh(
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.5, s, s]],
            vec![[0, 1, 2]],
        ))
        .unwrap();

        for (a, b) in flat.connections.iter().zip(&tilted.connections) {
            assert!((a.area_coefficient - b.area_coefficient).abs() < 1e-14);
        }
        for (a, b) in flat.areas.iter().zip(&tilted.areas) {
            assert!((a - b).abs() < 1e-14);
        }
    }

    #[test]
    fn sum_of_areas_is_triangle_area_for_an_obtuse_triangle() {
        let control_volumes = calculate_control_volumes(&mesh(
            vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [1.0, 0.5, 0.0]],
            vec![[0, 1, 2]],
        ))
        .unwrap();

        let total: f64 = control_volumes.areas.iter().sum();
        assert!((total - 1.0).abs() < 1e-14);
        assert_eq!(control_volumes.negative_connection_count(), 1);
        assert_eq!(control_volumes.negative_volume_count(), 2);
    }

    #[test]
    fn triangle_with_no_area_is_an_error() {
        let error = calculate_control_volumes(&mesh(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0, 1, 2]],
        ))
        .unwrap_err();
        assert!(
            matches!(error, FehmError::DegenerateTriangle { cell: 0 }),
            "{error}"
        );
    }
}
