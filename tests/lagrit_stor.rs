//! The control volumes and coefficients of fehm-rs against the `.stor` file
//! which LaGriT writes for the same DFN mesh (LaGriT test level01/recon1,
//! `output_dfn_recon1.inp` converted by `scripts/avs_to_vtu.py`).

use fehm::coefficients::calculate_control_volumes;
use fehm::io::vtu::parse_mesh;

const MESH: &str = include_str!("data/lagrit_dfn/output_dfn_recon1.vtu");
const STOR: &str = include_str!("data/lagrit_dfn/output_dfn_recon1.stor");

/// The data of an ASCII LaGriT `.stor` file with one area coefficient for each connection
struct Stor {
    volumes: Vec<f64>,
    /// (node, connected node, area coefficient), nodes from 0, without the diagonal
    connections: Vec<(usize, usize, f64)>,
}

fn parse_stor(text: &str) -> Stor {
    let mut lines = text.lines();
    lines.next();
    lines.next();
    let header: Vec<usize> = lines
        .next()
        .unwrap()
        .split_whitespace()
        .map(|value| value.parse().unwrap())
        .collect();
    let (entry_count, node_count, area_coefficient_count) = (header[0], header[1], header[3]);
    assert_eq!(area_coefficient_count, 1);

    let mut values = lines.flat_map(str::split_whitespace);
    let mut next_block = |count: usize| -> Vec<f64> {
        (&mut values)
            .take(count)
            .map(|value| value.parse().unwrap())
            .collect()
    };

    let volumes = next_block(node_count);
    let row_pointers: Vec<usize> = next_block(node_count + 1)
        .iter()
        .map(|&v| v as usize)
        .collect();
    let columns: Vec<usize> = next_block(entry_count)
        .iter()
        .map(|&v| v as usize - 1)
        .collect();
    let coefficient_pointers: Vec<usize> = next_block(entry_count)
        .iter()
        .map(|&v| v as usize)
        .collect();
    next_block(node_count + 1);
    next_block(node_count);
    let coefficients: Vec<f64> = values.map(|value| value.parse().unwrap()).collect();

    let mut connections = Vec::new();
    for node in 0..node_count {
        let first = row_pointers[node] - (node_count + 1);
        let last = row_pointers[node + 1] - (node_count + 1);
        for entry in first..last {
            if columns[entry] != node {
                // LaGriT writes −(area/length)
                let coefficient = -coefficients[coefficient_pointers[entry] - 1];
                connections.push((node, columns[entry], coefficient));
            }
        }
    }

    Stor {
        volumes,
        connections,
    }
}

#[test]
fn volumes_and_coefficients_match_lagrit() {
    let mesh = parse_mesh(MESH).unwrap();
    let control_volumes = calculate_control_volumes(&mesh).unwrap();
    let stor = parse_stor(STOR);

    assert_eq!(stor.volumes.len(), mesh.points.len());

    // The .stor file has 13 significant digits
    let tolerance = 1e-11;
    let largest_volume = stor
        .volumes
        .iter()
        .fold(0.0, |max: f64, v| max.max(v.abs()));
    // LaGriT writes the Voronoi areas; the volumes of fehm-rs also have the aperture
    for (node, (volume, expected)) in control_volumes.areas.iter().zip(&stor.volumes).enumerate() {
        assert!(
            (volume - expected).abs() <= tolerance * largest_volume,
            "node {node}: area {volume}, LaGriT {expected}"
        );
    }

    // Each connection is two times in the .stor file, once in each row
    assert_eq!(
        stor.connections.len(),
        2 * control_volumes.connections.len()
    );
    let largest_coefficient = stor
        .connections
        .iter()
        .fold(0.0, |max: f64, connection| max.max(connection.2.abs()));
    let mut largest_difference = (0.0, [0, 0]);
    for (node, other, expected) in stor.connections {
        let nodes = [node.min(other), node.max(other)];
        let index = control_volumes
            .connections
            .binary_search_by_key(&nodes, |connection| connection.nodes)
            .unwrap_or_else(|_| panic!("LaGriT has the connection {nodes:?}, fehm-rs does not"));
        let difference = (control_volumes.connections[index].area_coefficient - expected).abs();
        if difference > largest_difference.0 {
            largest_difference = (difference, nodes);
        }
    }
    // The coordinates of the .inp file have 13 significant digits, and the
    // cotangent of an angle near 0 or 180° magnifies their rounding error:
    // the largest difference is about 2e-11 of the largest coefficient.
    assert!(
        largest_difference.0 <= 1e-10 * largest_coefficient,
        "largest difference {:e} (connection {:?}), largest coefficient {largest_coefficient}",
        largest_difference.0,
        largest_difference.1
    );
}
