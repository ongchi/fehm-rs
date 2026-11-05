//! Control volumes of a regular voxel grid: one node at the center of each
//! voxel with a permeability, and one connection for each face which two
//! such voxels share (two-point flux).

use crate::FehmError;
use crate::coefficients::{Connection, ControlVolumes};

/// A regular grid of voxels, x fastest
#[derive(Clone, Debug)]
pub struct VoxelGrid {
    /// Number of voxels in x, y and z
    pub shape: [usize; 3],
    /// Size of one voxel in x, y and z (m)
    pub voxel_size: [f64; 3],
    /// The outer corner of the voxel (0, 0, 0) (m)
    pub origin: [f64; 3],
    /// For each voxel, index `x + nx · (y + ny · z)`: permeability (m²);
    /// `None` for a voxel with no flow (air, rock), which gets no node
    pub permeabilities: Vec<Option<f64>>,
}

/// Control volumes of a voxel grid, and the map between voxels and nodes
#[derive(Clone, Debug)]
pub struct VoxelControlVolumes {
    /// `areas` and the `area_coefficient` of each connection are 0: they
    /// are values of a triangle mesh. `aperture_coefficient` is face area /
    /// center distance (m), the coefficient of a 3-D `.stor` file.
    pub control_volumes: ControlVolumes,
    /// For each node: z of the voxel center (m)
    pub elevations: Vec<f64>,
    /// For each node: its voxel
    pub voxel_of_node: Vec<usize>,
    /// For each voxel: its node; `None` for a voxel with no permeability
    pub node_of_voxel: Vec<Option<usize>>,
}

/// Calculates the control volumes of the grid. The transmissibility of a
/// face is face area / center distance · the harmonic mean of the two
/// permeabilities (m³).
pub fn calculate_voxel_control_volumes(grid: &VoxelGrid) -> Result<VoxelControlVolumes, FehmError> {
    check_grid(grid)?;
    let [nx, ny, nz] = grid.shape;
    let [dx, dy, dz] = grid.voxel_size;

    // For each node: its voxel and the permeability of the voxel
    let nodes: Vec<(usize, f64)> = grid
        .permeabilities
        .iter()
        .enumerate()
        .filter_map(|(voxel, permeability)| permeability.map(|value| (voxel, value)))
        .collect();
    let mut node_of_voxel = vec![None; grid.permeabilities.len()];
    for (node, &(voxel, _)) in nodes.iter().enumerate() {
        node_of_voxel[voxel] = Some(node);
    }
    let voxel_of_node: Vec<usize> = nodes.iter().map(|&(voxel, _)| voxel).collect();
    let elevations = voxel_of_node
        .iter()
        .map(|&voxel| grid.origin[2] + (voxel / (nx * ny)) as f64 * dz + 0.5 * dz)
        .collect();

    // For each axis: voxel index step, face area / center distance, and the count along it
    let neighbors = [
        (1, dy * dz / dx, nx),
        (nx, dx * dz / dy, ny),
        (nx * ny, dx * dy / dz, nz),
    ];
    // The neighbors come in increasing voxel order, thus the connections are sorted by nodes
    let mut connections = Vec::new();
    for (node, &(voxel, permeability)) in nodes.iter().enumerate() {
        let position = [voxel % nx, (voxel / nx) % ny, voxel / (nx * ny)];
        for (axis, &(step, coefficient, count)) in neighbors.iter().enumerate() {
            if position[axis] + 1 == count {
                continue;
            }
            let Some(neighbor) = node_of_voxel[voxel + step] else {
                continue;
            };
            let neighbor_permeability = nodes[neighbor].1;
            let mean_permeability =
                2.0 * permeability * neighbor_permeability / (permeability + neighbor_permeability);
            connections.push(Connection {
                nodes: [node, neighbor],
                area_coefficient: 0.0,
                aperture_coefficient: coefficient,
                transmissibility: coefficient * mean_permeability,
            });
        }
    }

    Ok(VoxelControlVolumes {
        control_volumes: ControlVolumes {
            areas: vec![0.0; nodes.len()],
            volumes: vec![dx * dy * dz; nodes.len()],
            connections,
        },
        elevations,
        voxel_of_node,
        node_of_voxel,
    })
}

fn check_grid(grid: &VoxelGrid) -> Result<(), FehmError> {
    let voxel_count = grid
        .shape
        .iter()
        .try_fold(1_usize, |count, &size| count.checked_mul(size));
    if voxel_count != Some(grid.permeabilities.len()) {
        return Err(FehmError::VoxelArrayLength {
            shape: grid.shape,
            length: grid.permeabilities.len(),
        });
    }
    if !grid
        .voxel_size
        .iter()
        .all(|size| size.is_finite() && *size > 0.0)
    {
        return Err(FehmError::VoxelSize {
            voxel_size: grid.voxel_size,
        });
    }
    for (voxel, permeability) in grid.permeabilities.iter().enumerate() {
        if let Some(permeability) = *permeability
            && !(permeability.is_finite() && permeability > 0.0)
        {
            return Err(FehmError::VoxelPermeability {
                voxel,
                permeability,
            });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn grid(shape: [usize; 3], permeabilities: Vec<Option<f64>>) -> VoxelGrid {
        VoxelGrid {
            shape,
            voxel_size: [1.0, 2.0, 3.0],
            origin: [0.0, 0.0, 10.0],
            permeabilities,
        }
    }

    #[test]
    fn grid_of_two_by_two_by_two_has_twelve_faces() {
        let voxels =
            calculate_voxel_control_volumes(&grid([2, 2, 2], vec![Some(1e-12); 8])).unwrap();
        let control_volumes = &voxels.control_volumes;

        assert_eq!(control_volumes.volumes, vec![6.0; 8]);
        assert_eq!(
            voxels.elevations,
            [11.5, 11.5, 11.5, 11.5, 14.5, 14.5, 14.5, 14.5]
        );
        assert_eq!(control_volumes.connections.len(), 12);
        assert!(
            control_volumes
                .connections
                .windows(2)
                .all(|pair| pair[0].nodes < pair[1].nodes)
        );

        // Face area / center distance: x 2·3/1, y 1·3/2, z 1·2/3
        let coefficient_of = |nodes: [usize; 2]| {
            control_volumes
                .connections
                .iter()
                .find(|connection| connection.nodes == nodes)
                .unwrap()
                .aperture_coefficient
        };
        assert_eq!(coefficient_of([0, 1]), 6.0);
        assert_eq!(coefficient_of([0, 2]), 1.5);
        assert!((coefficient_of([0, 4]) - 2.0 / 3.0).abs() < 1e-15);
    }

    #[test]
    fn transmissibility_uses_the_harmonic_mean_permeability() {
        let voxels =
            calculate_voxel_control_volumes(&grid([2, 1, 1], vec![Some(1.0), Some(3.0)])).unwrap();
        let connection = &voxels.control_volumes.connections[0];
        assert_eq!(connection.transmissibility, 6.0 * 1.5);
    }

    #[test]
    fn voxel_without_permeability_has_no_node_and_no_connection() {
        let voxels =
            calculate_voxel_control_volumes(&grid([3, 1, 1], vec![Some(1.0), None, Some(1.0)]))
                .unwrap();
        assert_eq!(voxels.node_of_voxel, [Some(0), None, Some(1)]);
        assert_eq!(voxels.voxel_of_node, [0, 2]);
        assert!(voxels.control_volumes.connections.is_empty());
    }

    #[test]
    fn permeability_array_of_wrong_length_is_an_error() {
        let error = calculate_voxel_control_volumes(&grid([2, 1, 1], vec![Some(1.0)])).unwrap_err();
        assert!(matches!(
            error,
            FehmError::VoxelArrayLength { length: 1, .. }
        ));
    }

    #[test]
    fn voxel_size_of_zero_is_an_error() {
        let mut grid = grid([1, 1, 1], vec![Some(1.0)]);
        grid.voxel_size[1] = 0.0;
        let error = calculate_voxel_control_volumes(&grid).unwrap_err();
        assert!(matches!(error, FehmError::VoxelSize { .. }));
    }

    #[test]
    fn permeability_of_zero_is_an_error() {
        let error = calculate_voxel_control_volumes(&grid([2, 1, 1], vec![Some(1.0), Some(0.0)]))
            .unwrap_err();
        assert!(matches!(
            error,
            FehmError::VoxelPermeability { voxel: 1, .. }
        ));
    }
}
