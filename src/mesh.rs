/// Triangle mesh of a discrete fracture network
///
/// Fractures which intersect share the points on their intersection line.
#[derive(Clone, Debug)]
pub struct Mesh {
    pub points: Vec<[f64; 3]>,
    /// The three point indices of each triangle
    pub triangles: Vec<[usize; 3]>,
    /// For each triangle: the fracture which it is a part of
    pub fracture_ids: Vec<i64>,
    /// For each triangle: fracture aperture (m)
    pub apertures: Vec<f64>,
    /// For each triangle: fracture permeability (m²)
    pub permeabilities: Vec<f64>,
    /// The point data arrays `boundary__<name>`, in the order of the file
    pub boundaries: Vec<BoundaryNodes>,
}

/// The points of one point data array `boundary__<name>`
#[derive(Clone, Debug)]
pub struct BoundaryNodes {
    /// `<name>` of the array name
    pub name: String,
    /// The points with a value above 0, in increasing order
    pub nodes: Vec<usize>,
}
