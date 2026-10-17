"""Writes a triangle mesh as an ASCII VTK UnstructuredGrid file (.vtu), in
the format which fehm-rs reads."""


def write_triangle_vtu(path, points, triangles, cell_arrays, point_arrays):
    """Writes the mesh.

    points: list of (x, y, z)
    triangles: list of (a, b, c), point indices from 0
    cell_arrays: list of (name, VTK type, values), one value for each triangle
    point_arrays: list of (name, VTK type, values), one value for each point
    """
    def data_array(name, vtk_type, values, components=1):
        text = " ".join(str(value) for value in values)
        return (f'<DataArray type="{vtk_type}" Name="{name}" '
                f'NumberOfComponents="{components}" format="ascii">\n{text}\n</DataArray>\n')

    with open(path, "w") as file:
        file.write('<?xml version="1.0"?>\n')
        file.write('<VTKFile type="UnstructuredGrid" version="1.0">\n<UnstructuredGrid>\n')
        file.write(f'<Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(triangles)}">\n')
        file.write("<Points>\n")
        file.write(data_array("Points", "Float64",
                              [repr(float(c)) for point in points for c in point], 3))
        file.write("</Points>\n<Cells>\n")
        file.write(data_array("connectivity", "Int64",
                              [index for triangle in triangles for index in triangle]))
        file.write(data_array("offsets", "Int64", range(3, 3 * len(triangles) + 1, 3)))
        file.write(data_array("types", "UInt8", [5] * len(triangles)))
        file.write("</Cells>\n<CellData>\n")
        for name, vtk_type, values in cell_arrays:
            file.write(data_array(name, vtk_type, [repr(value) for value in values]))
        file.write("</CellData>\n<PointData>\n")
        for name, vtk_type, values in point_arrays:
            file.write(data_array(name, vtk_type, values))
        file.write("</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n")
