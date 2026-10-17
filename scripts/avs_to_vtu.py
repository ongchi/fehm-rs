"""Converts a LaGriT AVS UCD triangle mesh (.inp) to the .vtu format of fehm-rs.

Usage: python3 scripts/avs_to_vtu.py <mesh.inp> <mesh.vtu>

The material of each triangle becomes its fracture_id. All triangles have
the aperture 1e-4 m and the permeability b²/12 (cubic law). The point data
arrays boundary__west and boundary__east mark the points in the lowest and
in the highest 10 % of the x range of the mesh.
"""

import sys

from vtu import write_triangle_vtu

APERTURE = 1e-4
BOUNDARY_FRACTION = 0.1


def read_avs_triangles(path):
    with open(path) as file:
        lines = file.read().splitlines()

    point_count, cell_count = (int(value) for value in lines[0].split()[:2])
    points = [tuple(float(value) for value in line.split()[1:4])
              for line in lines[1:1 + point_count]]
    triangles = []
    materials = []
    for line in lines[1 + point_count:1 + point_count + cell_count]:
        fields = line.split()
        if fields[2] != "tri":
            raise ValueError(f"cell {fields[0]} is a {fields[2]}, not a tri")
        materials.append(int(fields[1]))
        triangles.append(tuple(int(value) - 1 for value in fields[3:6]))
    return points, triangles, materials


def main(inp_path, vtu_path):
    points, triangles, materials = read_avs_triangles(inp_path)
    x_values = [point[0] for point in points]
    low, high = min(x_values), max(x_values)
    width = BOUNDARY_FRACTION * (high - low)
    write_triangle_vtu(vtu_path, points, triangles,
                       [("fracture_id", "Int32", materials),
                        ("aperture", "Float64", [APERTURE] * len(triangles)),
                        ("permeability", "Float64", [APERTURE ** 2 / 12.0] * len(triangles))],
                       [("boundary__west", "Int32", [1 if x <= low + width else 0 for x in x_values]),
                        ("boundary__east", "Int32", [1 if x >= high - width else 0 for x in x_values])])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
