"""Makes the small test meshes of tests/data.

Usage: python3 scripts/make_test_meshes.py tests/data

Each fracture is a rectangle in a coordinate plane, meshed with a grid of
1 m squares, each split into two right triangles. Fractures which cross
share the points on their intersection line.

- single_fracture.vtu: one vertical fracture in the plane y = 0,
  x = 0..10, z = 0..5. Boundaries: west (x = 0), east (x = 10),
  bottom (z = 0), top (z = 5).
- crossing_fractures.vtu: fracture 1 in the plane y = 0 (x, z = -5..5) and
  fracture 2 in the plane x = 0 (y, z = -5..5), with different apertures.
  Boundaries: west (x = -5), east (x = 5), south (y = -5), north (y = 5).
- isolated_fracture.vtu: crossing_fractures.vtu with fracture 3, a
  horizontal fracture at z = 10 which touches no other fracture and no
  boundary.
"""

import os
import sys

from vtu import write_triangle_vtu


class MeshBuilder:
    def __init__(self):
        self.points = []
        self.point_index = {}
        self.triangles = []
        self.fracture_ids = []
        self.apertures = []

    def add_point(self, point):
        key = tuple(round(c, 9) for c in point)
        if key not in self.point_index:
            self.point_index[key] = len(self.points)
            self.points.append(point)
        return self.point_index[key]

    def add_rectangle(self, fracture_id, aperture, origin, u_axis, v_axis, u_count, v_count):
        """Adds a fracture: origin + i * u_axis + j * v_axis, i = 0..u_count, j = 0..v_count"""
        def point(i, j):
            return tuple(origin[k] + i * u_axis[k] + j * v_axis[k] for k in range(3))

        for i in range(u_count):
            for j in range(v_count):
                corners = [self.add_point(point(i + di, j + dj))
                           for di, dj in [(0, 0), (1, 0), (1, 1), (0, 1)]]
                for triangle in [(corners[0], corners[1], corners[2]),
                                 (corners[0], corners[2], corners[3])]:
                    self.triangles.append(triangle)
                    self.fracture_ids.append(fracture_id)
                    self.apertures.append(aperture)

    def write(self, path, boundaries):
        """boundaries: list of (name, function of a point which is true on the boundary)"""
        # Cubic law: k = b² / 12
        permeabilities = [aperture ** 2 / 12.0 for aperture in self.apertures]
        point_arrays = [(f"boundary__{name}", "Int32",
                         [1 if is_on(point) else 0 for point in self.points])
                        for name, is_on in boundaries]
        write_triangle_vtu(path, self.points, self.triangles,
                           [("fracture_id", "Int32", self.fracture_ids),
                            ("aperture", "Float64", self.apertures),
                            ("permeability", "Float64", permeabilities)],
                           point_arrays)


def near(a, b):
    return abs(a - b) < 1e-9


def add_crossing_fractures(builder):
    builder.add_rectangle(1, 1e-4, (-5.0, 0.0, -5.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0), 10, 10)
    builder.add_rectangle(2, 2e-4, (0.0, -5.0, -5.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), 10, 10)


CROSSING_BOUNDARIES = [
    ("west", lambda p: near(p[0], -5.0)),
    ("east", lambda p: near(p[0], 5.0)),
    ("south", lambda p: near(p[1], -5.0)),
    ("north", lambda p: near(p[1], 5.0)),
]


def main(folder):
    os.makedirs(folder, exist_ok=True)

    single = MeshBuilder()
    single.add_rectangle(1, 1e-4, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0), 10, 5)
    single.write(os.path.join(folder, "single_fracture.vtu"), [
        ("west", lambda p: near(p[0], 0.0)),
        ("east", lambda p: near(p[0], 10.0)),
        ("bottom", lambda p: near(p[2], 0.0)),
        ("top", lambda p: near(p[2], 5.0)),
    ])

    crossing = MeshBuilder()
    add_crossing_fractures(crossing)
    crossing.write(os.path.join(folder, "crossing_fractures.vtu"), CROSSING_BOUNDARIES)

    isolated = MeshBuilder()
    add_crossing_fractures(isolated)
    isolated.add_rectangle(3, 1e-4, (1.0, 1.0, 10.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), 2, 2)
    isolated.write(os.path.join(folder, "isolated_fracture.vtu"), CROSSING_BOUNDARIES)


if __name__ == "__main__":
    main(sys.argv[1])
