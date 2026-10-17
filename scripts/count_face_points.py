"""Counts the points of an AVS UCD mesh (.inp) on each face of its bounding box.

Usage: python3 scripts/count_face_points.py <mesh.inp>
"""

import sys


def main(path):
    with open(path) as file:
        lines = file.read().splitlines()
    point_count = int(lines[0].split()[0])
    points = [[float(value) for value in line.split()[1:4]] for line in lines[1:1 + point_count]]

    for axis, name in enumerate("xyz"):
        values = [point[axis] for point in points]
        low, high = min(values), max(values)
        tolerance = 1e-6 * (high - low)
        low_count = sum(1 for value in values if value - low <= tolerance)
        high_count = sum(1 for value in values if high - value <= tolerance)
        print(f"{name}: {low:.6g} .. {high:.6g}; points on min face {low_count}, on max face {high_count}")


if __name__ == "__main__":
    main(sys.argv[1])
