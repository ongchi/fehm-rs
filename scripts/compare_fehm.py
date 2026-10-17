"""Compares the pressures of fehm-rs with those of the Fortran FEHM.

Usage: python3 scripts/compare_fehm.py <fehm run.fin> <fehm-rs pressure.vtu>

Prints the largest pressure difference, relative to the range of the
fehm-rs pressures, and the temperature range of FEHM.
"""

import sys


def read_fin_block(path, name):
    """The values of a block (pressure, temperature, ...) of a FEHM .fin file"""
    with open(path) as file:
        lines = file.read().splitlines()
    node_count = int(lines[3].split()[0])
    start = next(index for index, line in enumerate(lines) if line.strip() == name) + 1
    values = []
    for line in lines[start:]:
        values.extend(float(value) for value in line.split())
        if len(values) >= node_count:
            return values[:node_count]
    raise ValueError(f"the block {name} of {path} is too short")


def read_vtu_array(path, name):
    """The values of a DataArray of a .vtu file of fehm-rs: the values are
    between the line of the start tag and the line </DataArray>"""
    with open(path) as file:
        lines = file.read().splitlines()
    start = next(index for index, line in enumerate(lines)
                 if line.startswith("<DataArray") and f'Name="{name}"' in line) + 1
    end = lines.index("</DataArray>", start)
    return [float(value) for line in lines[start:end] for value in line.split()]


def main(fin_path, vtu_path):
    fehm = read_fin_block(fin_path, "pressure")
    temperatures = read_fin_block(fin_path, "temperature")
    fehm_rs = read_vtu_array(vtu_path, "pressure")
    if len(fehm) != len(fehm_rs):
        raise ValueError(f"FEHM has {len(fehm)} nodes, fehm-rs {len(fehm_rs)}")

    pressure_range = max(fehm_rs) - min(fehm_rs)
    differences = [abs(a - b) for a, b in zip(fehm, fehm_rs)]
    largest = max(range(len(differences)), key=differences.__getitem__)
    print(f"nodes: {len(fehm)}")
    print(f"pressure range (fehm-rs): {pressure_range:.6e} MPa")
    print(f"largest difference: {differences[largest]:.6e} MPa at node {largest + 1} "
          f"({differences[largest] / pressure_range:.3e} of the range)")
    print(f"FEHM temperature range: {min(temperatures):.10f} .. {max(temperatures):.10f} °C")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
