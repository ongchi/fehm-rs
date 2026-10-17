#!/bin/sh
# Runs one case with fehm-rs and with the Fortran FEHM, and compares the pressures.
#
# Usage: scripts/compare_with_fehm.sh <input.toml> <work_folder> <xfehm executable>
#
# fehm-rs writes to <work_folder>/fehm-rs, the FEHM input and output are in
# <work_folder>/fehm.
set -eu

input_file=$1
work_folder=$2
xfehm=$3

cargo run --quiet --release -- "$input_file" "$work_folder/fehm-rs" > "$work_folder.fehm-rs.log"
cargo run --quiet --release --example export_fehm -- "$input_file" "$work_folder/fehm"

(
    cd "$work_folder/fehm"
    rm -f run.fin run.out
    "$xfehm" fehmn.files > xfehm.log 2>&1
)

if [ ! -f "$work_folder/fehm/run.fin" ]; then
    echo "FEHM wrote no run.fin; see $work_folder/fehm/run.out" >&2
    exit 1
fi
grep "timestep less than daymin" "$work_folder/fehm/run.out" && exit 1

python3 scripts/compare_fehm.py "$work_folder/fehm/run.fin" "$work_folder/fehm-rs/pressure.vtu"
