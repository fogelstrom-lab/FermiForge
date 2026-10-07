#!/usr/bin/env bash
# Continue the completed AA20/pmax=5 run; preserve its data and join histories.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
parent="runs/261002-091852-radial-anderson-aa20-diagnostics-0plus-T030-Fs54-radial-reference"
bash "$project_dir/tools/run_radial_symmetry.sh" \
  --input "$parent/input.nml" \
  --restart "$parent/final_fields_2d.dat" \
  --history-parent-run "$parent" \
  --build-dir work/radial-anderson-diagnostics-build \
  --case-name radial-anderson-aa20-0plus-pmax5-tight-continuation \
  --iteration-method anderson --anderson-cycle 20 --simple-cycle 0 \
  --anderson-history 10 --anderson-pmax 5 --anderson-progress 0.1 \
  --max-iterations 600 --tolerance 2e-6 --checkpoint-interval 1 \
  --save-iteration-diagnostics "$@"
