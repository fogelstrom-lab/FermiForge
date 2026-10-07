#!/usr/bin/env bash
# Restart the preserved 0+ run with two complete AA20/simple3 cycles.
# Additional runner options may be supplied, e.g. --max-iterations 69.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
reference_run="runs/260925-221559-single-harmonic-0plus-T030-Fs0-localized-0plus"
export OMPI_MCA_pml="${OMPI_MCA_pml:-ob1}"
export OMPI_MCA_btl="${OMPI_MCA_btl:-self,sm}"
exec caffeinate -i python3 tools/run_axisymmetric_core_benchmark.py \
  --input "$reference_run/input.nml" \
  --restart "$reference_run/final_fields_2d.dat" \
  --build-dir work/cycled-iteration-build \
  --case-name single-harmonic-0plus-cycled \
  --initialization aa20-simple3 \
  --ranks 10 --max-iterations 46 --checkpoint-interval 1 \
  --anderson-cycle 20 --simple-cycle 3 --simple-mixing 0.1 "$@"
