#!/usr/bin/env bash
# Same initial field as the AA20/simple3 experiment: isolate the schedule.
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
  --case-name single-harmonic-0plus-aa-restart \
  --initialization aa20-no-simple \
  --ranks 10 --max-iterations 60 --checkpoint-interval 1 \
  --anderson-cycle 20 --simple-cycle 0 "$@"
