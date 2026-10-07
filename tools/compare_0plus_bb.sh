#!/usr/bin/env bash
# Radial-style BB settings; same preserved starting field as the AA comparison.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
reference_run="runs/260925-221559-single-harmonic-0plus-T030-Fs0-localized-0plus"
export OMPI_MCA_pml="${OMPI_MCA_pml:-ob1}"
export OMPI_MCA_btl="${OMPI_MCA_btl:-self,sm}"
exec caffeinate -i python3 tools/run_axisymmetric_core_benchmark.py \
  --input "$reference_run/input.nml" \
  --restart "$reference_run/final_fields_2d.dat" \
  --build-dir work/bb-iteration-build \
  --case-name single-harmonic-0plus-bb \
  --initialization bb-radial-style \
  --ranks 10 --max-iterations 100 --checkpoint-interval 1 \
  --iteration-method bb --anderson-cycle 0 \
  --bb-curvature absolute --bb-initial-mixing 1 \
  --bb-minimum-mixing 0.001 --bb-maximum-mixing 100 --bb-growth-limit 0 "$@"
