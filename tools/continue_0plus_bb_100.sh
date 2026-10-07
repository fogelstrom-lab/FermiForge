#!/usr/bin/env bash
# Preserve the first 100 updates; archive and plot the appended continuation.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
parent="runs/260928-232830-single-harmonic-0plus-bb-bb-radial-style"
export OMPI_MCA_pml="${OMPI_MCA_pml:-ob1}"
export OMPI_MCA_btl="${OMPI_MCA_btl:-self,sm}"
exec caffeinate -i python3 tools/run_axisymmetric_core_benchmark.py \
  --input "$parent/input.nml" \
  --restart "$parent/final_fields_2d.dat" \
  --history-parent-run "$parent" \
  --build-dir work/bb-continuation-build \
  --case-name single-harmonic-0plus-bb-continued \
  --initialization bb-radial-style-101-200 \
  --ranks 10 --max-iterations 100 --checkpoint-interval 1 "$@"
