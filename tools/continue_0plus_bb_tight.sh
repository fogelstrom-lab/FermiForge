#!/usr/bin/env bash
# Tighten nonlinear convergence without changing the discretized problem.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
parent="runs/260929-161751-single-harmonic-0plus-bb-continued-bb-radial-style-101-200"
export OMPI_MCA_pml="${OMPI_MCA_pml:-ob1}"
export OMPI_MCA_btl="${OMPI_MCA_btl:-self,sm}"
exec caffeinate -i python3 tools/run_axisymmetric_core_benchmark.py \
  --input "$parent/input.nml" \
  --restart "$parent/final_fields_2d.dat" \
  --history-parent-run "$parent" \
  --build-dir work/bb-tight-continuation-build \
  --case-name single-harmonic-0plus-bb-tight \
  --initialization tolerance-2e-8 \
  --ranks 10 --max-iterations 100 --checkpoint-interval 1 \
  --tolerance 2e-8 "$@"
