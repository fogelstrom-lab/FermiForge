#!/usr/bin/env bash
# Same map and accelerators as the 2D runner; only the +x ray is independent.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
export OMPI_MCA_pml="${OMPI_MCA_pml:-ob1}"
export OMPI_MCA_btl="${OMPI_MCA_btl:-self,sm}"
exec caffeinate -i python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/radial_symmetry_a_phase.nml \
  --build-dir work/radial-symmetry-build \
  --case-name a-phase-radial-symmetry --ranks 10 "$@"
