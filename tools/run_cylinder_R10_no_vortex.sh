#!/usr/bin/env bash
# Isolated, timestamped cylinder tests; no changes to new_src/qcv.inp or data.
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
mode="${1:-relax}"
ranks="${2:-10}"
case "$mode" in
  smoke) input=radial_cylinder_R10_no_vortex_smoke.inp ;;
  relax) input=radial_cylinder_R10_no_vortex.inp ;;
  angular48) input=radial_cylinder_R10_no_vortex_angular48.inp ;;
  *) echo 'Usage: bash tools/run_cylinder_R10_no_vortex.sh [smoke|relax|angular48] [ranks]' >&2; exit 2 ;;
esac
if [[ "$(uname -s)" == Darwin ]]; then
  # Local Mac Open MPI 5: use shared memory; leave explicit user settings intact.
  export OMPI_MCA_pml="${OMPI_MCA_pml:-ob1}"
  export OMPI_MCA_btl="${OMPI_MCA_btl:-self,sm}"
fi
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
cd "$project_dir"
exec python3 tools/run_and_plot.py --solver new-src \
  --input "examples/$input" --ranks "$ranks" \
  --case-name "cylinder-R10-n0-T030-Fs0-$mode"
