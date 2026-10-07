#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
seed="${1:-0plus}"
case "$seed" in
  0plus|plus0) ;;
  *) echo "Usage: bash $0 [0plus|plus0]" >&2; exit 2 ;;
esac
bash "$project_dir/tools/run_radial_symmetry.sh" \
  --input "$project_dir/examples/radial_bb_localized_${seed}.nml" \
  --build-dir "$project_dir/work/radial-bb-diagnostics-build" \
  --case-name "radial-bb-diagnostics-${seed}-T030-Fs54" \
  --save-iteration-diagnostics
