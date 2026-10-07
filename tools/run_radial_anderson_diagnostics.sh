#!/usr/bin/env bash
# Matched scratch seed and physics; AA20 restarts, no simple steps.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
seed="${1:-0plus}"
if [[ $# -gt 0 ]]; then shift; fi
case "$seed" in
  0plus|plus0) ;;
  *) echo "Usage: bash $0 [0plus|plus0] [runner options]" >&2; exit 2 ;;
esac
bash "$project_dir/tools/run_radial_symmetry.sh" \
  --input "$project_dir/examples/radial_bb_localized_${seed}.nml" \
  --build-dir "$project_dir/work/radial-anderson-diagnostics-build" \
  --case-name "radial-anderson-aa20-diagnostics-${seed}-T030-Fs54" \
  --iteration-method anderson --anderson-cycle 20 --simple-cycle 0 \
  --anderson-history 10 --anderson-pmax 5 --anderson-progress 0.1 \
  --max-iterations 600 --tolerance 2e-5 --checkpoint-interval 1 \
  --save-iteration-diagnostics "$@"
