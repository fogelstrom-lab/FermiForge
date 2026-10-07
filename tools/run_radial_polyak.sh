#!/usr/bin/env bash
# Matched scratch comparison with radial BB; separate build and dated outputs.
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
  --build-dir "$project_dir/work/radial-polyak-build" \
  --case-name "radial-polyak-scratch-${seed}-T030-Fs54" \
  --iteration-method polyak --anderson-cycle 0 \
  --polyak-step-size 2.0 --polyak-drag 0.5 \
  --no-save-iteration-diagnostics "$@"
