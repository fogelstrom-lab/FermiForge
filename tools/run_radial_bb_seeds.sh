#!/usr/bin/env bash
# Run the matched scratch seeds sequentially: never two ten-rank jobs together.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case "${1:-both}" in
  both) seeds=(0plus plus0) ;;
  0plus|plus0) seeds=("$1") ;;
  *) echo "Usage: bash $0 [both|0plus|plus0]" >&2; exit 2 ;;
esac
for seed in "${seeds[@]}"; do
  echo "Starting scratch radial BB: $seed, T/Tc=0.30, Fs1=5.4, 60 radial points"
  bash "$project_dir/tools/run_radial_symmetry.sh" \
    --input "$project_dir/examples/radial_bb_localized_${seed}.nml" \
    --build-dir "$project_dir/work/radial-bb-seeds-build" \
    --case-name "radial-bb-scratch-${seed}-T030-Fs54"
done
