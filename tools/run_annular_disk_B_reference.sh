#!/usr/bin/env bash
# Prescribed graded rings; independent 2D fields, not radial symmetry.
set -euo pipefail
cd "$(dirname "$0")/.."
reference="${1:-runs/261006-192650-696944-modern-cylinder-R10-n0-bulk-b}"
if (($#)); then shift; fi
python3 tools/run_modern_cylinder.py \
  --spatial-mode full_2d --disk-grid annular --disk-rings 16 \
  --disk-stretch 2 --disk-tangent-spacing 1.25 \
  --initialization -1 --initial-run "$reference" \
  --radius 10 --radial-points 100 --fs1 0 --azimuths 32 \
  --ranks 10 --iterations 5 --tolerance 2e-7 "$@"
