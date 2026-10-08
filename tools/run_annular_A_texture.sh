#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
texture="${1:-a_mermin_ho}"
if (($#)); then shift; fi
python3 tools/run_modern_cylinder.py \
  --initialization -2 --texture "$texture" --spatial-mode full_2d \
  --disk-grid annular --disk-radial-layout core_wall \
  --disk-rings 24 --disk-stretch 1.5 --disk-tangent-spacing 1 \
  --radius 10 --fs1 0 --azimuths 32 --ranks 10 --iterations 20 \
  --tolerance 2e-7 "$@"
