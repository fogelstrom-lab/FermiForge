#!/usr/bin/env bash
# Independent 2D normal-core cylinder, initialized from the converged radial run.
set -euo pipefail
cd "$(dirname "$0")/.."
reference="${1:-runs/261007-103929-261905-modern-cylinder-R20-n1-normal-core-Fs5.4}"
if (($#)); then shift; fi
python3 tools/run_modern_cylinder.py \
  --spatial-mode full_2d --disk-grid annular --disk-radial-layout core_wall \
  --disk-rings 32 --disk-stretch 2 --disk-tangent-spacing 1.25 \
  --initialization 0 --initial-run "$reference" --radius 20 --fs1 5.4 \
  --radial-points 201 --step 0.1 --half-length 160 \
  --azimuths 32 --ranks 10 --iterations 5 --tolerance 2e-7 "$@"
