#!/usr/bin/env bash
# Radially constrained A-phase-core vortex in a specular B-phase cylinder.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 tools/run_modern_cylinder.py \
  --initialization 1 --radius 20 --fs1 5.4 \
  --radial-points 201 --step 0.1 --half-length 160 \
  --azimuths 32 --ranks 10 --iterations 300 --tolerance 2e-7 "$@"
