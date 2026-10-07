#!/usr/bin/env bash
# Same settings as the R20 A-core experiment; only the initial seed changes.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 tools/run_modern_cylinder.py \
  --initialization 0 --radius 20 --fs1 5.4 \
  --radial-points 201 --step 0.1 --half-length 160 \
  --azimuths 32 --ranks 10 --iterations 300 --tolerance 2e-7 "$@"
