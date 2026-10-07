#!/usr/bin/env bash
# Matched pmax=3 scratch comparison with the pmax=5 and pmax=10 runs.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bash "$project_dir/tools/run_radial_anderson_diagnostics.sh" 0plus \
  --case-name radial-anderson-aa20-0plus-pmax3-scratch-T030-Fs54 \
  --anderson-pmax 3 --tolerance 2e-6 --max-iterations 600 "$@"
