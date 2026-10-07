#!/usr/bin/env bash
# Matched scratch comparison: only AA cap and tight stopping criterion differ.
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bash "$project_dir/tools/run_radial_anderson_diagnostics.sh" 0plus \
  --case-name radial-anderson-aa20-0plus-pmax10-scratch-T030-Fs54 \
  --anderson-pmax 10 --tolerance 2e-6 --max-iterations 600 "$@"
