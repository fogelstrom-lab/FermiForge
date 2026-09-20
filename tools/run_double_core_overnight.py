#!/usr/bin/env python3
"""Run the larger, finer FermiForge double-core overnight calculation."""

from __future__ import annotations

import sys
from pathlib import Path

from run_double_core_from_scratch import main


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OVERNIGHT_INPUT = PROJECT_ROOT / "examples" / "2d_double_core_overnight.nml"


if __name__ == "__main__":
    raise SystemExit(
        main(
            [
                "--input",
                str(OVERNIGHT_INPUT),
                "--run-name",
                "double-core-overnight",
                *sys.argv[1:],
            ]
        )
    )
