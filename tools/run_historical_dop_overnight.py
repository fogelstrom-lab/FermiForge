#!/usr/bin/env python3
"""Run the historical dop seed on the large overnight 2D mesh."""

from __future__ import annotations

import sys
from pathlib import Path

from run_double_core_from_scratch import main


PROJECT_ROOT = Path(__file__).resolve().parents[1]
HISTORICAL_DOP_INPUT = (
    PROJECT_ROOT / "examples" / "2d_historical_dop_overnight.nml"
)


if __name__ == "__main__":
    raise SystemExit(
        main(
            [
                "--input",
                str(HISTORICAL_DOP_INPUT),
                "--run-name",
                "historical-dop-overnight",
                *sys.argv[1:],
            ]
        )
    )
