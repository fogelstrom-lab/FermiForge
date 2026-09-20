#!/usr/bin/env python3
"""Run the first testable projected, dependent-halo double-core workflow."""

from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def restart_with_plot_capable_python() -> None:
    """Prefer a local Python with Matplotlib without making it mandatory."""

    if "--no-plot" in sys.argv or "--help" in sys.argv or "-h" in sys.argv:
        return
    if importlib.util.find_spec("matplotlib") is not None:
        return

    candidates = (
        shutil.which("python3"),
        "/opt/homebrew/bin/python3",
        "/usr/local/bin/python3",
        "/opt/local/bin/python3",
    )
    current = Path(sys.executable).resolve()
    environment = os.environ.copy()
    environment.setdefault(
        "MPLCONFIGDIR",
        str(Path(tempfile.gettempdir()) / "fermiforge-matplotlib-cache"),
    )
    for candidate_name in candidates:
        if not candidate_name:
            continue
        candidate = Path(candidate_name).expanduser()
        if not candidate.is_file() or candidate.resolve() == current:
            continue
        check = subprocess.run(
            [str(candidate), "-c", "import matplotlib"],
            env=environment,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if check.returncode == 0:
            print(f"using {candidate} for Matplotlib output", file=sys.stderr)
            os.execve(
                str(candidate),
                [str(candidate), str(Path(__file__).resolve()), *sys.argv[1:]],
                environment,
            )

from run_double_core_radius_ladder import main


if __name__ == "__main__":
    restart_with_plot_capable_python()
    raise SystemExit(
        main(
            [
                "--zero-mode-projection",
                "gauge_translation_orientation",
                "--inactive-halo-policy",
                "state_asymptotic",
                *sys.argv[1:],
            ]
        )
    )
