#!/usr/bin/env python3
"""Plot whether the two hard cores remain separated during iteration."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MATPLOTLIB_CACHE = PROJECT_ROOT / "work" / "matplotlib"
MATPLOTLIB_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MATPLOTLIB_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot double-core positions and amplitudes versus iteration."
    )
    parser.add_argument("history", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=220)
    arguments = parser.parse_args()

    lines = arguments.history.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("#"):
        raise SystemExit(f"missing history header in {arguments.history}")
    names = lines[0][1:].split()
    data = np.loadtxt(arguments.history, comments="#", ndmin=2)
    if data.shape[1] != len(names):
        raise SystemExit(
            f"history declares {len(names)} columns but contains {data.shape[1]}"
        )
    column = {name: data[:, index] for index, name in enumerate(names)}
    required = {
        "iteration",
        "half_core_y_minus",
        "half_core_y_plus",
        "half_core_separation",
        "center_pair_amplitude",
        "half_core_amplitude_minus",
        "half_core_amplitude_plus",
    }
    missing = sorted(required.difference(column))
    if missing:
        raise SystemExit(f"history is missing columns: {', '.join(missing)}")

    iteration = column["iteration"]
    figure, axes = plt.subplots(2, 1, figsize=(8.6, 7.0), sharex=True,
                                constrained_layout=True)
    axes[0].plot(
        iteration,
        column["half_core_y_plus"],
        "o-",
        color="#0072B2",
        label=r"positive minimum $y_+$",
    )
    axes[0].plot(
        iteration,
        column["half_core_y_minus"],
        "s-",
        color="#D55E00",
        label=r"negative minimum $y_-$",
    )
    axes[0].plot(
        iteration,
        0.5 * column["half_core_separation"],
        "^-",
        color="#000000",
        label=r"measured $a/2$",
    )
    axes[0].axhline(0.0, color="0.65", linewidth=0.7)
    axes[0].set_ylabel(r"position / $\xi_0$")
    axes[0].set_title("Double-core separation diagnostic")
    axes[0].legend(loc="best")

    axes[1].plot(
        iteration,
        column["center_pair_amplitude"],
        "^-",
        color="#000000",
        label="center",
    )
    axes[1].plot(
        iteration,
        column["half_core_amplitude_plus"],
        "o-",
        color="#0072B2",
        label="positive minimum",
    )
    axes[1].plot(
        iteration,
        column["half_core_amplitude_minus"],
        "s--",
        color="#D55E00",
        label="negative minimum",
    )
    axes[1].set_xlabel("completed Anderson update")
    axes[1].set_ylabel("pair amplitude")
    axes[1].legend(loc="best")
    for axis in axes:
        axis.grid(True, color="0.87", linewidth=0.7)
        axis.tick_params(direction="in", top=True, right=True)

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output, dpi=arguments.dpi)
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
