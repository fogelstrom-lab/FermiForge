#!/usr/bin/env python3
"""Plot FermiForge 2D nonlinear convergence and Anderson mixing history."""

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
        description="Plot a FermiForge 2D self-consistency iteration history."
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
    iteration = column["iteration"]
    if "update_rms" in column and "update_max" in column:
        rms = np.maximum(column["update_rms"], np.finfo(float).tiny)
        maximum = np.maximum(column["update_max"], np.finfo(float).tiny)
        residual_label = "projected update residual"
    else:
        rms = np.maximum(column["map_rms"], np.finfo(float).tiny)
        maximum = np.maximum(column["map_max"], np.finfo(float).tiny)
        residual_label = "map residual"
    mixing = column["p"]
    stochastic = column.get("stochastic_fallback", np.zeros_like(iteration)) > 0.5

    figure, axis = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    axis.axvspan(
        iteration[0] - 0.5,
        iteration[-1] + 0.5,
        color="#efedf5",
        alpha=0.65,
        label="legacy Anderson engine",
        zorder=0,
    )
    axis.semilogy(
        iteration,
        rms,
        "o-",
        color="#2166ac",
        label=f"RMS {residual_label}",
    )
    axis.semilogy(
        iteration,
        maximum,
        "s-",
        color="#b2182b",
        label=f"maximum {residual_label}",
    )
    if np.any(stochastic):
        axis.scatter(
            iteration[stochastic],
            maximum[stochastic],
            marker="X",
            s=75,
            color="#000000",
            label="stochastic fallback",
            zorder=4,
        )
    axis.set_xlabel("nonlinear map evaluation")
    axis.set_ylabel("residual")
    axis.grid(True, which="both", color="0.86", linewidth=0.7)

    mixing_axis = axis.twinx()
    mixing_axis.plot(
        iteration, mixing, "^-", color="#1b7837", linewidth=1.2, label="mixing p"
    )
    mixing_axis.set_ylabel("Anderson mixing p", color="#1b7837")
    mixing_axis.tick_params(axis="y", colors="#1b7837")

    handles, labels = axis.get_legend_handles_labels()
    extra_handles, extra_labels = mixing_axis.get_legend_handles_labels()
    axis.legend(handles + extra_handles, labels + extra_labels, loc="best")
    axis.set_title("2D self-consistency convergence")

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output, dpi=arguments.dpi)
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
