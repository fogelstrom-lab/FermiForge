#!/usr/bin/env python3
"""Plot double-core order-parameter profiles in the style of dcvlong Fig. 4."""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MATPLOTLIB_CACHE = PROJECT_ROOT / "work" / "matplotlib"
MATPLOTLIB_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MATPLOTLIB_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_2d_fields import (COMPONENTS, HARMONIC_COMPONENTS, FieldMap, PlotError,
                           load_field_map, cartesian_to_harmonics)


ROW_COLOURS = ("#000000", "#D55E00", "#0072B2")
COLUMN_LINESTYLES = ("-", "--", "-.")


@dataclass(frozen=True)
class AxisProfile:
    coordinate: np.ndarray
    values: np.ndarray


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot normalized real and imaginary A_mn profiles on the positive "
            "x and y symmetry axes, corresponding to Fig. 4 of dcvlong.pdf."
        )
    )
    parser.add_argument("field_map", type=Path, help="field map to plot")
    parser.add_argument('--basis', choices=('cartesian', 'harmonic'), default='cartesian',
                        help='local components, without removing vortex/angular phase')
    parser.add_argument(
        "--reference-field",
        type=Path,
        help="optional field map shown as a first comparison row",
    )
    parser.add_argument("--label", default="relaxed state")
    parser.add_argument("--title", default="Double-core symmetry-axis order-parameter profiles")
    parser.add_argument("--reference-label", default="initial seed")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--bulk-gap",
        type=float,
        help="normalization amplitude (default: robust outer-boundary estimate)",
    )
    parser.add_argument(
        "--visibility-threshold",
        type=float,
        default=2.0e-3,
        help=(
            "omit Re/Im channels whose maximum is below this fraction of "
            "the normalization amplitude (default: 2e-3)"
        ),
    )
    parser.add_argument("--dpi", type=int, default=240)
    return parser.parse_args(argv)


def zero_index(coordinate: np.ndarray, name: str) -> int:
    index = int(np.argmin(np.abs(coordinate)))
    spacing = np.diff(coordinate)
    scale = float(np.min(np.abs(spacing))) if spacing.size else 1.0
    tolerance = max(512.0 * np.finfo(float).eps, 1.0e-8 * scale)
    if abs(float(coordinate[index])) > tolerance:
        raise PlotError(f"the field mesh has no {name}=0 symmetry axis")
    return index


def positive_axis_profiles(field: FieldMap, basis: str = 'cartesian') -> tuple[AxisProfile, AxisProfile]:
    values = (cartesian_to_harmonics(field.order_parameter)
              if basis == 'harmonic' else field.order_parameter)
    x_zero = zero_index(field.x, "x")
    y_zero = zero_index(field.y, "y")
    x_mask = field.x >= -512.0 * np.finfo(float).eps
    y_mask = field.y >= -512.0 * np.finfo(float).eps
    x_profile = AxisProfile(
        coordinate=field.x[x_mask],
        values=values[:, :, y_zero, :][:, :, x_mask],
    )
    y_profile = AxisProfile(
        coordinate=field.y[y_mask],
        values=values[:, :, :, x_zero][:, :, y_mask],
    )
    return x_profile, y_profile


def outer_bulk_amplitude(field: FieldMap) -> float:
    pair_amplitude = np.sqrt(
        np.sum(np.abs(field.order_parameter) ** 2, axis=(0, 1)) / 3.0
    )
    boundary = np.concatenate(
        (
            pair_amplitude[0, :],
            pair_amplitude[-1, :],
            pair_amplitude[1:-1, 0],
            pair_amplitude[1:-1, -1],
        )
    )
    estimate = float(np.median(boundary[np.isfinite(boundary)]))
    if not np.isfinite(estimate) or estimate <= 0.0:
        raise PlotError("could not estimate a positive outer bulk amplitude")
    return estimate


def plot_profile_panel(
    axis: plt.Axes,
    profile: AxisProfile,
    normalization: float,
    threshold: float,
    reverse_axis: bool,
    basis: str = 'cartesian',
) -> float:
    largest = 0.0
    line_count = 0
    marker_stride = max(1, profile.coordinate.size // 10)
    components = HARMONIC_COMPONENTS if basis == 'harmonic' else COMPONENTS
    symbol = 'C' if basis == 'harmonic' else 'A'
    for name, spin, orbital in components:
        normalized = profile.values[spin, orbital] / normalization
        for part_name, part, marker in (
            ("Re", normalized.real, None),
            ("Im", normalized.imag, "o"),
        ):
            peak = float(np.max(np.abs(part)))
            largest = max(largest, peak)
            if peak < threshold:
                continue
            axis.plot(
                profile.coordinate,
                part,
                color=ROW_COLOURS[spin],
                linestyle=COLUMN_LINESTYLES[orbital],
                linewidth=1.65,
                marker=marker,
                markersize=2.6 if marker else 0.0,
                markevery=marker_stride if marker else None,
                markerfacecolor="white" if marker else None,
                markeredgewidth=0.7 if marker else None,
                label=rf"{part_name} ${symbol}_{{{name}}}$",
            )
            line_count += 1
    if line_count == 0:
        raise PlotError("no order-parameter channel exceeded the visibility threshold")
    axis.axhline(0.0, color="0.72", linewidth=0.7, zorder=0)
    axis.grid(True, axis="y", color="0.90", linewidth=0.6)
    if reverse_axis:
        axis.set_xlim(float(profile.coordinate[-1]), 0.0)
    else:
        axis.set_xlim(0.0, float(profile.coordinate[-1]))
    axis.tick_params(direction="in", top=True, right=True)
    axis.legend(
        loc="best",
        fontsize=7.0,
        frameon=True,
        framealpha=0.90,
        ncol=2,
        handlelength=2.7,
        columnspacing=0.8,
    )
    return largest


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    if arguments.visibility_threshold < 0.0:
        raise SystemExit("visibility-threshold cannot be negative")
    fields: list[tuple[str, FieldMap]] = []
    if arguments.reference_field is not None:
        fields.append(
            (arguments.reference_label, load_field_map(arguments.reference_field))
        )
    fields.append((arguments.label, load_field_map(arguments.field_map)))
    normalization = (
        arguments.bulk_gap
        if arguments.bulk_gap is not None
        else outer_bulk_amplitude(fields[-1][1])
    )
    if normalization <= 0.0:
        raise SystemExit("bulk-gap must be positive")

    figure, axes = plt.subplots(
        len(fields),
        2,
        figsize=(12.0, 4.25 * len(fields)),
        squeeze=False,
        sharey=True,
        constrained_layout=True,
    )
    largest = 1.0
    for row, (label, field) in enumerate(fields):
        x_profile, y_profile = positive_axis_profiles(field, arguments.basis)
        largest = max(
            largest,
            plot_profile_panel(
                axes[row, 0],
                x_profile,
                normalization,
                arguments.visibility_threshold,
                reverse_axis=True,
                basis=arguments.basis,
            ),
            plot_profile_panel(
                axes[row, 1],
                y_profile,
                normalization,
                arguments.visibility_threshold,
                reverse_axis=False,
                basis=arguments.basis,
            ),
        )
        axes[row, 0].set_ylabel(r"$C_{s k}/\Delta_{\rm bulk}$" if arguments.basis == 'harmonic'
                                else r"$A_{mn}/\Delta_{\rm bulk}$")
        axes[row, 0].set_title(f"{label}: positive x axis")
        axes[row, 1].set_title(f"{label}: positive y axis")
        axes[row, 0].set_xlabel(r"$x/\xi_0$")
        axes[row, 1].set_xlabel(r"$y/\xi_0$")

    limit = 1.08 * largest
    for axis_row in axes:
        for axis in axis_row:
            axis.set_ylim(-limit, limit)
    figure.suptitle(
        arguments.title + "\n" +
        (r"Local harmonics: spin, orbital = $(+,0,-)$; angular phase retained; "
         if arguments.basis == 'harmonic' else 'FermiForge counterpart of dcvlong Fig. 4; ') +
        rf"$\Delta_{{\rm bulk}}={normalization:.6g}$",
        fontsize=12,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output, dpi=arguments.dpi)
    plt.close(figure)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
