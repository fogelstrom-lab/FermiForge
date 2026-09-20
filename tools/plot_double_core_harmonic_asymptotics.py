#!/usr/bin/env python3
"""Plot double-core harmonic tails and their 1/r + 1/r^2 continuation."""

from __future__ import annotations

import argparse
import json
import os
import re
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

from plot_2d_fields import (
    HARMONIC_COMPONENTS,
    FieldMap,
    PlotError,
    cartesian_to_harmonics,
    load_field_map,
)


REAL_COLOUR = "#000000"
IMAGINARY_COLOUR = "#0072B2"
REAL_FIT_COLOUR = "#D55E00"
IMAGINARY_FIT_COLOUR = "#56B4E9"
X_COLOUR = "#332288"
Y_COLOUR = "#CC6677"


@dataclass(frozen=True)
class TailParameters:
    bulk_gap: float
    bulk_phase: float
    bulk_rotation: np.ndarray
    fit_inner_radius: float
    matching_radius: float
    active_radius_x: float
    active_radius_y: float


@dataclass(frozen=True)
class HarmonicTail:
    axis_name: str
    radius: np.ndarray
    departure: np.ndarray
    coefficient_one: np.ndarray
    coefficient_two: np.ndarray
    fitted_departure: np.ndarray
    active_radius: float


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("field_map", type=Path)
    parser.add_argument(
        "--input",
        type=Path,
        help="run input namelist (default: input.nml beside the field map)",
    )
    parser.add_argument(
        "--metrics",
        type=Path,
        help="run metrics file (default: metrics.txt beside the field map)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="output directory (default: plots beside the field map)",
    )
    parser.add_argument(
        "--minimum-tail-radius",
        type=float,
        help=(
            "smallest radius shown in scaled-tail panels "
            "(default: half the smaller active radius)"
        ),
    )
    parser.add_argument("--dpi", type=int, default=240)
    return parser.parse_args(argv)


def namelist_float(text: str, key: str) -> float:
    match = re.search(
        rf"(?im)^\s*{re.escape(key)}\s*=\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)",
        text,
    )
    if match is None:
        raise PlotError(f"input namelist does not define '{key}'")
    return float(match.group(1).replace("d", "e").replace("D", "E"))


def read_metrics(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.strip()] = value.strip()
    return result


def metric_float(metrics: dict[str, str], key: str) -> float:
    try:
        return float(metrics[key].split()[0].replace("d", "e").replace("D", "E"))
    except (KeyError, ValueError, IndexError) as error:
        raise PlotError(f"metrics do not define numeric '{key}'") from error


def metric_row(metrics: dict[str, str], key: str) -> np.ndarray:
    try:
        values = np.asarray(
            [float(value.replace("d", "e").replace("D", "E")) for value in metrics[key].split()],
            dtype=float,
        )
    except (KeyError, ValueError) as error:
        raise PlotError(f"metrics do not define numeric row '{key}'") from error
    if values.shape != (3,):
        raise PlotError(f"metrics row '{key}' does not contain three values")
    return values


def load_parameters(input_path: Path, metrics_path: Path) -> TailParameters:
    input_text = input_path.read_text(encoding="utf-8")
    metrics = read_metrics(metrics_path)
    if metrics.get("update_region") == "centered_ellipse":
        active_radius_x = metric_float(metrics, "update_radius_x")
        active_radius_y = metric_float(metrics, "update_radius_y")
    else:
        active_radius_x = metric_float(metrics, "update_radius")
        active_radius_y = active_radius_x
    rotation = np.vstack(
        [metric_row(metrics, f"bulk_rotation_row_{row}") for row in range(1, 4)]
    )
    return TailParameters(
        bulk_gap=namelist_float(input_text, "asymptotic_bulk_gap"),
        bulk_phase=metric_float(metrics, "bulk_reference_phase_offset"),
        bulk_rotation=rotation,
        fit_inner_radius=namelist_float(input_text, "asymptotic_fit_inner_radius"),
        matching_radius=namelist_float(input_text, "asymptotic_matching_radius"),
        active_radius_x=active_radius_x,
        active_radius_y=active_radius_y,
    )


def zero_index(coordinates: np.ndarray, name: str) -> int:
    index = int(np.argmin(np.abs(coordinates)))
    scale = max(float(np.max(np.abs(coordinates))), 1.0)
    if abs(float(coordinates[index])) > 512.0 * np.finfo(float).eps * scale:
        raise PlotError(f"field mesh has no {name}=0 axis")
    return index


def harmonic_matrix(matrix: np.ndarray) -> np.ndarray:
    expanded = matrix[:, :, np.newaxis, np.newaxis]
    return cartesian_to_harmonics(expanded)[:, :, 0, 0]


def interpolate_matrix(radius: np.ndarray, values: np.ndarray, point: float) -> np.ndarray:
    if point < float(radius[0]) or point > float(radius[-1]):
        raise PlotError(f"fit radius {point:g} is outside the field axis")
    result = np.empty((3, 3), dtype=np.complex128)
    for spin in range(3):
        for orbital in range(3):
            result[spin, orbital] = np.interp(
                point, radius, values[spin, orbital].real
            ) + 1j * np.interp(point, radius, values[spin, orbital].imag)
    return result


def harmonic_decay_power(spin: int, orbital: int) -> int:
    """Return the symmetry-allowed power from dcvlong Eqs. (19)-(20)."""
    projection_zero = 1
    mixed_zero_channel = (spin == projection_zero) != (
        orbital == projection_zero
    )
    return 1 if mixed_zero_channel else 2


def make_tail(
    field: FieldMap,
    harmonics: np.ndarray,
    parameters: TailParameters,
    axis_name: str,
) -> HarmonicTail:
    x_zero = zero_index(field.x, "x")
    y_zero = zero_index(field.y, "y")
    if axis_name == "x":
        coordinates = field.x
        values = harmonics[:, :, y_zero, :]
        angle = 0.0
        active_radius = parameters.active_radius_x
    elif axis_name == "y":
        coordinates = field.y
        values = harmonics[:, :, :, x_zero]
        angle = 0.5 * np.pi
        active_radius = parameters.active_radius_y
    else:
        raise PlotError(f"unknown axis '{axis_name}'")

    positive = coordinates > 512.0 * np.finfo(float).eps
    radius = coordinates[positive]
    values = values[:, :, positive]
    bulk_cartesian = (
        parameters.bulk_gap
        * np.exp(1j * (angle + parameters.bulk_phase))
        * parameters.bulk_rotation
    )
    bulk_harmonic = harmonic_matrix(bulk_cartesian)
    departure = (values - bulk_harmonic[:, :, np.newaxis]) / parameters.bulk_gap

    outer = interpolate_matrix(radius, departure, parameters.matching_radius)
    coefficient_one = np.zeros((3, 3), dtype=np.complex128)
    coefficient_two = np.zeros((3, 3), dtype=np.complex128)
    for spin in range(3):
        for orbital in range(3):
            if harmonic_decay_power(spin, orbital) == 1:
                coefficient_one[spin, orbital] = outer[spin, orbital]
            else:
                coefficient_two[spin, orbital] = outer[spin, orbital]
    radial_ratio = parameters.matching_radius / radius
    fitted = (
        coefficient_one[:, :, np.newaxis] * radial_ratio
        + coefficient_two[:, :, np.newaxis] * radial_ratio**2
    )
    return HarmonicTail(
        axis_name=axis_name,
        radius=radius,
        departure=departure,
        coefficient_one=coefficient_one,
        coefficient_two=coefficient_two,
        fitted_departure=fitted,
        active_radius=active_radius,
    )


def plot_departure_profiles(
    x_tail: HarmonicTail, y_tail: HarmonicTail, parameters: TailParameters
) -> plt.Figure:
    figure, axes = plt.subplots(
        3, 3, figsize=(12.0, 9.8), sharex=True, sharey=True, constrained_layout=True
    )
    positive_values: list[np.ndarray] = []
    for tail in (x_tail, y_tail):
        positive_values.append(np.abs(tail.departure[np.abs(tail.departure) > 0.0]))
    concatenated = np.concatenate([values for values in positive_values if values.size])
    lower = max(float(np.nanpercentile(concatenated, 1.0)) * 0.2, 1.0e-8)
    upper = max(float(np.nanmax(concatenated)) * 1.2, 10.0 * lower)

    for name, spin, orbital in HARMONIC_COMPONENTS:
        axis = axes[spin, orbital]
        axis.plot(
            x_tail.radius,
            np.abs(x_tail.departure[spin, orbital]),
            color=X_COLOUR,
            linewidth=1.5,
            marker="o",
            markersize=2.4,
            markevery=max(1, x_tail.radius.size // 18),
            label="positive x axis",
        )
        axis.plot(
            y_tail.radius,
            np.abs(y_tail.departure[spin, orbital]),
            color=Y_COLOUR,
            linewidth=1.5,
            linestyle="--",
            marker="s",
            markersize=2.4,
            markevery=max(1, y_tail.radius.size // 18),
            label="positive y axis",
        )
        axis.axvline(x_tail.active_radius, color=X_COLOUR, linewidth=0.8, alpha=0.7)
        axis.axvline(y_tail.active_radius, color=Y_COLOUR, linewidth=0.8, alpha=0.7)
        axis.axvspan(
            parameters.fit_inner_radius,
            parameters.matching_radius,
            color="#E69F00",
            alpha=0.10,
        )
        axis.axvline(parameters.matching_radius, color="0.35", linewidth=0.8)
        axis.set_yscale("log")
        axis.set_ylim(lower, upper)
        axis.set_xlim(0.0, max(float(x_tail.radius[-1]), float(y_tail.radius[-1])))
        axis.grid(True, which="major", color="0.90", linewidth=0.55)
        axis.tick_params(direction="in", top=True, right=True)
        axis.set_title(rf"$|\delta A_{{{name}}}|/\Delta_{{\rm bulk}}$")
        if spin == 2:
            axis.set_xlabel(r"$r/\xi_0$")
        if orbital == 0:
            axis.set_ylabel("harmonic departure")
    axes[0, 0].legend(loc="best", fontsize=7.4, framealpha=0.92)
    figure.suptitle(
        "Bulk-subtracted harmonic profiles on the positive symmetry axes\n"
        "vertical coloured lines: active radii; orange band: current solver fit window",
        fontsize=12,
    )
    return figure


def tail_error(tail: HarmonicTail, spin: int, orbital: int, minimum: float) -> float:
    lower = max(minimum, 0.60 * tail.active_radius)
    selected = (tail.radius >= lower) & (tail.radius <= tail.active_radius)
    if not np.any(selected):
        return float("nan")
    difference = tail.departure[spin, orbital, selected] - tail.fitted_departure[
        spin, orbital, selected
    ]
    reference = tail.departure[spin, orbital, selected]
    denominator = max(float(np.sqrt(np.mean(np.abs(reference) ** 2))), 1.0e-14)
    return float(np.sqrt(np.mean(np.abs(difference) ** 2)) / denominator)


def plot_scaled_tail(
    tail: HarmonicTail,
    parameters: TailParameters,
    minimum_radius: float,
) -> plt.Figure:
    selected = tail.radius >= minimum_radius
    if np.count_nonzero(selected) < 4:
        raise PlotError("too few axis points beyond the requested minimum tail radius")
    radius = tail.radius[selected]
    inverse_coordinate = parameters.matching_radius / radius

    figure, axes = plt.subplots(3, 3, figsize=(12.0, 9.8), constrained_layout=True)
    for name, spin, orbital in HARMONIC_COMPONENTS:
        axis = axes[spin, orbital]
        power = harmonic_decay_power(spin, orbital)
        scaled = tail.departure[spin, orbital, selected] / (
            inverse_coordinate**power
        )
        coefficient = (
            tail.coefficient_one[spin, orbital]
            if power == 1
            else tail.coefficient_two[spin, orbital]
        )
        fit_scaled = np.full_like(scaled, coefficient)
        axis.plot(
            inverse_coordinate,
            scaled.real,
            color=REAL_COLOUR,
            linewidth=1.25,
            marker="o",
            markersize=3.0,
            markerfacecolor="white",
            label="Re data",
        )
        axis.plot(
            inverse_coordinate,
            scaled.imag,
            color=IMAGINARY_COLOUR,
            linewidth=1.25,
            linestyle="--",
            marker="s",
            markersize=2.8,
            markerfacecolor="white",
            label="Im data",
        )
        axis.plot(
            inverse_coordinate,
            fit_scaled.real,
            color=REAL_FIT_COLOUR,
            linewidth=1.5,
            label=rf"Re $C_{power}$",
        )
        axis.plot(
            inverse_coordinate,
            fit_scaled.imag,
            color=IMAGINARY_FIT_COLOUR,
            linewidth=1.5,
            linestyle="--",
            label=rf"Im $C_{power}$",
        )
        active_coordinate = parameters.matching_radius / tail.active_radius
        inner_coordinate = parameters.matching_radius / parameters.fit_inner_radius
        axis.axvline(active_coordinate, color="#CC6677", linewidth=0.9)
        axis.axvline(inner_coordinate, color="#E69F00", linewidth=0.9)
        axis.axvline(1.0, color="0.35", linewidth=0.9)
        axis.axhline(0.0, color="0.78", linewidth=0.6, zorder=0)
        axis.grid(True, color="0.91", linewidth=0.55)
        axis.tick_params(direction="in", top=True, right=True)
        axis.set_xlim(float(np.max(inverse_coordinate)), float(np.min(inverse_coordinate)))
        mismatch = tail_error(tail, spin, orbital, minimum_radius)
        axis.set_title(
            rf"$A_{{{name}}}\sim r^{{-{power}}}$; $\epsilon={mismatch:.2g}$"
        )
        if spin == 2:
            axis.set_xlabel(r"$R_c/r$  (outward $\longrightarrow$)")
        if orbital == 0:
            axis.set_ylabel(
                rf"$(r/R_c)^{{{power}}}\,\delta A/\Delta_{{\rm bulk}}$"
            )
    axes[0, 0].legend(loc="best", fontsize=6.5, framealpha=0.92)
    figure.suptitle(
        f"{tail.axis_name}-axis harmonic asymptotic diagnostic\n"
        r"dcvlong Eqs. (19)-(20): the correctly scaled channel is constant; "
        r"$\epsilon$ tests the active outer tail",
        fontsize=12,
    )
    return figure


def save_figure(figure: plt.Figure, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def diagnostic_report(
    x_tail: HarmonicTail,
    y_tail: HarmonicTail,
    parameters: TailParameters,
    minimum_radius: float,
) -> dict[str, object]:
    axes: dict[str, object] = {}
    for tail in (x_tail, y_tail):
        components: dict[str, object] = {}
        for name, spin, orbital in HARMONIC_COMPONENTS:
            components[name] = {
                "decay_power": harmonic_decay_power(spin, orbital),
                "coefficient_1_re": float(tail.coefficient_one[spin, orbital].real),
                "coefficient_1_im": float(tail.coefficient_one[spin, orbital].imag),
                "coefficient_2_re": float(tail.coefficient_two[spin, orbital].real),
                "coefficient_2_im": float(tail.coefficient_two[spin, orbital].imag),
                "active_outer_tail_relative_rms_mismatch": tail_error(
                    tail, spin, orbital, minimum_radius
                ),
            }
        axes[tail.axis_name] = {
            "active_radius": tail.active_radius,
            "components": components,
        }
    return {
        "kind": "fermiforge_double_core_harmonic_asymptotic_diagnostic",
        "definition": (
            "dcvlong Eqs. (19)-(20): +0, 0+, 0-, -0 use "
            "delta A=C1*(Rc/r); ++, +-, 00, -+, -- use "
            "delta A=C2*(Rc/r)^2. Coefficients are dimensionless harmonic "
            "departures normalized by the bulk gap"
        ),
        "bulk_gap": parameters.bulk_gap,
        "fit_inner_radius": parameters.fit_inner_radius,
        "matching_radius": parameters.matching_radius,
        "minimum_tail_radius": minimum_radius,
        "axes": axes,
    }


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        field_path = arguments.field_map.expanduser().resolve()
        run_directory = field_path.parent
        input_path = (arguments.input or run_directory / "input.nml").expanduser().resolve()
        metrics_path = (arguments.metrics or run_directory / "metrics.txt").expanduser().resolve()
        output_directory = (
            arguments.output_dir or run_directory / "plots"
        ).expanduser().resolve()
        parameters = load_parameters(input_path, metrics_path)
        minimum_tail_radius = (
            arguments.minimum_tail_radius
            if arguments.minimum_tail_radius is not None
            else 0.5
            * min(parameters.active_radius_x, parameters.active_radius_y)
        )
        if minimum_tail_radius <= 0.0:
            raise PlotError("minimum-tail-radius must be positive")
        field = load_field_map(field_path)
        harmonics = cartesian_to_harmonics(field.order_parameter)
        x_tail = make_tail(field, harmonics, parameters, "x")
        y_tail = make_tail(field, harmonics, parameters, "y")

        outputs = (
            output_directory / "double_core_harmonic_departure_profiles.png",
            output_directory / "double_core_harmonic_asymptotics_x.png",
            output_directory / "double_core_harmonic_asymptotics_y.png",
        )
        save_figure(
            plot_departure_profiles(x_tail, y_tail, parameters),
            outputs[0],
            arguments.dpi,
        )
        save_figure(
            plot_scaled_tail(
                x_tail, parameters, minimum_tail_radius
            ),
            outputs[1],
            arguments.dpi,
        )
        save_figure(
            plot_scaled_tail(
                y_tail, parameters, minimum_tail_radius
            ),
            outputs[2],
            arguments.dpi,
        )
        report_path = output_directory / "harmonic_tail_diagnostics.json"
        report_path.write_text(
            json.dumps(
                diagnostic_report(
                    x_tail, y_tail, parameters, minimum_tail_radius
                ),
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    except (OSError, ValueError, PlotError) as error:
        print(f"error: {error}", file=__import__("sys").stderr)
        return 2

    for path in outputs:
        print(f"wrote {path}")
    print(f"wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
