#!/usr/bin/env python3
"""Compare passive off-domain quasiclassical maps with the imposed tails."""

from __future__ import annotations

import argparse
import json
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

from plot_2d_fields import HARMONIC_COMPONENTS, PlotError, cartesian_to_harmonics
from plot_double_core_harmonic_asymptotics import load_parameters, read_metrics


POSITIVE_COLOUR = "#332288"
NEGATIVE_COLOUR = "#D55E00"
FLOOR = 1.0e-14


@dataclass(frozen=True)
class ProbeData:
    coordinate: np.ndarray
    ray_id: np.ndarray
    input_order_parameter: np.ndarray
    mapped_order_parameter: np.ndarray
    input_mean_field: np.ndarray
    mapped_mean_field: np.ndarray


@dataclass(frozen=True)
class ShadowProbeData:
    coordinate: np.ndarray
    ray_id: np.ndarray
    baseline_order_parameter: np.ndarray
    first_mapped_order_parameter: np.ndarray
    shadow_order_parameter: np.ndarray
    final_mapped_order_parameter: np.ndarray
    baseline_mean_field: np.ndarray
    first_mapped_mean_field: np.ndarray
    shadow_mean_field: np.ndarray
    final_mapped_mean_field: np.ndarray


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe_file", type=Path)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shadow-file", type=Path)
    parser.add_argument("--shadow-history", type=Path)
    parser.add_argument("--dpi", type=int, default=240)
    return parser.parse_args(argv)


def complex_matrices(rows: np.ndarray, start: int) -> np.ndarray:
    pairs = rows[:, start : start + 18].reshape((-1, 3, 3, 2))
    return pairs[..., 0] + 1j * pairs[..., 1]


def load_probes(path: Path) -> ProbeData:
    rows = np.loadtxt(path, comments="#", ndmin=2)
    if rows.shape[1] != 45:
        raise PlotError(
            f"asymptotic probe file has {rows.shape[1]} columns; expected 45"
        )
    ray_id = np.rint(rows[:, 2]).astype(int)
    if not set(ray_id).issubset({-2, -1, 1, 2}):
        raise PlotError("asymptotic probe file contains an unknown ray id")
    return ProbeData(
        coordinate=rows[:, 0:2],
        ray_id=ray_id,
        input_order_parameter=complex_matrices(rows, 3),
        mapped_order_parameter=complex_matrices(rows, 21),
        input_mean_field=rows[:, 39:42],
        mapped_mean_field=rows[:, 42:45],
    )


def load_shadow_probes(path: Path) -> ShadowProbeData:
    rows = np.loadtxt(path, comments="#", ndmin=2)
    if rows.shape[1] != 87:
        raise PlotError(
            f"asymptotic shadow file has {rows.shape[1]} columns; expected 87"
        )
    ray_id = np.rint(rows[:, 2]).astype(int)
    if not set(ray_id).issubset({-2, -1, 1, 2}):
        raise PlotError("asymptotic shadow file contains an unknown ray id")

    def state(block: int) -> tuple[np.ndarray, np.ndarray]:
        start = 3 + 21 * block
        return complex_matrices(rows, start), rows[:, start + 18 : start + 21]

    baseline_order, baseline_mean = state(0)
    first_order, first_mean = state(1)
    shadow_order, shadow_mean = state(2)
    final_order, final_mean = state(3)
    return ShadowProbeData(
        coordinate=rows[:, 0:2],
        ray_id=ray_id,
        baseline_order_parameter=baseline_order,
        first_mapped_order_parameter=first_order,
        shadow_order_parameter=shadow_order,
        final_mapped_order_parameter=final_order,
        baseline_mean_field=baseline_mean,
        first_mapped_mean_field=first_mean,
        shadow_mean_field=shadow_mean,
        final_mapped_mean_field=final_mean,
    )


def harmonic_matrices(matrices: np.ndarray) -> np.ndarray:
    expanded = np.transpose(matrices, (1, 2, 0))[:, :, np.newaxis, :]
    harmonics = cartesian_to_harmonics(expanded)[:, :, 0, :]
    return np.transpose(harmonics, (2, 0, 1))


def bulk_harmonics(
    coordinates: np.ndarray,
    bulk_gap: float,
    bulk_phase: float,
    bulk_rotation: np.ndarray,
) -> np.ndarray:
    angle = np.arctan2(coordinates[:, 1], coordinates[:, 0])
    phase = np.exp(1j * (angle + bulk_phase))
    matrices = bulk_gap * phase[:, np.newaxis, np.newaxis] * bulk_rotation
    return harmonic_matrices(matrices)


def tail_metrics(
    residual: np.ndarray, reference: np.ndarray, bulk_gap: float
) -> dict[str, float | None]:
    residual_rms = float(np.sqrt(np.mean(np.abs(residual) ** 2))) / bulk_gap
    reference_rms = float(np.sqrt(np.mean(np.abs(reference) ** 2))) / bulk_gap
    relative = residual_rms / reference_rms if reference_rms > FLOOR else None
    return {
        "tail_rms_over_bulk": reference_rms,
        "map_minus_tail_rms_over_bulk": residual_rms,
        "relative_rms_map_minus_tail": relative,
    }


def mismatch_label(metrics: dict[str, float | None]) -> str:
    relative = metrics["relative_rms_map_minus_tail"]
    if relative is None:
        return f"abs {metrics['map_minus_tail_rms_over_bulk']:.1e}"
    return f"rel {relative:.2g}"


def ray_selection(
    data: ProbeData | ShadowProbeData, ray_id: int
) -> tuple[np.ndarray, np.ndarray]:
    selected = np.flatnonzero(data.ray_id == ray_id)
    if selected.size == 0:
        raise PlotError(f"no samples were written for ray {ray_id}")
    radius = np.linalg.norm(data.coordinate[selected], axis=1)
    order = np.argsort(radius)
    return selected[order], radius[order]


def positive_limits(values: Sequence[np.ndarray]) -> tuple[float, float]:
    finite = np.concatenate(
        [np.asarray(value)[np.isfinite(value) & (np.asarray(value) > 0.0)] for value in values]
    )
    if finite.size == 0:
        return 1.0e-14, 1.0
    lower = max(float(np.min(finite)) * 0.5, 1.0e-14)
    upper = max(float(np.max(finite)) * 2.0, 10.0 * lower)
    return lower, upper


def plot_harmonic_axis(
    data: ProbeData,
    input_harmonic: np.ndarray,
    mapped_harmonic: np.ndarray,
    bulk_harmonic: np.ndarray,
    bulk_gap: float,
    axis_name: str,
    active_radius: float,
) -> plt.Figure:
    ray_number = 1 if axis_name == "x" else 2
    sides = (
        (ray_number, "+", POSITIVE_COLOUR, "o", "-"),
        (-ray_number, "-", NEGATIVE_COLOUR, "s", "--"),
    )
    departures = input_harmonic - bulk_harmonic
    residual = mapped_harmonic - input_harmonic
    all_values: list[np.ndarray] = []
    for ray_id, _, _, _, _ in sides:
        selected, _ = ray_selection(data, ray_id)
        all_values.extend(
            [np.abs(departures[selected]) / bulk_gap, np.abs(residual[selected]) / bulk_gap]
        )
    lower, upper = positive_limits(all_values)

    figure, axes = plt.subplots(
        3, 3, figsize=(12.2, 9.8), sharex=True, sharey=True,
        constrained_layout=True,
    )
    for name, spin, orbital in HARMONIC_COMPONENTS:
        axis = axes[spin, orbital]
        mismatch: list[str] = []
        for ray_id, sign, colour, marker, linestyle in sides:
            selected, radius = ray_selection(data, ray_id)
            expected = np.abs(departures[selected, spin, orbital]) / bulk_gap
            difference = np.abs(residual[selected, spin, orbital]) / bulk_gap
            axis.plot(
                radius,
                np.maximum(expected, FLOOR),
                color=colour,
                linewidth=1.0,
                linestyle=":",
                alpha=0.65,
                label=f"{sign}{axis_name} imposed tail",
            )
            axis.plot(
                radius,
                np.maximum(difference, FLOOR),
                color=colour,
                linewidth=1.35,
                linestyle=linestyle,
                marker=marker,
                markersize=3.2,
                markerfacecolor="white",
                label=f"{sign}{axis_name} map - tail",
            )
            metrics = tail_metrics(
                residual[selected, spin, orbital],
                departures[selected, spin, orbital],
                bulk_gap,
            )
            mismatch.append(
                f"{sign}: {mismatch_label(metrics)}"
            )
        axis.axvline(active_radius, color="0.35", linewidth=0.8)
        axis.set_yscale("log")
        axis.set_ylim(lower, upper)
        axis.grid(True, which="major", color="0.90", linewidth=0.55)
        axis.tick_params(direction="in", top=True, right=True)
        axis.set_title(rf"$A_{{{name}}}$; RMS mismatch " + ", ".join(mismatch))
        if spin == 2:
            axis.set_xlabel(r"$r/\xi_0$")
        if orbital == 0:
            axis.set_ylabel(r"magnitude$/\Delta_{\rm bulk}$")
    axes[0, 0].legend(loc="best", fontsize=6.7, framealpha=0.92)
    figure.suptitle(
        f"Independent quasiclassical check outside the active {axis_name}-radius\n"
        "dotted: imposed Eq. (19)-(20) tail; marked: freshly mapped minus imposed",
        fontsize=12,
    )
    return figure


def plot_mean_field(
    data: ProbeData, bulk_gap: float, active_radius_x: float, active_radius_y: float
) -> plt.Figure:
    residual = data.mapped_mean_field - data.input_mean_field
    figure, axes = plt.subplots(
        3, 2, figsize=(10.4, 9.0), sharey=True, constrained_layout=True
    )
    values = [np.abs(data.input_mean_field) / bulk_gap, np.abs(residual) / bulk_gap]
    lower, upper = positive_limits(values)
    for column, (axis_name, ray_number, active_radius) in enumerate(
        (("x", 1, active_radius_x), ("y", 2, active_radius_y))
    ):
        for component in range(3):
            axis = axes[component, column]
            for ray_id, sign, colour, marker, linestyle in (
                (ray_number, "+", POSITIVE_COLOUR, "o", "-"),
                (-ray_number, "-", NEGATIVE_COLOUR, "s", "--"),
            ):
                selected, radius = ray_selection(data, ray_id)
                expected = np.abs(data.input_mean_field[selected, component]) / bulk_gap
                difference = np.abs(residual[selected, component]) / bulk_gap
                axis.plot(radius, np.maximum(expected, FLOOR), color=colour,
                          linestyle=":", linewidth=1.0, alpha=0.65,
                          label=f"{sign}{axis_name} imposed")
                axis.plot(radius, np.maximum(difference, FLOOR), color=colour,
                          linestyle=linestyle, linewidth=1.35, marker=marker,
                          markersize=3.2, markerfacecolor="white",
                          label=f"{sign}{axis_name} map - imposed")
            axis.axvline(active_radius, color="0.35", linewidth=0.8)
            axis.set_yscale("log")
            axis.set_ylim(lower, upper)
            axis.grid(True, which="major", color="0.90", linewidth=0.55)
            axis.tick_params(direction="in", top=True, right=True)
            axis.set_title(rf"${axis_name}$ ray, $\nu_{{{component + 1}}}$")
            if component == 2:
                axis.set_xlabel(r"$r/\xi_0$")
            if column == 0:
                axis.set_ylabel(r"magnitude$/\Delta_{\rm bulk}$")
    axes[0, 0].legend(loc="best", fontsize=7.0, framealpha=0.92)
    figure.suptitle(
        "Independent Fermi-liquid mean-field tail check\n"
        "dotted: imposed $1/r+1/r^2$ tail; marked: freshly mapped minus imposed",
        fontsize=12,
    )
    return figure


def plot_shadow_harmonic_axis(
    data: ShadowProbeData,
    baseline_harmonic: np.ndarray,
    first_harmonic: np.ndarray,
    shadow_harmonic: np.ndarray,
    final_harmonic: np.ndarray,
    bulk_harmonic: np.ndarray,
    bulk_gap: float,
    axis_name: str,
    active_radius: float,
) -> plt.Figure:
    ray_number = 1 if axis_name == "x" else 2
    sides = (
        (ray_number, "+", POSITIVE_COLOUR, "o"),
        (-ray_number, "-", NEGATIVE_COLOUR, "s"),
    )
    departure = baseline_harmonic - bulk_harmonic
    first_defect = first_harmonic - baseline_harmonic
    correction = shadow_harmonic - baseline_harmonic
    terminal_defect = final_harmonic - shadow_harmonic
    all_values: list[np.ndarray] = []
    for ray_id, _, _, _ in sides:
        selected, _ = ray_selection(data, ray_id)
        all_values.extend(
            [
                np.abs(departure[selected]) / bulk_gap,
                np.abs(first_defect[selected]) / bulk_gap,
                np.abs(correction[selected]) / bulk_gap,
                np.abs(terminal_defect[selected]) / bulk_gap,
            ]
        )
    lower, upper = positive_limits(all_values)

    figure, axes = plt.subplots(
        3, 3, figsize=(12.4, 10.0), sharex=True, sharey=True,
        constrained_layout=True,
    )
    for name, spin, orbital in HARMONIC_COMPONENTS:
        axis = axes[spin, orbital]
        for ray_id, sign, colour, marker in sides:
            selected, radius = ray_selection(data, ray_id)
            expected = np.abs(departure[selected, spin, orbital]) / bulk_gap
            initial = np.abs(first_defect[selected, spin, orbital]) / bulk_gap
            change = np.abs(correction[selected, spin, orbital]) / bulk_gap
            terminal = np.abs(terminal_defect[selected, spin, orbital]) / bulk_gap
            axis.plot(
                radius, np.maximum(expected, FLOOR), color=colour,
                linestyle=":", linewidth=1.0, alpha=0.55,
                label=f"{sign}{axis_name} imposed tail",
            )
            axis.plot(
                radius, np.maximum(initial, FLOOR), color=colour,
                linestyle="-.", linewidth=0.9, marker=marker,
                markersize=3.0, markerfacecolor="white", alpha=0.75,
                label=f"{sign}{axis_name} first defect",
            )
            axis.plot(
                radius, np.maximum(change, FLOOR), color=colour,
                linestyle="-", linewidth=1.45, marker=marker,
                markersize=3.4, markerfacecolor=colour,
                label=f"{sign}{axis_name} shadow - tail",
            )
            axis.plot(
                radius, np.maximum(terminal, FLOOR), color=colour,
                linestyle="--", linewidth=1.1, marker="x", markersize=3.8,
                label=f"{sign}{axis_name} terminal defect",
            )
        axis.axvline(active_radius, color="0.35", linewidth=0.8)
        axis.set_yscale("log")
        axis.set_ylim(lower, upper)
        axis.grid(True, which="major", color="0.90", linewidth=0.55)
        axis.tick_params(direction="in", top=True, right=True)
        axis.set_title(rf"$A_{{{name}}}$", fontsize=10.5)
        if spin == 2:
            axis.set_xlabel(r"$r/\xi_0$")
        if orbital == 0:
            axis.set_ylabel(r"magnitude$/\Delta_{\rm bulk}$")
    axes[0, 0].legend(loc="best", fontsize=5.9, framealpha=0.92)
    figure.suptitle(
        f"Held-out shadow relaxation outside the active {axis_name}-radius\n"
        "correction tests extrapolation; terminal defect tests shadow convergence",
        fontsize=12,
    )
    return figure


def plot_shadow_mean_field(
    data: ShadowProbeData,
    bulk_gap: float,
    active_radius_x: float,
    active_radius_y: float,
) -> plt.Figure:
    first_defect = data.first_mapped_mean_field - data.baseline_mean_field
    correction = data.shadow_mean_field - data.baseline_mean_field
    terminal_defect = data.final_mapped_mean_field - data.shadow_mean_field
    values = [
        np.abs(data.baseline_mean_field) / bulk_gap,
        np.abs(first_defect) / bulk_gap,
        np.abs(correction) / bulk_gap,
        np.abs(terminal_defect) / bulk_gap,
    ]
    lower, upper = positive_limits(values)
    figure, axes = plt.subplots(
        3, 2, figsize=(10.6, 9.2), sharey=True, constrained_layout=True
    )
    for column, (axis_name, ray_number, active_radius) in enumerate(
        (("x", 1, active_radius_x), ("y", 2, active_radius_y))
    ):
        for component in range(3):
            axis = axes[component, column]
            for ray_id, sign, colour, marker in (
                (ray_number, "+", POSITIVE_COLOUR, "o"),
                (-ray_number, "-", NEGATIVE_COLOUR, "s"),
            ):
                selected, radius = ray_selection(data, ray_id)
                curves = (
                    (data.baseline_mean_field[selected, component], ":", None,
                     f"{sign}{axis_name} imposed tail"),
                    (first_defect[selected, component], "-.", marker,
                     f"{sign}{axis_name} first defect"),
                    (correction[selected, component], "-", marker,
                     f"{sign}{axis_name} shadow - tail"),
                    (terminal_defect[selected, component], "--", "x",
                     f"{sign}{axis_name} terminal defect"),
                )
                for curve, linestyle, curve_marker, label in curves:
                    axis.plot(
                        radius, np.maximum(np.abs(curve) / bulk_gap, FLOOR),
                        color=colour, linestyle=linestyle,
                        linewidth=1.35 if linestyle == "-" else 1.0,
                        marker=curve_marker, markersize=3.2,
                        markerfacecolor=(
                            colour if linestyle == "-" else "white"
                        ),
                        alpha=0.6 if linestyle == ":" else 1.0,
                        label=label,
                    )
            axis.axvline(active_radius, color="0.35", linewidth=0.8)
            axis.set_yscale("log")
            axis.set_ylim(lower, upper)
            axis.grid(True, which="major", color="0.90", linewidth=0.55)
            axis.tick_params(direction="in", top=True, right=True)
            axis.set_title(rf"${axis_name}$ ray, $\nu_{{{component + 1}}}$")
            if component == 2:
                axis.set_xlabel(r"$r/\xi_0$")
            if column == 0:
                axis.set_ylabel(r"magnitude$/\Delta_{\rm bulk}$")
    axes[0, 0].legend(loc="best", fontsize=6.2, framealpha=0.92)
    figure.suptitle(
        "Held-out Fermi-liquid mean-field shadow relaxation",
        fontsize=12,
    )
    return figure


def plot_shadow_history(
    path: Path, convergence_tolerance: float | None = None
) -> plt.Figure:
    rows = np.loadtxt(path, comments="#", ndmin=2)
    if rows.shape[1] != 9:
        raise PlotError(f"shadow history has {rows.shape[1]} columns; expected 9")
    figure, axis = plt.subplots(figsize=(7.8, 4.8), constrained_layout=True)
    axis.semilogy(rows[:, 0], np.maximum(rows[:, 2], FLOOR), color=POSITIVE_COLOUR,
                  marker="o", markerfacecolor="white", label="RMS defect")
    axis.semilogy(rows[:, 0], np.maximum(rows[:, 3], FLOOR), color=NEGATIVE_COLOUR,
                  marker="s", markerfacecolor="white", linestyle="--",
                  label="maximum-component defect")
    if convergence_tolerance is not None and convergence_tolerance > 0.0:
        axis.axhline(
            convergence_tolerance,
            color="0.25",
            linewidth=1.0,
            linestyle=":",
            label="convergence tolerance",
        )
    axis.set_xlabel("shadow map evaluation")
    axis.set_ylabel("self-consistency defect")
    axis.grid(True, which="both", color="0.90", linewidth=0.55)
    axis.tick_params(direction="in", top=True, right=True)
    axis.legend(loc="best")
    axis.set_title("Held-out exterior shadow-iteration convergence")
    return figure


def diagnostic_report(
    data: ProbeData,
    input_harmonic: np.ndarray,
    mapped_harmonic: np.ndarray,
    bulk_harmonic: np.ndarray,
    bulk_gap: float,
) -> dict[str, object]:
    departure = input_harmonic - bulk_harmonic
    residual = mapped_harmonic - input_harmonic
    rays: dict[str, object] = {}
    for ray_id, label in ((1, "+x"), (-1, "-x"), (2, "+y"), (-2, "-y")):
        selected, radius = ray_selection(data, ray_id)
        components: dict[str, object] = {}
        for name, spin, orbital in HARMONIC_COMPONENTS:
            components[name] = tail_metrics(
                residual[selected, spin, orbital],
                departure[selected, spin, orbital],
                bulk_gap,
            )
        mean_field_residual = data.mapped_mean_field[selected] - data.input_mean_field[selected]
        rays[label] = {
            "point_count": int(selected.size),
            "minimum_radius": float(radius[0]),
            "maximum_radius": float(radius[-1]),
            "harmonic_mismatch": components,
            "mean_field_mismatch": [
                tail_metrics(
                    mean_field_residual[:, component],
                    data.input_mean_field[selected, component],
                    bulk_gap,
                )
                for component in range(3)
            ],
        }
    return {
        "kind": "fermiforge_asymptotic_ray_map_check",
        "definition": (
            "The input is the imposed asymptotic state at passive points outside "
            "the active domain. The mapped state is a fresh full quasiclassical "
            "self-consistency evaluation at the same points."
        ),
        "rays": rays,
    }


def shadow_diagnostic_report(
    data: ShadowProbeData,
    baseline_harmonic: np.ndarray,
    first_harmonic: np.ndarray,
    shadow_harmonic: np.ndarray,
    final_harmonic: np.ndarray,
    bulk_harmonic: np.ndarray,
    bulk_gap: float,
    metrics: dict[str, str],
) -> dict[str, object]:
    departure = baseline_harmonic - bulk_harmonic
    first_defect = first_harmonic - baseline_harmonic
    correction = shadow_harmonic - baseline_harmonic
    terminal_defect = final_harmonic - shadow_harmonic
    rays: dict[str, object] = {}
    for ray_id, label in ((1, "+x"), (-1, "-x"), (2, "+y"), (-2, "-y")):
        selected, radius = ray_selection(data, ray_id)
        components: dict[str, object] = {}
        for name, spin, orbital in HARMONIC_COMPONENTS:
            reference = departure[selected, spin, orbital]
            components[name] = {
                "first_map_defect": tail_metrics(
                    first_defect[selected, spin, orbital], reference, bulk_gap
                ),
                "shadow_correction": tail_metrics(
                    correction[selected, spin, orbital], reference, bulk_gap
                ),
                "terminal_defect": tail_metrics(
                    terminal_defect[selected, spin, orbital], reference, bulk_gap
                ),
            }
        mean_components: list[object] = []
        first_mean = data.first_mapped_mean_field - data.baseline_mean_field
        correction_mean = data.shadow_mean_field - data.baseline_mean_field
        terminal_mean = data.final_mapped_mean_field - data.shadow_mean_field
        for component in range(3):
            reference = data.baseline_mean_field[selected, component]
            mean_components.append(
                {
                    "first_map_defect": tail_metrics(
                        first_mean[selected, component], reference, bulk_gap
                    ),
                    "shadow_correction": tail_metrics(
                        correction_mean[selected, component], reference, bulk_gap
                    ),
                    "terminal_defect": tail_metrics(
                        terminal_mean[selected, component], reference, bulk_gap
                    ),
                }
            )
        rays[label] = {
            "point_count": int(selected.size),
            "minimum_radius": float(radius[0]),
            "maximum_radius": float(radius[-1]),
            "harmonics": components,
            "mean_field": mean_components,
        }
    production_terminal_status = metrics.get("terminal_status", "unknown")
    shadow_terminal_status = metrics.get(
        "asymptotic_probe_terminal_status", "unknown"
    )

    def metric_value(key: str) -> float | None:
        try:
            return float(metrics[key])
        except (KeyError, ValueError):
            return None

    if shadow_terminal_status != "converged":
        validation_status = "shadow_not_converged"
    elif production_terminal_status != "converged":
        validation_status = "provisional_production_not_converged"
    else:
        validation_status = "passed"
    return {
        "kind": "fermiforge_asymptotic_shadow_iteration_check",
        "definition": {
            "first_map_defect": "F(B)-B",
            "shadow_correction": "S*-B; the correction to the imposed tail",
            "terminal_defect": "F(S*)-S*; convergence of the held-out solve",
        },
        "production_terminal_status": production_terminal_status,
        "shadow_terminal_status": shadow_terminal_status,
        "validation_status": validation_status,
        "validation_passed": validation_status == "passed",
        "global_shadow_stencil": {
            "terminal_rms": metric_value("asymptotic_probe_terminal_rms"),
            "terminal_maximum": metric_value(
                "asymptotic_probe_terminal_maximum"
            ),
            "terminal_relative_l2": metric_value(
                "asymptotic_probe_terminal_relative_l2"
            ),
            "convergence_tolerance": metric_value(
                "asymptotic_probe_convergence_tolerance"
            ),
            "iteration_points": metric_value(
                "asymptotic_probe_iteration_points"
            ),
        },
        "interpretation": (
            "The extrapolation is supported only when the terminal defect is "
            "small and the shadow correction is small relative to the imposed "
            "tail. A result is provisional unless the production solve itself "
            "has converged."
        ),
        "rays": rays,
    }


def save(figure: plt.Figure, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        probe_path = arguments.probe_file.expanduser().resolve()
        input_path = arguments.input.expanduser().resolve()
        metrics_path = arguments.metrics.expanduser().resolve()
        output_directory = arguments.output_dir.expanduser().resolve()
        parameters = load_parameters(input_path, metrics_path)
        metrics = read_metrics(metrics_path)
        data = load_probes(probe_path)
        input_harmonic = harmonic_matrices(data.input_order_parameter)
        mapped_harmonic = harmonic_matrices(data.mapped_order_parameter)
        bulk_harmonic = bulk_harmonics(
            data.coordinate,
            parameters.bulk_gap,
            parameters.bulk_phase,
            parameters.bulk_rotation,
        )
        outputs = [
            output_directory / "double_core_asymptotic_probe_harmonics_x.png",
            output_directory / "double_core_asymptotic_probe_harmonics_y.png",
            output_directory / "double_core_asymptotic_probe_mean_field.png",
        ]
        save(
            plot_harmonic_axis(
                data, input_harmonic, mapped_harmonic, bulk_harmonic,
                parameters.bulk_gap, "x", parameters.active_radius_x,
            ),
            outputs[0], arguments.dpi,
        )
        save(
            plot_harmonic_axis(
                data, input_harmonic, mapped_harmonic, bulk_harmonic,
                parameters.bulk_gap, "y", parameters.active_radius_y,
            ),
            outputs[1], arguments.dpi,
        )
        save(
            plot_mean_field(
                data, parameters.bulk_gap,
                parameters.active_radius_x, parameters.active_radius_y,
            ),
            outputs[2], arguments.dpi,
        )
        report_path = output_directory / "asymptotic_probe_diagnostics.json"
        report_path.write_text(
            json.dumps(
                diagnostic_report(
                    data, input_harmonic, mapped_harmonic, bulk_harmonic,
                    parameters.bulk_gap,
                ),
                indent=2,
            ) + "\n",
            encoding="utf-8",
        )
        report_paths = [report_path]

        if arguments.shadow_file is not None:
            shadow_path = arguments.shadow_file.expanduser().resolve()
            shadow_data = load_shadow_probes(shadow_path)
            baseline_harmonic = harmonic_matrices(
                shadow_data.baseline_order_parameter
            )
            first_harmonic = harmonic_matrices(
                shadow_data.first_mapped_order_parameter
            )
            shadow_harmonic = harmonic_matrices(
                shadow_data.shadow_order_parameter
            )
            final_harmonic = harmonic_matrices(
                shadow_data.final_mapped_order_parameter
            )
            shadow_bulk_harmonic = bulk_harmonics(
                shadow_data.coordinate,
                parameters.bulk_gap,
                parameters.bulk_phase,
                parameters.bulk_rotation,
            )
            shadow_outputs = [
                output_directory
                / "double_core_asymptotic_shadow_harmonics_x.png",
                output_directory
                / "double_core_asymptotic_shadow_harmonics_y.png",
                output_directory
                / "double_core_asymptotic_shadow_mean_field.png",
            ]
            save(
                plot_shadow_harmonic_axis(
                    shadow_data,
                    baseline_harmonic,
                    first_harmonic,
                    shadow_harmonic,
                    final_harmonic,
                    shadow_bulk_harmonic,
                    parameters.bulk_gap,
                    "x",
                    parameters.active_radius_x,
                ),
                shadow_outputs[0],
                arguments.dpi,
            )
            save(
                plot_shadow_harmonic_axis(
                    shadow_data,
                    baseline_harmonic,
                    first_harmonic,
                    shadow_harmonic,
                    final_harmonic,
                    shadow_bulk_harmonic,
                    parameters.bulk_gap,
                    "y",
                    parameters.active_radius_y,
                ),
                shadow_outputs[1],
                arguments.dpi,
            )
            save(
                plot_shadow_mean_field(
                    shadow_data,
                    parameters.bulk_gap,
                    parameters.active_radius_x,
                    parameters.active_radius_y,
                ),
                shadow_outputs[2],
                arguments.dpi,
            )
            outputs.extend(shadow_outputs)

            shadow_report_path = (
                output_directory / "asymptotic_shadow_diagnostics.json"
            )
            shadow_report_path.write_text(
                json.dumps(
                    shadow_diagnostic_report(
                        shadow_data,
                        baseline_harmonic,
                        first_harmonic,
                        shadow_harmonic,
                        final_harmonic,
                        shadow_bulk_harmonic,
                        parameters.bulk_gap,
                        metrics,
                    ),
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
            report_paths.append(shadow_report_path)

        if arguments.shadow_history is not None:
            shadow_history_path = arguments.shadow_history.expanduser().resolve()
            shadow_history_output = (
                output_directory
                / "double_core_asymptotic_shadow_convergence.png"
            )
            try:
                shadow_tolerance = float(
                    metrics["asymptotic_probe_convergence_tolerance"]
                )
            except (KeyError, ValueError):
                shadow_tolerance = None
            save(
                plot_shadow_history(shadow_history_path, shadow_tolerance),
                shadow_history_output,
                arguments.dpi,
            )
            outputs.append(shadow_history_output)
    except (OSError, ValueError, PlotError) as error:
        print(f"error: {error}", file=__import__("sys").stderr)
        return 2

    for path in outputs:
        print(f"wrote {path}")
    for path in report_paths:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
