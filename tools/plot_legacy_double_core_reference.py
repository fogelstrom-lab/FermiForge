#!/usr/bin/env python3
"""Create modern reference figures directly from a legacy 2D vortex archive.

The legacy state is split over ``op_x``, ``op_y``, ``op_z``, and ``curr``.
This tool reconstructs the same :class:`FieldMap` used by the modern plotting
tools, validates the redundant legacy columns, and writes a defensible subset
of the modern double-core figure set.  It deliberately does not invent an
initial state or iteration-dependent probe/shadow histories that are absent
from the archive.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MATPLOTLIB_CACHE = PROJECT_ROOT / "work" / "matplotlib"
MATPLOTLIB_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MATPLOTLIB_CACHE))


def ensure_plotting_python() -> None:
    """Re-execute with an installed NumPy/Matplotlib Python when necessary."""
    try:
        __import__("numpy")
        __import__("matplotlib")
        return
    except ModuleNotFoundError:
        pass
    if os.environ.get("FERMIFORGE_PLOTTING_REEXEC") == "1":
        raise SystemExit("no Python installation with NumPy and Matplotlib was found")
    candidates = (
        "/opt/homebrew/opt/python@3.14/bin/python3.14",
        shutil.which("python3"),
        "/opt/homebrew/bin/python3",
        "/usr/local/bin/python3",
        "/opt/local/bin/python3",
    )
    current = Path(sys.executable).resolve()
    for candidate_name in candidates:
        if not candidate_name:
            continue
        candidate = Path(candidate_name).expanduser().resolve()
        if candidate == current or not candidate.is_file():
            continue
        check = subprocess.run(
            [str(candidate), "-c", "import matplotlib, numpy"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            env=os.environ,
        )
        if check.returncode != 0:
            continue
        environment = os.environ.copy()
        environment["FERMIFORGE_PLOTTING_REEXEC"] = "1"
        os.execve(
            str(candidate),
            [str(candidate), str(Path(__file__).resolve()), *sys.argv[1:]],
            environment,
        )
    raise SystemExit("no Python installation with NumPy and Matplotlib was found")


ensure_plotting_python()

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_2d_fields import (
    HARMONIC_COMPONENTS,
    FieldMap,
    PlotError,
    cartesian_to_harmonics,
    plot_amplitudes,
    plot_basis_amplitudes,
    plot_basis_phase_cosines,
    plot_density_and_current,
    plot_phase_cosines,
    save_figure,
)
from plot_double_core_axis_profiles import (
    plot_profile_panel,
    positive_axis_profiles,
)


FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"


@dataclass(frozen=True)
class LegacyValidation:
    temperature_over_tc: float
    input_f1s: float
    bulk_gap: float
    maximum_gap_norm_difference: float
    maximum_in_plane_magnitude_difference: float
    maximum_legacy_density_difference: float
    half_core_y_minus: float
    half_core_y_plus: float


@dataclass(frozen=True)
class ConvergenceBlock:
    fixed_iteration: np.ndarray
    fixed_average: np.ndarray
    fixed_maximum: np.ndarray
    aa_iteration: np.ndarray
    aa_depth: np.ndarray
    aa_average: np.ndarray
    aa_maximum: np.ndarray
    aa_mixing: np.ndarray


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "archive",
        type=Path,
        help="legacy directory containing op_x, op_y, op_z, curr, and qcv.log",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="output directory (default: ARCHIVE/plots)",
    )
    parser.add_argument(
        "--core-half-width",
        type=float,
        default=30.0,
        help="half width of the additional core-window plots (default: 30 xi0)",
    )
    parser.add_argument(
        "--phase-mask-fraction",
        type=float,
        default=1.0e-6,
        help="mask phases below this fraction of the common amplitude maximum",
    )
    parser.add_argument(
        "--quiver-target",
        type=int,
        default=24,
        help="approximate number of mean-field arrows along each axis",
    )
    parser.add_argument("--dpi", type=int, default=220)
    return parser.parse_args(argv)


def require_numeric_table(path: Path, columns: int) -> np.ndarray:
    if not path.is_file():
        raise PlotError(f"legacy archive is missing {path.name}: {path}")
    table = np.loadtxt(path, ndmin=2)
    if table.shape[1] != columns:
        raise PlotError(
            f"{path} has {table.shape[1]} columns; expected {columns}"
        )
    if not np.all(np.isfinite(table)):
        raise PlotError(f"{path} contains non-finite values")
    return table


def infer_grid(coordinates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.unique(coordinates[:, 0])
    y = np.unique(coordinates[:, 1])
    if x.size * y.size != coordinates.shape[0]:
        raise PlotError("legacy coordinates do not form a complete rectangle")
    expected_x = np.tile(x, y.size)
    expected_y = np.repeat(y, x.size)
    scale = max(float(np.max(np.abs(coordinates))), 1.0)
    tolerance = 256.0 * np.finfo(float).eps * scale
    if not np.allclose(coordinates[:, 0], expected_x, rtol=0.0, atol=tolerance):
        raise PlotError("legacy records are not ordered with x varying fastest")
    if not np.allclose(coordinates[:, 1], expected_y, rtol=0.0, atol=tolerance):
        raise PlotError("legacy records do not contain complete y rows")
    return x, y


def require_matching_coordinates(
    reference: np.ndarray, actual: np.ndarray, path: Path
) -> None:
    if reference.shape != actual.shape:
        raise PlotError(f"coordinate count differs in {path}")
    scale = max(float(np.max(np.abs(reference))), 1.0)
    tolerance = 256.0 * np.finfo(float).eps * scale
    if not np.allclose(reference, actual, rtol=0.0, atol=tolerance):
        raise PlotError(f"coordinates differ in {path}")


def last_bulk_gap(log_path: Path) -> float:
    if not log_path.is_file():
        raise PlotError(f"legacy archive is missing {log_path.name}: {log_path}")
    text = log_path.read_text(encoding="utf-8", errors="replace")
    run_starts = [
        match.start() for match in re.finditer(r"(?m)^\s*got the input", text)
    ]
    final_run = text[run_starts[-1] :] if run_starts else text
    input_record = re.search(
        rf"(?m)^\s*({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
        rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s+"
        rf"({FLOAT_PATTERN})\s*$",
        final_run,
    )
    if input_record is not None:
        value = float(
            input_record.group(2).replace("d", "e").replace("D", "E")
        )
    else:
        matches = re.findall(
            rf"(?im)^\s*Delta\(T\)\s*:\s*({FLOAT_PATTERN})",
            final_run,
        )
        if not matches:
            raise PlotError(f"could not find Delta(T) in {log_path}")
        value = float(matches[-1].replace("d", "e").replace("D", "E"))
    if value <= 0.0:
        raise PlotError("legacy Delta(T) must be positive")
    return value


def legacy_input_parameters(input_path: Path) -> tuple[float, float]:
    if not input_path.is_file():
        raise PlotError(f"legacy archive is missing {input_path.name}: {input_path}")
    records: list[float] = []
    for line in input_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(rf"^\s*({FLOAT_PATTERN})", line)
        if match is None:
            continue
        records.append(float(match.group(1).replace("d", "e").replace("D", "E")))
        if len(records) == 2:
            break
    if len(records) < 2:
        raise PlotError(f"could not read temperature and F1s from {input_path}")
    return records[0], records[1]


def half_core_positions(field: FieldMap) -> tuple[float, float]:
    x_zero = int(np.argmin(np.abs(field.x)))
    amplitude = np.sqrt(field.pair_density[:, x_zero])
    negative = np.flatnonzero(field.y < 0.0)
    positive = np.flatnonzero(field.y > 0.0)
    if negative.size == 0 or positive.size == 0:
        raise PlotError("legacy mesh must straddle y=0")
    y_minus = float(field.y[negative[np.argmin(amplitude[negative])]])
    y_plus = float(field.y[positive[np.argmin(amplitude[positive])]])
    return y_minus, y_plus


def load_legacy_state(archive: Path) -> tuple[FieldMap, LegacyValidation]:
    reference_coordinates: np.ndarray | None = None
    x: np.ndarray | None = None
    y: np.ndarray | None = None
    order_parameter: np.ndarray | None = None

    for spin, leaf in enumerate(("op_x", "op_y", "op_z")):
        table = require_numeric_table(archive / leaf, 8)
        coordinates = table[:, :2]
        if reference_coordinates is None:
            reference_coordinates = coordinates.copy()
            x, y = infer_grid(reference_coordinates)
            order_parameter = np.empty(
                (3, 3, y.size, x.size), dtype=np.complex128
            )
        else:
            require_matching_coordinates(
                reference_coordinates, coordinates, archive / leaf
            )
        assert x is not None and y is not None and order_parameter is not None
        for orbital in range(3):
            column = 2 + 2 * orbital
            order_parameter[spin, orbital] = (
                table[:, column] + 1j * table[:, column + 1]
            ).reshape(y.size, x.size)
        del table

    assert reference_coordinates is not None
    assert x is not None and y is not None and order_parameter is not None

    current_table = require_numeric_table(archive / "curr", 7)
    require_matching_coordinates(
        reference_coordinates, current_table[:, :2], archive / "curr"
    )
    current = np.stack(
        [
            current_table[:, component].reshape(y.size, x.size)
            for component in (4, 5, 6)
        ],
        axis=0,
    )
    pair_density = np.sum(np.abs(order_parameter) ** 2, axis=(0, 1)) / 3.0
    gap_norm = np.sqrt(pair_density)
    maximum_gap_difference = float(
        np.max(np.abs(gap_norm - current_table[:, 2].reshape(y.size, x.size)))
    )
    maximum_in_plane_difference = float(
        np.max(
            np.abs(
                np.hypot(current[0], current[1])
                - current_table[:, 3].reshape(y.size, x.size)
            )
        )
    )

    field = FieldMap(
        x=x,
        y=y,
        order_parameter=order_parameter,
        pair_density=pair_density,
        current=current,
        metadata={
            "source": "legacy_split_2d_archive",
            "density_kind": "pair_density",
            "current_kind": "current_related_mean_field",
        },
    )

    temperature, input_f1s = legacy_input_parameters(archive / "qcv.inp")
    bulk_gap = last_bulk_gap(archive / "qcv.log")
    density_table = require_numeric_table(archive / "dens", 5)
    require_matching_coordinates(
        reference_coordinates, density_table[:, :2], archive / "dens"
    )
    normalized_pair_density = pair_density / bulk_gap**2
    maximum_density_difference = float(
        np.max(
            np.abs(
                normalized_pair_density
                - density_table[:, 2].reshape(y.size, x.size)
            )
        )
    )
    y_minus, y_plus = half_core_positions(field)
    validation = LegacyValidation(
        temperature_over_tc=temperature,
        input_f1s=input_f1s,
        bulk_gap=bulk_gap,
        maximum_gap_norm_difference=maximum_gap_difference,
        maximum_in_plane_magnitude_difference=maximum_in_plane_difference,
        maximum_legacy_density_difference=maximum_density_difference,
        half_core_y_minus=y_minus,
        half_core_y_plus=y_plus,
    )
    return field, validation


def crop_field(field: FieldMap, half_width: float) -> FieldMap:
    if half_width <= 0.0:
        raise PlotError("core-half-width must be positive")
    x_mask = np.abs(field.x) <= half_width + 64.0 * np.finfo(float).eps
    y_mask = np.abs(field.y) <= half_width + 64.0 * np.finfo(float).eps
    if np.count_nonzero(x_mask) < 2 or np.count_nonzero(y_mask) < 2:
        raise PlotError("core-half-width selects fewer than two grid points per axis")
    return FieldMap(
        x=field.x[x_mask],
        y=field.y[y_mask],
        order_parameter=field.order_parameter[:, :, y_mask, :][:, :, :, x_mask],
        pair_density=field.pair_density[y_mask, :][:, x_mask],
        current=field.current[:, y_mask, :][:, :, x_mask],
        metadata=dict(field.metadata),
    )


def add_provenance_title(figure: plt.Figure, suffix: str) -> None:
    title = figure._suptitle.get_text() if figure._suptitle is not None else ""
    figure.suptitle(f"{title}\n{suffix}" if title else suffix, fontsize=12)


def write_state_figures(
    field: FieldMap,
    output_directory: Path,
    prefix: str,
    title_suffix: str,
    phase_mask_fraction: float,
    quiver_target: int,
    dpi: int,
) -> list[Path]:
    harmonics = cartesian_to_harmonics(field.order_parameter)
    figures = (
        (plot_amplitudes(field), "order_parameter_amplitude"),
        (
            plot_phase_cosines(field, phase_mask_fraction),
            "order_parameter_phase_cosine",
        ),
        (
            plot_basis_amplitudes(
                field.x,
                field.y,
                harmonics,
                HARMONIC_COMPONENTS,
                "Order-parameter amplitude - spherical/harmonic basis",
            ),
            "order_parameter_harmonic_amplitude",
        ),
        (
            plot_basis_phase_cosines(
                field.x,
                field.y,
                harmonics,
                HARMONIC_COMPONENTS,
                phase_mask_fraction,
                "Order-parameter phase cosine - spherical/harmonic basis",
            ),
            "order_parameter_harmonic_phase_cosine",
        ),
        (
            plot_density_and_current(field, quiver_target, 1.0e-12),
            "density_and_current_related_mean_field",
        ),
    )
    paths: list[Path] = []
    for figure, kind in figures:
        add_provenance_title(figure, title_suffix)
        paths.extend(
            save_figure(
                figure,
                output_directory,
                prefix,
                kind,
                ("png",),
                dpi,
            )
        )
    return paths


def make_axis_profile_figure(field: FieldMap, bulk_gap: float) -> plt.Figure:
    x_profile, y_profile = positive_axis_profiles(field)
    figure, axes = plt.subplots(
        1, 2, figsize=(12.0, 4.8), squeeze=False, sharey=True,
        constrained_layout=True,
    )
    largest = max(
        1.0,
        plot_profile_panel(axes[0, 0], x_profile, bulk_gap, 2.0e-3, True),
        plot_profile_panel(axes[0, 1], y_profile, bulk_gap, 2.0e-3, False),
    )
    axes[0, 0].set_ylabel(r"$A_{mn}/\Delta_{\rm bulk}$")
    axes[0, 0].set_title("legacy reference: positive x axis")
    axes[0, 1].set_title("legacy reference: positive y axis")
    axes[0, 0].set_xlabel(r"$x/\xi_0$")
    axes[0, 1].set_xlabel(r"$y/\xi_0$")
    limit = 1.08 * largest
    axes[0, 0].set_ylim(-limit, limit)
    axes[0, 1].set_ylim(-limit, limit)
    figure.suptitle(
        "Legacy double-core symmetry-axis order-parameter profiles\n"
        rf"dcvlong Fig. 4 convention; $\Delta_{{\rm bulk}}={bulk_gap:.7g}$",
        fontsize=12,
    )
    return figure


def last_run_text(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    starts = [match.start() for match in re.finditer(r"(?m)^\s*got the input", text)]
    return text[starts[-1] :] if starts else text


def parse_convergence(log_path: Path) -> ConvergenceBlock:
    text = last_run_text(log_path)
    fixed_pattern = re.compile(
        rf"(?im)^\s*it\s+nr\s+(\d+)\s+errav\s*=\s*({FLOAT_PATTERN})"
        rf"\s+errmax\s*=\s*({FLOAT_PATTERN})"
    )
    aa_pattern = re.compile(
        rf"(?m)^\s*(\d+)\s+(\d+)\s+({FLOAT_PATTERN})\s+"
        rf"({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"
    )
    fixed_matches = list(fixed_pattern.finditer(text))
    if not fixed_matches:
        raise PlotError(f"no legacy fixed-iteration table found in {log_path}")
    fixed_rows = [match.groups() for match in fixed_matches]
    aa_rows = aa_pattern.findall(text[fixed_matches[-1].end() :])
    if not aa_rows:
        raise PlotError(f"no legacy Anderson table found in {log_path}")

    def real(value: str) -> float:
        return float(value.replace("d", "e").replace("D", "E"))

    return ConvergenceBlock(
        fixed_iteration=np.asarray([int(row[0]) for row in fixed_rows]),
        fixed_average=np.asarray([real(row[1]) for row in fixed_rows]),
        fixed_maximum=np.asarray([real(row[2]) for row in fixed_rows]),
        aa_iteration=np.asarray([int(row[0]) for row in aa_rows]),
        aa_depth=np.asarray([int(row[1]) for row in aa_rows]),
        aa_average=np.asarray([real(row[2]) for row in aa_rows]),
        aa_maximum=np.asarray([real(row[3]) for row in aa_rows]),
        aa_mixing=np.asarray([real(row[4]) for row in aa_rows]),
    )


def make_convergence_figure(history: ConvergenceBlock) -> plt.Figure:
    fixed_x = np.arange(1, history.fixed_iteration.size + 1, dtype=float)
    aa_x = np.arange(
        history.fixed_iteration.size + 1,
        history.fixed_iteration.size + history.aa_iteration.size + 1,
        dtype=float,
    )
    figure, axis = plt.subplots(figsize=(8.7, 5.4), constrained_layout=True)
    axis.axvspan(
        fixed_x[0] - 0.5,
        fixed_x[-1] + 0.5,
        color="#f0f0f0",
        label="legacy preliminary iterator",
        zorder=0,
    )
    axis.axvspan(
        aa_x[0] - 0.5,
        aa_x[-1] + 0.5,
        color="#efedf5",
        label="legacy Anderson engine",
        zorder=0,
    )
    axis.semilogy(
        fixed_x,
        np.maximum(history.fixed_average, np.finfo(float).tiny),
        "o--",
        color="#2166ac",
        label="reported average error",
    )
    axis.semilogy(
        fixed_x,
        np.maximum(history.fixed_maximum, np.finfo(float).tiny),
        "s--",
        color="#b2182b",
        label="reported maximum error",
    )
    axis.semilogy(
        aa_x,
        np.maximum(history.aa_average, np.finfo(float).tiny),
        "o-",
        color="#2166ac",
    )
    axis.semilogy(
        aa_x,
        np.maximum(history.aa_maximum, np.finfo(float).tiny),
        "s-",
        color="#b2182b",
    )
    axis.axvline(aa_x[0] - 0.5, color="0.40", linewidth=0.8)
    axis.set_xlabel("sequential update record in final qcv.log block")
    axis.set_ylabel("legacy reported error")
    axis.grid(True, which="both", color="0.86", linewidth=0.7)

    mixing_axis = axis.twinx()
    mixing_axis.plot(
        aa_x,
        history.aa_mixing,
        "^-",
        color="#1b7837",
        linewidth=1.2,
        label="Anderson mixing p",
    )
    mixing_axis.set_ylabel("Anderson mixing p", color="#1b7837")
    mixing_axis.tick_params(axis="y", colors="#1b7837")
    handles, labels = axis.get_legend_handles_labels()
    extra_handles, extra_labels = mixing_axis.get_legend_handles_labels()
    axis.legend(handles + extra_handles, labels + extra_labels, loc="best")
    axis.set_title(
        "Legacy double-core convergence\n"
        "preliminary and Anderson counters are separate in the source log"
    )
    return figure


def write_manifest(
    output_path: Path,
    archive: Path,
    field: FieldMap,
    validation: LegacyValidation,
    generated_paths: Sequence[Path],
    core_half_width: float,
) -> None:
    harmonics = cartesian_to_harmonics(field.order_parameter)
    in_plane_magnitude = np.hypot(field.current[0], field.current[1])
    payload = {
        "source_archive": str(archive.resolve()),
        "grid": {
            "nx": int(field.x.size),
            "ny": int(field.y.size),
            "x_min": float(field.x[0]),
            "x_max": float(field.x[-1]),
            "y_min": float(field.y[0]),
            "y_max": float(field.y[-1]),
            "core_plot_half_width": core_half_width,
        },
        "physics": {
            "temperature_over_tc": validation.temperature_over_tc,
            "input_F1s": validation.input_f1s,
            "bulk_gap": validation.bulk_gap,
            "half_core_y_minus": validation.half_core_y_minus,
            "half_core_y_plus": validation.half_core_y_plus,
            "hard_core_separation": (
                validation.half_core_y_plus - validation.half_core_y_minus
            ),
        },
        "validation": {
            "maximum_stored_gap_norm_difference": (
                validation.maximum_gap_norm_difference
            ),
            "maximum_stored_in_plane_magnitude_difference": (
                validation.maximum_in_plane_magnitude_difference
            ),
            "maximum_legacy_normalized_density_difference": (
                validation.maximum_legacy_density_difference
            ),
        },
        "comparison_extrema": {
            "cartesian_amplitude_maximum": float(
                np.max(np.abs(field.order_parameter))
            ),
            "harmonic_amplitude_maximum": float(np.max(np.abs(harmonics))),
            "pair_density_minimum": float(np.min(field.pair_density)),
            "pair_density_maximum": float(np.max(field.pair_density)),
            "in_plane_mean_field_magnitude_maximum": float(
                np.max(in_plane_magnitude)
            ),
            "axial_mean_field_absolute_maximum": float(
                np.max(np.abs(field.current[2]))
            ),
        },
        "conventions": {
            "harmonics": (
                "recomputed from op_x/op_y/op_z in the FermiForge orthonormal "
                "spherical basis; auxiliary Cpp...Cmm files use a factor-two "
                "mixed-component normalization and a different mixed-zero "
                "filename ordering"
            ),
            "density": "sum(abs(A_alpha_i)**2)/3",
            "vector_field": (
                "curr columns 5:7, labelled current-related Fermi-liquid mean "
                "field rather than physical mass current"
            ),
            "convergence": (
                "final qcv.log block only; preliminary and Anderson counters "
                "are shown as sequential but separately shaded records"
            ),
        },
        "unavailable_from_archive": [
            "initial-state field plots",
            "iteration-dependent double-core structure history",
            "modern asymptotic probe and shadow histories",
        ],
        "png_files": [path.name for path in generated_paths],
    }
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    archive = arguments.archive.resolve()
    output_directory = (arguments.output_dir or archive / "plots").resolve()
    if arguments.phase_mask_fraction < 0.0:
        raise SystemExit("phase-mask-fraction cannot be negative")
    if arguments.core_half_width <= 0.0:
        raise SystemExit("core-half-width must be positive")

    try:
        field, validation = load_legacy_state(archive)
        core_field = crop_field(field, arguments.core_half_width)
        output_directory.mkdir(parents=True, exist_ok=True)
        generated: list[Path] = []
        generated.extend(
            write_state_figures(
                field,
                output_directory,
                "legacy_reference_full",
                (
                    f"legacy T/Tc={validation.temperature_over_tc:.2f}, "
                    f"F1s={validation.input_f1s:g}; full "
                    f"[{field.x[0]:g},{field.x[-1]:g}] x "
                    f"[{field.y[0]:g},{field.y[-1]:g}] xi0 archive"
                ),
                arguments.phase_mask_fraction,
                arguments.quiver_target,
                arguments.dpi,
            )
        )
        generated.extend(
            write_state_figures(
                core_field,
                output_directory,
                "legacy_reference_core",
                (
                    f"legacy T/Tc={validation.temperature_over_tc:.2f}, "
                    f"F1s={validation.input_f1s:g}; "
                    f"core window +/-{arguments.core_half_width:g} xi0"
                ),
                arguments.phase_mask_fraction,
                arguments.quiver_target,
                arguments.dpi,
            )
        )
        axis_figure = make_axis_profile_figure(field, validation.bulk_gap)
        generated.extend(
            save_figure(
                axis_figure,
                output_directory,
                "legacy_reference",
                "double_core_axis_profiles",
                ("png",),
                arguments.dpi,
            )
        )
        history = parse_convergence(archive / "qcv.log")
        convergence_figure = make_convergence_figure(history)
        generated.extend(
            save_figure(
                convergence_figure,
                output_directory,
                "legacy_reference",
                "convergence",
                ("png",),
                arguments.dpi,
            )
        )
        manifest_path = output_directory / "legacy_reference_manifest.json"
        write_manifest(
            manifest_path,
            archive,
            field,
            validation,
            generated,
            arguments.core_half_width,
        )
    except (OSError, ValueError, PlotError) as error:
        print(f"error: {error}", file=__import__("sys").stderr)
        return 2

    print(f"loaded legacy archive: {archive}")
    print(f"grid: {field.x.size} x {field.y.size}")
    print(
        "hard-core minima y-/y+: "
        f"{validation.half_core_y_minus:.3f} "
        f"{validation.half_core_y_plus:.3f} xi0"
    )
    print(
        "validation maxima [gap norm, in-plane field, normalized density]: "
        f"{validation.maximum_gap_norm_difference:.3e} "
        f"{validation.maximum_in_plane_magnitude_difference:.3e} "
        f"{validation.maximum_legacy_density_difference:.3e}"
    )
    for path in generated:
        print(f"wrote {path}")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
