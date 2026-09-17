#!/usr/bin/env python3
"""Plot a FermiForge structured 2D field map with accessible colour scales."""

from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[1]
MATPLOTLIB_CACHE = PROJECT_ROOT / "work" / "matplotlib"
MATPLOTLIB_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MATPLOTLIB_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


COMPONENTS = (
    ("xx", 0, 0), ("xy", 0, 1), ("xz", 0, 2),
    ("yx", 1, 0), ("yy", 1, 1), ("yz", 1, 2),
    ("zx", 2, 0), ("zy", 2, 1), ("zz", 2, 2),
)

WHITE_INDIGO = LinearSegmentedColormap.from_list(
    "fermiforge_white_indigo",
    ("#ffffff", "#dadaeb", "#9e9ac8", "#6a51a3", "#3f007d"),
)
BLUE_WHITE_RED = LinearSegmentedColormap.from_list(
    "fermiforge_blue_white_red",
    ("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"),
)
BLUE_WHITE_RED.set_bad("#bdbdbd")


class PlotError(RuntimeError):
    """Expected input or plotting failure."""


@dataclass(frozen=True)
class FieldMap:
    x: np.ndarray
    y: np.ndarray
    order_parameter: np.ndarray
    pair_density: np.ndarray
    current: np.ndarray
    metadata: Dict[str, str]


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot nine order-parameter amplitudes, cosines of their phases, "
            "pair density, and current components from a FermiForge 2D map."
        )
    )
    parser.add_argument("field_map", type=Path, help="FermiForge 2D ASCII map")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="plot directory (default: a plots directory beside the map)",
    )
    parser.add_argument(
        "--prefix",
        help="output filename prefix (default: input filename stem)",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "pdf"),
        default=("png", "pdf"),
        help="figure formats (default: png pdf)",
    )
    parser.add_argument("--dpi", type=int, default=220, help="PNG resolution")
    parser.add_argument(
        "--phase-mask-fraction",
        type=float,
        default=1.0e-6,
        help=(
            "mask phase where |A| is below this fraction of the common "
            "amplitude maximum (default: 1e-6)"
        ),
    )
    parser.add_argument(
        "--quiver-target",
        type=int,
        default=24,
        help="approximate number of current arrows along each axis",
    )
    parser.add_argument(
        "--relative-zero-threshold",
        type=float,
        default=1.0e-12,
        help=(
            "treat a current component as numerical zero when its maximum is "
            "below this fraction of the in-plane maximum (default: 1e-12)"
        ),
    )
    return parser.parse_args(argv)


def read_header(path: Path) -> tuple[Dict[str, str], list[str]]:
    metadata: Dict[str, str] = {}
    columns: list[str] = []
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped.startswith("#"):
                break
            payload = stripped[1:].strip()
            if payload.startswith("columns="):
                columns = payload.split("=", 1)[1].split()
            elif "=" in payload:
                key, value = payload.split("=", 1)
                metadata[key.strip()] = value.strip()
    if not columns:
        raise PlotError(f"{path} has no '# columns=' header")
    return metadata, columns


def load_field_map(path: Path) -> FieldMap:
    if not path.is_file():
        raise PlotError(f"field map does not exist: {path}")
    metadata, columns = read_header(path)
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != len(columns):
        raise PlotError(
            f"{path} contains {data.shape[1]} values per row but declares "
            f"{len(columns)} columns"
        )
    column = {name: data[:, index] for index, name in enumerate(columns)}
    required = {"x", "y", "pair_density", "j_x", "j_y", "j_z"}
    for name, _, _ in COMPONENTS:
        required.update((f"A_{name}_re", f"A_{name}_im"))
    missing = sorted(required.difference(column))
    if missing:
        raise PlotError(f"{path} is missing columns: {', '.join(missing)}")

    x_values = np.unique(column["x"])
    y_values = np.unique(column["y"])
    expected_rows = x_values.size * y_values.size
    if data.shape[0] != expected_rows:
        raise PlotError(
            f"map is not a complete rectangular grid: {data.shape[0]} rows "
            f"for {x_values.size} by {y_values.size} coordinates"
        )

    shape = (y_values.size, x_values.size)

    def to_grid(values: np.ndarray) -> np.ndarray:
        grid = np.full(shape, np.nan, dtype=values.dtype)
        x_index = np.searchsorted(x_values, column["x"])
        y_index = np.searchsorted(y_values, column["y"])
        if np.any(np.isfinite(grid[y_index, x_index])):
            raise PlotError("field map contains duplicate coordinate pairs")
        grid[y_index, x_index] = values
        if np.any(~np.isfinite(grid)):
            raise PlotError("field map contains missing or non-finite grid values")
        return grid

    order_parameter = np.empty((3, 3, *shape), dtype=np.complex128)
    for name, spin, orbital in COMPONENTS:
        order_parameter[spin, orbital] = to_grid(column[f"A_{name}_re"]) + 1j * to_grid(
            column[f"A_{name}_im"]
        )
    current = np.stack(
        (to_grid(column["j_x"]), to_grid(column["j_y"]), to_grid(column["j_z"])),
        axis=0,
    )
    return FieldMap(
        x=x_values,
        y=y_values,
        order_parameter=order_parameter,
        pair_density=to_grid(column["pair_density"]),
        current=current,
        metadata=metadata,
    )


def style_axes(axis: plt.Axes, x_label: bool, y_label: bool) -> None:
    axis.set_aspect("equal")
    if x_label:
        axis.set_xlabel("x")
    if y_label:
        axis.set_ylabel("y")
    axis.tick_params(direction="out", length=3, width=0.8)


def plot_amplitudes(field_map: FieldMap) -> plt.Figure:
    amplitudes = np.abs(field_map.order_parameter)
    maximum = float(np.nanmax(amplitudes))
    if not np.isfinite(maximum) or maximum <= 0.0:
        maximum = 1.0

    figure, axes = plt.subplots(3, 3, figsize=(10.6, 9.4), constrained_layout=True)
    image = None
    for name, spin, orbital in COMPONENTS:
        axis = axes[spin, orbital]
        image = axis.pcolormesh(
            field_map.x,
            field_map.y,
            amplitudes[spin, orbital],
            shading="auto",
            cmap=WHITE_INDIGO,
            vmin=0.0,
            vmax=maximum,
            rasterized=True,
        )
        axis.set_title(rf"$|A_{{{name}}}|$")
        style_axes(axis, spin == 2, orbital == 0)
    assert image is not None
    colorbar = figure.colorbar(image, ax=axes, shrink=0.92, pad=0.02)
    colorbar.set_label("amplitude (common scale)")
    figure.suptitle("Order-parameter amplitude")
    return figure


def plot_phase_cosines(field_map: FieldMap, mask_fraction: float) -> plt.Figure:
    amplitudes = np.abs(field_map.order_parameter)
    common_maximum = float(np.nanmax(amplitudes))
    threshold = max(0.0, mask_fraction) * common_maximum

    figure, axes = plt.subplots(3, 3, figsize=(10.6, 9.4), constrained_layout=True)
    image = None
    for name, spin, orbital in COMPONENTS:
        amplitude = amplitudes[spin, orbital]
        phase_cosine = np.full_like(amplitude, np.nan, dtype=float)
        np.divide(
            field_map.order_parameter[spin, orbital].real,
            amplitude,
            out=phase_cosine,
            where=amplitude > threshold,
        )
        axis = axes[spin, orbital]
        image = axis.pcolormesh(
            field_map.x,
            field_map.y,
            np.ma.masked_invalid(phase_cosine),
            shading="auto",
            cmap=BLUE_WHITE_RED,
            vmin=-1.0,
            vmax=1.0,
            rasterized=True,
        )
        axis.set_title(rf"$\cos[\arg(A_{{{name}}})]$")
        style_axes(axis, spin == 2, orbital == 0)
    assert image is not None
    colorbar = figure.colorbar(image, ax=axes, shrink=0.92, pad=0.02)
    colorbar.set_ticks((-1.0, -0.5, 0.0, 0.5, 1.0))
    colorbar.set_label(r"$\cos(\mathrm{phase})$")
    figure.suptitle("Order-parameter phase cosine (grey: phase undefined)")
    return figure


def plot_density_and_current(
    field_map: FieldMap, quiver_target: int, relative_zero_threshold: float
) -> plt.Figure:
    density = field_map.pair_density
    current_x, current_y, current_z = field_map.current
    in_plane_magnitude = np.hypot(current_x, current_y)
    current_kind = field_map.metadata.get("current_kind", "current")
    current_label = current_kind.replace("_", " ")

    figure, axes = plt.subplots(1, 3, figsize=(15.2, 4.8), constrained_layout=True)

    density_maximum = max(float(np.nanmax(density)), np.finfo(float).tiny)
    density_image = axes[0].pcolormesh(
        field_map.x,
        field_map.y,
        density,
        shading="auto",
        cmap=WHITE_INDIGO,
        vmin=0.0,
        vmax=density_maximum,
        rasterized=True,
    )
    axes[0].set_title(r"pair density $\sum_{\alpha i}|A_{\alpha i}|^2/3$")
    style_axes(axes[0], True, True)
    figure.colorbar(density_image, ax=axes[0], shrink=0.86, pad=0.02)

    current_maximum = max(float(np.nanmax(in_plane_magnitude)), np.finfo(float).tiny)
    current_image = axes[1].pcolormesh(
        field_map.x,
        field_map.y,
        in_plane_magnitude,
        shading="auto",
        cmap=WHITE_INDIGO,
        vmin=0.0,
        vmax=current_maximum,
        rasterized=True,
    )
    x_stride = max(1, field_map.x.size // max(1, quiver_target))
    y_stride = max(1, field_map.y.size // max(1, quiver_target))
    grid_x, grid_y = np.meshgrid(field_map.x, field_map.y)
    axes[1].quiver(
        grid_x[::y_stride, ::x_stride],
        grid_y[::y_stride, ::x_stride],
        current_x[::y_stride, ::x_stride],
        current_y[::y_stride, ::x_stride],
        color="white",
        edgecolor="black",
        linewidth=0.45,
        pivot="mid",
        angles="xy",
        scale_units="xy",
    )
    axes[1].set_title(rf"in-plane {current_label}: $|j_\perp|$ and direction")
    style_axes(axes[1], True, True)
    figure.colorbar(current_image, ax=axes[1], shrink=0.86, pad=0.02)

    axial_limit = float(np.nanmax(np.abs(current_z)))
    numerical_zero = (
        np.isfinite(axial_limit)
        and current_maximum > np.finfo(float).tiny
        and axial_limit <= max(0.0, relative_zero_threshold) * current_maximum
    )
    axial_title = rf"axial {current_label}: $j_z$"
    if numerical_zero:
        actual_axial_maximum = axial_limit
        axial_limit = max(
            max(0.0, relative_zero_threshold) * current_maximum,
            np.finfo(float).tiny,
        )
        axial_title += "\n" + rf"numerical zero; $|j_z|_{{max}}={actual_axial_maximum:.2e}$"
    elif not np.isfinite(axial_limit) or axial_limit <= 0.0:
        axial_limit = 1.0
    axial_image = axes[2].pcolormesh(
        field_map.x,
        field_map.y,
        current_z,
        shading="auto",
        cmap=BLUE_WHITE_RED,
        vmin=-axial_limit,
        vmax=axial_limit,
        rasterized=True,
    )
    axes[2].set_title(axial_title)
    style_axes(axes[2], True, True)
    figure.colorbar(axial_image, ax=axes[2], shrink=0.86, pad=0.02)

    figure.suptitle("Density and current map")
    return figure


def safe_prefix(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", value.strip()).strip("-.")
    return cleaned or "fermiforge-2d"


def save_figure(
    figure: plt.Figure,
    output_directory: Path,
    prefix: str,
    kind: str,
    formats: Iterable[str],
    dpi: int,
) -> list[Path]:
    paths: list[Path] = []
    for format_name in formats:
        path = output_directory / f"{prefix}_{kind}.{format_name}"
        save_options = {"bbox_inches": "tight"}
        if format_name == "png":
            save_options["dpi"] = dpi
        figure.savefig(path, **save_options)
        paths.append(path)
    plt.close(figure)
    return paths


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        field_map = load_field_map(arguments.field_map)
        output_directory = arguments.output_dir or arguments.field_map.parent / "plots"
        output_directory.mkdir(parents=True, exist_ok=True)
        prefix = safe_prefix(arguments.prefix or arguments.field_map.stem)

        paths: list[Path] = []
        paths.extend(
            save_figure(
                plot_amplitudes(field_map),
                output_directory,
                prefix,
                "order_parameter_amplitude",
                arguments.formats,
                arguments.dpi,
            )
        )
        paths.extend(
            save_figure(
                plot_phase_cosines(field_map, arguments.phase_mask_fraction),
                output_directory,
                prefix,
                "order_parameter_phase_cosine",
                arguments.formats,
                arguments.dpi,
            )
        )
        paths.extend(
            save_figure(
                plot_density_and_current(
                    field_map,
                    arguments.quiver_target,
                    arguments.relative_zero_threshold,
                ),
                output_directory,
                prefix,
                "density_and_current",
                arguments.formats,
                arguments.dpi,
            )
        )
    except (OSError, ValueError, PlotError) as error:
        print(f"error: {error}", file=__import__("sys").stderr)
        return 2

    print(f"loaded {arguments.field_map}")
    for path in paths:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
