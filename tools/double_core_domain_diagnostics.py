#!/usr/bin/env python3
"""Measure whether a double-core active ellipse contains the relaxed texture."""

from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path
from typing import Sequence


class DomainDiagnosticError(RuntimeError):
    """The field or residual data cannot support a domain diagnostic."""


def read_numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            rows.append([float(value) for value in stripped.split()])
        except ValueError as error:
            raise DomainDiagnosticError(f"nonnumeric row in {path}") from error
    if not rows:
        raise DomainDiagnosticError(f"no numeric records in {path}")
    return rows


def relative_step(left: Sequence[float], right: Sequence[float]) -> float:
    numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))
    denominator = max(math.sqrt(sum(value * value for value in left)), 1.0e-30)
    return numerator / denominator


def one_sided_axis_diagnostic(
    points: list[tuple[float, list[float]]], radius: float, interior_steps: int
) -> dict[str, float]:
    points.sort(key=lambda item: item[0])
    inside_index = max(
        (index for index, item in enumerate(points) if item[0] <= radius),
        default=-1,
    )
    if inside_index < interior_steps or inside_index + 1 >= len(points):
        raise DomainDiagnosticError(
            f"active radius {radius:g} is not bracketed by enough axis points"
        )
    inside_radius, inside_field = points[inside_index]
    outside_radius, outside_field = points[inside_index + 1]
    interface_step = relative_step(inside_field, outside_field)
    interior = [
        relative_step(points[index - 1][1], points[index][1])
        for index in range(inside_index - interior_steps + 1, inside_index + 1)
    ]
    interior_reference = max(statistics.median(interior), 1.0e-14)
    return {
        "inside_radius": inside_radius,
        "outside_radius": outside_radius,
        "relative_interface_step": interface_step,
        "interior_relative_step": interior_reference,
        "interface_amplification": interface_step / interior_reference,
    }


def axis_diagnostic(
    field_rows: list[list[float]], axis: str, radius: float, interior_steps: int
) -> dict[str, object]:
    coordinate_index = 0 if axis == "x" else 1
    transverse_index = 1 - coordinate_index
    scale = max(max(abs(row[coordinate_index]) for row in field_rows), 1.0)
    tolerance = 128.0 * math.ulp(scale)
    sides: dict[str, dict[str, float]] = {}
    for sign, label in ((-1.0, "negative"), (1.0, "positive")):
        points = [
            (abs(row[coordinate_index]), row[2:])
            for row in field_rows
            if abs(row[transverse_index]) <= tolerance
            and sign * row[coordinate_index] >= -tolerance
        ]
        sides[label] = one_sided_axis_diagnostic(points, radius, interior_steps)
    return {
        "negative": sides["negative"],
        "positive": sides["positive"],
        "relative_interface_step": max(
            float(sides[label]["relative_interface_step"]) for label in sides
        ),
        "interface_amplification": max(
            float(sides[label]["interface_amplification"]) for label in sides
        ),
    }


def residual_diagnostics(
    residual_rows: list[list[float]],
    radius_x: float,
    radius_y: float,
    collar_fraction: float,
    core_fraction: float,
) -> dict[str, object]:
    core_values: list[float] = []
    directional: dict[str, list[tuple[float, float]]] = {"x": [], "y": []}
    for row in residual_rows:
        if len(row) < 10:
            raise DomainDiagnosticError("sampled residual row has too few columns")
        x, y = row[0], row[1]
        normalized_x, normalized_y = x / radius_x, y / radius_y
        rho = math.hypot(normalized_x, normalized_y)
        point_rms, maximum = row[-2], row[-1]
        if rho <= core_fraction:
            core_values.append(point_rms)
        if 1.0 - collar_fraction <= rho <= 1.0 + 1.0e-12:
            axis = "x" if abs(normalized_x) >= abs(normalized_y) else "y"
            directional[axis].append((point_rms, maximum))
    if not core_values or any(not values for values in directional.values()):
        raise DomainDiagnosticError("residual data do not cover core and boundary sectors")
    core_rms = math.sqrt(sum(value * value for value in core_values) / len(core_values))
    result: dict[str, object] = {"core_rms": core_rms}
    for axis, values in directional.items():
        boundary_rms = math.sqrt(
            sum(value[0] * value[0] for value in values) / len(values)
        )
        result[axis] = {
            "point_count": len(values),
            "boundary_rms": boundary_rms,
            "boundary_maximum": max(value[1] for value in values),
            "boundary_to_core_rms": boundary_rms / max(core_rms, 1.0e-30),
        }
    return result


def diagnose_domain(
    field_path: Path,
    residual_path: Path,
    radius_x: float,
    radius_y: float,
    *,
    collar_fraction: float = 0.15,
    core_fraction: float = 0.50,
    interface_amplification_limit: float = 3.0,
    interface_step_floor: float = 1.0e-3,
    boundary_to_core_limit: float = 0.25,
    interior_steps: int = 3,
) -> dict[str, object]:
    if radius_x <= 0.0 or radius_y <= 0.0:
        raise DomainDiagnosticError("active ellipse radii must be positive")
    if not 0.0 < collar_fraction < 1.0 or not 0.0 < core_fraction < 1.0:
        raise DomainDiagnosticError("collar and core fractions must lie in (0,1)")
    field_rows = read_numeric_rows(field_path)
    residual_rows = read_numeric_rows(residual_path)
    residual = residual_diagnostics(
        residual_rows, radius_x, radius_y, collar_fraction, core_fraction
    )
    axes: dict[str, dict[str, object]] = {}
    for axis, radius in (("x", radius_x), ("y", radius_y)):
        interface = axis_diagnostic(field_rows, axis, radius, interior_steps)
        collar = dict(residual[axis])
        reasons: list[str] = []
        if (
            float(interface["relative_interface_step"]) > interface_step_floor
            and float(interface["interface_amplification"])
            > interface_amplification_limit
        ):
            reasons.append("interface_step")
        if float(collar["boundary_to_core_rms"]) > boundary_to_core_limit:
            reasons.append("boundary_residual")
        axes[axis] = {
            "radius": radius,
            "interface": interface,
            "residual_collar": collar,
            "needs_growth": bool(reasons),
            "reasons": reasons,
        }
    return {
        "field_file": str(field_path),
        "residual_file": str(residual_path),
        "radius_x": radius_x,
        "radius_y": radius_y,
        "criteria": {
            "collar_fraction": collar_fraction,
            "core_fraction": core_fraction,
            "interface_amplification_limit": interface_amplification_limit,
            "interface_step_floor": interface_step_floor,
            "boundary_to_core_limit": boundary_to_core_limit,
            "interior_steps": interior_steps,
        },
        "core_residual_rms": residual["core_rms"],
        "axes": axes,
        "adequate": not any(bool(value["needs_growth"]) for value in axes.values()),
    }


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("field", type=Path)
    parser.add_argument("residuals", type=Path)
    parser.add_argument("--radius-x", type=float, required=True)
    parser.add_argument("--radius-y", type=float, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--collar-fraction", type=float, default=0.15)
    parser.add_argument("--core-fraction", type=float, default=0.50)
    parser.add_argument("--interface-amplification-limit", type=float, default=3.0)
    parser.add_argument("--interface-step-floor", type=float, default=1.0e-3)
    parser.add_argument("--boundary-to-core-limit", type=float, default=0.25)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        report = diagnose_domain(
            arguments.field.expanduser().resolve(),
            arguments.residuals.expanduser().resolve(),
            arguments.radius_x,
            arguments.radius_y,
            collar_fraction=arguments.collar_fraction,
            core_fraction=arguments.core_fraction,
            interface_amplification_limit=arguments.interface_amplification_limit,
            interface_step_floor=arguments.interface_step_floor,
            boundary_to_core_limit=arguments.boundary_to_core_limit,
        )
        serialized = json.dumps(report, indent=2) + "\n"
        if arguments.output is None:
            print(serialized, end="")
        else:
            arguments.output.expanduser().resolve().write_text(
                serialized, encoding="utf-8"
            )
        return 0 if bool(report["adequate"]) else 1
    except (OSError, ValueError, DomainDiagnosticError) as error:
        print(f"error: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
