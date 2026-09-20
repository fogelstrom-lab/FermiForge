#!/usr/bin/env python3
"""Adapt the double-core active ellipse from solution-based diagnostics."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import re
import sys
from pathlib import Path
from typing import Sequence

from double_core_domain_diagnostics import DomainDiagnosticError, diagnose_domain
from run_double_core_from_scratch import build_solver
from run_double_core_radius_ladder import (
    LadderError,
    read_metrics,
    replace_namelist_value,
    stream_command,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "examples" / "2d_double_core_overnight.nml"
DEFAULT_BUILD = PROJECT_ROOT / "work" / "double-core-scratch-build"


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--restart",
        type=Path,
        help="source checkpoint; defaults to the newest overnight or adaptive checkpoint",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD)
    parser.add_argument("--run-root", type=Path, default=PROJECT_ROOT / "runs")
    parser.add_argument("--block-iterations", type=int, default=8)
    parser.add_argument("--maximum-stages", type=int, default=4)
    parser.add_argument("--growth-factor", type=float, default=1.25)
    parser.add_argument("--guard-cells", type=float, default=2.0)
    parser.add_argument("--initial-radius-x", type=float)
    parser.add_argument("--initial-radius-y", type=float)
    parser.add_argument("--maximum-radius-x", type=float)
    parser.add_argument("--maximum-radius-y", type=float)
    parser.add_argument("--collar-fraction", type=float, default=0.15)
    parser.add_argument("--core-fraction", type=float, default=0.50)
    parser.add_argument("--interface-amplification-limit", type=float, default=3.0)
    parser.add_argument("--interface-step-floor", type=float, default=1.0e-3)
    parser.add_argument("--boundary-to-core-limit", type=float, default=0.25)
    parser.add_argument("--ranks", type=int, default=10)
    parser.add_argument("--jobs", type=int, default=10)
    parser.add_argument("--diagnostics-only", action="store_true")
    parser.add_argument(
        "--hold-domain",
        action="store_true",
        help="relax one block at the current radii without domain growth",
    )
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args(argv)


def namelist_float(text: str, key: str) -> float:
    match = re.search(
        rf"(?im)^\s*{re.escape(key)}\s*=\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)",
        text,
    )
    if match is None:
        raise LadderError(f"input namelist does not define numeric '{key}'")
    return float(match.group(1).replace("d", "e").replace("D", "E"))


def latest_source_checkpoint(run_root: Path) -> Path:
    candidates = []
    for pattern in (
        "*-double-core-overnight/sparse_checkpoint_state.dat",
        "*-double-core-adaptive-domain/stages/*/sparse_checkpoint_state.dat",
    ):
        candidates.extend(path for path in run_root.glob(pattern) if path.is_file())
    if not candidates:
        raise LadderError("no overnight or adaptive checkpoint found; pass --restart")

    def production_score(path: Path) -> tuple[int, float]:
        metrics_path = path.parent / "metrics.txt"
        point_count = 0
        if metrics_path.is_file():
            metrics = read_metrics(metrics_path)
            try:
                point_count = int(metrics["mesh_x_points"]) * int(
                    metrics["mesh_y_points"]
                )
            except (KeyError, ValueError):
                point_count = 0
        return point_count, path.stat().st_mtime

    # This is the production adaptive-domain driver.  Prefer the largest
    # available field, then its newest checkpoint, so a later compact smoke
    # test cannot silently replace the 119x119 production source.
    return max(candidates, key=production_score).resolve()


def infer_active_radii(
    checkpoint: Path, radius_x: float | None, radius_y: float | None
) -> tuple[float, float]:
    metrics_path = checkpoint.parent / "metrics.txt"
    if not metrics_path.is_file() and (radius_x is None or radius_y is None):
        raise LadderError("cannot infer active radii without sibling metrics.txt")
    metrics = read_metrics(metrics_path) if metrics_path.is_file() else {}
    region = metrics.get("update_region", "")
    if radius_x is None:
        key = "update_radius_x" if region == "centered_ellipse" else "update_radius"
        try:
            radius_x = float(metrics[key])
        except (KeyError, ValueError) as error:
            raise LadderError("cannot infer source x radius") from error
    if radius_y is None:
        key = "update_radius_y" if region == "centered_ellipse" else "update_radius"
        try:
            radius_y = float(metrics[key])
        except (KeyError, ValueError) as error:
            raise LadderError("cannot infer source y radius") from error
    if radius_x <= 0.0 or radius_y <= 0.0:
        raise LadderError("source active radii must be positive")
    return radius_x, radius_y


def source_products(checkpoint: Path) -> tuple[Path, Path]:
    field = checkpoint.parent / "final_fields_2d.dat"
    residuals = checkpoint.parent / "sampled_residuals.dat"
    if not field.is_file() or not residuals.is_file():
        raise LadderError(
            "checkpoint needs sibling final_fields_2d.dat and sampled_residuals.dat"
        )
    return field, residuals


def diagnostic_options(arguments: argparse.Namespace) -> dict[str, float]:
    return {
        "collar_fraction": arguments.collar_fraction,
        "core_fraction": arguments.core_fraction,
        "interface_amplification_limit": arguments.interface_amplification_limit,
        "interface_step_floor": arguments.interface_step_floor,
        "boundary_to_core_limit": arguments.boundary_to_core_limit,
    }


def mesh_spacing_at(template: str, radius: float) -> float:
    fine_limit = namelist_float(template, "scratch_fine_region_half_width")
    medium_limit = namelist_float(template, "scratch_medium_region_half_width")
    if radius <= fine_limit:
        return namelist_float(template, "scratch_fine_spacing")
    if radius <= medium_limit:
        return namelist_float(template, "scratch_medium_spacing")
    return namelist_float(template, "scratch_coarse_spacing")


def grow_radius(radius: float, template: str, factor: float) -> float:
    spacing = mesh_spacing_at(template, radius)
    return max(radius * factor, radius + spacing)


def stage_geometry(
    template: str,
    radius_x: float,
    radius_y: float,
    guard_cells: float,
) -> tuple[float, float, float, float]:
    def axis_surfaces(active_radius: float) -> tuple[float, float]:
        matching = active_radius - guard_cells * mesh_spacing_at(
            template, active_radius
        )
        fit_inner = matching - guard_cells * mesh_spacing_at(template, matching)
        if fit_inner <= 0.0 or matching <= fit_inner:
            raise LadderError(
                "active radius is too small for two guarded fit surfaces"
            )
        return fit_inner, matching

    fit_inner_x, matching_x = axis_surfaces(radius_x)
    fit_inner_y, matching_y = axis_surfaces(radius_y)
    return fit_inner_x, fit_inner_y, matching_x, matching_y


def make_stage_input(
    template: str,
    radius_x: float,
    radius_y: float,
    fit_inner_x: float,
    fit_inner_y: float,
    matching_x: float,
    matching_y: float,
) -> str:
    replacements = {
        "update_region": "'centered_ellipse'",
        "update_radius": format(max(radius_x, radius_y), ".17g"),
        "update_radius_x": format(radius_x, ".17g"),
        "update_radius_y": format(radius_y, ".17g"),
        "checkpoint_scope": "'full_state'",
        "bulk_reference_policy": "'theoretical'",
        "asymptotic_fit_inner_radius": format(
            min(fit_inner_x, fit_inner_y), ".17g"
        ),
        "asymptotic_matching_radius": format(max(matching_x, matching_y), ".17g"),
        "asymptotic_fit_inner_radius_x": format(fit_inner_x, ".17g"),
        "asymptotic_fit_inner_radius_y": format(fit_inner_y, ".17g"),
        "asymptotic_matching_radius_x": format(matching_x, ".17g"),
        "asymptotic_matching_radius_y": format(matching_y, ".17g"),
        "bulk_reference_fit_radius": format(min(matching_x, matching_y), ".17g"),
    }
    result = template
    for key, value in replacements.items():
        result = replace_namelist_value(result, key, value)
    return result


def radius_label(value: float) -> str:
    return format(value, ".5g").replace("-", "m").replace(".", "p")


def sparse_state(path: Path) -> dict[tuple[float, float], list[float]]:
    values: dict[tuple[float, float], list[float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        columns = [float(value) for value in stripped.split()]
        values[(round(columns[0], 12), round(columns[1], 12))] = columns[2:]
    if not values:
        raise LadderError(f"sparse state has no values: {path}")
    return values


def common_checkpoint_change(previous: Path, current: Path) -> dict[str, float | int]:
    left, right = sparse_state(previous), sparse_state(current)
    common = left.keys() & right.keys()
    if not common:
        raise LadderError("successive sparse checkpoints have no common points")
    difference_squared = 0.0
    reference_squared = 0.0
    maximum_point_relative = 0.0
    for coordinate in common:
        a, b = left[coordinate], right[coordinate]
        local_difference = sum((x - y) ** 2 for x, y in zip(a, b))
        local_reference = sum(x * x for x in a)
        difference_squared += local_difference
        reference_squared += local_reference
        maximum_point_relative = max(
            maximum_point_relative,
            math.sqrt(local_difference / max(local_reference, 1.0e-30)),
        )
    return {
        "common_points": len(common),
        "relative_l2": math.sqrt(
            difference_squared / max(reference_squared, 1.0e-30)
        ),
        "maximum_point_relative": maximum_point_relative,
    }


def convert_full_field_to_sparse(field: Path, destination: Path) -> int:
    """Preserve the complete prior field for cells admitted by a later stage."""
    records: list[str] = []
    for line_number, line in enumerate(
        field.read_text(encoding="utf-8").splitlines(), start=1
    ):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        columns = stripped.split()
        if len(columns) != 24:
            raise LadderError(
                f"full field record {line_number} has {len(columns)} columns, "
                "expected 24"
            )
        # The full map includes pair density between the 18 order-parameter
        # values and the three mean fields.  Sparse restart records omit that
        # derived column.
        records.append(" ".join([*columns[:20], *columns[21:24]]))
    if not records:
        raise LadderError(f"full field has no records: {field}")
    destination.write_text(
        "# FermiForge sparse spinful state version 1\n"
        "# label=adaptive_domain_complete_source_state\n"
        f"# point_count={len(records)}\n"
        "# x y Re/Im[A(1,1)..A(3,3)] current_mean_field(1:3)\n"
        + "\n".join(records)
        + "\n",
        encoding="utf-8",
    )
    return len(records)


def find_stage_report(stage_root: Path, run_name: str) -> Path:
    reports = list(stage_root.glob(f"*-{run_name}/scratch_run_report.json"))
    if len(reports) != 1:
        raise LadderError(f"expected one report for {run_name}, found {len(reports)}")
    return reports[0].resolve()


def write_report(path: Path, report: dict[str, object]) -> None:
    report["updated_at"] = dt.datetime.now().astimezone().isoformat()
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        if (
            arguments.block_iterations < 1
            or arguments.maximum_stages < 1
            or arguments.ranks < 1
            or arguments.jobs < 1
            or arguments.growth_factor <= 1.0
            or arguments.guard_cells <= 0.0
            or (
                arguments.maximum_radius_x is not None
                and arguments.maximum_radius_x <= 0.0
            )
            or (
                arguments.maximum_radius_y is not None
                and arguments.maximum_radius_y <= 0.0
            )
        ):
            raise LadderError("invalid adaptive-domain controls")

        run_root = arguments.run_root.expanduser().resolve()
        checkpoint = (
            arguments.restart.expanduser().resolve()
            if arguments.restart is not None
            else latest_source_checkpoint(run_root)
        )
        if not checkpoint.is_file():
            raise LadderError(f"source checkpoint was not found: {checkpoint}")
        template_path = arguments.input.expanduser().resolve()
        template = template_path.read_text(encoding="utf-8")
        half_width = namelist_float(template, "scratch_half_width")
        radius_x, radius_y = infer_active_radii(
            checkpoint, arguments.initial_radius_x, arguments.initial_radius_y
        )
        field, residuals = source_products(checkpoint)
        diagnostics = diagnose_domain(
            field,
            residuals,
            radius_x,
            radius_y,
            **diagnostic_options(arguments),
        )

        timestamp = dt.datetime.now().strftime("%y%m%d-%H%M%S")
        adaptive_directory = run_root / f"{timestamp}-double-core-adaptive-domain"
        stage_root = adaptive_directory / "stages"
        input_root = adaptive_directory / "inputs"
        adaptive_directory.mkdir(parents=True)
        stage_root.mkdir()
        input_root.mkdir()
        report_path = adaptive_directory / "adaptive_domain_report.json"
        report: dict[str, object] = {
            "kind": "fermiforge_double_core_adaptive_domain",
            "source_checkpoint": str(checkpoint),
            "source_full_field": str(field),
            "input_template": str(template_path),
            "controls": {
                "block_iterations": arguments.block_iterations,
                "maximum_stages": arguments.maximum_stages,
                "growth_factor": arguments.growth_factor,
                "guard_cells": arguments.guard_cells,
                "hold_domain": arguments.hold_domain,
                **diagnostic_options(arguments),
            },
            "source": {
                "radius_x": radius_x,
                "radius_y": radius_y,
                "diagnostics": diagnostics,
            },
            "stages": [],
            "terminal_status": "diagnostics_only" if arguments.diagnostics_only else "running",
            "adequate": bool(diagnostics["adequate"]),
        }
        write_report(report_path, report)
        print(f"source checkpoint: {checkpoint}")
        for axis in ("x", "y"):
            result = diagnostics["axes"][axis]
            print(
                f"{axis} radius {result['radius']:g}: "
                f"grow={result['needs_growth']} reasons={result['reasons']}"
            )
        if arguments.diagnostics_only or (
            bool(diagnostics["adequate"]) and not arguments.hold_domain
        ):
            report["terminal_status"] = (
                "adequate" if bool(diagnostics["adequate"]) else "diagnostics_only"
            )
            write_report(report_path, report)
            print(f"report: {report_path}")
            return 0 if bool(diagnostics["adequate"]) else 1

        build_directory = arguments.build_dir.expanduser().resolve()
        executable = build_directory / "fermiforge_double_core_2d"
        if not arguments.no_build:
            executable = build_solver(build_directory, arguments.jobs)
        elif not executable.is_file():
            raise LadderError(f"solver executable was not found: {executable}")

        complete_restart = adaptive_directory / "source_full_state_restart.dat"
        complete_restart_points = convert_full_field_to_sparse(
            field, complete_restart
        )
        report["complete_source_restart"] = str(complete_restart)
        report["complete_source_restart_points"] = complete_restart_points
        checkpoint = complete_restart
        write_report(report_path, report)

        terminal_status = "maximum_stages"
        for stage_number in range(1, arguments.maximum_stages + 1):
            if arguments.hold_domain:
                grow_x = False
                grow_y = False
                requested_x = radius_x
                requested_y = radius_y
            else:
                grow_x = bool(diagnostics["axes"]["x"]["needs_growth"])
                grow_y = bool(diagnostics["axes"]["y"]["needs_growth"])
                requested_x = (
                    grow_radius(radius_x, template, arguments.growth_factor)
                    if grow_x
                    else radius_x
                )
                requested_y = (
                    grow_radius(radius_y, template, arguments.growth_factor)
                    if grow_y
                    else radius_y
                )
            guard_x = arguments.guard_cells * mesh_spacing_at(template, requested_x)
            guard_y = arguments.guard_cells * mesh_spacing_at(template, requested_y)
            limit_x = arguments.maximum_radius_x or (half_width - guard_x)
            limit_y = arguments.maximum_radius_y or (half_width - guard_y)
            new_x = min(requested_x, limit_x)
            new_y = min(requested_y, limit_y)
            if not arguments.hold_domain and (
                (grow_x and new_x <= radius_x) or (grow_y and new_y <= radius_y)
            ):
                terminal_status = "cell_limit"
                break
            fit_inner_x, fit_inner_y, matching_x, matching_y = stage_geometry(
                template, new_x, new_y, arguments.guard_cells
            )
            if max(matching_x, matching_y) >= half_width:
                terminal_status = "cell_limit"
                break

            stage_kind = "hold" if arguments.hold_domain else "adaptive"
            run_name = (
                f"{stage_kind}-{stage_number:02d}-rx{radius_label(new_x)}-"
                f"ry{radius_label(new_y)}"
            )
            archived_input = input_root / f"{run_name}.nml"
            archived_input.write_text(
                make_stage_input(
                    template,
                    new_x,
                    new_y,
                    fit_inner_x,
                    fit_inner_y,
                    matching_x,
                    matching_y,
                ),
                encoding="utf-8",
            )
            command = [
                sys.executable,
                str(PROJECT_ROOT / "tools" / "run_double_core_from_scratch.py"),
                "--input",
                str(archived_input),
                "--run-name",
                run_name,
                "--run-root",
                str(stage_root),
                "--iterations",
                str(arguments.block_iterations),
                "--checkpoint-every",
                "1",
                "--ranks",
                str(arguments.ranks),
                "--jobs",
                str(arguments.jobs),
                "--build-dir",
                str(build_directory),
                "--no-build",
                "--sparse-restart",
                str(checkpoint),
            ]
            if arguments.no_plot:
                command.append("--no-plot")
            print(
                f"\n=== adaptive stage {stage_number}: "
                f"{new_x:g} x {new_y:g} xi0 ===",
                flush=True,
            )
            stream_command(command, adaptive_directory / f"stage-{stage_number:02d}.log")
            stage_report_path = find_stage_report(stage_root, run_name)
            stage_report = json.loads(stage_report_path.read_text(encoding="utf-8"))
            if not bool(stage_report.get("passed")):
                raise LadderError(f"adaptive stage {stage_number} failed")
            next_checkpoint = Path(str(stage_report["checkpoint"])).resolve()
            next_field, next_residuals = source_products(next_checkpoint)
            next_diagnostics = diagnose_domain(
                next_field,
                next_residuals,
                new_x,
                new_y,
                **diagnostic_options(arguments),
            )
            stage_record = {
                "stage": stage_number,
                "radius_x": new_x,
                "radius_y": new_y,
                "fit_inner_radius_x": fit_inner_x,
                "fit_inner_radius_y": fit_inner_y,
                "matching_radius_x": matching_x,
                "matching_radius_y": matching_y,
                "fit_inner_radius": min(fit_inner_x, fit_inner_y),
                "matching_radius": max(matching_x, matching_y),
                "grew_x": grow_x,
                "grew_y": grow_y,
                "checkpoint_change": common_checkpoint_change(
                    checkpoint, next_checkpoint
                ),
                "diagnostics": next_diagnostics,
                "run_report": str(stage_report_path),
                "checkpoint": str(next_checkpoint),
                "axis_profile_plot": str(
                    stage_report_path.parent
                    / "plots"
                    / "double_core_axis_profiles.png"
                ),
            }
            report["stages"].append(stage_record)
            checkpoint = next_checkpoint
            radius_x, radius_y = new_x, new_y
            diagnostics = next_diagnostics
            report["adequate"] = bool(diagnostics["adequate"])
            write_report(report_path, report)
            if arguments.hold_domain:
                terminal_status = "hold_domain_complete"
                break
            if bool(diagnostics["adequate"]):
                terminal_status = "adequate"
                break

        report["terminal_status"] = terminal_status
        report["adequate"] = bool(diagnostics["adequate"])
        report["final_checkpoint"] = str(checkpoint)
        report["final_radius_x"] = radius_x
        report["final_radius_y"] = radius_y
        if terminal_status == "cell_limit":
            report["recommendation"] = (
                "enlarge or remesh the physical cell before further active-domain growth"
            )
        write_report(report_path, report)
        print(f"\nterminal status: {terminal_status}")
        print(f"report:          {report_path}")
        print(f"checkpoint:      {checkpoint}")
        return 0 if terminal_status in ("adequate", "hold_domain_complete") else 1
    except (OSError, ValueError, DomainDiagnosticError, LadderError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
