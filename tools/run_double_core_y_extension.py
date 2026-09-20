#!/usr/bin/env python3
"""Extend a converged double-core state along y through nested checkpoints."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path
from typing import Sequence

from run_double_core_from_scratch import build_solver
from run_double_core_radius_ladder import (
    LadderError,
    replace_namelist_value,
    stream_command,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "examples" / "2d_double_core_y_extension.nml"
DEFAULT_BUILD = PROJECT_ROOT / "work" / "double-core-scratch-build"

STAGES = (
    {
        "name": "y28",
        "radius_x": 22.0,
        "radius_y": 28.0,
        "fit_inner_radius": 26.0,
        "matching_radius": 32.0,
    },
    {
        "name": "y34",
        "radius_x": 22.0,
        "radius_y": 34.0,
        "fit_inner_radius": 32.0,
        "matching_radius": 36.0,
    },
)


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Continue the large double-core calculation through nested "
            "22x28 and 22x34 xi0 self-consistent ellipses."
        )
    )
    parser.add_argument(
        "--restart",
        type=Path,
        help=(
            "sparse checkpoint used by the first selected stage; by default "
            "the newest double-core-overnight checkpoint is selected"
        ),
    )
    parser.add_argument(
        "--start-stage",
        choices=tuple(stage["name"] for stage in STAGES),
        default="y28",
        help="first stage to run (use y34 when resuming from a y28 checkpoint)",
    )
    parser.add_argument("--iterations-per-stage", type=int, default=8)
    parser.add_argument("--ranks", type=int, default=10)
    parser.add_argument("--jobs", type=int, default=10)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD)
    parser.add_argument("--run-root", type=Path, default=PROJECT_ROOT / "runs")
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args(argv)


def latest_overnight_checkpoint(run_root: Path) -> Path:
    candidates = [
        path
        for path in run_root.glob("*-double-core-overnight/sparse_checkpoint_state.dat")
        if path.is_file()
    ]
    if not candidates:
        raise LadderError(
            "no double-core overnight checkpoint was found; pass --restart explicitly"
        )
    return max(candidates, key=lambda path: path.stat().st_mtime).resolve()


def selected_stages(first_name: str) -> list[dict[str, float | str]]:
    first_index = next(
        index for index, stage in enumerate(STAGES) if stage["name"] == first_name
    )
    return [dict(stage) for stage in STAGES[first_index:]]


def stage_input(template: str, stage: dict[str, float | str]) -> str:
    replacements = {
        "update_region": "'centered_ellipse'",
        "update_radius_x": format(float(stage["radius_x"]), ".17g"),
        "update_radius_y": format(float(stage["radius_y"]), ".17g"),
        "asymptotic_fit_inner_radius": format(
            float(stage["fit_inner_radius"]), ".17g"
        ),
        "asymptotic_matching_radius": format(
            float(stage["matching_radius"]), ".17g"
        ),
        "bulk_reference_fit_radius": format(
            float(stage["matching_radius"]), ".17g"
        ),
    }
    result = template
    for key, value in replacements.items():
        result = replace_namelist_value(result, key, value)
    return result


def find_stage_report(stage_root: Path, stage_name: str) -> Path:
    reports = list(
        stage_root.glob(f"*-double-core-{stage_name}/scratch_run_report.json")
    )
    if len(reports) != 1:
        raise LadderError(
            f"expected one report for stage {stage_name}, found {len(reports)}"
        )
    return reports[0].resolve()


def write_ladder_report(
    path: Path,
    source_checkpoint: Path,
    arguments: argparse.Namespace,
    results: list[dict[str, object]],
    complete: bool,
) -> None:
    report = {
        "kind": "fermiforge_double_core_y_extension",
        "updated_at": dt.datetime.now().astimezone().isoformat(),
        "source_checkpoint": str(source_checkpoint),
        "start_stage": arguments.start_stage,
        "iterations_per_stage": arguments.iterations_per_stage,
        "ranks": arguments.ranks,
        "results": results,
        "complete": complete,
        "passed": complete and all(bool(item["report"]["passed"]) for item in results),
    }
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        if (
            arguments.iterations_per_stage < 1
            or arguments.ranks < 1
            or arguments.jobs < 1
        ):
            raise LadderError("iterations, ranks, and jobs must be positive")

        run_root = arguments.run_root.expanduser().resolve()
        restart = (
            arguments.restart.expanduser().resolve()
            if arguments.restart is not None
            else latest_overnight_checkpoint(run_root)
        )
        if not restart.is_file():
            raise LadderError(f"sparse restart was not found: {restart}")

        input_path = arguments.input.expanduser().resolve()
        template = input_path.read_text(encoding="utf-8")
        build_directory = arguments.build_dir.expanduser().resolve()
        executable = build_directory / "fermiforge_double_core_2d"
        if not arguments.no_build:
            executable = build_solver(build_directory, arguments.jobs)
        elif not executable.is_file():
            raise LadderError(f"solver executable was not found: {executable}")

        timestamp = dt.datetime.now().strftime("%y%m%d-%H%M%S")
        ladder_directory = run_root / f"{timestamp}-double-core-y-extension"
        stage_root = ladder_directory / "stages"
        template_root = ladder_directory / "inputs"
        stage_root.mkdir(parents=True)
        template_root.mkdir()
        report_path = ladder_directory / "y_extension_report.json"
        source_checkpoint = restart
        results: list[dict[str, object]] = []
        write_ladder_report(
            report_path, source_checkpoint, arguments, results, complete=False
        )

        print(f"source checkpoint: {source_checkpoint}")
        print(f"ladder directory:  {ladder_directory}")
        print(
            "stages:            "
            + ", ".join(str(stage["name"]) for stage in selected_stages(arguments.start_stage))
        )

        for stage in selected_stages(arguments.start_stage):
            name = str(stage["name"])
            archived_template = template_root / f"{name}.nml"
            archived_template.write_text(
                stage_input(template, stage), encoding="utf-8"
            )
            command = [
                sys.executable,
                str(PROJECT_ROOT / "tools" / "run_double_core_from_scratch.py"),
                "--input",
                str(archived_template),
                "--run-name",
                f"double-core-{name}",
                "--run-root",
                str(stage_root),
                "--iterations",
                str(arguments.iterations_per_stage),
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
                str(restart),
            ]
            if arguments.no_plot:
                command.append("--no-plot")
            print(
                f"\n=== {name}: active ellipse "
                f"{stage['radius_x']:g} x {stage['radius_y']:g} xi0 ===",
                flush=True,
            )
            stream_command(command, ladder_directory / f"{name}.log")
            stage_report_path = find_stage_report(stage_root, name)
            stage_report = json.loads(stage_report_path.read_text(encoding="utf-8"))
            if not bool(stage_report.get("passed")):
                raise LadderError(f"stage {name} did not pass its acceptance checks")
            checkpoint_value = stage_report.get("checkpoint")
            if not checkpoint_value:
                raise LadderError(f"stage {name} did not report a checkpoint")
            restart = Path(str(checkpoint_value)).resolve()
            if not restart.is_file():
                raise LadderError(f"stage {name} checkpoint was not found: {restart}")
            results.append(
                {
                    "stage": stage,
                    "input": str(archived_template),
                    "report_path": str(stage_report_path),
                    "checkpoint": str(restart),
                    "axis_profile_plot": str(
                        stage_report_path.parent
                        / "plots"
                        / "double_core_axis_profiles.png"
                    ),
                    "report": stage_report,
                }
            )
            write_ladder_report(
                report_path, source_checkpoint, arguments, results, complete=False
            )

        write_ladder_report(
            report_path, source_checkpoint, arguments, results, complete=True
        )
        print(f"\nreport:           {report_path}")
        print(f"final checkpoint: {restart}")
        return 0
    except (OSError, ValueError, LadderError, StopIteration) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
