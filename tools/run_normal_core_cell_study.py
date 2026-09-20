#!/usr/bin/env python3
"""Run converged normal-core solves in several cells and compare the results."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "benchmarks" / "normal_core_2d" / "cell_study.json"
RUNNER = PROJECT_ROOT / "tools" / "run_axisymmetric_core_benchmark.py"


class StudyError(RuntimeError):
    """Expected configuration, run, or comparison failure."""


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the converged normal-core 2D solver for several computational "
            "cells and compare every result with the radial qcv reference."
        )
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--ranks", type=int, help="override MPI ranks")
    parser.add_argument("--jobs", type=int, help="override build jobs")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--run-root", type=Path, default=PROJECT_ROOT / "runs")
    return parser.parse_args(argv)


def read_metrics(path: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        metrics[key.strip()] = value.strip()
    return metrics


def metric_float(metrics: dict[str, str], key: str) -> float:
    try:
        return float(metrics[key])
    except (KeyError, ValueError) as error:
        raise StudyError(f"metrics do not contain a numeric '{key}'") from error


def read_field_values(path: Path, radius: float) -> dict[tuple[float, float], tuple[list[complex], list[float]]]:
    values: dict[tuple[float, float], tuple[list[complex], list[float]]] = {}
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            row = [float(value) for value in line.split()]
            if len(row) != 24:
                raise StudyError(f"unexpected field-map record in {path}")
            x, y = row[0:2]
            if math.hypot(x, y) > radius + 1.0e-12:
                continue
            gap = [complex(row[index], row[index + 1]) for index in range(2, 20, 2)]
            current = row[21:24]
            values[(round(x, 12), round(y, 12))] = (gap, current)
    return values


def compare_field_maps(
    baseline_path: Path, candidate_path: Path, radius: float
) -> dict[str, float | int]:
    baseline = read_field_values(baseline_path, radius)
    candidate = read_field_values(candidate_path, radius)
    common = sorted(baseline.keys() & candidate.keys())
    if not common:
        raise StudyError("cell comparison found no common Cartesian nodes")

    gap_difference_squared = 0.0
    gap_reference_squared = 0.0
    current_difference_squared = 0.0
    current_reference_squared = 0.0
    gap_maximum_difference = 0.0
    gap_maximum_reference = 0.0
    current_maximum_difference = 0.0
    current_maximum_reference = 0.0
    for coordinate in common:
        baseline_gap, baseline_current = baseline[coordinate]
        candidate_gap, candidate_current = candidate[coordinate]
        for reference, value in zip(baseline_gap, candidate_gap):
            difference = abs(value - reference)
            gap_difference_squared += difference * difference
            gap_reference_squared += abs(reference) ** 2
            gap_maximum_difference = max(gap_maximum_difference, difference)
            gap_maximum_reference = max(gap_maximum_reference, abs(reference))
        for reference, value in zip(baseline_current, candidate_current):
            difference = abs(value - reference)
            current_difference_squared += difference * difference
            current_reference_squared += reference * reference
            current_maximum_difference = max(current_maximum_difference, difference)
            current_maximum_reference = max(current_maximum_reference, abs(reference))

    tiny = sys.float_info.min
    return {
        "common_points": len(common),
        "gap_max_relative_to_baseline": gap_maximum_difference
        / max(gap_maximum_reference, tiny),
        "gap_relative_l2_to_baseline": math.sqrt(
            gap_difference_squared / max(gap_reference_squared, tiny)
        ),
        "current_max_relative_to_baseline": current_maximum_difference
        / max(current_maximum_reference, tiny),
        "current_relative_l2_to_baseline": math.sqrt(
            current_difference_squared / max(current_reference_squared, tiny)
        ),
    }


def append_option(command: list[str], option: str, value: Any) -> None:
    if value is not None:
        command.extend((option, str(value)))


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        config_path = arguments.config.expanduser().resolve()
        config = json.loads(config_path.read_text(encoding="utf-8"))
        cases = config.get("cases")
        if not isinstance(cases, list) or len(cases) < 2:
            raise StudyError("cell study needs at least two cases")

        timestamp = dt.datetime.now().strftime("%y%m%d-%H%M%S")
        study_directory = arguments.run_root.expanduser().resolve() / (
            f"{timestamp}-normal-core-cell-study"
        )
        study_directory.mkdir(parents=True)
        archived_config = study_directory / "cell_study.json"
        archived_config.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

        shared = config.get("shared", {})
        ranks = arguments.ranks or int(shared.get("ranks", 10))
        jobs = arguments.jobs or int(shared.get("jobs", ranks))
        input_path = (PROJECT_ROOT / shared.get("input", "examples/2d_normal_core_multiscale.nml")).resolve()
        build_directory = study_directory / "build"
        results: list[dict[str, Any]] = []

        for index, case in enumerate(cases):
            name = str(case["name"])
            case_root = study_directory / "cases" / name
            command = [
                sys.executable,
                str(RUNNER),
                "--input",
                str(input_path),
                "--case-name",
                f"normal-core-{name}",
                "--initialization",
                "radial-qcv-reference",
                "--ranks",
                str(ranks),
                "--jobs",
                str(jobs),
                "--build-dir",
                str(build_directory),
                "--run-root",
                str(case_root),
            ]
            if index > 0:
                command.append("--no-build")
            if arguments.no_plot or not bool(shared.get("plot", True)):
                command.append("--no-plot")
            options = {
                "--half-width": case.get("half_width"),
                "--active-radius": case.get("active_radius"),
                "--fine-region": case.get("fine_region_half_width"),
                "--medium-region": case.get("medium_region_half_width"),
                "--fine-spacing": case.get("fine_spacing"),
                "--medium-spacing": case.get("medium_spacing"),
                "--coarse-spacing": case.get("coarse_spacing"),
                "--outer-radius": case.get("asymptotic_outer_radius"),
                "--trajectory-step": case.get("trajectory_maximum_step"),
                "--max-iterations": shared.get("maximum_iterations"),
                "--tolerance": shared.get("convergence_tolerance"),
                "--anderson-history": shared.get("anderson_history_limit"),
                "--anderson-progress": shared.get("anderson_progress_threshold"),
                "--anderson-pmax": shared.get("anderson_maximum_mixing"),
                "--checkpoint-interval": shared.get("checkpoint_interval"),
                "--perturbation": shared.get("initial_perturbation_amplitude"),
                "--perturbation-radius": shared.get("initial_perturbation_radius"),
            }
            for option, value in options.items():
                append_option(command, option, value)

            print(f"\n=== normal-core cell case: {name} ===", flush=True)
            completed = subprocess.run(command, cwd=PROJECT_ROOT, check=False)
            if completed.returncode != 0:
                raise StudyError(f"cell case '{name}' failed")
            run_directories = sorted(path for path in case_root.iterdir() if path.is_dir())
            if len(run_directories) != 1:
                raise StudyError(f"could not identify output for cell case '{name}'")
            run_directory = run_directories[0]
            metrics = read_metrics(run_directory / "metrics.txt")
            results.append(
                {
                    "name": name,
                    "run_directory": str(run_directory),
                    "final_field_map": str(run_directory / "final_fields_2d.dat"),
                    "metrics": metrics,
                }
            )

        comparison_radius = float(config.get("comparison_radius", 7.0))
        baseline_map = Path(results[0]["final_field_map"])
        results[0]["cell_comparison"] = {
            "common_points": 0,
            "gap_max_relative_to_baseline": 0.0,
            "gap_relative_l2_to_baseline": 0.0,
            "current_max_relative_to_baseline": 0.0,
            "current_relative_l2_to_baseline": 0.0,
        }
        for result in results[1:]:
            result["cell_comparison"] = compare_field_maps(
                baseline_map, Path(result["final_field_map"]), comparison_radius
            )

        acceptance = config.get("acceptance", {})
        checks: list[dict[str, Any]] = []
        for result in results:
            metrics = result["metrics"]
            comparison = result["cell_comparison"]
            case_checks = {
                "converged": metrics.get("terminal_status") == "converged",
                "solver_passed": metrics.get("passed", "F").upper() == "T",
                "normalization": metric_float(
                    metrics, "maximum_normalization_error_over_solve"
                )
                <= float(acceptance.get("maximum_normalization_error", 1.0e-10)),
                "qcv_gap_l2": metric_float(metrics, "final_gap_radial_relative_l2")
                <= float(acceptance.get("maximum_qcv_gap_relative_l2", 5.0e-2)),
                "qcv_current_l2": metric_float(
                    metrics, "final_current_radial_relative_l2"
                )
                <= float(acceptance.get("maximum_qcv_current_relative_l2", 5.0e-2)),
                "cell_gap_l2": float(comparison["gap_relative_l2_to_baseline"])
                <= float(acceptance.get("maximum_cell_gap_relative_l2", 1.0e-2)),
                "cell_current_l2": float(comparison["current_relative_l2_to_baseline"])
                <= float(acceptance.get("maximum_cell_current_relative_l2", 1.0e-2)),
            }
            checks.append({"name": result["name"], **case_checks})

        summary_rows: list[dict[str, Any]] = []
        for result, case_checks in zip(results, checks):
            metrics = result["metrics"]
            comparison = result["cell_comparison"]
            summary_rows.append(
                {
                    "case": result["name"],
                    "status": metrics.get("terminal_status", "missing"),
                    "passed": all(value for key, value in case_checks.items() if key != "name"),
                    "iterations": metrics.get("iterations_completed", ""),
                    "points": metrics.get("cartesian_points", ""),
                    "active_points": metrics.get("active_points", ""),
                    "minimum_spacing": metrics.get("minimum_spacing", ""),
                    "maximum_spacing": metrics.get("maximum_spacing", ""),
                    "qcv_gap_relative_l2": metrics.get("final_gap_radial_relative_l2", ""),
                    "qcv_current_relative_l2": metrics.get("final_current_radial_relative_l2", ""),
                    "cell_gap_relative_l2": comparison["gap_relative_l2_to_baseline"],
                    "cell_current_relative_l2": comparison[
                        "current_relative_l2_to_baseline"
                    ],
                    "normalization_error": metrics.get(
                        "maximum_normalization_error_over_solve", ""
                    ),
                    "map_seconds": metrics.get("cumulative_map_seconds", ""),
                    "run_directory": result["run_directory"],
                }
            )

        csv_path = study_directory / "cell_study_summary.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
            writer.writeheader()
            writer.writerows(summary_rows)
        report = {
            "kind": "fermiforge_normal_core_cell_study",
            "created_at": dt.datetime.now().astimezone().isoformat(),
            "configuration": str(archived_config),
            "comparison_radius": comparison_radius,
            "acceptance": acceptance,
            "checks": checks,
            "results": results,
            "passed": all(row["passed"] for row in summary_rows),
        }
        report_path = study_directory / "cell_study_report.json"
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

        print("\ncase                         status       qcv gap L2   cell gap L2  pass")
        for row in summary_rows:
            print(
                f"{row['case']:<28} {row['status']:<12} "
                f"{float(row['qcv_gap_relative_l2']):11.3e} "
                f"{float(row['cell_gap_relative_l2']):11.3e} "
                f"{str(row['passed']):>5}"
            )
        print(f"\nstudy directory: {study_directory}")
        print(f"summary:         {csv_path}")
        print(f"report:          {report_path}")
        return 0 if report["passed"] else 1
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError, StudyError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
