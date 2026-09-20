#!/usr/bin/env python3
"""Run an expanding, sparse-checkpointed double-core active-radius ladder."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "examples" / "2d_double_core_patch_smoke.nml"
DEFAULT_BUILD = PROJECT_ROOT / "work" / "double-core-ladder-build"


class LadderError(RuntimeError):
    """Expected build, input, launch, or output failure."""


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Expand the paired double-core update disks through a sequence of "
            "radii, carrying only the sparse relaxed state between cases."
        )
    )
    parser.add_argument("--radii", nargs="+", type=float, default=(0.4, 0.6))
    parser.add_argument("--iterations", type=int, default=2)
    parser.add_argument("--ranks", type=int, default=10)
    parser.add_argument("--jobs", type=int, default=10)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD)
    parser.add_argument("--run-root", type=Path, default=PROJECT_ROOT / "runs")
    parser.add_argument(
        "--initial-sparse-restart",
        type=Path,
        help="optional sparse checkpoint carried into the first radius",
    )
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument(
        "--zero-mode-projection",
        choices=(
            "none",
            "translations",
            "translation_orientation",
            "gauge_translation_orientation",
        ),
        default="none",
        help="collective-coordinate components removed before Anderson mixing",
    )
    parser.add_argument(
        "--inactive-halo-policy",
        choices=("frozen", "state_asymptotic"),
        default="frozen",
        help="treatment of stored nodes outside the asymptotic matching circle",
    )
    parser.add_argument(
        "--independent",
        action="store_true",
        help="start every radius from the archived field instead of continuing",
    )
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args(argv)


def find_program(name: str, fallbacks: Sequence[Path] = ()) -> Path:
    found = shutil.which(name)
    if found:
        return Path(found).resolve()
    for candidate in fallbacks:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise LadderError(f"required program '{name}' was not found")


def replace_namelist_value(text: str, key: str, value: str) -> str:
    pattern = re.compile(
        rf'''(?im)^(\s*{re.escape(key)}\s*=\s*)'''
        r'''(?:'[^']*'|"[^"]*"|[^,\n/]+)'''
    )
    updated, count = pattern.subn(rf"\g<1>{value}", text, count=1)
    if count != 1:
        raise LadderError(f"input namelist does not define '{key}'")
    return updated


def fortran_string(path: Path | str) -> str:
    return "'" + str(path).replace("'", "''") + "'"


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
        raise LadderError(f"metrics do not contain numeric '{key}'") from error


def metric_int(metrics: dict[str, str], key: str) -> int:
    try:
        return int(metrics[key])
    except (KeyError, ValueError) as error:
        raise LadderError(f"metrics do not contain integer '{key}'") from error


def read_history(path: Path) -> list[dict[str, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("#"):
        raise LadderError(f"iteration history has no header: {path}")
    names = lines[0][1:].split()
    rows: list[dict[str, float]] = []
    for line in lines[1:]:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = [float(value) for value in line.split()]
        if len(values) != len(names):
            raise LadderError(f"malformed iteration history row in {path}")
        rows.append(dict(zip(names, values)))
    if not rows:
        raise LadderError(f"iteration history is empty: {path}")
    return rows


def sparse_checkpoint_point_count(path: Path) -> int:
    count = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            count += 1
    if count < 1:
        raise LadderError(f"sparse checkpoint has no field records: {path}")
    return count


def stream_command(command: Sequence[str], log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as log:
        log.write("command: " + " ".join(command) + "\n")
        log.flush()
        process = subprocess.Popen(
            command,
            cwd=PROJECT_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert process.stdout is not None
        for line in process.stdout:
            log.write(line)
            log.flush()
            print(line, end="")
        return_code = process.wait()
    if return_code != 0:
        raise LadderError(
            f"command failed with status {return_code}; see {log_path}"
        )


def build_solver(build_directory: Path, jobs: int) -> Path:
    cmake = find_program(
        "cmake", (Path("/opt/local/bin/cmake"), Path("/opt/homebrew/bin/cmake"))
    )
    compiler = find_program("gfortran", (Path("/opt/homebrew/bin/gfortran"),))
    build_directory.mkdir(parents=True, exist_ok=True)
    configure = [
        str(cmake),
        "-S",
        str(PROJECT_ROOT),
        "-B",
        str(build_directory),
        f"-DCMAKE_Fortran_COMPILER={compiler}",
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    build = [
        str(cmake),
        "--build",
        str(build_directory),
        "--parallel",
        str(jobs),
        "--target",
        "benchmark_legacy_double_core_map_2d",
    ]
    stream_command(configure, build_directory / "configure.log")
    stream_command(build, build_directory / "build.log")
    executable = build_directory / "benchmark_legacy_double_core_map_2d"
    if not executable.is_file():
        raise LadderError(f"build did not create {executable}")
    return executable


def mpi_prefix(ranks: int) -> list[str]:
    mpiexec = find_program(
        "mpiexec",
        (Path("/opt/homebrew/bin/mpiexec"), Path("/opt/local/bin/mpiexec")),
    )
    command = [str(mpiexec)]
    version = subprocess.run(
        [str(mpiexec), "--version"],
        capture_output=True,
        text=True,
        check=False,
    )
    if platform.system() == "Darwin" and "Open MPI" in (
        version.stdout + version.stderr
    ):
        command.extend(
            (
                "--host",
                f"localhost:{ranks}",
                "--map-by",
                f"ppr:{ranks}:node",
                "--bind-to",
                "none",
            )
        )
    command.extend(("-n", str(ranks)))
    return command


def radius_name(radius: float) -> str:
    return "r" + format(radius, ".6g").replace("-", "m").replace(".", "p")


def write_plot(rows: list[dict[str, object]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    radii = [float(row["radius"]) for row in rows]
    maxima = [float(row["update_maximum"]) for row in rows]
    relative = [float(row["sample_relative_l2"]) for row in rows]
    points = [int(row["updated_points"]) for row in rows]
    seconds = [float(row["cumulative_map_seconds"]) for row in rows]
    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    axes[0].semilogy(radii, maxima, "o-", color="#0072B2", label="update max")
    axes[0].semilogy(
        radii, relative, "s-", color="#D55E00", label="sample relative L2"
    )
    axes[0].set_xlabel(r"update radius [$\xi_0$]")
    axes[0].set_ylabel("residual")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend()
    axes[1].plot(radii, points, "o-", color="#009E73", label="updated points")
    axes[1].set_xlabel(r"update radius [$\xi_0$]")
    axes[1].set_ylabel("updated points", color="#009E73")
    axes[1].tick_params(axis="y", labelcolor="#009E73")
    timing_axis = axes[1].twinx()
    timing_axis.plot(radii, seconds, "s--", color="#CC79A7", label="map time")
    timing_axis.set_ylabel("cumulative map time [s]", color="#CC79A7")
    timing_axis.tick_params(axis="y", labelcolor="#CC79A7")
    axes[1].grid(True, alpha=0.25)
    figure.suptitle("FermiForge double-core active-radius ladder")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        radii = list(arguments.radii)
        if not radii or any(radius <= 0.0 for radius in radii):
            raise LadderError("all update radii must be positive")
        if any(right <= left for left, right in zip(radii, radii[1:])):
            raise LadderError("update radii must be strictly increasing")
        if arguments.iterations < 1 or arguments.ranks < 1 or arguments.jobs < 1:
            raise LadderError("iterations, ranks, and jobs must be positive")

        input_path = arguments.input.expanduser().resolve()
        template = input_path.read_text(encoding="utf-8")
        build_directory = arguments.build_dir.expanduser().resolve()
        executable = build_directory / "benchmark_legacy_double_core_map_2d"
        if not arguments.no_build:
            executable = build_solver(build_directory, arguments.jobs)
        elif not executable.is_file():
            raise LadderError(f"solver executable was not found: {executable}")

        timestamp = dt.datetime.now().strftime("%y%m%d-%H%M%S")
        ladder_directory = arguments.run_root.expanduser().resolve() / (
            f"{timestamp}-double-core-radius-ladder"
        )
        ladder_directory.mkdir(parents=True)
        rows: list[dict[str, object]] = []
        previous_sparse: Path | None = None
        if arguments.initial_sparse_restart is not None:
            previous_sparse = arguments.initial_sparse_restart.expanduser().resolve()
            if not previous_sparse.is_file():
                raise LadderError(
                    f"initial sparse checkpoint was not found: {previous_sparse}"
                )
        mpi = mpi_prefix(arguments.ranks)

        for radius in radii:
            case_directory = ladder_directory / radius_name(radius)
            case_directory.mkdir()
            metrics_path = case_directory / "metrics.txt"
            history_path = case_directory / "iteration_history.dat"
            samples_path = case_directory / "sampled_residuals.dat"
            sparse_path = case_directory / "sparse_final_state.dat"
            restart_path = (
                previous_sparse
                if previous_sparse is not None and not arguments.independent
                else None
            )
            expected_restart_points = (
                sparse_checkpoint_point_count(restart_path)
                if restart_path is not None
                else 0
            )
            case_input = template
            replacements = {
                "maximum_iterations": str(arguments.iterations),
                "update_region": "'half_core_disks'",
                "update_radius": format(radius, ".17g"),
                "zero_mode_projection": fortran_string(
                    arguments.zero_mode_projection
                ),
                "inactive_halo_policy": fortran_string(
                    arguments.inactive_halo_policy
                ),
                "sparse_restart_file": fortran_string(restart_path or ""),
                "sparse_final_file": fortran_string(sparse_path),
                "final_field_file": "''",
                "history_file": fortran_string(history_path),
                "metrics_file": fortran_string(metrics_path),
                "sample_file": fortran_string(samples_path),
            }
            for key, value in replacements.items():
                case_input = replace_namelist_value(case_input, key, value)
            archived_input = case_directory / "input.nml"
            archived_input.write_text(case_input, encoding="utf-8")

            print(f"\n=== double-core update radius {radius:g} xi0 ===", flush=True)
            stream_command(
                [*mpi, str(executable), str(archived_input)],
                case_directory / "run.log",
            )
            metrics = read_metrics(metrics_path)
            history = read_history(history_path)
            first, last = history[0], history[-1]
            row: dict[str, object] = {
                "radius": radius,
                "status": metrics.get("terminal_status", "missing"),
                "updated_points": metric_int(metrics, "updated_points"),
                "anderson_values": metric_int(metrics, "anderson_vector_size"),
                "restart_points": metric_int(
                    metrics, "sparse_restart_applied_points"
                ),
                "expected_restart_points": expected_restart_points,
                "zero_mode_projection": metrics.get(
                    "zero_mode_projection", "none"
                ),
                "zero_mode_removed_fraction": float(
                    metrics.get("zero_mode_removed_fraction", "0")
                ),
                "inactive_halo_policy": metrics.get(
                    "inactive_halo_policy", "frozen"
                ),
                "dependent_halo_points": metric_int(
                    metrics, "dependent_asymptotic_halo_points"
                ),
                "last_halo_maximum_change": metric_float(
                    metrics, "last_halo_maximum_change"
                ),
                "sample_relative_l2": metric_float(
                    metrics, "sample_relative_l2"
                ),
                "core_relative_l2": metric_float(metrics, "core_relative_l2"),
                "update_initial_maximum": first["update_max"],
                "update_maximum": last["update_max"],
                "final_mixing": last["p"],
                "normalization_error": metric_float(
                    metrics, "maximum_normalization_error"
                ),
                "cumulative_map_seconds": metric_float(
                    metrics, "cumulative_map_seconds"
                ),
                "case_directory": str(case_directory),
                "sparse_state": str(sparse_path),
            }
            row["passed"] = (
                row["status"] in ("iteration_limit", "converged")
                and float(row["normalization_error"]) <= 1.0e-10
                and sparse_path.is_file()
                and (
                    restart_path is None
                    or int(row["restart_points"]) == expected_restart_points
                )
            )
            rows.append(row)
            previous_sparse = sparse_path

        csv_path = ladder_directory / "radius_ladder_summary.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        report = {
            "kind": "fermiforge_double_core_radius_ladder",
            "created_at": dt.datetime.now().astimezone().isoformat(),
            "input_template": str(input_path),
            "initial_sparse_restart": (
                str(arguments.initial_sparse_restart.expanduser().resolve())
                if arguments.initial_sparse_restart is not None
                else None
            ),
            "continuation": not arguments.independent,
            "ranks": arguments.ranks,
            "iterations_per_radius": arguments.iterations,
            "zero_mode_projection": arguments.zero_mode_projection,
            "inactive_halo_policy": arguments.inactive_halo_policy,
            "results": rows,
            "passed": all(bool(row["passed"]) for row in rows),
        }
        json_path = ladder_directory / "radius_ladder_report.json"
        json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        plot_path = ladder_directory / "radius_ladder.png"
        if not arguments.no_plot:
            try:
                write_plot(rows, plot_path)
            except ImportError as error:
                print(f"plot skipped: {error}", file=sys.stderr)

        print(f"\nladder directory: {ladder_directory}")
        print(f"summary:          {csv_path}")
        print(f"report:           {json_path}")
        if plot_path.is_file():
            print(f"plot:             {plot_path}")
        return 0 if report["passed"] else 1
    except (OSError, ValueError, LadderError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
