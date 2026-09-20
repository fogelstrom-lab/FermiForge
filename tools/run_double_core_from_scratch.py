#!/usr/bin/env python3
"""Build and run the first analytic-from-scratch double-core iteration."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

from run_double_core_radius_ladder import (
    LadderError,
    find_program,
    fortran_string,
    metric_float,
    metric_int,
    mpi_prefix,
    read_metrics,
    replace_namelist_value,
    sparse_checkpoint_point_count,
    stream_command,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "examples" / "2d_double_core_from_scratch_smoke.nml"
DEFAULT_BUILD = PROJECT_ROOT / "work" / "double-core-scratch-build"


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Construct a regularized London double-core seed on a new 2D "
            "mesh and begin MPI/Anderson self-consistency iteration."
        )
    )
    parser.add_argument(
        "--iterations",
        type=int,
        help=(
            "maximum Anderson updates (default: maximum_iterations from "
            "the selected input)"
        ),
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=1,
        help=(
            "write a restart checkpoint after this many completed updates; "
            "zero disables checkpoints (default: 1)"
        ),
    )
    parser.add_argument("--ranks", type=int, default=10)
    parser.add_argument("--jobs", type=int, default=10)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD)
    parser.add_argument("--run-root", type=Path, default=PROJECT_ROOT / "runs")
    parser.add_argument(
        "--run-name",
        default="double-core-from-scratch",
        help="suffix used for the timestamped run directory",
    )
    parser.add_argument(
        "--sparse-restart",
        type=Path,
        help="overlay a previous scratch-run sparse checkpoint on the analytic seed",
    )
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args(argv)


def namelist_integer(text: str, key: str) -> int:
    match = re.search(
        rf"(?im)^\s*{re.escape(key)}\s*=\s*([+-]?\d+)\s*(?:,|$)", text
    )
    if match is None:
        raise LadderError(f"input namelist does not define integer '{key}'")
    return int(match.group(1))


def optional_namelist_integer(text: str, key: str, default: int = 0) -> int:
    match = re.search(
        rf"(?im)^\s*{re.escape(key)}\s*=\s*([+-]?\d+)\s*(?:,|$)", text
    )
    return default if match is None else int(match.group(1))


def namelist_string(text: str, key: str) -> str:
    match = re.search(
        rf"(?im)^\s*{re.escape(key)}\s*=\s*['\"]([^'\"]+)['\"]", text
    )
    if match is None:
        raise LadderError(f"input namelist does not define string '{key}'")
    return match.group(1).strip()


def set_namelist_value(text: str, key: str, value: str) -> str:
    """Replace a value, or add it before the terminating namelist slash."""
    try:
        return replace_namelist_value(text, key, value)
    except LadderError:
        terminators = list(re.finditer(r"(?m)^\s*/\s*$", text))
        if len(terminators) != 1:
            raise LadderError(
                f"cannot add '{key}' to an input without one namelist terminator"
            )
        position = terminators[0].start()
        return text[:position] + f"  {key} = {value}\n" + text[position:]


def optional_metric_float(metrics: dict[str, str], key: str) -> float | None:
    value = metrics.get(key)
    if value is None:
        return None
    try:
        return float(value)
    except ValueError as error:
        raise LadderError(f"metrics contain nonnumeric '{key}'") from error


def optional_metric_int(metrics: dict[str, str], key: str, default: int = 0) -> int:
    value = metrics.get(key)
    if value is None:
        return default
    try:
        return int(value)
    except ValueError as error:
        raise LadderError(f"metrics contain nonnumeric '{key}'") from error


def build_solver(build_directory: Path, jobs: int) -> Path:
    cmake = find_program(
        "cmake", (Path("/opt/local/bin/cmake"), Path("/opt/homebrew/bin/cmake"))
    )
    compiler = find_program("gfortran", (Path("/opt/homebrew/bin/gfortran"),))
    build_directory.mkdir(parents=True, exist_ok=True)
    stream_command(
        [
            str(cmake),
            "-S",
            str(PROJECT_ROOT),
            "-B",
            str(build_directory),
            f"-DCMAKE_Fortran_COMPILER={compiler}",
            "-DCMAKE_BUILD_TYPE=Release",
        ],
        build_directory / "configure.log",
    )
    stream_command(
        [
            str(cmake),
            "--build",
            str(build_directory),
            "--parallel",
            str(jobs),
            "--target",
            "fermiforge_double_core_2d",
        ],
        build_directory / "build.log",
    )
    executable = build_directory / "fermiforge_double_core_2d"
    if not executable.is_file():
        raise LadderError(f"build did not create {executable}")
    return executable


def plotting_python() -> Path | None:
    candidates = (
        sys.executable,
        shutil.which("python3"),
        "/opt/homebrew/bin/python3",
        "/usr/local/bin/python3",
        "/opt/local/bin/python3",
    )
    seen: set[Path] = set()
    environment = os.environ.copy()
    environment.setdefault(
        "MPLCONFIGDIR",
        str(Path(tempfile.gettempdir()) / "fermiforge-matplotlib-cache"),
    )
    for candidate_name in candidates:
        if not candidate_name:
            continue
        candidate = Path(candidate_name).expanduser().resolve()
        if candidate in seen or not candidate.is_file():
            continue
        seen.add(candidate)
        check = subprocess.run(
            [str(candidate), "-c", "import matplotlib, numpy"],
            env=environment,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if check.returncode == 0:
            return candidate
    return None


def make_plots(
    run_directory: Path,
    completed_iterations: int,
    initial_label: str = "regularized London seed",
) -> list[Path]:
    python = plotting_python()
    if python is None:
        raise LadderError("no Python installation with Matplotlib and NumPy was found")
    environment = os.environ.copy()
    environment.setdefault(
        "MPLCONFIGDIR", str(run_directory / "matplotlib-cache")
    )
    plot_directory = run_directory / "plots"
    commands = [
        [
            str(python),
            str(PROJECT_ROOT / "tools" / "plot_2d_fields.py"),
            str(run_directory / "initial_fields_2d.dat"),
            "--output-dir",
            str(plot_directory),
            "--prefix",
            "initial",
            "--formats",
            "png",
        ],
        [
            str(python),
            str(PROJECT_ROOT / "tools" / "plot_2d_fields.py"),
            str(run_directory / "final_fields_2d.dat"),
            "--output-dir",
            str(plot_directory),
            "--prefix",
            "final",
            "--formats",
            "png",
        ],
        [
            str(python),
            str(PROJECT_ROOT / "tools" / "plot_2d_iteration_history.py"),
            str(run_directory / "iteration_history.dat"),
            "--output",
            str(plot_directory / "convergence.png"),
        ],
        [
            str(python),
            str(PROJECT_ROOT / "tools" / "plot_double_core_axis_profiles.py"),
            str(run_directory / "final_fields_2d.dat"),
            "--reference-field",
            str(run_directory / "initial_fields_2d.dat"),
            "--reference-label",
            initial_label,
            "--label",
            f"state after {completed_iterations} updates",
            "--output",
            str(plot_directory / "double_core_axis_profiles.png"),
        ],
        [
            str(python),
            str(
                PROJECT_ROOT
                / "tools"
                / "plot_double_core_harmonic_asymptotics.py"
            ),
            str(run_directory / "final_fields_2d.dat"),
            "--input",
            str(run_directory / "input.nml"),
            "--metrics",
            str(run_directory / "metrics.txt"),
            "--output-dir",
            str(plot_directory),
        ],
        [
            str(python),
            str(PROJECT_ROOT / "tools" / "plot_double_core_structure_history.py"),
            str(run_directory / "iteration_history.dat"),
            "--output",
            str(plot_directory / "double_core_structure_history.png"),
        ],
    ]
    probe_file = run_directory / "asymptotic_probe_states.dat"
    if probe_file.is_file():
        probe_command = [
            str(python),
            str(
                PROJECT_ROOT
                / "tools"
                / "plot_double_core_asymptotic_probes.py"
            ),
            str(probe_file),
            "--input",
            str(run_directory / "input.nml"),
            "--metrics",
            str(run_directory / "metrics.txt"),
            "--output-dir",
            str(plot_directory),
        ]
        shadow_file = run_directory / "asymptotic_shadow_states.dat"
        shadow_history = run_directory / "asymptotic_probe_history.dat"
        if shadow_file.is_file():
            probe_command.extend(["--shadow-file", str(shadow_file)])
        if shadow_history.is_file():
            probe_command.extend(["--shadow-history", str(shadow_history)])
        commands.append(probe_command)
    for index, command in enumerate(commands, start=1):
        with (run_directory / f"plot_{index}.log").open(
            "w", encoding="utf-8"
        ) as log:
            result = subprocess.run(
                command,
                cwd=PROJECT_ROOT,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        if result.returncode != 0:
            raise LadderError(
                f"plot command {index} failed; see {run_directory / f'plot_{index}.log'}"
            )
    return sorted(plot_directory.glob("*.png"))


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        if arguments.ranks < 1 or arguments.jobs < 1:
            raise LadderError("ranks and jobs must be positive")
        if arguments.checkpoint_every < 0:
            raise LadderError("checkpoint-every cannot be negative")
        input_path = arguments.input.expanduser().resolve()
        template = input_path.read_text(encoding="utf-8")
        if re.search(r"(?im)^\s*asymptotic_ray_probe_count\s*=", template) is None:
            update_region = namelist_string(template, "update_region")
            template = set_namelist_value(
                template,
                "asymptotic_ray_probe_count",
                "8" if update_region in ("centered_disk", "centered_ellipse") else "0",
            )
        requested_iterations = (
            arguments.iterations
            if arguments.iterations is not None
            else namelist_integer(template, "maximum_iterations")
        )
        if requested_iterations < 1:
            raise LadderError("iterations must be positive")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", arguments.run_name):
            raise LadderError("run-name contains unsafe filename characters")
        timestamp = dt.datetime.now().strftime("%y%m%d-%H%M%S")
        run_directory = arguments.run_root.expanduser().resolve() / (
            f"{timestamp}-{arguments.run_name}"
        )
        run_directory.mkdir(parents=True)

        build_directory = arguments.build_dir.expanduser().resolve()
        executable = build_directory / "fermiforge_double_core_2d"
        if not arguments.no_build:
            executable = build_solver(build_directory, arguments.jobs)
        elif not executable.is_file():
            raise LadderError(f"solver executable was not found: {executable}")

        paths = {
            "initial_field_file": run_directory / "initial_fields_2d.dat",
            "final_field_file": run_directory / "final_fields_2d.dat",
            "sparse_final_file": run_directory / "sparse_final_state.dat",
            "history_file": run_directory / "iteration_history.dat",
            "metrics_file": run_directory / "metrics.txt",
            "sample_file": run_directory / "sampled_residuals.dat",
        }
        ray_probe_count = namelist_integer(template, "asymptotic_ray_probe_count")
        shadow_probe_iterations = optional_namelist_integer(
            template, "asymptotic_probe_iterations"
        )
        if ray_probe_count > 0:
            paths["asymptotic_probe_file"] = (
                run_directory / "asymptotic_probe_states.dat"
            )
        if shadow_probe_iterations > 0:
            paths["asymptotic_probe_history_file"] = (
                run_directory / "asymptotic_probe_history.dat"
            )
            paths["asymptotic_shadow_file"] = (
                run_directory / "asymptotic_shadow_states.dat"
            )
        checkpoint_path = run_directory / "sparse_checkpoint_state.dat"
        input_text = replace_namelist_value(
            template, "maximum_iterations", str(requested_iterations)
        )
        input_text = replace_namelist_value(
            input_text, "checkpoint_interval", str(arguments.checkpoint_every)
        )
        input_text = replace_namelist_value(
            input_text,
            "checkpoint_file",
            fortran_string(checkpoint_path) if arguments.checkpoint_every else "''",
        )
        for key, path in paths.items():
            input_text = set_namelist_value(
                input_text, key, fortran_string(path)
            )
        archived_input = run_directory / "input.nml"
        expected_restart_points = 0
        if arguments.sparse_restart is not None:
            restart_path = arguments.sparse_restart.expanduser().resolve()
            if not restart_path.is_file():
                raise LadderError(f"sparse restart was not found: {restart_path}")
            expected_restart_points = sparse_checkpoint_point_count(restart_path)
            input_text = replace_namelist_value(
                input_text, "sparse_restart_file", fortran_string(restart_path)
            )
        archived_input.write_text(input_text, encoding="utf-8")

        print(f"run directory:         {run_directory}")
        print(f"input:                 {input_path}")
        print(f"MPI ranks:             {arguments.ranks}")
        print(f"requested iterations:  {requested_iterations}")
        print(f"checkpoint interval:   {arguments.checkpoint_every}")
        print()
        launcher_warning: str | None = None
        try:
            stream_command(
                [*mpi_prefix(arguments.ranks), str(executable), str(archived_input)],
                run_directory / "run.log",
            )
        except LadderError as error:
            # Open MPI/PRRTE can occasionally report a nonzero launcher exit
            # after a long macOS sleep/resume cycle even though every numerical
            # product and the restart checkpoint have already been flushed.
            # A newly-created run directory cannot contain stale products, so
            # allow the normal validation below to decide whether this is a
            # recoverable launcher-only failure.  Partial runs remain fatal.
            required_products = list(paths.values())
            if arguments.checkpoint_every:
                required_products.append(checkpoint_path)
            if not all(path.is_file() for path in required_products):
                raise
            launcher_warning = str(error)
            print(
                "warning: MPI launcher returned nonzero after all expected "
                "products were written; validating and recovering the run",
                file=sys.stderr,
            )
        metrics = read_metrics(paths["metrics_file"])
        result: dict[str, object] = {
            "terminal_status": metrics.get("terminal_status", "missing"),
            "initialization_mode": metrics.get("initialization_mode", "missing"),
            "mesh_x_points": metric_int(metrics, "mesh_x_points"),
            "mesh_y_points": metric_int(metrics, "mesh_y_points"),
            "updated_points": metric_int(metrics, "updated_points"),
            "dependent_halo_points": metric_int(
                metrics, "dependent_asymptotic_halo_points"
            ),
            "iterations_completed": metric_int(metrics, "iterations_completed"),
            "sparse_restart_applied_points": metric_int(
                metrics, "sparse_restart_applied_points"
            ),
            "sample_relative_l2": metric_float(metrics, "sample_relative_l2"),
            "update_region": metrics.get("update_region", "missing"),
            "update_radius": metric_float(metrics, "update_radius"),
            "update_radius_x": optional_metric_float(metrics, "update_radius_x"),
            "update_radius_y": optional_metric_float(metrics, "update_radius_y"),
            "normalization_error": metric_float(
                metrics, "maximum_normalization_error"
            ),
            "asymptotic_probe_terminal_status": metrics.get(
                "asymptotic_probe_terminal_status", "disabled"
            ),
            "asymptotic_probe_iterations_completed": optional_metric_int(
                metrics, "asymptotic_probe_iterations_completed"
            ),
            "asymptotic_probe_convergence_tolerance": optional_metric_float(
                metrics, "asymptotic_probe_convergence_tolerance"
            ),
            "asymptotic_probe_terminal_rms": optional_metric_float(
                metrics, "asymptotic_probe_terminal_rms"
            ),
            "asymptotic_probe_terminal_maximum": optional_metric_float(
                metrics, "asymptotic_probe_terminal_maximum"
            ),
            "asymptotic_probe_terminal_relative_l2": optional_metric_float(
                metrics, "asymptotic_probe_terminal_relative_l2"
            ),
            "asymptotic_probe_maximum_normalization_error": (
                optional_metric_float(
                    metrics, "asymptotic_probe_maximum_normalization_error"
                )
            ),
            "historical_seed_kind": metrics.get("historical_seed_kind"),
            "historical_seed_core_length_scale": optional_metric_float(
                metrics, "historical_seed_core_length_scale"
            ),
            "seed_center_pair_amplitude": optional_metric_float(
                metrics, "seed_center_pair_amplitude"
            ),
            "seed_positive_half_core_pair_amplitude": optional_metric_float(
                metrics, "seed_positive_half_core_pair_amplitude"
            ),
            "seed_negative_half_core_pair_amplitude": optional_metric_float(
                metrics, "seed_negative_half_core_pair_amplitude"
            ),
            "measured_half_core_separation": metric_float(
                metrics, "measured_half_core_separation"
            ),
            "measured_center_pair_amplitude": metric_float(
                metrics, "measured_center_pair_amplitude"
            ),
        }
        if shadow_probe_iterations == 0:
            asymptotic_validation_status = "not_requested"
        elif result["asymptotic_probe_terminal_status"] != "converged":
            asymptotic_validation_status = "shadow_not_converged"
        elif result["terminal_status"] != "converged":
            asymptotic_validation_status = "provisional_production_not_converged"
        else:
            asymptotic_validation_status = "passed"
        result["asymptotic_validation_status"] = asymptotic_validation_status
        result["asymptotic_validation_passed"] = (
            asymptotic_validation_status == "passed"
        )
        iteration_count_is_valid = (
            result["terminal_status"] == "converged"
            and 1 <= int(result["iterations_completed"]) <= requested_iterations
        ) or (
            result["terminal_status"] == "iteration_limit"
            and int(result["iterations_completed"]) == requested_iterations
        )
        calculation_passed = (
            result["initialization_mode"] in (
                "regularized_london",
                "historical_nop",
                "historical_aop",
                "historical_dop",
            )
            and result["terminal_status"] in ("iteration_limit", "converged")
            and iteration_count_is_valid
            and int(result["sparse_restart_applied_points"])
            == expected_restart_points
            and float(result["normalization_error"]) <= 1.0e-10
            and all(path.is_file() for path in paths.values())
            and (
                arguments.checkpoint_every == 0 or checkpoint_path.is_file()
            )
            and (
                shadow_probe_iterations == 0
                or (
                    result["asymptotic_probe_terminal_status"]
                    in ("iteration_limit", "converged")
                    and int(result["asymptotic_probe_iterations_completed"])
                    >= 1
                    and result["asymptotic_probe_maximum_normalization_error"]
                    is not None
                    and float(
                        result[
                            "asymptotic_probe_maximum_normalization_error"
                        ]
                    )
                    <= 1.0e-10
                )
            )
        )
        plot_paths: list[Path] = []
        if not arguments.no_plot:
            plot_paths = make_plots(
                run_directory,
                int(result["iterations_completed"]),
                (
                    "inherited checkpoint state"
                    if expected_restart_points > 0
                    else "regularized London seed"
                ),
            )
        report = {
            "kind": "fermiforge_double_core_from_scratch",
            "created_at": dt.datetime.now().astimezone().isoformat(),
            "input_template": str(input_path),
            "ranks": arguments.ranks,
            "requested_iterations": requested_iterations,
            "checkpoint_every": arguments.checkpoint_every,
            "checkpoint": (
                str(checkpoint_path) if arguments.checkpoint_every else None
            ),
            "sparse_restart": (
                str(arguments.sparse_restart.expanduser().resolve())
                if arguments.sparse_restart is not None
                else None
            ),
            "expected_restart_points": expected_restart_points,
            "launcher_warning": launcher_warning,
            "result": result,
            "asymptotic_validation_passed": result[
                "asymptotic_validation_passed"
            ],
            "plots": [str(path) for path in plot_paths],
            "passed": calculation_passed and (
                arguments.no_plot
                or len(plot_paths)
                == (
                    16
                    + (3 if paths.get("asymptotic_probe_file") else 0)
                    + (4 if paths.get("asymptotic_shadow_file") else 0)
                )
            ),
        }
        report_path = run_directory / "scratch_run_report.json"
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

        print(f"\nrun directory: {run_directory}")
        print(f"report:        {report_path}")
        print(f"passed:        {report['passed']}")
        return 0 if report["passed"] else 1
    except (OSError, ValueError, LadderError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
