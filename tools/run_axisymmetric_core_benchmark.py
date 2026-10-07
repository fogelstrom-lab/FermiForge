#!/usr/bin/env python3
"""Build, run, archive, and plot an axisymmetric-core 2D benchmark."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from contextlib import ExitStack
from iteration_history_chain import concatenate_history, read_history
from pathlib import Path
from typing import Iterable, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "examples" / "2d_normal_core_benchmark.nml"
DEFAULT_BUILD_DIRECTORY = PROJECT_ROOT / "work" / "normal-core-build"
DEFAULT_RUN_ROOT = PROJECT_ROOT / "runs"


class RunnerError(RuntimeError):
    """Expected build, execution, or plotting failure."""


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run and archive a FermiForge axisymmetric-core 2D benchmark."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--case-name", default="normal-core")
    parser.add_argument("--initialization", default="radial-reference")
    parser.add_argument(
        "--initialization-mode",
        choices=("radial_reference", "restart", "restart_core", "historical_nop", "historical_aop", "historical_dop", "historical_qop",
                 "localized_0plus", "localized_plus0", "localized_0minus", "localized_minus0"),
        help="override the solver initialization mode",
    )
    parser.add_argument(
        "--restart", type=Path, help="restart from a saved 2D field map"
    )
    parser.add_argument("--ranks", type=int, default=10)
    parser.add_argument("--jobs", type=int, default=10)
    parser.add_argument("--cells", type=int, help="override number_of_cells")
    parser.add_argument("--half-width", type=float, help="override half_width")
    parser.add_argument(
        "--mesh-kind",
        choices=("uniform", "multiscale", "smooth", "extended_smooth"),
        help="override mesh_kind",
    )
    parser.add_argument(
        "--fine-region", type=float, help="override fine_region_half_width"
    )
    parser.add_argument(
        "--medium-region", type=float, help="override medium_region_half_width"
    )
    parser.add_argument("--fine-spacing", type=float, help="override fine_spacing")
    parser.add_argument(
        "--medium-spacing", type=float, help="override medium_spacing"
    )
    parser.add_argument(
        "--coarse-spacing", type=float, help="override coarse_spacing"
    )
    parser.add_argument(
        "--active-radius", type=float, help="override circular active_radius"
    )
    parser.add_argument(
        "--outer-radius", type=float, help="override asymptotic_outer_radius"
    )
    parser.add_argument(
        "--fit-inner-radius",
        type=float,
        help="override asymptotic_fit_inner_radius",
    )
    parser.add_argument(
        "--matching-radius",
        type=float,
        help="override asymptotic_matching_radius",
    )
    parser.add_argument(
        "--trajectory-step", type=float, help="override trajectory_maximum_step"
    )
    parser.add_argument(
        "--endpoint-policy",
        choices=(
            "local",
            "free_vortex",
            "radial_reference",
            "radial_asymptotic",
            "state_asymptotic",
        ),
        help="override endpoint_policy",
    )
    parser.add_argument(
        "--max-iterations", type=int, help="override maximum_iterations"
    )
    parser.add_argument("--anderson-cycle", type=int, help="AA updates per cycle; 0 disables cycling")
    parser.add_argument("--simple-cycle", type=int, help="simple updates per cycle; 0 gives AA-only restarts (default 3)")
    parser.add_argument("--simple-mixing", type=float, help="simple damping factor in (0,1] (default 0.1)")
    parser.add_argument("--iteration-method", choices=("anderson", "bb", "polyak"))
    parser.add_argument("--polyak-step-size", type=float, help="Polyak residual multiplier (default 2)")
    parser.add_argument("--polyak-drag", type=float, help="Polyak drag in (0,1]; retained momentum is 1-drag")
    parser.add_argument("--temperature", type=float, help="T/Tc; enables startup Ozaki generation")
    parser.add_argument("--fs1", type=float, help="Landau F1s; overrides legacy feedback_scale")
    parser.add_argument("--save-iteration-diagnostics", action=argparse.BooleanOptionalAction, default=None,
                        help="save signed radial state/residual/update and regional errors; BB also records secants")
    parser.add_argument("--ozaki-cutoff", type=float, help="Ozaki selection cutoff (default 50)")
    parser.add_argument("--bulk-gap-mode", choices=("auto", "ozaki", "legacy", "manual"),
                        help="bulk gap: matching Ozaki quadrature, legacy Matsubara sum, or explicit input")
    parser.add_argument("--spatial-mode", choices=("full_2d", "radial_symmetry"),
                        help="iterate the full plane or only the independent +x radial ray")
    parser.add_argument("--bb-curvature", choices=("positive", "absolute"))
    parser.add_argument("--bb-initial-mixing", type=float)
    parser.add_argument("--bb-minimum-mixing", type=float)
    parser.add_argument("--bb-maximum-mixing", type=float)
    parser.add_argument("--bb-growth-limit", type=float, help="growth damping threshold >1; 0 disables")
    parser.add_argument(
        "--tolerance", type=float, help="override convergence_tolerance"
    )
    parser.add_argument(
        "--anderson-history", type=int, help="override Anderson history limit"
    )
    parser.add_argument(
        "--anderson-progress",
        type=float,
        help="override Anderson progress threshold",
    )
    parser.add_argument(
        "--anderson-pmax", type=float, help="override Anderson maximum mixing"
    )
    parser.add_argument(
        "--checkpoint-interval", type=int, help="checkpoint every N iterations"
    )
    parser.add_argument(
        "--perturbation", type=float, help="initial radial perturbation amplitude"
    )
    parser.add_argument(
        "--perturbation-radius", type=float, help="initial perturbation radius"
    )
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIRECTORY)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--history-parent-run", type=Path,
                        help="archive and concatenate the parent run history; requires its saved restart field")
    parser.add_argument(
        "--formats", nargs="+", choices=("png", "pdf"), default=("png", "pdf")
    )
    return parser.parse_args(argv)


def find_program(name: str, fallbacks: Iterable[Path] = ()) -> Path:
    found = shutil.which(name)
    if found:
        return Path(found).absolute()
    for path in fallbacks:
        if path.is_file() and os.access(path, os.X_OK):
            return path.absolute()
    raise RunnerError(f"required program '{name}' was not found")


def find_plotting_python(environment: dict[str, str]) -> Path:
    """Return a Python interpreter that can import the plotting dependencies."""
    candidates = [
        Path(sys.executable),
        Path("/opt/homebrew/bin/python3"),
        Path("/opt/homebrew/opt/python@3.14/bin/python3.14"),
        Path("/opt/local/bin/python3"),
    ]
    discovered = shutil.which("python3")
    if discovered:
        candidates.append(Path(discovered))

    checked: set[Path] = set()
    failures: list[str] = []
    for candidate in candidates:
        path = candidate.expanduser().resolve()
        if path in checked or not path.is_file() or not os.access(path, os.X_OK):
            continue
        checked.add(path)
        result = subprocess.run(
            [str(path), "-c", "import numpy, matplotlib"],
            env=environment,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if result.returncode == 0:
            return path
        detail = result.stderr.strip().splitlines()
        failures.append(f"{path}: {detail[-1] if detail else 'import failed'}")
    raise RunnerError(
        "no Python interpreter with NumPy and Matplotlib was found; checked "
        + "; ".join(failures)
    )


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", value.strip()).strip("-.")
    return cleaned or "case"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_logged(
    command: Sequence[str],
    log_path: Path,
    environment: dict[str, str],
    append: bool = False,
    echo: bool = False,
    stderr_log_path: Path | None = None,
) -> None:
    mode = "a" if append else "w"
    with ExitStack() as stack:
        log = stack.enter_context(log_path.open(mode, encoding="utf-8"))
        error_output = subprocess.STDOUT
        if stderr_log_path is not None:
            error_output = stack.enter_context(stderr_log_path.open(mode, encoding="utf-8"))
            log.write(f"stderr diagnostics: {stderr_log_path}\n")
            if echo:
                print(f"MPI/runtime diagnostics: {stderr_log_path}", flush=True)
        log.write("command: " + " ".join(command) + "\n")
        log.flush()
        if echo:
            process = subprocess.Popen(
                command,
                cwd=PROJECT_ROOT,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=error_output,
                text=True,
            )
            assert process.stdout is not None
            stack.callback(process.stdout.close)
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end="", flush=True)
            return_code = process.wait()
        else:
            result = subprocess.run(
                command,
                cwd=PROJECT_ROOT,
                env=environment,
                stdout=log,
                stderr=error_output,
                check=False,
            )
            return_code = result.returncode
    if return_code != 0:
        raise RunnerError(
            f"command failed with status {return_code}; see {log_path}"
            + (f" and {stderr_log_path}" if stderr_log_path is not None else "")
        )


def configure_mpi_shutdown(environment: dict[str, str], platform: str = sys.platform) -> None:
    """Local Open MPI mitigation, not warning suppression; preserve user overrides.

    PMIx 5.0.9 can return success when its finalize ACK timer expires.
    A longer bounded allowance tolerates temporary launcher delays, while
    client diagnostics identify timeout versus acknowledgement in mpi.log.
    """
    if platform == "darwin":
        environment.setdefault("PMIX_MCA_pmix_finalize_timeout", "60")
        environment.setdefault("PMIX_MCA_pmix_client_base_verbose", "2")


def replace_namelist_value(text: str, key: str, value: str) -> str:
    # Slashes, commas and exclamation marks inside quoted paths are data,
    # not namelist terminators. Fortran escapes a quote by doubling it.
    pattern = re.compile(
        rf"(?im)^([ \t]*{re.escape(key)}[ \t]*=[ \t]*)"
        r"(?:'(?:[^'\n]|'')*'|\"(?:[^\"\n]|\"\")*\"|[^,!\n/]+)"
    )
    updated, count = pattern.subn(lambda match: match.group(1) + value, text, count=1)
    if count == 0 and key in {
        "anderson_cycle_iterations", "simple_cycle_iterations", "simple_mixing",
        "anderson_history_limit", "anderson_progress_threshold", "anderson_maximum_mixing",
        "temperature", "ozaki_mode", "ozaki_cutoff", "bulk_gap_mode", "fs1", "save_iteration_diagnostics",
        "iteration_method", "spatial_mode", "bb_curvature", "bb_initial_mixing",
        "bb_minimum_mixing", "bb_maximum_mixing", "bb_growth_limit", "polyak_step_size", "polyak_drag"
    }:
        # New optional controls must also work with archived, pre-cycle inputs.
        updated, count = re.subn(r"(?m)^([ \t]*)/([ \t]*(?:!.*)?)$",
                                 lambda m: f"  {key} = {value}\n" + m.group(0),
                                 text, count=1)
    if count != 1:
        raise RunnerError(f"input namelist does not define '{key}'")
    return updated


def build_benchmark(
    build_directory: Path, jobs: int, environment: dict[str, str]
) -> Path:
    cmake = find_program(
        "cmake", (Path("/opt/local/bin/cmake"), Path("/opt/homebrew/bin/cmake"))
    )
    gfortran = find_program("gfortran", (Path("/opt/homebrew/bin/gfortran"),))
    build_directory.mkdir(parents=True, exist_ok=True)
    log = build_directory / "axisymmetric_core_build.log"
    configure = [
        str(cmake),
        "-S",
        str(PROJECT_ROOT),
        "-B",
        str(build_directory),
        f"-DCMAKE_Fortran_COMPILER={gfortran}",
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    build = [
        str(cmake),
        "--build",
        str(build_directory),
        "--target",
        "benchmark_axisymmetric_core_2d",
        "--parallel",
        str(max(1, jobs)),
    ]
    run_logged(configure, log, environment)
    run_logged(build, log, environment, append=True)
    executable = build_directory / "benchmark_axisymmetric_core_2d"
    if not executable.is_file():
        raise RunnerError(f"build did not produce {executable}")
    return executable


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        if arguments.ranks < 1:
            raise RunnerError("MPI rank count must be positive")
        input_path = arguments.input.expanduser().resolve()
        if not input_path.is_file():
            raise RunnerError(f"input file does not exist: {input_path}")

        label = f"{safe_name(arguments.case_name)}-{safe_name(arguments.initialization)}"
        timestamp = dt.datetime.now().strftime("%y%m%d-%H%M%S")
        run_directory = arguments.run_root.expanduser().resolve() / f"{timestamp}-{label}"
        serial = 1
        while run_directory.exists():
            serial += 1
            run_directory = run_directory.with_name(
                f"{timestamp}-{label}-{serial:02d}"
            )
        run_directory.mkdir(parents=True)
        archived_input = run_directory / "input.nml"
        parent_history = None
        if arguments.history_parent_run is not None:
            parent_run = arguments.history_parent_run.expanduser().resolve()
            if arguments.restart is None or arguments.restart.expanduser().resolve() not in (
                    parent_run / "final_fields_2d.dat", parent_run / "checkpoint_fields_2d.dat"):
                raise RunnerError("history parent requires restarting its final or checkpoint field")
            parent_metrics = dict(line.split("=",1) for line in
                                  (parent_run / "metrics.txt").read_text().splitlines() if "=" in line)
            if parent_metrics.get("terminal_status","").strip() not in ("iteration_limit","converged"):
                raise RunnerError("history parent must have completed its iteration segment")
            source_history = parent_run / "combined_iteration_history.dat"
            if not source_history.is_file():
                source_history = parent_run / "iteration_history.dat"
            try:
                read_history(source_history)
                _, local_parent_rows = read_history(parent_run / "iteration_history.dat")
                if len(local_parent_rows) != int(parent_metrics["iterations_completed"]):
                    raise ValueError("parent history does not match completed iteration count")
                if arguments.restart.name == "checkpoint_fields_2d.dat":
                    with arguments.restart.open() as checkpoint_input:
                        checkpoint_header = "".join(next(checkpoint_input, "") for _ in range(3))
                    match = re.search(r"checkpoint_iteration_(\d+)", checkpoint_header)
                    if match is None or int(match[1]) != len(local_parent_rows):
                        raise ValueError("checkpoint does not match the end of the parent history")
            except ValueError as error:
                raise RunnerError(str(error)) from error
            parent_history = run_directory / "parent_iteration_history.dat"
            shutil.copyfile(source_history, parent_history)
        input_text = input_path.read_text(encoding="utf-8")
        overrides = {
            "temperature": arguments.temperature,
            "fs1": arguments.fs1,
            "save_iteration_diagnostics": arguments.save_iteration_diagnostics,
            "ozaki_mode": "generate" if arguments.temperature is not None else None,
            "ozaki_cutoff": arguments.ozaki_cutoff,
            "bulk_gap_mode": arguments.bulk_gap_mode,
            "number_of_cells": arguments.cells,
            "half_width": arguments.half_width,
            "fine_region_half_width": arguments.fine_region,
            "medium_region_half_width": arguments.medium_region,
            "fine_spacing": arguments.fine_spacing,
            "medium_spacing": arguments.medium_spacing,
            "coarse_spacing": arguments.coarse_spacing,
            "active_radius": arguments.active_radius,
            "asymptotic_outer_radius": arguments.outer_radius,
            "asymptotic_fit_inner_radius": arguments.fit_inner_radius,
            "asymptotic_matching_radius": arguments.matching_radius,
            "trajectory_maximum_step": arguments.trajectory_step,
            "maximum_iterations": arguments.max_iterations,
            "anderson_cycle_iterations": arguments.anderson_cycle,
            "iteration_method": arguments.iteration_method,
            "polyak_step_size": arguments.polyak_step_size,
            "polyak_drag": arguments.polyak_drag,
            "spatial_mode": arguments.spatial_mode,
            "bb_curvature": arguments.bb_curvature,
            "bb_initial_mixing": arguments.bb_initial_mixing,
            "bb_minimum_mixing": arguments.bb_minimum_mixing,
            "bb_maximum_mixing": arguments.bb_maximum_mixing,
            "bb_growth_limit": arguments.bb_growth_limit,
            "simple_cycle_iterations": arguments.simple_cycle,
            "simple_mixing": arguments.simple_mixing,
            "convergence_tolerance": arguments.tolerance,
            "anderson_history_limit": arguments.anderson_history,
            "anderson_progress_threshold": arguments.anderson_progress,
            "anderson_maximum_mixing": arguments.anderson_pmax,
            "checkpoint_interval": arguments.checkpoint_interval,
            "initial_perturbation_amplitude": arguments.perturbation,
            "initial_perturbation_radius": arguments.perturbation_radius,
        }
        for key, value in overrides.items():
            if value is not None:
                input_text = replace_namelist_value(
                    input_text, key, (".true." if value else ".false.") if isinstance(value, bool)
                    else repr(value) if isinstance(value, str) else str(value))
        if arguments.endpoint_policy is not None:
            input_text = replace_namelist_value(
                input_text, "endpoint_policy", f"'{arguments.endpoint_policy}'"
            )
        if arguments.mesh_kind is not None:
            input_text = replace_namelist_value(
                input_text, "mesh_kind", f"'{arguments.mesh_kind}'"
            )
        if arguments.initialization_mode is not None:
            input_text = replace_namelist_value(
                input_text,
                "initialization_mode",
                f"'{arguments.initialization_mode}'",
            )
        if arguments.restart is not None:
            restart_path = arguments.restart.expanduser().resolve()
            if not restart_path.is_file():
                raise RunnerError(f"restart field map does not exist: {restart_path}")
            input_text = replace_namelist_value(
                input_text, "initialization_mode", repr(arguments.initialization_mode or 'restart')
            )
            input_text = replace_namelist_value(
                input_text, "restart_field_file", f"'{restart_path}'"
            )
            input_text = replace_namelist_value(
                input_text, "initial_perturbation_amplitude", "0.0"
            )
        archived_input.write_text(input_text, encoding="utf-8")

        environment = os.environ.copy()
        configure_mpi_shutdown(environment)
        environment["PATH"] = ":".join(
            ("/opt/local/bin", "/opt/homebrew/bin", environment.get("PATH", ""))
        )
        matplotlib_cache = PROJECT_ROOT / "work" / "matplotlib"
        matplotlib_cache.mkdir(parents=True, exist_ok=True)
        environment.setdefault("MPLCONFIGDIR", str(matplotlib_cache))

        build_directory = arguments.build_dir.expanduser().resolve()
        executable = build_directory / "benchmark_axisymmetric_core_2d"
        if not arguments.no_build:
            executable = build_benchmark(build_directory, arguments.jobs, environment)
        if not executable.is_file():
            raise RunnerError(f"benchmark executable does not exist: {executable}")

        mpiexec = find_program("mpiexec", (Path("/opt/homebrew/bin/mpiexec"),))
        input_map = run_directory / "input_fields_2d.dat"
        mapped_map = run_directory / "mapped_fields_2d.dat"
        final_map = run_directory / "final_fields_2d.dat"
        metrics = run_directory / "metrics.txt"
        history = run_directory / "iteration_history.dat"
        checkpoint = run_directory / "checkpoint_fields_2d.dat"
        mpi_command = [str(mpiexec)]
        if sys.platform == "darwin":
            mpi_command.extend(
                [
                    "--host",
                    f"localhost:{arguments.ranks}",
                    "--map-by",
                    f"ppr:{arguments.ranks}:node",
                    "--bind-to",
                    "none",
                ]
            )
        mpi_command.extend(
            [
                "-n",
                str(arguments.ranks),
                str(executable),
                str(archived_input),
                str(input_map),
                str(mapped_map),
                str(metrics),
                str(final_map),
                str(history),
                str(checkpoint),
            ]
        )
        launcher_error = None
        try:
            run_logged(mpi_command, run_directory / "run.log", environment, echo=True,
                       stderr_log_path=run_directory / "mpi.log")
        except RunnerError as error:
            launcher_error = str(error)
            print(f"warning: {error}; attempting saved-output postprocessing", file=sys.stderr)
        for output in (input_map, mapped_map, final_map, metrics, history):
            if not output.is_file():
                raise RunnerError(f"benchmark did not produce {output}")

        # Check what Fortran actually read, not merely the requested CLI value.
        expected_match = re.search(r"(?im)^\s*maximum_iterations\s*=\s*(\d+)", input_text)
        metric_values = dict(
            line.split('=', 1) for line in metrics.read_text().splitlines() if '=' in line
        )
        if expected_match:
            expected_iterations = int(expected_match.group(1))
            if int(metric_values.get('maximum_iterations', '-1')) != expected_iterations:
                raise RunnerError('solver iteration limit differs from the prepared input; inspect input.nml')
            if expected_iterations > 0 and metric_values.get('terminal_status', '').strip() == 'single_map':
                raise RunnerError('solver performed a single map despite requested iterations')

        combined_history = None
        if parent_history is not None:
            combined_history = run_directory / "combined_iteration_history.dat"
            try:
                concatenate_history(parent_history, history, combined_history)
            except ValueError as error:
                raise RunnerError(str(error)) from error
        plots: list[str] = []
        if not arguments.no_plot:
            python = find_plotting_python(environment)
            plot_script = PROJECT_ROOT / "tools" / "plot_2d_fields.py"
            plot_directory = run_directory / "plots"
            for role, field_map in (("input", input_map), ("final", final_map)):
                command = [
                    str(python),
                    str(plot_script),
                    str(field_map),
                    "--output-dir",
                    str(plot_directory),
                    "--prefix",
                    f"{timestamp[:6]}_{label}_{role}",
                    "--formats",
                    *arguments.formats,
                ]
                run_logged(
                    command,
                    run_directory / "plot.log",
                    environment,
                    append=role != "input",
                )
            convergence_plot = plot_directory / f"{timestamp[:6]}_{label}_convergence.png"
            for basis, filename in (('cartesian', 'axis_profiles.png'),
                                    ('harmonic', 'harmonic_axis_profiles.png')):
                run_logged([
                    str(python), str(PROJECT_ROOT / 'tools/plot_double_core_axis_profiles.py'),
                    str(final_map), '--reference-field', str(input_map),
                    '--reference-label', 'starting state', '--label', 'saved final state',
                    '--title', 'Vortex symmetry-axis order-parameter profiles',
                    '--basis', basis, '--output', str(plot_directory / filename),
                ], run_directory / 'plot.log', environment, append=True)
            radial_reference = Path(str(metrics) + '.radial_reference.dat')
            if radial_reference.is_file():
                run_logged([
                    str(python), str(PROJECT_ROOT / 'tools/recover_axisymmetric_plots.py'),
                    str(run_directory), '--radial-reference', str(radial_reference),
                ], run_directory / 'plot.log', environment, append=True)
            run_logged(
                [
                    str(python),
                    str(PROJECT_ROOT / "tools" / "plot_2d_iteration_history.py"),
                    str(history),
                    "--output",
                    str(convergence_plot),
                ],
                run_directory / "plot.log",
                environment,
                append=True,
            )
            run_logged([
                str(python), str(PROJECT_ROOT / 'tools/plot_residual_locations.py'),
                str(run_directory), '--output', str(plot_directory),
            ], run_directory / 'plot.log', environment, append=True)
            plots = [str(path) for path in sorted(plot_directory.glob("*"))]
            if combined_history is not None:
                run_logged([
                    str(python), str(PROJECT_ROOT / "tools/plot_2d_iteration_history.py"),
                    str(combined_history), "--output", str(plot_directory / "combined_convergence.png"),
                ], run_directory / "plot.log", environment, append=True)
                plots.append(str(plot_directory / "combined_convergence.png"))

        manifest = {
            "history_parent_run": str(arguments.history_parent_run.resolve()) if parent_history else None,
            "parent_history_sha256": sha256(parent_history) if parent_history else None,
            "combined_history": str(combined_history) if combined_history else None,
            "iterator_history_restarted": parent_history is not None,
            "launcher_error": launcher_error,
            "kind": "fermiforge_axisymmetric_core_2d_benchmark",
            "created_at": dt.datetime.now().astimezone().isoformat(),
            "case_name": safe_name(arguments.case_name),
            "initialization": safe_name(arguments.initialization),
            "mpi_ranks": arguments.ranks,
            "mpi_shutdown_environment": {
                key: environment[key] for key in (
                    "PMIX_MCA_pmix_finalize_timeout", "PMIX_MCA_pmix_client_base_verbose")
                if key in environment
            },
            "overrides": {
                key: value
                for key, value in {
                    "number_of_cells": arguments.cells,
                    "half_width": arguments.half_width,
                    "mesh_kind": arguments.mesh_kind,
                    "fine_region_half_width": arguments.fine_region,
                    "medium_region_half_width": arguments.medium_region,
                    "fine_spacing": arguments.fine_spacing,
                    "medium_spacing": arguments.medium_spacing,
                    "coarse_spacing": arguments.coarse_spacing,
                    "active_radius": arguments.active_radius,
                    "asymptotic_outer_radius": arguments.outer_radius,
                    "asymptotic_fit_inner_radius": arguments.fit_inner_radius,
                    "asymptotic_matching_radius": arguments.matching_radius,
                    "trajectory_maximum_step": arguments.trajectory_step,
                    "endpoint_policy": arguments.endpoint_policy,
                    "maximum_iterations": arguments.max_iterations,
                    "anderson_cycle_iterations": arguments.anderson_cycle,
                    "iteration_method": arguments.iteration_method,
                    "polyak_step_size": arguments.polyak_step_size,
                    "polyak_drag": arguments.polyak_drag,
                    "bb_curvature": arguments.bb_curvature,
                    "bb_initial_mixing": arguments.bb_initial_mixing,
                    "bb_minimum_mixing": arguments.bb_minimum_mixing,
                    "bb_maximum_mixing": arguments.bb_maximum_mixing,
                    "bb_growth_limit": arguments.bb_growth_limit,
                    "simple_cycle_iterations": arguments.simple_cycle,
                    "simple_mixing": arguments.simple_mixing,
                    "convergence_tolerance": arguments.tolerance,
                    "anderson_history_limit": arguments.anderson_history,
                    "anderson_progress_threshold": arguments.anderson_progress,
                    "anderson_maximum_mixing": arguments.anderson_pmax,
                    "checkpoint_interval": arguments.checkpoint_interval,
                    "initial_perturbation_amplitude": arguments.perturbation,
                    "initial_perturbation_radius": arguments.perturbation_radius,
                    "initialization_mode": arguments.initialization_mode,
                    "restart_field_file": (
                        str(arguments.restart.expanduser().resolve())
                        if arguments.restart is not None
                        else None
                    ),
                }.items()
                if value is not None
            },
            "input_source": str(input_path),
            "input_sha256": sha256(archived_input),
            "executable": str(executable),
            "mpi_log": str(run_directory / "mpi.log"),
            "executable_sha256": sha256(executable),
            "input_field_map_sha256": sha256(input_map),
            "mapped_field_map_sha256": sha256(mapped_map),
            "final_field_map_sha256": sha256(final_map),
            "metrics_sha256": sha256(metrics),
            "iteration_history_sha256": sha256(history),
            "checkpoint_field_map_sha256": (
                sha256(checkpoint) if checkpoint.is_file() else None
            ),
            "plots": plots,
        }
        (run_directory / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
    except (OSError, RunnerError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    print(f"run directory: {run_directory}")
    print(f"metrics:       {metrics}")
    print(f"history:       {history}")
    if plots:
        print(f"plots:         {run_directory / 'plots'}")
    if launcher_error:
        print(f"Saved plots recovered, but MPI launch failed: {launcher_error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
