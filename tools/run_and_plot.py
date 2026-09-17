#!/usr/bin/env python3
"""Build, run, archive, and plot the cylindrically symmetric 3He solvers.

The legacy fixed-form and transitional ``new_src`` solvers are built in
isolated work directories, so a test run never overwrites reference output or
an interactive build.  Each MPI run gets a timestamped archive.  Plotting uses
NumPy and Matplotlib's noninteractive Agg backend and therefore needs neither
Qt nor a display server.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LEGACY_SOURCE = PROJECT_ROOT / "incoming" / "legacy-f77"
NEW_SOURCE = PROJECT_ROOT / "new_src"
SHARED_SOURCE = PROJECT_ROOT / "src"
DEFAULT_RUN_ROOT = PROJECT_ROOT / "runs"
DEFAULT_LEGACY_BUILD_DIR = PROJECT_ROOT / "work" / "legacy-runner-build"
DEFAULT_NEW_BUILD_DIR = PROJECT_ROOT / "work" / "new-src-runner" / "new_src"

LEGACY_BUILD_FILES = (
    "Makefile",
    "qcv.f",
    "iter.f",
    "gap.f",
    "riccati.f",
    "linpacks.f",
    "trajectories.f",
    "init_calc.f",
    "mpicalls.f",
    "interpol.f",
    "qcv.dat",
    "iter.dat",
    "gap.dat",
)

NEW_BUILD_FILES = (
    "Makefile",
    "global_dec.f90",
    "bulkgap.f90",
    "init_calc.f90",
    "mpicalls.f90",
    "interpol.f90",
    "riccati.f90",
    "getnewses.f90",
    "packses.f90",
    "iter_AA.f90",
    "iter.f90",
    "iter_NN.f90",
    "Iter_BR.f90",
    "iter_BB.f90",
    "qcv.f90",
)

SHARED_BUILD_FILES = (
    "he3_kinds.f90",
    "legacy_anderson_mixing.f90",
)

CARTESIAN_COMPONENTS = (
    "Axx", "Axy", "Axz",
    "Ayx", "Ayy", "Ayz",
    "Azx", "Azy", "Azz",
)

HARMONIC_COMPONENTS = (
    "A++", "A+0", "A+-",
    "A0+", "A00", "A0-",
    "A-+", "A-0", "A--",
)


class RunnerError(RuntimeError):
    """An expected build, run, data, or plotting failure."""


@dataclass(frozen=True)
class ConvergenceRecord:
    """One recorded nonlinear residual evaluation."""

    global_iteration: int
    phase_iteration: int
    average_error: float
    maximum_error: float
    engine: str
    stage: int


def parse_arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the cylindrical Fortran/MPI vortex solver in an isolated "
            "directory and plot its output, or plot an existing run."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        help=(
            "solver input file (default: the selected solver's qcv.inp)"
        ),
    )
    parser.add_argument(
        "--solver",
        choices=("legacy", "new-src"),
        default="legacy",
        help=(
            "solver source to build and input syntax to use "
            "(default: legacy)"
        ),
    )
    parser.add_argument(
        "--ranks",
        type=int,
        default=10,
        help="number of MPI processes (default: 10)",
    )
    parser.add_argument(
        "--case-name",
        default="cylindrical",
        help="short name used in the new run-directory name",
    )
    parser.add_argument(
        "--run-root",
        type=Path,
        default=DEFAULT_RUN_ROOT,
        help="directory under which unique run directories are created",
    )
    parser.add_argument(
        "--executable",
        type=Path,
        help="use this solver executable instead of building the selected source",
    )
    parser.add_argument(
        "--plot-only",
        type=Path,
        metavar="RUN_DIR",
        help="plot files in an existing run or reference-data directory",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        help="plot output directory (default: RUN_DIR/plots)",
    )
    parser.add_argument(
        "--plot-backend",
        choices=("matplotlib", "gnuplot"),
        default="matplotlib",
        help="plotting implementation (default: matplotlib)",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="run the solver without producing plots",
    )
    return parser.parse_args(argv)


def find_program(name: str, fallbacks: Iterable[Path] = ()) -> Path:
    found = shutil.which(name)
    if found:
        # MPI compiler wrappers can dispatch according to argv[0]. Preserve a
        # wrapper symlink such as mpifort instead of resolving it to the generic
        # Open MPI backend executable.
        return Path(found).absolute()
    for candidate in fallbacks:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.absolute()
    choices = ", ".join(str(path) for path in fallbacks)
    suffix = f"; also checked {choices}" if choices else ""
    raise RunnerError(f"required program '{name}' was not found{suffix}")


def copy_legacy_build_inputs(build_dir: Path) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    for name in LEGACY_BUILD_FILES:
        source = LEGACY_SOURCE / name
        if not source.is_file():
            raise RunnerError(f"missing legacy build input: {source}")
        shutil.copy2(source, build_dir / name)
    for source in LEGACY_SOURCE.glob("mpif*.h"):
        shutil.copy2(source, build_dir / source.name)


def build_legacy_solver(build_dir: Path, parallel_jobs: int) -> Tuple[Path, Path]:
    copy_legacy_build_inputs(build_dir)
    make = find_program("make", (Path("/usr/bin/make"),))
    gfortran = find_program("gfortran", (Path("/opt/homebrew/bin/gfortran"),))
    mpifort = find_program("mpifort", (Path("/opt/homebrew/bin/mpifort"),))
    log_path = build_dir / "build.log"
    command = [
        str(make),
        "-B",
        "-j",
        str(max(1, parallel_jobs)),
        "qcv",
        f"FC={gfortran}",
        f"MPIcomp={mpifort}",
    ]
    build_environment = os.environ.copy()
    build_environment["PATH"] = (
        str(gfortran.parent)
        + os.pathsep
        + str(mpifort.parent)
        + os.pathsep
        + build_environment.get("PATH", "")
    )
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command,
            cwd=build_dir,
            env=build_environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        raise RunnerError(f"legacy build failed; see {log_path}")
    executable = build_dir / "qcv"
    if not executable.is_file():
        raise RunnerError(f"build completed without producing {executable}")
    return executable, log_path


def copy_new_build_inputs(build_dir: Path) -> None:
    """Copy the active free-form solver into an isolated build tree."""

    build_dir.mkdir(parents=True, exist_ok=True)
    shared_build_dir = build_dir.parent / "src"
    shared_build_dir.mkdir(parents=True, exist_ok=True)
    for name in NEW_BUILD_FILES:
        source = NEW_SOURCE / name
        if not source.is_file():
            raise RunnerError(f"missing new_src build input: {source}")
        shutil.copy2(source, build_dir / name)
    for name in SHARED_BUILD_FILES:
        source = SHARED_SOURCE / name
        if not source.is_file():
            raise RunnerError(f"missing shared build input: {source}")
        shutil.copy2(source, shared_build_dir / name)


def build_new_solver(build_dir: Path, parallel_jobs: int) -> Tuple[Path, Path]:
    copy_new_build_inputs(build_dir)
    make = find_program("make", (Path("/usr/bin/make"),))
    mpifort = find_program("mpifort", (Path("/opt/homebrew/bin/mpifort"),))
    log_path = build_dir / "build.log"
    command = [
        str(make),
        "-B",
        "-j",
        str(max(1, parallel_jobs)),
        "qcv",
        f"MPIFC={mpifort}",
    ]
    build_environment = os.environ.copy()
    build_environment["PATH"] = (
        str(mpifort.parent)
        + os.pathsep
        + build_environment.get("PATH", "")
    )
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command,
            cwd=build_dir,
            env=build_environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        raise RunnerError(f"new_src build failed; see {log_path}")
    executable = build_dir / "qcv"
    if not executable.is_file():
        raise RunnerError(f"build completed without producing {executable}")
    return executable, log_path


def first_numeric_token(line: str) -> float:
    token = line.split(":", 1)[0].strip().split()[0]
    return float(token.replace("D", "E").replace("d", "e"))


def read_legacy_input_summary(input_path: Path) -> Dict[str, object]:
    lines = input_path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 11:
        raise RunnerError(f"legacy input has fewer than 11 records: {input_path}")
    return {
        "temperature_over_tc": first_numeric_token(lines[0]),
        "vorticity": first_numeric_token(lines[1]),
        "cylindrical": int(first_numeric_token(lines[2])),
        "fermi_liquid_f1s": first_numeric_token(lines[3]),
        "radius_legacy_units": first_numeric_token(lines[4]),
        "azimuthal_directions": int(first_numeric_token(lines[5])),
        "requested_error": first_numeric_token(lines[6]),
        "initial_state": int(first_numeric_token(lines[9])),
        "maximum_iterations": int(first_numeric_token(lines[10])),
    }


ITERATION_ENGINES = {
    0: "NN",
    1: "BR",
    2: "BB",
    3: "AA",
}


def read_new_input_summary(input_path: Path) -> Dict[str, object]:
    lines = input_path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 9:
        raise RunnerError(f"new_src input has fewer than 9 records: {input_path}")
    engine_number = int(first_numeric_token(lines[7]))
    if engine_number not in ITERATION_ENGINES:
        raise RunnerError(
            f"unknown new_src iteration-engine number {engine_number} in {input_path}"
        )
    summary: Dict[str, object] = {
        "temperature_over_tc": first_numeric_token(lines[0]),
        "vorticity": first_numeric_token(lines[1]),
        "cylindrical": int(first_numeric_token(lines[2])),
        "fermi_liquid_f1s": first_numeric_token(lines[3]),
        "azimuthal_directions": int(first_numeric_token(lines[4])),
        "requested_error": first_numeric_token(lines[5]),
        "initial_state": int(first_numeric_token(lines[6])),
        "iteration_engine_number": engine_number,
        "iteration_engine": ITERATION_ENGINES[engine_number],
        "maximum_iterations": int(first_numeric_token(lines[8])),
    }
    if len(lines) >= 10:
        summary["anderson_p_max"] = first_numeric_token(lines[9])
    return summary


def read_input_summary(input_path: Path, solver: str) -> Dict[str, object]:
    if solver == "legacy":
        return read_legacy_input_summary(input_path)
    if solver == "new-src":
        return read_new_input_summary(input_path)
    raise RunnerError(f"unsupported solver selection: {solver}")


def ozaki_table_temperature(path: Path) -> float:
    first_line = path.read_text(encoding="utf-8").splitlines()[0]
    values = first_line.split()
    if len(values) < 2:
        raise RunnerError(f"cannot read temperature from Ozaki table: {path}")
    return float(values[1].replace("D", "E").replace("d", "e"))


def select_ozaki_table(temperature: float, source_dir: Path) -> Path:
    candidates: List[Tuple[float, Path]] = []
    for path in sorted(source_dir.glob("ozaki_T=*.dat")):
        candidates.append((abs(ozaki_table_temperature(path) - temperature), path))
    if not candidates:
        raise RunnerError("no temperature-labelled Ozaki tables were found")
    difference, path = min(candidates, key=lambda item: item[0])
    if difference > 1.0e-8:
        available = ", ".join(
            f"{ozaki_table_temperature(candidate):g}" for _, candidate in candidates
        )
        raise RunnerError(
            f"no Ozaki table matches T/Tc={temperature:g}; available: {available}"
        )
    return path


def safe_case_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", name.strip()).strip("-.")
    return cleaned or "run"


def create_unique_run_directory(run_root: Path, case_name: str) -> Path:
    run_root.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    stem = f"{stamp}-{safe_case_name(case_name)}"
    candidate = run_root / stem
    serial = 1
    while candidate.exists():
        candidate = run_root / f"{stem}-{serial:02d}"
        serial += 1
    candidate.mkdir()
    return candidate.resolve()


def write_json(path: Path, value: Dict[str, object]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_solver(args: argparse.Namespace) -> Path:
    if args.ranks < 1:
        raise RunnerError("--ranks must be positive")
    source_dir = LEGACY_SOURCE if args.solver == "legacy" else NEW_SOURCE
    default_input = source_dir / "qcv.inp"
    input_path = (args.input or default_input).expanduser().resolve()
    if not input_path.is_file():
        raise RunnerError(f"input file does not exist: {input_path}")

    parameters = read_input_summary(input_path, args.solver)
    temperature = float(parameters["temperature_over_tc"])
    ozaki = select_ozaki_table(temperature, source_dir)

    if args.executable:
        executable = args.executable.expanduser().resolve()
        build_log: Optional[Path] = None
        if not executable.is_file():
            raise RunnerError(f"solver executable does not exist: {executable}")
    else:
        parallel_jobs = min(args.ranks, os.cpu_count() or 1)
        if args.solver == "legacy":
            executable, build_log = build_legacy_solver(
                DEFAULT_LEGACY_BUILD_DIR, parallel_jobs=parallel_jobs
            )
        else:
            executable, build_log = build_new_solver(
                DEFAULT_NEW_BUILD_DIR, parallel_jobs=parallel_jobs
            )

    run_dir = create_unique_run_directory(args.run_root.expanduser().resolve(), args.case_name)
    run_executable = run_dir / "qcv"
    shutil.copy2(executable, run_executable)
    shutil.copy2(input_path, run_dir / "qcv.inp")
    shutil.copy2(source_dir / "gauss11.dat", run_dir / "gauss11.dat")
    shutil.copy2(ozaki, run_dir / "ozaki.dat")
    archived_build_log: Optional[Path] = None
    if build_log and build_log.is_file():
        archived_build_log = run_dir / "build.log"
        shutil.copy2(build_log, archived_build_log)

    # Desktop Open MPI installations can advertise fewer scheduler slots than
    # the user intends to run, especially inside an application sandbox.
    # --oversubscribe permits the requested local process count without
    # changing solver arithmetic.
    if args.ranks == 1:
        # Open MPI supports singleton execution.  Avoiding the launcher also
        # makes one-rank smoke tests usable in restricted CI environments.
        command = [str(run_executable)]
    else:
        mpirun = find_program("mpirun", (Path("/opt/homebrew/bin/mpirun"),))
        command = [
            str(mpirun),
            "--map-by",
            "slot:OVERSUBSCRIBE",
            "--bind-to",
            "none",
            "-n",
            str(args.ranks),
            str(run_executable),
        ]
    metadata: Dict[str, object] = {
        "status": "running",
        "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version,
        "runner_version": 3,
        "solver": args.solver,
        "case_name": safe_case_name(args.case_name),
        "mpi_processes": args.ranks,
        "command": command,
        "input_source": str(input_path),
        "input_sha256": sha256_file(run_dir / "qcv.inp"),
        "ozaki_source": str(ozaki),
        "ozaki_sha256": sha256_file(run_dir / "ozaki.dat"),
        "executable_source": str(executable),
        "executable_sha256": sha256_file(run_executable),
        "build_log": archived_build_log.name if archived_build_log else None,
        "parameters": parameters,
    }
    if args.solver == "legacy":
        legacy_nx = 49
        radial_spacing = float(parameters["radius_legacy_units"]) / legacy_nx
        metadata.update(
            {
                "mesh": "fixed uniform radial grid",
                "radial_cells": legacy_nx,
                "uniform_radial_spacing": radial_spacing,
                "warning": (
                    "This executable uses the legacy fixed 49-cell uniform radial "
                    "grid; it does not yet exercise adaptive refinement."
                ),
            }
        )
    else:
        metadata.update(
            {
                "mesh": "static tangent-mapped radial grid",
                "warning": (
                    "The new_src executable uses a graded but static radial grid; "
                    "it does not yet adapt from an error indicator."
                ),
            }
        )
    metadata_path = run_dir / "run_metadata.json"
    write_json(metadata_path, metadata)

    print(f"Run directory: {run_dir}")
    print(f"MPI processes: {args.ranks}")
    print(f"Solver: {args.solver}")
    if args.solver == "legacy":
        print(f"Legacy radial spacing: {metadata['uniform_radial_spacing']:.8g}")
    print(metadata["warning"])

    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    start = time.monotonic()
    try:
        with (run_dir / "qcv.inp").open("r", encoding="utf-8") as input_stream, \
                stdout_path.open("w", encoding="utf-8") as stdout_stream, \
                stderr_path.open("w", encoding="utf-8") as stderr_stream:
            result = subprocess.run(
                command,
                cwd=run_dir,
                stdin=input_stream,
                stdout=stdout_stream,
                stderr=stderr_stream,
                check=False,
            )
    except KeyboardInterrupt:
        metadata.update(
            {
                "status": "interrupted",
                "return_code": 130,
                "elapsed_seconds": time.monotonic() - start,
                "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            }
        )
        write_json(metadata_path, metadata)
        raise
    elapsed = time.monotonic() - start
    metadata.update(
        {
            "status": "completed" if result.returncode == 0 else "failed",
            "return_code": result.returncode,
            "elapsed_seconds": elapsed,
            "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
    )
    xgrid_path = run_dir / "xgrid.dat"
    if xgrid_path.is_file():
        try:
            xgrid = read_numeric_table(xgrid_path, 1)
            coordinate_column = 1 if len(xgrid[0]) >= 2 else 0
            metadata["radial_grid_points"] = len(xgrid)
            metadata["radial_extent"] = [
                xgrid[0][coordinate_column], xgrid[-1][coordinate_column]
            ]
        except RunnerError as exc:
            metadata["radial_grid_note"] = str(exc)
    write_json(metadata_path, metadata)
    if result.returncode != 0:
        raise RunnerError(
            f"MPI run failed with status {result.returncode}; see {stderr_path}"
        )
    return run_dir


def read_numeric_table(path: Path, minimum_columns: int) -> List[List[float]]:
    rows: List[List[float]] = []
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                row = [float(token.replace("D", "E").replace("d", "e"))
                       for token in stripped.split()]
            except ValueError as exc:
                raise RunnerError(f"non-numeric row in {path}:{line_number}") from exc
            if len(row) < minimum_columns:
                raise RunnerError(
                    f"expected at least {minimum_columns} columns in "
                    f"{path}:{line_number}, found {len(row)}"
                )
            if not all(math.isfinite(value) for value in row):
                raise RunnerError(f"non-finite value in {path}:{line_number}")
            rows.append(row)
    if not rows:
        raise RunnerError(f"no numeric data found in {path}")
    return rows


def load_matplotlib():
    cache_directory = PROJECT_ROOT / "work" / "matplotlib-cache"
    cache_directory.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_directory))
    try:
        import numpy as np
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ImportError as exc:
        homebrew_python = Path("/opt/homebrew/bin/python3")
        suggestion = (
            f" Run the command with {homebrew_python} instead."
            if homebrew_python.is_file()
            else " Install the packages listed in requirements.txt."
        )
        raise RunnerError(
            "the Matplotlib backend requires NumPy and Matplotlib." + suggestion
        ) from exc
    return np, plt


def write_complex_grid_matplotlib(
    data_path: Path,
    output_path: Path,
    labels: Sequence[str],
    title: str,
) -> None:
    np, plt = load_matplotlib()
    data = np.asarray(read_numeric_table(data_path, 1 + 2 * len(labels)), dtype=float)
    radius = data[:, 0]
    component_magnitudes = np.hypot(data[:, 1::2], data[:, 2::2])
    common_limit = 1.05 * float(np.max(component_magnitudes))
    if common_limit <= 0.0:
        common_limit = 1.0
    figure, axes = plt.subplots(3, 3, figsize=(18, 14), constrained_layout=True)
    figure.suptitle(title, fontsize=16)
    for index, (axis, label) in enumerate(zip(axes.flat, labels)):
        real_part = data[:, 1 + 2 * index]
        imaginary_part = data[:, 2 + 2 * index]
        magnitude = np.hypot(real_part, imaginary_part)
        axis.plot(radius, real_part, linewidth=1.8, label="Re")
        axis.plot(radius, imaginary_part, linewidth=1.8, label="Im")
        axis.plot(radius, magnitude, linewidth=1.8, linestyle="--", label="|A|")
        axis.set_title(label)
        axis.set_xlabel("radius")
        axis.set_ylabel("component")
        axis.set_ylim(-common_limit, common_limit)
        axis.grid(True, color="0.88", linewidth=0.8)
        axis.legend(loc="best", framealpha=0.9)
    figure.savefig(output_path, dpi=160)
    plt.close(figure)


def write_current_matplotlib(data_path: Path, output_path: Path) -> None:
    np, plt = load_matplotlib()
    data = np.asarray(read_numeric_table(data_path, 5), dtype=float)
    figure, axes = plt.subplots(1, 2, figsize=(16, 7), constrained_layout=True)
    figure.suptitle("Pair amplitude and current-related mean field", fontsize=16)

    axes[0].plot(data[:, 0], data[:, 1], linewidth=2.0, label="pair amplitude")
    axes[0].set_xlabel("radius")
    axes[0].set_ylabel("amplitude")
    axes[0].grid(True, color="0.88", linewidth=0.8)
    axes[0].legend(loc="best")

    if data.shape[1] >= 6:
        components = data[:, 3:6]
    else:
        components = data[:, 2:5]
    magnitude = np.linalg.norm(components, axis=1)
    axes[1].plot(
        data[:, 0], magnitude, linewidth=2.0, linestyle="--", label="|v|"
    )
    field_labels = ("vx", "vy", "vz")
    for column, label in enumerate(field_labels):
        axes[1].plot(
            data[:, 0], components[:, column], linewidth=2.0, label=label
        )
    axes[1].set_xlabel("radius")
    axes[1].set_ylabel("field")
    axes[1].grid(True, color="0.88", linewidth=0.8)
    axes[1].legend(loc="best")

    figure.savefig(output_path, dpi=160)
    plt.close(figure)


def convergence_groups(
    records: Sequence[ConvergenceRecord],
) -> List[List[ConvergenceRecord]]:
    groups: List[List[ConvergenceRecord]] = []
    for record in records:
        key = (record.engine, record.stage)
        if not groups or (groups[-1][0].engine, groups[-1][0].stage) != key:
            groups.append([record])
        else:
            groups[-1].append(record)
    return groups


def plot_convergence_axis(axis, records: Sequence[ConvergenceRecord]) -> None:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    engine_colors = {
        "NN": "#6b7280",
        "AA": "#2563eb",
        "BR": "#e76f51",
        "BB": "#7c3aed",
    }
    groups = convergence_groups(records)
    engine_counts: Dict[str, int] = {}
    for group in groups:
        engine_counts[group[0].engine] = engine_counts.get(group[0].engine, 0) + 1

    for group in groups:
        engine = group[0].engine
        color = engine_colors.get(engine, "#0f766e")
        iterations = [record.global_iteration for record in group]
        averages = [max(record.average_error, sys.float_info.min) for record in group]
        maxima = [max(record.maximum_error, sys.float_info.min) for record in group]
        axis.semilogy(
            iterations, averages, color=color, marker="o", markersize=3.0,
            linewidth=1.7
        )
        axis.semilogy(
            iterations, maxima, color=color, marker="x", markersize=3.5,
            linestyle="--", linewidth=1.5
        )
        left = iterations[0] - 0.45
        right = iterations[-1] + 0.45
        axis.axvspan(left, right, color=color, alpha=0.055, linewidth=0)
        label = engine
        if engine_counts[engine] > 1:
            label += f"-{group[0].stage}"
        midpoint = 0.5 * (iterations[0] + iterations[-1])
        axis.text(
            midpoint, 0.985, label, color=color, fontsize=8.5,
            ha="center", va="top", transform=axis.get_xaxis_transform()
        )
    for left_group, right_group in zip(groups, groups[1:]):
        boundary = 0.5 * (
            left_group[-1].global_iteration + right_group[0].global_iteration
        )
        axis.axvline(boundary, color="0.72", linewidth=0.8)

    style_handles = [
        Line2D([0], [0], color="0.15", linewidth=1.8, marker="o", label="average"),
        Line2D(
            [0], [0], color="0.15", linewidth=1.6, linestyle="--",
            marker="x", label="maximum"
        ),
    ]
    engines_in_order = list(dict.fromkeys(record.engine for record in records))
    engine_handles = [
        Patch(
            facecolor=engine_colors.get(engine, "#0f766e"), alpha=0.22,
            edgecolor=engine_colors.get(engine, "#0f766e"), label=engine
        )
        for engine in engines_in_order
    ]
    axis.set_title("Self-consistency convergence")
    axis.set_xlabel("nonlinear evaluation")
    axis.set_ylabel("residual error")
    axis.grid(True, which="both", color="0.88", linewidth=0.8)
    axis.legend(handles=style_handles + engine_handles, loc="best", fontsize=8.5)


def write_convergence_matplotlib(
    records: Sequence[ConvergenceRecord], output_path: Path
) -> None:
    _, plt = load_matplotlib()
    figure, axis = plt.subplots(figsize=(11, 7.5), constrained_layout=True)
    plot_convergence_axis(axis, records)
    figure.savefig(output_path, dpi=160)
    plt.close(figure)


def read_run_metadata(run_dir: Path) -> Dict[str, object]:
    metadata_path = run_dir / "run_metadata.json"
    if not metadata_path.is_file():
        return {}
    try:
        value = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    return value if isinstance(value, dict) else {}


def summary_title(run_dir: Path, metadata: Dict[str, object]) -> str:
    details: List[str] = [run_dir.name]
    parameters = metadata.get("parameters", {})
    if isinstance(parameters, dict):
        temperature = parameters.get("temperature_over_tc")
        if isinstance(temperature, (int, float)):
            details.append(f"T/Tc={temperature:g}")
        cylindrical = parameters.get("cylindrical")
        if isinstance(cylindrical, (int, float)):
            details.append("container" if int(cylindrical) else "free vortex")
    solver = metadata.get("solver")
    if isinstance(solver, str):
        details.append(solver)
    ranks = metadata.get("mpi_processes")
    if isinstance(ranks, int):
        details.append(f"MPI ranks={ranks}")
    return "  |  ".join(details)


INITIALIZATION_TARGETS = {
    -2: "bulk-a",
    -1: "bulk-b",
    0: "n-core",
    1: "a-core",
    2: "d-core",
    3: "restart",
}


def initialization_target_name(value: object) -> str:
    try:
        number = int(value)
    except (TypeError, ValueError):
        return "init-unknown"
    return INITIALIZATION_TARGETS.get(number, f"init-{number}")


def infer_initial_state(run_dir: Path, metadata: Dict[str, object]) -> object:
    parameters = metadata.get("parameters", {})
    if isinstance(parameters, dict) and "initial_state" in parameters:
        return parameters["initial_state"]
    input_path = run_dir / "qcv.inp"
    if not input_path.is_file():
        return None
    lines = input_path.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        if len(lines) >= 11:
            return int(first_numeric_token(lines[9]))
        if len(lines) >= 9:
            return int(first_numeric_token(lines[6]))
    except (IndexError, ValueError):
        return None
    return None


def plot_name_stem(run_dir: Path, metadata: Dict[str, object]) -> str:
    archive_match = re.match(r"^(\d{8})-(\d{6})-(.+)$", run_dir.name)
    if archive_match:
        date_stamp = archive_match.group(1)[2:]
        inferred_case_name = archive_match.group(3)
    else:
        date_stamp = ""
        inferred_case_name = run_dir.name

    if not date_stamp:
        started = metadata.get("started_utc")
        if isinstance(started, str):
            try:
                date_stamp = dt.datetime.fromisoformat(
                    started.replace("Z", "+00:00")
                ).strftime("%y%m%d")
            except ValueError:
                date_stamp = ""
    if not date_stamp:
        date_stamp = dt.datetime.now().strftime("%y%m%d")

    case_name = metadata.get("case_name", inferred_case_name)
    if not isinstance(case_name, str):
        case_name = inferred_case_name
    initialization = initialization_target_name(
        infer_initial_state(run_dir, metadata)
    )
    return f"{date_stamp}_{safe_case_name(case_name)}_{initialization}"


def unique_plot_name_stem(plot_dir: Path, requested_stem: str) -> str:
    candidate = requested_stem
    serial = 1
    while any(plot_dir.glob(f"{candidate}_*")):
        candidate = f"{requested_stem}_{serial:02d}"
        serial += 1
    return candidate


def write_summary_matplotlib(
    run_dir: Path,
    png_path: Path,
    pdf_path: Path,
    records: Sequence[ConvergenceRecord],
) -> None:
    np, plt = load_matplotlib()
    op_path = run_dir / "op_xyz"
    if not op_path.is_file():
        raise RunnerError("the panel summary requires op_xyz")
    op_data = np.asarray(
        read_numeric_table(op_path, 1 + 2 * len(CARTESIAN_COMPONENTS)),
        dtype=float,
    )
    radius = op_data[:, 0]
    order_parameter = (
        op_data[:, 1::2] + 1j * op_data[:, 2::2]
    )
    common_order_parameter_limit = 1.05 * float(np.max(np.abs(order_parameter)))
    if common_order_parameter_limit <= 0.0:
        common_order_parameter_limit = 1.0

    figure, axes = plt.subplots(2, 3, figsize=(17, 10.5), constrained_layout=True)
    orbital_colors = ("#0072b2", "#d55e00", "#009e73")
    orbital_labels = ("x", "y", "z")
    spin_labels = ("x", "y", "z")
    panel_letters = ("a", "b", "c", "d", "e", "f")

    for spin_index, axis in enumerate(axes[0, :]):
        for orbital_index, (color, orbital_label) in enumerate(
            zip(orbital_colors, orbital_labels)
        ):
            component = 3 * spin_index + orbital_index
            axis.plot(
                radius,
                np.abs(order_parameter[:, component]),
                color=color,
                linewidth=1.8,
                label=rf"$|A_{{{spin_labels[spin_index]}{orbital_label}}}|$",
            )
        axis.set_title(
            rf"({panel_letters[spin_index]}) spin row $\alpha={spin_labels[spin_index]}$"
        )
        axis.set_xlabel(r"radius $r$")
        axis.set_ylabel(r"order-parameter magnitude $|A_{\alpha i}|$")
        axis.set_ylim(0.0, common_order_parameter_limit)
        axis.grid(True, color="0.88", linewidth=0.8)
        axis.legend(loc="best", fontsize=8.5)

    density_axis = axes[1, 0]
    pair_density = np.sum(np.abs(order_parameter) ** 2, axis=1) / 3.0
    density_scale = float(np.max(pair_density))
    if density_scale <= 0.0:
        density_scale = 1.0
    density_axis.plot(
        radius, pair_density / density_scale, color="#cc79a7", linewidth=2.1
    )
    density_axis.set_title(r"(d) normalized pair-density proxy")
    density_axis.set_xlabel(r"radius $r$")
    density_axis.set_ylabel(r"$\mathrm{Tr}(AA^\dagger)/(3\,\max)$")
    density_axis.grid(True, color="0.88", linewidth=0.8)
    density_axis.text(
        0.02,
        0.04,
        r"order-parameter measure; not response-defined $\rho_s$",
        transform=density_axis.transAxes,
        fontsize=8.5,
        color="0.30",
    )

    field_axis = axes[1, 1]
    current_path = run_dir / "curr"
    if current_path.is_file():
        current_data = np.asarray(read_numeric_table(current_path, 5), dtype=float)
        if current_data.shape[1] >= 6:
            components = current_data[:, 3:6]
        else:
            components = current_data[:, 2:5]
        field_axis.plot(
            current_data[:, 0],
            np.linalg.norm(components, axis=1),
            color="0.15",
            linewidth=1.7,
            linestyle="--",
            label=r"$|\boldsymbol{v}|$",
        )
        for index, (color, label) in enumerate(zip(orbital_colors, orbital_labels)):
            field_axis.plot(
                current_data[:, 0], components[:, index], color=color,
                linewidth=1.7, label=rf"$v_{label}$"
            )
        field_axis.legend(loc="best", fontsize=8.5)
    else:
        field_axis.text(
            0.5, 0.5, "curr not available", ha="center", va="center",
            transform=field_axis.transAxes
        )
    field_axis.set_title("(e) current-related mean field")
    field_axis.set_xlabel(r"radius $r$")
    field_axis.set_ylabel(r"solver field $v_i$")
    field_axis.grid(True, color="0.88", linewidth=0.8)

    convergence_axis = axes[1, 2]
    if records:
        plot_convergence_axis(convergence_axis, records)
        convergence_axis.set_title("(f) convergence and iteration engine")
    else:
        convergence_axis.text(
            0.5, 0.5, "convergence history not available",
            ha="center", va="center", transform=convergence_axis.transAxes
        )
        convergence_axis.set_title("(f) convergence")
        convergence_axis.set_xlabel("nonlinear evaluation")
        convergence_axis.set_ylabel("residual error")
        convergence_axis.grid(True, color="0.88", linewidth=0.8)

    figure.suptitle(summary_title(run_dir, read_run_metadata(run_dir)), fontsize=15)
    figure.savefig(png_path, dpi=180)
    figure.savefig(pdf_path)
    plt.close(figure)


def gnuplot_quote(path: Path) -> str:
    return "'" + str(path.resolve()).replace("\\", "\\\\").replace("'", "\\'") + "'"


def run_gnuplot(script_path: Path) -> None:
    gnuplot = find_program("gnuplot", (Path("/opt/local/bin/gnuplot"),))
    plot_environment = os.environ.copy()
    # A user-level GNUTERM=qt setting can fail before the script's explicit
    # PNG terminal is read when gnuplot was built without Qt support.
    plot_environment.pop("GNUTERM", None)
    result = subprocess.run(
        [str(gnuplot), str(script_path)],
        cwd=script_path.parent,
        env=plot_environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RunnerError(
            f"gnuplot failed for {script_path}:\n{result.stdout.rstrip()}"
        )


def write_complex_grid_plot(
    data_path: Path,
    output_path: Path,
    script_path: Path,
    labels: Sequence[str],
    title: str,
) -> None:
    rows = read_numeric_table(data_path, 1 + 2 * len(labels))
    common_limit = 1.05 * max(
        math.hypot(row[1 + 2 * index], row[2 + 2 * index])
        for row in rows
        for index in range(len(labels))
    )
    if common_limit <= 0.0:
        common_limit = 1.0
    commands = [
        "set terminal pngcairo size 1800,1400 enhanced font ',12'",
        f"set output {gnuplot_quote(output_path)}",
        "set datafile separator whitespace",
        "set grid back lc rgb '#dddddd'",
        "set key top right opaque",
        "set xlabel 'radius'",
        "set ylabel 'component'",
        f"set yrange [{-common_limit:.17g}:{common_limit:.17g}]",
        f"set multiplot layout 3,3 rowsfirst title '{title}'",
    ]
    source = gnuplot_quote(data_path)
    for index, label in enumerate(labels):
        real_column = 2 + 2 * index
        imaginary_column = real_column + 1
        commands.extend(
            [
                f"set title '{label}'",
                (
                    f"plot {source} using 1:{real_column} with lines lw 2 "
                    "title 'Re', \\\n"
                    f"     {source} using 1:{imaginary_column} with lines lw 2 "
                    "title 'Im', \\\n"
                    f"     {source} using 1:(sqrt(column({real_column})**2+"
                    f"column({imaginary_column})**2)) with lines lw 2 dt 2 title '|A|'"
                ),
            ]
        )
    commands.extend(["unset multiplot", "unset output"])
    script_path.write_text("\n".join(commands) + "\n", encoding="utf-8")
    run_gnuplot(script_path)


def write_current_plot(data_path: Path, output_path: Path, script_path: Path) -> None:
    rows = read_numeric_table(data_path, 5)
    component_start = 4 if len(rows[0]) >= 6 else 3
    source = gnuplot_quote(data_path)
    commands = [
        "set terminal pngcairo size 1600,700 enhanced font ',13'",
        f"set output {gnuplot_quote(output_path)}",
        "set datafile separator whitespace",
        "set grid back lc rgb '#dddddd'",
        "set key top right opaque",
        "set multiplot layout 1,2 rowsfirst title 'Pair amplitude and current-related mean field'",
        "set xlabel 'radius'",
        "set ylabel 'amplitude'",
        f"plot {source} using 1:2 with lines lw 2 title 'gap amplitude'",
        "set ylabel 'field'",
        (
            f"plot {source} using 1:{component_start} with lines lw 2 title 'vx', \\\n"
            f"     {source} using 1:{component_start + 1} with lines lw 2 title 'vy', \\\n"
            f"     {source} using 1:{component_start + 2} with lines lw 2 title 'vz'"
        ),
        "unset multiplot",
        "unset output",
    ]
    script_path.write_text("\n".join(commands) + "\n", encoding="utf-8")
    run_gnuplot(script_path)


def fortran_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def parse_convergence_text(text_value: str) -> List[ConvergenceRecord]:
    number = r"[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+\-]?\d+)?"
    modern_pattern = re.compile(
        rf"solved\s*\(\s*([A-Za-z]+)\s*\)\s*:\s*(\d+).*?"
        rf"error(?:\s*\(\s*avg\s*,\s*max\s*\))?\s*=\s*"
        rf"\(\s*({number})\s+({number})\s*\)",
        flags=re.IGNORECASE,
    )
    initial_pattern = re.compile(
        rf"it\s+nr\s+(\d+)\s+errav\s*=\s*({number})\s+"
        rf"errmax\s*=\s*({number})",
        flags=re.IGNORECASE,
    )
    accelerated_pattern = re.compile(
        rf"^\s*(\d+)\s+(\d+)\s+({number})\s+({number})\s+({number})\s*$"
    )

    records: List[ConvergenceRecord] = []
    stage_counts: Dict[str, int] = {}
    active_engine: Optional[str] = None
    previous_phase_iteration: Optional[int] = None
    legacy_aa_stage = 0
    legacy_aa_iteration: Optional[int] = None

    for line in text_value.splitlines():
        match = modern_pattern.search(line)
        if match:
            engine = match.group(1).upper()
            phase_iteration = int(match.group(2))
            if engine != active_engine or (
                previous_phase_iteration is not None
                and phase_iteration <= previous_phase_iteration
            ):
                stage_counts[engine] = stage_counts.get(engine, 0) + 1
            active_engine = engine
            previous_phase_iteration = phase_iteration
            records.append(
                ConvergenceRecord(
                    len(records) + 1,
                    phase_iteration,
                    fortran_float(match.group(3)),
                    fortran_float(match.group(4)),
                    engine,
                    stage_counts[engine],
                )
            )
            continue

        match = initial_pattern.search(line)
        if match:
            phase_iteration = int(match.group(1))
            if active_engine != "NN":
                stage_counts["NN"] = stage_counts.get("NN", 0) + 1
            active_engine = "NN"
            previous_phase_iteration = phase_iteration
            records.append(
                ConvergenceRecord(
                    len(records) + 1,
                    phase_iteration,
                    fortran_float(match.group(2)),
                    fortran_float(match.group(3)),
                    "NN",
                    stage_counts["NN"],
                )
            )
            continue

        match = accelerated_pattern.match(line)
        if match:
            phase_iteration = int(match.group(1))
            if legacy_aa_iteration is None or phase_iteration <= legacy_aa_iteration:
                legacy_aa_stage += 1
            legacy_aa_iteration = phase_iteration
            active_engine = "AA"
            previous_phase_iteration = phase_iteration
            records.append(
                ConvergenceRecord(
                    len(records) + 1,
                    phase_iteration,
                    fortran_float(match.group(3)),
                    fortran_float(match.group(4)),
                    "AA",
                    legacy_aa_stage,
                )
            )
    return records


def selected_engine_from_run(run_dir: Path) -> Optional[str]:
    metadata = read_run_metadata(run_dir)
    parameters = metadata.get("parameters", {})
    if isinstance(parameters, dict):
        engine = parameters.get("iteration_engine")
        if isinstance(engine, str):
            return engine.upper()
    input_path = run_dir / "qcv.inp"
    if input_path.is_file():
        lines = input_path.read_text(encoding="utf-8", errors="replace").splitlines()
        if 9 <= len(lines) <= 10:
            try:
                number = int(first_numeric_token(lines[7]))
            except (ValueError, IndexError):
                number = -1
            if number in ITERATION_ENGINES:
                return ITERATION_ENGINES[number]
    logs = list(run_dir.glob("error_log_*.dat"))
    if logs:
        newest = max(logs, key=lambda path: path.stat().st_mtime)
        return newest.stem.removeprefix("error_log_").upper()
    return None


def records_from_error_log(run_dir: Path) -> List[ConvergenceRecord]:
    engine = selected_engine_from_run(run_dir)
    if not engine:
        return []
    path = run_dir / f"error_log_{engine}.dat"
    if not path.is_file():
        return []
    rows = read_numeric_table(path, 3)
    return [
        ConvergenceRecord(
            global_iteration=index,
            phase_iteration=int(row[0]),
            average_error=row[1],
            maximum_error=row[2],
            engine=engine,
            stage=1,
        )
        for index, row in enumerate(rows, start=1)
    ]


def write_convergence_table(
    records: Sequence[ConvergenceRecord], table_path: Path
) -> None:
    with table_path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# global_iteration phase_iteration average_error maximum_error "
            "engine stage\n"
        )
        for record in records:
            stream.write(
                f"{record.global_iteration} {record.phase_iteration} "
                f"{record.average_error:.17e} {record.maximum_error:.17e} "
                f"{record.engine} {record.stage}\n"
            )


def extract_convergence(run_dir: Path, table_path: Path) -> List[ConvergenceRecord]:
    records: List[ConvergenceRecord] = []
    for source_path in (run_dir / "stdout.log", run_dir / "fort.30"):
        if source_path.is_file():
            records = parse_convergence_text(
                source_path.read_text(encoding="utf-8", errors="replace")
            )
            if records:
                break
    if not records:
        records = records_from_error_log(run_dir)
    if records:
        write_convergence_table(records, table_path)
    return records


def write_convergence_plot(data_path: Path, output_path: Path, script_path: Path) -> None:
    source = gnuplot_quote(data_path)
    commands = [
        "set terminal pngcairo size 1100,750 enhanced font ',13'",
        f"set output {gnuplot_quote(output_path)}",
        "set datafile separator whitespace",
        "set grid back lc rgb '#dddddd'",
        "set key top right opaque",
        "set xlabel 'nonlinear iteration'",
        "set ylabel 'error'",
        "set logscale y",
        "set title 'Self-consistency convergence'",
        (
            f"plot {source} using 1:3 with linespoints lw 2 title 'average', \\\n"
            f"     {source} using 1:4 with linespoints lw 2 title 'maximum'"
        ),
        "unset output",
    ]
    script_path.write_text("\n".join(commands) + "\n", encoding="utf-8")
    run_gnuplot(script_path)


def plot_results(
    run_dir: Path,
    requested_plot_dir: Optional[Path] = None,
    backend: str = "matplotlib",
) -> Path:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.is_dir():
        raise RunnerError(f"result directory does not exist: {run_dir}")
    plot_dir = (
        requested_plot_dir.expanduser().resolve()
        if requested_plot_dir
        else run_dir / "plots"
    )
    plot_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_run_metadata(run_dir)
    requested_stem = plot_name_stem(run_dir, metadata)
    output_stem = unique_plot_name_stem(plot_dir, requested_stem)

    produced: List[Path] = []
    op_xyz = run_dir / "op_xyz"
    if op_xyz.is_file():
        output = plot_dir / f"{output_stem}_order_parameter_cartesian.png"
        if backend == "matplotlib":
            write_complex_grid_matplotlib(
                op_xyz, output, CARTESIAN_COMPONENTS,
                "Cartesian order-parameter components"
            )
        else:
            write_complex_grid_plot(
                op_xyz,
                output,
                plot_dir / f"{output_stem}_order_parameter_cartesian.gnuplot",
                CARTESIAN_COMPONENTS,
                "Cartesian order-parameter components",
            )
        produced.append(output)

    op_harm = run_dir / "op_harm"
    if op_harm.is_file():
        output = plot_dir / f"{output_stem}_order_parameter_harmonic.png"
        if backend == "matplotlib":
            write_complex_grid_matplotlib(
                op_harm, output, HARMONIC_COMPONENTS,
                "Spherical/harmonic order-parameter components"
            )
        else:
            write_complex_grid_plot(
                op_harm,
                output,
                plot_dir / f"{output_stem}_order_parameter_harmonic.gnuplot",
                HARMONIC_COMPONENTS,
                "Spherical/harmonic order-parameter components",
            )
        produced.append(output)

    current = run_dir / "curr"
    if current.is_file():
        output = plot_dir / f"{output_stem}_amplitude_and_fields.png"
        if backend == "matplotlib":
            write_current_matplotlib(current, output)
        else:
            write_current_plot(
                current,
                output,
                plot_dir / f"{output_stem}_amplitude_and_fields.gnuplot",
            )
        produced.append(output)

    convergence_table = plot_dir / f"{output_stem}_convergence.tsv"
    convergence_records = extract_convergence(run_dir, convergence_table)
    if convergence_records:
        output = plot_dir / f"{output_stem}_convergence.png"
        if backend == "matplotlib":
            write_convergence_matplotlib(convergence_records, output)
        else:
            write_convergence_plot(
                convergence_table,
                output,
                plot_dir / f"{output_stem}_convergence.gnuplot",
            )
        produced.append(output)

    if backend == "matplotlib" and op_xyz.is_file():
        summary_png = plot_dir / f"{output_stem}_run_summary.png"
        summary_pdf = plot_dir / f"{output_stem}_run_summary.pdf"
        write_summary_matplotlib(
            run_dir, summary_png, summary_pdf, convergence_records
        )
        produced.extend((summary_png, summary_pdf))

    if not produced:
        raise RunnerError(
            f"no recognized result files found in {run_dir}; expected op_xyz, "
            "op_harm, curr, stdout.log, fort.30, or error_log_*.dat"
        )

    manifest = {
        "source_directory": str(run_dir),
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "backend": backend,
        "filename_stem": output_stem,
        "initialization_target": initialization_target_name(
            infer_initial_state(run_dir, metadata)
        ),
        "plots": [str(path) for path in produced],
        "derived_quantities": {
            "pair_density_proxy": (
                "Tr(A A^dagger)/3, normalized by its maximum in the plotted "
                "profile; this is not the response-defined superfluid density"
            ),
            "current_related_field": (
                "the solver output fields vx, vy, vz; these combine current and "
                "diagonal/Fermi-liquid self-energy information"
            ),
        },
    }
    write_json(plot_dir / f"{output_stem}_plot_manifest.json", manifest)
    print("Plots:")
    for path in produced:
        print(f"  {path}")
    return plot_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_arguments(argv)
    try:
        if args.plot_only:
            plot_results(args.plot_only, args.plot_dir, args.plot_backend)
            return 0

        run_dir = run_solver(args)
        if not args.no_plot:
            plot_results(run_dir, args.plot_dir, args.plot_backend)
        print(f"Completed run: {run_dir}")
        return 0
    except RunnerError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
