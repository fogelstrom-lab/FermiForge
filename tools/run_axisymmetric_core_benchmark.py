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
    parser.add_argument("--ranks", type=int, default=10)
    parser.add_argument("--jobs", type=int, default=10)
    parser.add_argument("--cells", type=int, help="override number_of_cells")
    parser.add_argument("--half-width", type=float, help="override half_width")
    parser.add_argument(
        "--mesh-kind",
        choices=("uniform", "multiscale"),
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
        "--trajectory-step", type=float, help="override trajectory_maximum_step"
    )
    parser.add_argument(
        "--endpoint-policy",
        choices=("local", "free_vortex", "radial_reference"),
        help="override endpoint_policy",
    )
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIRECTORY)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
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
    command: Sequence[str], log_path: Path, environment: dict[str, str], append: bool = False
) -> None:
    mode = "a" if append else "w"
    with log_path.open(mode, encoding="utf-8") as log:
        log.write("command: " + " ".join(command) + "\n")
        log.flush()
        result = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if result.returncode != 0:
        raise RunnerError(
            f"command failed with status {result.returncode}; see {log_path}"
        )


def replace_namelist_value(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?im)^(\s*{re.escape(key)}\s*=\s*)[^,\n/]+")
    updated, count = pattern.subn(rf"\g<1>{value}", text, count=1)
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
        input_text = input_path.read_text(encoding="utf-8")
        overrides = {
            "number_of_cells": arguments.cells,
            "half_width": arguments.half_width,
            "fine_region_half_width": arguments.fine_region,
            "medium_region_half_width": arguments.medium_region,
            "fine_spacing": arguments.fine_spacing,
            "medium_spacing": arguments.medium_spacing,
            "coarse_spacing": arguments.coarse_spacing,
            "active_radius": arguments.active_radius,
            "asymptotic_outer_radius": arguments.outer_radius,
            "trajectory_maximum_step": arguments.trajectory_step,
        }
        for key, value in overrides.items():
            if value is not None:
                input_text = replace_namelist_value(input_text, key, str(value))
        if arguments.endpoint_policy is not None:
            input_text = replace_namelist_value(
                input_text, "endpoint_policy", f"'{arguments.endpoint_policy}'"
            )
        if arguments.mesh_kind is not None:
            input_text = replace_namelist_value(
                input_text, "mesh_kind", f"'{arguments.mesh_kind}'"
            )
        archived_input.write_text(input_text, encoding="utf-8")

        environment = os.environ.copy()
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
        metrics = run_directory / "metrics.txt"
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
            ]
        )
        run_logged(mpi_command, run_directory / "run.log", environment)
        for output in (input_map, mapped_map, metrics):
            if not output.is_file():
                raise RunnerError(f"benchmark did not produce {output}")

        plots: list[str] = []
        if not arguments.no_plot:
            python = find_plotting_python(environment)
            plot_script = PROJECT_ROOT / "tools" / "plot_2d_fields.py"
            plot_directory = run_directory / "plots"
            for role, field_map in (("input", input_map), ("mapped", mapped_map)):
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
            plots = [str(path) for path in sorted(plot_directory.glob("*"))]

        manifest = {
            "kind": "fermiforge_axisymmetric_core_2d_benchmark",
            "created_at": dt.datetime.now().astimezone().isoformat(),
            "case_name": safe_name(arguments.case_name),
            "initialization": safe_name(arguments.initialization),
            "mpi_ranks": arguments.ranks,
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
                    "trajectory_maximum_step": arguments.trajectory_step,
                    "endpoint_policy": arguments.endpoint_policy,
                }.items()
                if value is not None
            },
            "input_source": str(input_path),
            "input_sha256": sha256(archived_input),
            "executable": str(executable),
            "executable_sha256": sha256(executable),
            "input_field_map_sha256": sha256(input_map),
            "mapped_field_map_sha256": sha256(mapped_map),
            "metrics_sha256": sha256(metrics),
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
    if plots:
        print(f"plots:         {run_directory / 'plots'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
