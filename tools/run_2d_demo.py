#!/usr/bin/env python3
"""Build, run, archive, and plot the isolated FermiForge 2D demonstration."""

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
DEFAULT_INPUT = PROJECT_ROOT / "examples" / "2d_sampler_demo.nml"
DEFAULT_BUILD_DIRECTORY = PROJECT_ROOT / "work" / "fermiforge-2d-build"
DEFAULT_RUN_ROOT = PROJECT_ROOT / "runs"


class RunnerError(RuntimeError):
    """Expected build, execution, or plotting failure."""


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the manufactured FermiForge 2D mesh/interpolation example, "
            "archive its input and output, and create accessible plots."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Fortran namelist input (default: examples/2d_sampler_demo.nml)",
    )
    parser.add_argument(
        "--case-name",
        default="sampler-demo",
        help="short case label used in the archive and plot filenames",
    )
    parser.add_argument(
        "--run-root",
        type=Path,
        default=DEFAULT_RUN_ROOT,
        help="parent directory for timestamped run archives",
    )
    parser.add_argument(
        "--build-dir",
        type=Path,
        default=DEFAULT_BUILD_DIRECTORY,
        help="out-of-source CMake build directory",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=10,
        help="parallel build jobs (default: 10)",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="reuse an existing fermiforge_2d_demo executable",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="generate the field map without figures",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "pdf"),
        default=("png", "pdf"),
        help="plot formats (default: png pdf)",
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


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", value.strip()).strip("-.")
    return cleaned or "sampler-demo"


def create_run_directory(run_root: Path, case_name: str) -> Path:
    run_root.mkdir(parents=True, exist_ok=True)
    timestamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    base = run_root / f"{timestamp}-{safe_name(case_name)}-2d"
    candidate = base
    serial = 1
    while candidate.exists():
        serial += 1
        candidate = Path(f"{base}-{serial:02d}")
    candidate.mkdir(parents=True)
    return candidate


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_logged(
    command: Sequence[str],
    log_path: Path,
    *,
    environment: dict[str, str],
    append: bool = False,
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
        raise RunnerError(f"command failed with status {result.returncode}; see {log_path}")


def build_demo(build_directory: Path, jobs: int, environment: dict[str, str]) -> Path:
    cmake = find_program(
        "cmake", (Path("/opt/local/bin/cmake"), Path("/opt/homebrew/bin/cmake"))
    )
    ninja = find_program("ninja", (Path("/opt/local/bin/ninja"),))
    gfortran = find_program("gfortran", (Path("/opt/homebrew/bin/gfortran"),))
    build_directory.mkdir(parents=True, exist_ok=True)
    log_path = build_directory / "fermiforge_2d_demo_build.log"

    configure = [
        str(cmake),
        "-S",
        str(PROJECT_ROOT),
        "-B",
        str(build_directory),
        "-G",
        "Ninja",
        f"-DCMAKE_MAKE_PROGRAM={ninja}",
        f"-DCMAKE_Fortran_COMPILER={gfortran}",
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    build = [
        str(cmake),
        "--build",
        str(build_directory),
        "--target",
        "fermiforge_2d_demo",
        "--parallel",
        str(max(1, jobs)),
    ]
    run_logged(configure, log_path, environment=environment)
    run_logged(build, log_path, environment=environment, append=True)

    executable = build_directory / "fermiforge_2d_demo"
    if not executable.is_file():
        raise RunnerError(f"build did not produce {executable}")
    return executable


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(argv)
    try:
        input_path = arguments.input.expanduser().resolve()
        if not input_path.is_file():
            raise RunnerError(f"input file does not exist: {input_path}")

        run_directory = create_run_directory(
            arguments.run_root.expanduser().resolve(), arguments.case_name
        )
        archived_input = run_directory / "input.nml"
        shutil.copy2(input_path, archived_input)
        field_map = run_directory / "fields_2d.dat"
        plot_directory = run_directory / "plots"

        environment = os.environ.copy()
        matplotlib_cache = PROJECT_ROOT / "work" / "matplotlib"
        matplotlib_cache.mkdir(parents=True, exist_ok=True)
        environment.setdefault("MPLCONFIGDIR", str(matplotlib_cache))

        executable = arguments.build_dir.expanduser().resolve() / "fermiforge_2d_demo"
        if not arguments.no_build:
            executable = build_demo(
                arguments.build_dir.expanduser().resolve(), arguments.jobs, environment
            )
        if not executable.is_file():
            raise RunnerError(f"2D demo executable does not exist: {executable}")

        run_log = run_directory / "run.log"
        run_logged(
            [str(executable), str(archived_input), str(field_map)],
            run_log,
            environment=environment,
        )
        if not field_map.is_file():
            raise RunnerError("2D demo completed without producing a field map")

        plot_paths: list[str] = []
        if not arguments.no_plot:
            python_candidates = (
                Path("/opt/homebrew/bin/python3"),
                Path("/opt/local/bin/python3"),
                Path(sys.executable),
            )
            python = next(
                (
                    path.absolute()
                    for path in python_candidates
                    if path.is_file() and os.access(path, os.X_OK)
                ),
                None,
            )
            if python is None:
                raise RunnerError("no Python interpreter was found for plotting")
            plot_script = PROJECT_ROOT / "tools" / "plot_2d_fields.py"
            plot_prefix = (
                dt.datetime.now().strftime("%y%m%d")
                + "_"
                + safe_name(arguments.case_name)
                + "_manufactured"
            )
            plot_log = run_directory / "plot.log"
            plot_command = [
                str(python),
                str(plot_script),
                str(field_map),
                "--output-dir",
                str(plot_directory),
                "--prefix",
                plot_prefix,
                "--formats",
                *arguments.formats,
            ]
            run_logged(plot_command, plot_log, environment=environment)
            plot_paths = [str(path) for path in sorted(plot_directory.glob("*"))]

        manifest = {
            "kind": "fermiforge_2d_manufactured_demo",
            "created_at": dt.datetime.now().astimezone().isoformat(),
            "case_name": safe_name(arguments.case_name),
            "input_source": str(input_path),
            "input_sha256": sha256(archived_input),
            "executable": str(executable),
            "executable_sha256": sha256(executable),
            "field_map": str(field_map),
            "field_map_sha256": sha256(field_map),
            "plots": plot_paths,
            "physical_status": (
                "manufactured interpolation/visualization test; "
                "not a self-consistent quasiclassical solution"
            ),
        }
        (run_directory / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
    except (OSError, RunnerError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    print(f"run directory: {run_directory}")
    print(f"field map:    {field_map}")
    if plot_paths:
        print(f"plots:        {plot_directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
