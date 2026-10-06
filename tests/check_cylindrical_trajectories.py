#!/usr/bin/env python3
"""Compile and compare the cylinder port with both frozen legacy sources."""
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "new_src" / "3DFS_MPICodes"


def run(command, **kwargs):
    subprocess.run(command, check=True, **kwargs)


def main():
    with tempfile.TemporaryDirectory(prefix="fermiforge-cylinder-") as directory:
        work = Path(directory)
        legacy = work / "legacy_trajectories.f"
        legacy.write_text(
            (SOURCES / "trajectories.f").read_text()
            .replace("trajectories(", "legacy_trajectories(")
            .replace("zigzag(", "legacy_zigzag(")
        )
        for mode, precision in (
            ("default", []),
            ("real8", ["-fdefault-real-8", "-fdefault-double-8"]),
        ):
            # Renamed fixed-form symbols need more than 72 columns.
            flags = ["gfortran", "-O0", "-g", "-fcheck=all",
                     "-ffpe-trap=invalid,zero,overflow", *precision]
            executable = work / f"regression-{mode}"
            run([*flags, "-ffixed-line-length-none", str(legacy),
                 str(SOURCES / "trajectories.f90"),
                 str(ROOT / "tests/test_cylindrical_trajectories.f90"),
                 "-o", str(executable)], cwd=work)
            run([str(executable)], cwd=work)
            modern = work / f"poretraj-{mode}"
            original = work / f"legacy-poretraj-{mode}"
            run([*flags, "-std=f2018", "-fimplicit-none",
                 str(SOURCES / "trajectories.f90"),
                 str(SOURCES / "poretraj.f90"), "-o", str(modern)], cwd=work)
            run([*flags, str(SOURCES / "poretraj.f"), "-o", str(original)], cwd=work)
            for direction in ("0.48 0.64 0.6", "0.48 0.64 -0.6", "1 0 0", "0 0 1"):
                data = f"2\n0.13\n120\n0.3 -0.2 0.7\n{direction}\n"
                run([str(original)], input=data, text=True,
                    stdout=subprocess.DEVNULL, cwd=work)
                expected = (work / "traj").read_bytes()
                run([str(modern)], input=data, text=True,
                    stdout=subprocess.DEVNULL, cwd=work)
                if (work / "traj").read_bytes() != expected:
                    raise AssertionError(f"Driver output mismatch: {mode}, {direction}")
            print(f"PASS: {mode} driver output matches original poretraj.f")


if __name__ == "__main__":
    main()
