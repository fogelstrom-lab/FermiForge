#!/usr/bin/env python3
"""Short standalone checks for new_src's optional specular cylinder path."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
NEW = ROOT / "new_src"


def main():
    with tempfile.TemporaryDirectory(prefix="fermiforge-radial-cylinder-") as directory:
        work = Path(directory)
        flags = ["mpifort", "-std=f2018", "-fimplicit-none", "-O0", "-g",
                 "-fdefault-real-8", "-fdefault-double-8", "-fcheck=all",
                 "-ffpe-trap=invalid,zero,overflow"]
        sampler = work / "test_sampling"
        geometry = NEW / "3DFS_MPICodes/trajectories.f90"
        subprocess.run([*flags, str(NEW / "global_dec.f90"), str(geometry),
                        str(NEW / "interpol.f90"),
                        str(ROOT / "tests/test_radial_cylinder_sampling.f90"),
                        "-o", str(sampler)], cwd=work, check=True)
        subprocess.run([str(sampler)], cwd=work, check=True)
        executable = work / "test_map"
        sources = [NEW / name for name in (
            "global_dec.f90", "bulkgap.f90", "init_calc.f90", "mpicalls.f90")]
        sources += [geometry, NEW / "interpol.f90", NEW / "riccati.f90",
                    NEW / "getnewses.f90", ROOT / "tests/test_radial_cylinder_map.f90"]
        # On macOS the MPI-linked probe SIGILLs with GNU FP traps enabled;
        # the same probe passes without them. Keep bounds and explicit finite
        # checks for this probe, and retain FP traps in the non-MPI sampler.
        map_flags = [flag for flag in flags
                     if sys.platform != "darwin" or not flag.startswith("-ffpe-trap=")]
        subprocess.run([*map_flags, *map(str, sources), "-o", str(executable)],
                       cwd=work, check=True)
        shutil.copy2(NEW / "gauss11.dat", work / "gauss11.dat")
        shutil.copy2(NEW / "ozaki_T=0.3.dat", work / "ozaki.dat")
        records = (ROOT / "examples/radial_specular_cylinder_smoke.inp").read_text().splitlines()
        for radius in (4.0, 6.0):
            records[10] = str(radius)
            result = subprocess.run([str(executable)], input="\n".join(records)+"\n",
                                    text=True, cwd=work, capture_output=True, timeout=60)
            if result.returncode:
                raise RuntimeError(result.stdout + result.stderr)
            rows = (work / "cylinder_map.dat").read_text().splitlines()
            assert len(rows) == 100
            assert float(rows[-1].split()[0]) == radius
            print(f"PASS: radius {radius:g} reduced-quadrature cylinder map")
            expected = (work / "cylinder_map.dat").read_bytes()
            extended = records[:3]+["40.0", "4.0"]+records[3:]
            result = subprocess.run([str(executable)], input="\n".join(extended)+"\n",
                                    text=True, cwd=work, capture_output=True, timeout=60)
            if result.returncode:
                raise RuntimeError(result.stdout + result.stderr)
            assert (work / "cylinder_map.dat").read_bytes() == expected
            print("PASS: legacy and explicit-grid cylinder inputs give identical maps")
        # Reject silent use of a restart on an unrelated radial grid.
        records[6] = "3"
        records[10] = "4.0"
        result = subprocess.run([str(executable)], input="\n".join(records)+"\n",
                                text=True, cwd=work, capture_output=True, timeout=10)
        assert result.returncode != 0 and "grid does not match radius" in result.stderr
        records[6] = "1"
        for radius in ("0", "-1", "NaN", None):
            invalid = records[:10] if radius is None else records[:10]+[radius]
            result = subprocess.run([str(executable)], input="\n".join(invalid)+"\n",
                                    text=True, cwd=work, capture_output=True, timeout=10)
            assert result.returncode != 0 and ("radius" in result.stderr)
        print("PASS: mismatched restart and invalid/missing radius are rejected")


if __name__ == "__main__":
    main()
