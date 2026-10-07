"""Bounded MPI regression: 20 AA, 3 simple, fresh AA, on a tiny grid.

Run from the repository with: python3 tests/check_iteration_cycle.py BUILD_DIR
Uses temporary outputs; this is an engine test, not a physical benchmark.
"""
import os
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from run_axisymmetric_core_benchmark import replace_namelist_value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("build_dir", type=Path)
    parser.add_argument("--simple-steps", type=int, choices=(0, 3), default=3)
    parser.add_argument("--method", choices=("anderson", "bb", "bb-radial-style", "polyak"), default="anderson")
    args = parser.parse_args()
    executable = args.build_dir.resolve() / "benchmark_axisymmetric_core_2d"
    text = (ROOT / "examples/2d_normal_core_convergence_smoke.nml").read_text()
    # Reaching the iteration limit is intentional for this control-flow test.
    text = text.replace("\n/\n", "\n  require_convergence = .false.\n/\n")
    changes = dict(number_of_cells=4, half_width=2.0, active_radius=1.5,
                   fine_region_half_width=0.5, medium_region_half_width=1.1,
                   maximum_iterations=24, convergence_tolerance=1e-30,
                   anderson_cycle_iterations=20, simple_cycle_iterations=args.simple_steps,
                   simple_mixing=0.1)
    if args.method.startswith("bb"):
        changes.update(iteration_method="'bb'", maximum_iterations=4,
                       anderson_cycle_iterations=0)
    if args.method == "bb-radial-style":
        changes.update(bb_initial_mixing=1.0, bb_maximum_mixing=100.0,
                       bb_growth_limit=0.0, bb_curvature="'absolute'")
    if args.method == 'polyak':
        changes.update(iteration_method="'polyak'", maximum_iterations=4,
                       anderson_cycle_iterations=0,polyak_step_size=2.0,polyak_drag=0.5)
    for key, value in changes.items():
        text = replace_namelist_value(text, key, str(value))
    env = dict(os.environ, OMPI_MCA_pml="ob1", OMPI_MCA_btl="self,sm")
    with tempfile.TemporaryDirectory(prefix="fermiforge-cycle-") as directory:
        folder = Path(directory)
        input_file = folder / "input.nml"
        input_file.write_text(text)
        outputs = [folder / name for name in
                   ("initial.dat", "mapped.dat", "metrics.txt", "final.dat",
                    "history.dat", "checkpoint.dat")]
        result = subprocess.run(["mpiexec", "--host", "localhost:2", "--map-by",
                                 "ppr:2:node", "--bind-to", "none", "-n", "2",
                                 str(executable), str(input_file),
                                 *map(str, outputs)], cwd=ROOT, env=env,
                                capture_output=True, text=True, timeout=240)
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        lines = outputs[4].read_text().splitlines()
        names = lines[0].lstrip("#").split()
        rows = [dict(zip(names, map(float, line.split()))) for line in lines[1:]]
        if args.method == 'polyak':
            assert len(rows)==4
            assert all(row['engine']==3 and row['history']==0 and row['p']==2 for row in rows)
        elif args.method.startswith("bb"):
            assert len(rows) == 4
            assert all(row["engine"] == 2 and row["history"] == 0 for row in rows)
            initial, upper = (1.0, 100.0) if args.method == "bb-radial-style" else (0.1, 5.0)
            assert abs(rows[0]["p"]-initial) < 1e-14
            assert all(0.001 <= row["p"] <= upper for row in rows)
            if args.method == "bb-radial-style":
                assert all(row["bb_status"] != 5 for row in rows)
            assert rows[1]["bb_formula"] == 2 or rows[1]["bb_status"] > 0
        else:
            assert len(rows) == 24, len(rows)
            period = 20 + args.simple_steps
            assert [row["engine"] for row in rows] == [
                int(i % period < 20) for i in range(24)]
            assert [row["cycle_position"] for row in rows] == [
                i % period + 1 for i in range(24)]
            assert all(row["history"] == 0 and abs(row["p"]-0.1) < 1e-14
                       for row in rows[20:period])
            assert rows[period]["history"] == 1 and abs(rows[period]["p"]-0.01) < 1e-14
        subprocess.run([sys.executable, str(ROOT / "tools/plot_2d_iteration_history.py"),
                        str(outputs[4]), "--output", str(folder / "cycle.png")], check=True)
        assert (folder / "cycle.png").stat().st_size > 1000
    print(f"PASS: two-rank {args.method}, simple count {args.simple_steps}, history, mixing factors, and plot")


if __name__ == "__main__":
    main()
