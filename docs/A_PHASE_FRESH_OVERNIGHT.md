# Fresh A-phase overnight run

From the project root:

```bash
caffeinate -i env OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm \
  python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_a_phase_fresh_overnight.nml \
  --build-dir work/double-core-scratch-build \
  --case-name a-phase-fresh-overnight --initialization historical-aop \
  --ranks 10 --formats png
```

No restart argument. The input selects historical_aop and up to 50 updates,
with a checkpoint every update. Convergence can stop the run earlier.
The historical seed has zero initial Fermi-liquid mean field; the old radial
data are used only for comparison, not initialization or trajectory tails.

Settings match run 260922-204123: T/Tc=0.30 (from the same Ozaki table),
F1s=5.4 (feedback scale 1.9285714285714286), bulk gap 0.2799, half-width 20,
65x65 smooth mesh, stretch 2.3, active radius 16, quadratic trajectory
interpolation, trajectory step 0.25, fit radii 10 and 14, tail radius 70,
Anderson maximum p=5 and history=10, tolerance 2e-6. Only the seed and
iteration limit change. Historical nop/aop/dop select initialization only;
all share the same 2D equations and evolving asymptotic implementation.

The runner creates a new timestamped run directory. Saved products include
input/final fields, checkpoints, residual snapshots, convergence history,
field and harmonic plots, and corrected radial-reference comparisons.
Agreement with the legacy radial solver is a test outcome, not guaranteed:
its historical transport reconstruction differs from the full tensor inverse.
The mesh and matching radii also require independent convergence checks.
