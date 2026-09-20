# Normal-core 2D benchmark

## Purpose

This is the first physics regression for removing cylindrical symmetry.  An
axisymmetric normal-core radial solution is embedded on a Cartesian 2D mesh,
the full spin-matrix quasiclassical self-consistency map is evaluated at every
mesh point, and the result is compared with the embedded radial solution.

The benchmark separates four questions:

1. Can both historical uniform and current tangent radial grids be read
   without silently changing their coordinates?
2. Does radial rotation reproduce the source solver's Cartesian reconstruction
   at every 2D grid point?
3. Does the propagated and momentum-integrated 2D map complete without loss of
   propagator normalization?
4. Does the full repeated 2D self-consistency solve converge under the legacy
   Anderson engine?
5. Does its final field reproduce the radial qcv state and remain stable when
   the computational cell is enlarged?

Questions 1--3 are acceptance checks. The comparison in question 2 evaluates
the radial profile directly at each `(x,y)` point; it deliberately does not
infer harmonics back from an off-axis Cartesian field because the legacy
transport reconstruction enforces `A_xy=-A_yx` and is lossy in that subspace.
The RMS, maximum, and relative L2 map residuals are recorded at every
evaluation. The independently tested one-point oracle supplies the
physics-kernel acceptance check. The static multiscale mesh now supports a
complete fixed-mesh nonlinear solve and restart; automatic error-driven block
refinement remains under development.

## Reference status

`examples/2d_normal_core_benchmark.nml` uses the converged current-solver
reference under `reference/current_new_src_T0.30_F1s5.4/`. Its manifest records
the complete controls, final residuals, executable hash, and file hashes. The
older `new_src/op_xyz_n` and `new_src/curr_n` pair remains a useful secondary
50-point uniform-grid check, but is not the primary physics reference.

The A-phase-core extension now uses the same executable and metrics through
`examples/2d_a_phase_core_benchmark.nml` and
`benchmarks/a_phase_core_2d/`, with no core-specific code path.

## Run

From the FermiForge project root:

```text
/opt/homebrew/bin/python3 tools/run_axisymmetric_core_benchmark.py \
  --case-name normal-core --initialization radial-reference --ranks 10
```

This is the preferred one-map command: it builds the benchmark, creates a timestamped
archive, preserves the input and logs, hashes the executable and data, and
plots both the embedded reference and mapped field with the accessible colour
scales. The lower-level equivalent is:

```text
cmake -S . -B work/normal-core-build \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_BUILD_TYPE=Release
cmake --build work/normal-core-build --parallel 10 \
  --target benchmark_axisymmetric_core_2d
mkdir -p work/normal-core-2d
mpiexec --host localhost:10 --map-by ppr:10:node --bind-to none -n 10 \
  work/normal-core-build/benchmark_axisymmetric_core_2d \
  examples/2d_normal_core_benchmark.nml
```

The run writes the embedded radial input, the quasiclassical mapped field, the
final state, an iteration history, and a machine-readable metrics file. The
full converged multiscale and computational-cell commands are documented in
`docs/RUNNING_NORMAL_CORE_2D_SOLVER.md`.

## Required convergence matrix

Before accepting a numerical residual as the normal-core baseline, repeat at
least:

- Cartesian spacings `1.0`, `0.5`, and `0.25` at fixed physical extent;
- rectangle-matched free-endpoint radii `24`, `36`, and `48` (or larger if the
  trend has not settled);
- radial-reference endpoints beyond the source matching radius
  `Rc=45.9993349989`, including radii `60`, `70`, and `100`;
- one, four, and ten MPI ranks at an identical discretization;
- local, free-vortex, and radial-reference endpoint policies. The last is an
  axisymmetric oracle and must not be used as a general 2D boundary law.

The longer-term mesh target is a circular active region with fine core blocks
and coarser outer blocks. It must match the asymptotic expansion on a circle of
fixed radius, as in the radial solver, instead of inheriting a
direction-dependent matching radius from the rectangular mesh edge.

The first static version of that improvement is available through
`examples/2d_normal_core_multiscale.nml`. It uses nonuniform tensor-product
zones and a circular active mask inside a fixed radial-reference halo. It is a
validated stepping stone to block AMR, not yet automatic adaptation.
