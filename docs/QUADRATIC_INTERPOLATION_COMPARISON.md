# Controlled quadratic trajectory comparison

The reference is the user's ten-update linear run:
`runs/260921-115501-double-core-Fs1-zero-reference`.
The new input copies that archived input's physical and numerical controls,
changes only `trajectory_interpolation_order=2`, and clears old output paths
so the runner creates fresh timestamped products. Run ten updates from the
same original legacy state, not the already relaxed linear checkpoint:

```sh
cd /Users/mikael/Documents/Codex/3he-vortex-modernization
caffeinate -i python3 tools/run_double_core_from_scratch.py \
  --input examples/2d_double_core_Fs1_zero_quadratic.nml \
  --run-name double-core-Fs1-zero-quadratic --ranks 10 --iterations 10
```

This automatically stops after at most ten updates and generates the same
plots/checkpoints. The production comparison has not been launched by the
assistant. The mesh, quadrature, asymptotic fit, halo, mixing parameters,
and initial resampling match the archived linear run. Initial resampling
and interpolation used by the asymptotic fit remain bilinear intentionally:
this first test isolates trajectory interpolation rather than changing both
the transport and the outer-boundary construction simultaneously.

## Numerical definition

Interpolation acts on the real/imaginary Cartesian matrix components and
the current-related vector, not amplitude/phase. Each axis uses actual
nonuniform coordinates. In an interior cell, average the Lagrange quadratic
through its left neighbour and endpoints with the quadratic through its
endpoints and right neighbour. At a domain boundary use the single available
one-sided quadratic. Take the tensor product of these weights in x and y.
The result is degree two in each coordinate and uses up to 16 grid nodes.
Both constituent polynomials reproduce the cell endpoints, preserving
continuity across nodes and avoiding a nearest-node stencil switch halfway
through a cell. It reproduces all tensor-product quadratics exactly.

This is not a spline: derivatives can still jump at cell boundaries and
negative weights can produce overshoot. No limiter is applied in this first
comparison. At least three nodes on each axis are required; interpolation
does not extrapolate beyond the mesh. The default remains bilinear (order 1).
The larger stencil can increase memory use and runtime.

## What to compare

Compare raw initial map defects, iteration histories, residual bands at
|x| or |y|=6 and 12, core separation, and initial-to-final field changes.
Residual maxima may switch between symmetry-related points. A smaller
residual is not by itself proof of smaller discretization error because
the map itself has changed. Compare wall time as well as iteration count.
Snapshots are pre-update; map 10 evaluates the state after nine updates.
Use a separate zero-update restart from each ten-update checkpoint if
comparing final-state residuals. See RUNNING_FS1_ZERO_REFERENCE.md.

While running, regenerate plots between completed map writes:

```sh
python3 tools/plot_residual_locations.py runs/YOUR-QUADRATIC-RUN \
  --output runs/YOUR-QUADRATIC-RUN/plots
```
