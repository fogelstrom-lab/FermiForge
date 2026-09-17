# A-phase-core 2D benchmark

## Purpose

This is the second axisymmetric physics regression for the 2D FermiForge
path. It applies the same source-preserving radial embedding, spin-matrix
Riccati propagation, momentum and energy quadrature, MPI point distribution,
and residual diagnostics used for the normal core to a converged A-phase-core
radial state.

The benchmark has two levels:

1. an independent original-`new_src` point map at radial index 24, used as a
   strict transport and convention oracle; and
2. a complete map on the current uniform Cartesian grid, used to measure the
   resolution and domain-size requirements of the hard and soft core.

The full-grid residual is a convergence-study metric, not yet a production
tolerance. In particular, the uniform 17 by 17 baseline cannot be expected to
resolve the coherence-length structure and the extended soft core at once.

## Reference

The preserved reference is under
`reference/current_new_src_T0.30_F1s5.4/`. The user generated it by iterating
the radial `new_src` solver to full convergence, then performed a zero-iteration
restart check with `istart=3`. The preserved `restart_check.inp` documents that
reload check; it is not represented as the original convergence input.

The order-parameter and current files both contain 100 points on the current
radius-70 tangent grid and have identical radial coordinates. The independent
point-map fixture was regenerated from these named files with
`tests/reference/legacy_getnewop_point_driver.f90`.

## Run

From the project root:

```text
/opt/homebrew/bin/python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_a_phase_core_benchmark.nml \
  --case-name a-phase-core \
  --initialization converged-radial-reference \
  --ranks 10
```

The runner creates a timestamped archive containing the effective input, build
and run logs, hashes, metrics, both 2D fields, and accessible plots.

The first multiscale comparison is run with:

```text
/opt/homebrew/bin/python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_a_phase_core_multiscale.nml \
  --case-name a-phase-core-multiscale \
  --initialization converged-radial-reference \
  --ranks 10
```

This mesh has a fine core region, an intermediate soft-core region, a coarse
trajectory halo, and a circular active update mask. Residuals are reported by
zone so the soft core cannot be hidden by the larger outer field.

## Next convergence gate

Refine the intermediate soft-core zone in a controlled spacing/radius study,
then introduce dense block refinement with field transfer and Anderson-history
reset between mesh stages. The A-phase-core comparison must retain the same
physics controls and source conventions as the normal-core regression.
