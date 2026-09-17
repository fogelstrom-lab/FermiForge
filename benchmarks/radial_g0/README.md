# G0 radial free/cylinder regression

## Purpose

This directory defines the first validation gate in `docs/PROJECT_PLAN.md`.
It keeps two questions separate:

1. Does a run workflow preserve the legacy radial solver's behavior across MPI
   rank counts?
2. Does a modern radial implementation reproduce converged free-vortex and
   specular-cylinder reference solutions?

The first question now has smoke-test evidence. The second is not yet passed.

## Inputs

`inputs/free_smoke.inp` and `inputs/contained_smoke.inp` exercise both legacy
trajectory branches cheaply. They deliberately use a loose stopping tolerance,
four azimuthal directions, and one iteration. They are workflow tests, not
physics calculations.

`inputs/free_a_core_reference.inp` and
`inputs/contained_a_core_reference.inp` define matched production candidates at
`T/Tc=0.30`, `R=20`, `F1s=5.4`, with the A-core seed. A result is not accepted
merely because the process exits normally: both average and maximum residuals
must satisfy the recorded tolerance and the result must pass resolution checks.

## Evidence from 2026-09-15

The one-rank and ten-rank smoke outputs are exactly equal at the printed-file
precision for `op_xyz`, `op_harm`, and the normalized `curr` columns:

- `free_smoke_rank_comparison.json`;
- `contained_smoke_rank_comparison.json`.

The free and contained smoke outputs are intentionally different. The report
`free_vs_contained_smoke.json` gives relative L2 differences of approximately
`0.1284` for the order parameter and `0.0733` for the pair-amplitude/vector
field output. This confirms that the legacy `icyl` switch selects a materially
different trajectory/boundary calculation.

The run directories are recorded in those JSON reports under `runs/g0-smoke/`.
The unsuccessful socket-restricted attempts earlier on the same date remain
archived with `status=failed`; they are provenance, not benchmark evidence.

## Rejected baseline

`runs/20260914-150620-cylindrical` is useful historical output but is not an
accepted G0 baseline. Although the executable returned zero, the final legacy
AA stage ended at average error `9.1833e-4` and maximum error `1.9316e-3`, far
above the requested scaled tolerance printed as approximately `5.6e-7`.
Normal process completion is therefore not treated as convergence.

## Comparison command

```text
python3 tools/compare_radial_runs.py REFERENCE_RUN CANDIDATE_RUN \
  --atol 1e-10 --rtol 1e-8 --output comparison.json
```

The comparison normalizes the two `curr` layouts to pair amplitude plus
`vx,vy,vz`, then interpolates a candidate with a different radial grid onto the
reference grid. It does not align global gauge or spin-orbit orientation; those
must be fixed by the benchmark setup.

## Remaining work before G0 passes

- Obtain converged full free and contained legacy runs with matched numerical
  controls, including a continuation/restart schedule if the current 40-step
  stages are insufficient.
- Repeat both full cases with at least one and ten ranks and declare tolerances
  below the output formatting precision.
- Make `new_src` execute an explicit free/specular-cylinder boundary strategy.
  At present it reads `icyl`, but only calls `intord_v`, so its cylinder flag has
  no transport effect.
- Compare the modern result against the accepted legacy fields and selected
  integrated observables after grid, trajectory-step, angle, and energy
  convergence.
- Add a machine-readable terminal status that distinguishes `converged`,
  `iteration_limit`, and numerical failure.
