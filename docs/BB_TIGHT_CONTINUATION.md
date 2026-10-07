# Tight-tolerance BB continuation

Run from any directory:

```sh
bash /Users/mikael/Documents/Codex/3he-vortex-modernization/tools/continue_0plus_bb_tight.sh
```

Restarts the converged 198-evaluation solution from the September 29
continuation. Requests up to 100 additional evaluations on 10 MPI ranks,
checkpointing every iteration, with absolute maximum-residual tolerance
2e-8 instead of 2e-6. All physical, mesh, trajectory, quadrature, asymptotic
and BB parameters are inherited unchanged. BB history starts fresh at alpha=1.
The executable has a separate build directory and outputs have a new dated
directory. Existing results are preserved. Combined history continues at 199
and records the restart; the first 198 rows used the older stopping tolerance.

Budget approximately 15 hours at the latest 547 seconds/evaluation. The
100-evaluation limit may be reached before the tighter tolerance: a further
hundredfold reduction could take around 160–200 evaluations if the recent
trend persists, and maximum residual need not follow relative RMS exactly.
Reaching iteration_limit is not evidence of a numerical floor.

Live inspection, from the project directory:

```sh
python3 tools/plot_running_axisymmetric.py --latest single-harmonic-0plus-bb-tight --no-radial-comparison
```

The standard final plots and combined history are produced on completion.
Inspect map_relative_l2 for the legacy-compatible relative RMS. map_rms is
absolute RMS; the terminal average is mean absolute residual.

## Interpretation and subsequent accuracy check

A fixed deterministic discretized map may have a residual far below its
error relative to the continuum equations. Interpolation/integration error
therefore does not by itself impose a nonlinear-residual floor. Roundoff,
unstable iteration or nonsmooth state-dependent operations can limit residual
reduction. Propagator normalization near 1e-15 tests an algebraic identity,
not integration or interpolation accuracy. No numerical floor is established
by the current data.

After this run, use the same saved field for controlled map evaluations:
halve trajectory integration steps, refine spatial sampling in a separate
test, increase angular/frequency quadrature separately, and test fitting/
endpoint radii separately. Compare changes in mapped fields at common physical
points using core/outer and maximum as well as global norms. Transfer between
meshes introduces its own interpolation error. If a refined map changes by
much more than the remaining residual, further iteration of the original map
does not improve continuum accuracy. Quantitative solution errors require
refined solves; map differences alone can understate slowly relaxing modes.
