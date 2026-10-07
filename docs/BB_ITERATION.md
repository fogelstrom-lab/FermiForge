# BB iteration for the full 2D solver

The general 2D driver used by run_axisymmetric_core_benchmark.py now supports
an independent Barzilai–Borwein engine. No changes were made to the legacy
radial BB implementation. Existing inputs still default to Anderson.

## Controlled comparison

From the project directory, run:

```sh
bash tools/compare_0plus_bb.sh
```

The launcher uses the same preserved starting field as the AA restart-only
comparison, a separate build directory (work/bb-iteration-build), a new
timestamped run directory, 10 MPI ranks, up to 60 updates, and checkpoints
every update. It does not modify the running Anderson executable or its files.
Do not launch it concurrently on the laptop if comparable timings are desired.
The physical parameters, grid and quadratures remain those of the reference input.

As of 2026-09-28, this launcher uses **radial-style BB settings**: initial
alpha=1, minimum=0.001, maximum=100, absolute curvature, and growth damping
disabled (bb_growth_limit=0). The generic solver defaults below remain
conservative. The run suffix is bb-radial-style, so earlier bb-positive
results remain separate. This matches the radial BB parameters, not the
whole radial workflow: no five-step NN warmup is added, and the initial field
is the preserved 2D restart rather than a fresh radial seed.

For other inputs, use the normal runner with `--iteration-method bb
--anderson-cycle 0`. To start from a different saved field, supply `--restart`.
As for AA, restarts save physical fields but start fresh iterator history.
The BB option is in the current general 2D driver; it is not added to the
separate historical double-core comparison driver.

## Update and safeguards

Let r = F(X)-X, s = X_current-X_previous, and y = r_previous-r_current.
The accepted active vector contains 18 real order-parameter components and
3 real Fermi-liquid components per point. BB uses the actual vector displacement.
It alternates alpha_BB1 = (s dot s)/(s dot y) and
alpha_BB2 = (s dot y)/(y dot y), beginning with BB2 on update 2.
Then X_next = X + alpha*r. Only root keeps the BB history; the updated field
and common convergence report are broadcast through the existing MPI adapter.
BB does not allocate or accumulate an Anderson history.

Generic solver defaults (overridden by the comparison launcher):

| Control | Default | Meaning |
|---|---:|---|
| bb_initial_mixing | 0.1 | First step and ordinary fallback |
| bb_minimum_mixing | 0.001 | Minimum proposed step |
| bb_maximum_mixing | 5 | Maximum proposed step |
| bb_growth_limit | 2 | Residual-norm growth triggering damping |
| bb_curvature | positive | Require sufficiently positive secant curvature |

Each control has a matching command-line flag with hyphens instead of underscores.
For example `--bb-maximum-mixing 10`.
Setting `--bb-growth-limit 0` disables the growth test entirely without
division by zero. Values other than zero must exceed one.

Finite and tiny-denominator checks protect the secant calculation. Degenerate
or unsuitable curvature causes a fallback. When the current residual norm is
more than twice the previous norm, the next alpha is reduced to the smaller
of the initial value and half the previous step, respecting the lower bound.
**This is not a line search:** the current iterate is not rejected or rolled
back, and no additional map is evaluated to accept a proposed step.
It therefore does not guarantee monotonic convergence or stability.

The usual raw maximum-component residual decides convergence. Both fields,
active/frozen masks, state-derived asymptotics, halo refreshes, and checkpoint
output retain their existing meanings. Safeguards do not change the map.

To reproduce the earlier conservative launcher settings, use:

```sh
bash tools/compare_0plus_bb.sh --bb-curvature positive --bb-initial-mixing 0.1 --bb-maximum-mixing 5 --bb-growth-limit 2 --initialization bb-positive
```

The radial-style launcher takes the absolute secant curvature as in
iter_BB.f90. Numerical protections for degenerate/nonfinite calculations
remain; it is not a bitwise reproduction of the legacy code.

## Diagnostics and inspection

Terminal updates identify BB and print the selected formula, safeguard status
and whether the step was limited. History uses engine=2 and records p=alpha.
The retained legacy column names anderson_vector_norm/anderson_relative also
hold BB residual diagnostics for compatibility; they do not imply AA was used.

bb_status is 0 for an ordinary startup/secant/converged update, 3 for degenerate
displacements/residual differences, 4 for unsuitable curvature, 5 for excessive
residual growth, and 6 for a nonfinite secant proposal/update fallback.
bb_formula is 0 when no quotient was selected, otherwise 1 or 2.
bb_limited marks clipping of a finite quotient to the configured bounds.
The convergence plot labels BB and marks safeguarded updates.
Completed runs also write harmonic_axis_profiles.png next to axis_profiles.png,
with the same axes and starting-state comparison. The harmonics are complex local
components (no angular phase unwinding). When radial-comparison plots are
generated, harmonic_axis_profiles_radial_comparison.png is included as well;
the supplied radial reference is not necessarily a matching-parameter target.

Inspect checkpoints while running:

```sh
python3 tools/plot_running_axisymmetric.py --latest single-harmonic-0plus-bb --no-radial-comparison
```

Compare RMS/max residuals, core-size evolution and wall time with AA from the
same initial field. Good propagator normalization is necessary but does not
establish nonlinear convergence or the correct equilibrium core size.

## Verification

Unit tests cover startup, BB1/BB2, actual displacements, positive/absolute
curvature, constant residuals, growth damping, both step bounds, zero-vector
convergence and reset/reinitialization. Field-adapter tests cover both active
fields and frozen points. A bounded two-rank smoke test is:

```sh
python3 tests/check_iteration_cycle.py work/bb-iteration-build --method bb
```

It exercises four BB updates and the convergence plot, not physical accuracy.
