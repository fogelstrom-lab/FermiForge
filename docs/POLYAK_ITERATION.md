# Polyak momentum iteration

FermiForge now offers the fixed-step Polyak heavy-ball method described in
Holmvall et al., *SuperConga*, Applied Physics Reviews 10, 011317 (2023),
Appendix E, equations E4, E6 and E7, DOI 10.1063/5.0100324.
This is not the alternative Polyak step-size rule based on an objective value.

For the packed real vector of independent order-parameter and Fermi-liquid
components, let `r = F(x) - x`. The update is

```
v_next = (1 - drag) * v + step_size * r
x_next = x + v_next
```

Velocity starts at zero. Defaults `step_size=2` and `drag=0.5` match Listing 14
in the paper. A larger drag means LESS retained momentum. FermiForge permits
`0 < drag <= 1`; drag exactly one extends the paper's open interval to the
memoryless Picard limit. Step size must be finite and positive. These defaults
are a comparison point, not parameters established to be optimal for 3He.

## Run a matched radial test

From the project base directory:

```sh
bash tools/run_radial_polyak.sh 0plus
# Or, separately:
bash tools/run_radial_polyak.sh plus0
```

The launcher reuses the BB scratch input but overrides its engine. It keeps
T/Tc=0.30, Fs1=5.4, 60 independent radial nodes, the same physical/quadrature/
asymptotic controls, a 600-iteration limit, tolerance 2e-6, ten MPI ranks and
checkpoints every iteration. It creates a separate `work/radial-polyak-build`
and fresh dated run directory. Do not start a second ten-rank job while another
is using the laptop. This launcher follows the existing macOS caffeinate wrapper.

Overrides can be appended, for example:

```sh
bash tools/run_radial_polyak.sh 0plus --max-iterations 100 --polyak-step-size 1 --polyak-drag 0.5
```

The general runner also accepts these flags:

```
--iteration-method polyak --polyak-step-size 2 --polyak-drag 0.5 --anderson-cycle 0
```

Or set these values in `&axisymmetric_core_benchmark`:

```fortran
iteration_method = 'polyak'
polyak_step_size = 2.0
polyak_drag = 0.5
anderson_cycle_iterations = 0
```

Both radial_symmetry and full 2D modes of that executable use the same routine.
The separate historical double-core benchmark launcher is not extended in this
step. The generic serial and MPI state-update interfaces accept optional
`polyak=...`, so other drivers can reuse it without duplicating the algorithm.

## Interpretation and limits

- One expensive map per iteration; no line search, random fallback, clipping,
  adaptive step size or automatic momentum reset is hidden in this method.
- MPI rank zero holds the momentum vector and broadcasts the updated fields.
  Only independent active components enter it; the existing symmetry/halo
  reconstruction and asymptotic fitting remain unchanged.
- The convergence check is on `F(x)-x`, not on the momentum step. Once the
  maximum residual meets tolerance, the current state is retained and velocity
  cleared. The finite last update at an iteration limit is not itself a newly
  evaluated residual, as with the existing engines.
- History `engine=3` denotes Polyak in the modern runner; plots label it as
  such. `p` is the fixed residual multiplier (not the total effective step),
  and `md=0` because no Anderson subspace is used. Drag is stored in the input
  and metrics. Legacy new_src uses its own engine numbering, unchanged.
- Field-file restarts begin with zero velocity, so they are not uninterrupted
  momentum continuations. Momentum is not saved in field checkpoints.
- Standard convergence histories, residual maps and plots remain available.
  Signed state/residual/update recording and regional norms are available for
  every engine in radial mode using `--save-iteration-diagnostics`. This
  launcher defaults to recording off; append that option to turn it on.
  BB-specific step-estimate metadata remains BB-only.
- Momentum can accelerate a slow mode but can also amplify oscillations.
  There is no general stability guarantee for the nonlinear self-consistency
  map. Reduce the step size or increase drag if the residual grows persistently;
  compare total maps and elapsed time to the same tolerance from the same seed.

Implementation: `src/polyak_mixing.f90`. Regression coverage includes the exact
recurrence, startup/reset, convergence without momentum drift, the Picard limit,
linear-map convergence, and bounded one-/two-rank radial integration.
