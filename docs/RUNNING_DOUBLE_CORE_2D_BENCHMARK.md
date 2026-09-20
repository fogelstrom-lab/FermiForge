# Running the imported double-core one-map benchmark

This benchmark evaluates the modern MPI-distributed quasiclassical map on a
converged, nonaxisymmetric double-core state produced by the legacy full-2D
solver. It does not update the state. Its purpose is to measure the residual
between the imported field and one modern transport/self-consistency map
before attempting double-core relaxation.

## Why the first map is sampled

The imported state contains 601 by 601 points. A complete map would require
361,201 points times 704 momentum directions times eight Ozaki poles and is
not an appropriate first laptop test. The full field is nevertheless retained
on every MPI rank and used for interpolation along every trajectory.

The supplied input evaluates 27 target points:

- the center and 3 by 3 stencils around the two half-core centers at
  `y = +/-11.8 xi0`;
- a 3 by 3 outer set at `x,y = -40,0,+40 xi0`.

Only those target values are replaced by the map when the residual is formed.
The center is counted only once where the two masks overlap, leaving 19 core
points and eight distinct outer points. The central and half-core stencil
radii are independently configurable with `center_probe_stencil_radius` and
`probe_stencil_radius`.

## Boundary policy

The default benchmark now uses `state_asymptotic`. Before propagating any
trajectory it samples the imported state around a fitting circle and removes
the imposed vortex phase. The angular average is projected onto the general
B-phase bulk form

```text
A_mu_i = Delta exp[i (phi + phi0)] R_mu_i,  R in SO(3).
```

The fitted `phi0` and all nine entries of `R` are written to the metrics file.
They are not assumed to be zero and the identity. On each trajectory ray, the
remaining full complex 3 by 3 matrix is continued with the angle-dependent
form of Eqs. (19)--(20) of `dcvlong`: mixed z/in-plane components have the
`1/r` tail, while the in-plane block and zz component have the `1/r^2` tail.
The three Fermi-liquid mean-field components retain a general two-coefficient
`1/r + 1/r^2` fit. Thus the finite-radius spin-orbit texture is retained even
when the constant infinite-radius rotation is close to the identity.

For a control comparison, setting `endpoint_policy = 'local_supplied_field'`
ends trajectories at the supplied Cartesian field boundary.

## Run

From the project directory:

```text
mkdir -p work/double-core-one-map
mpirun --host localhost:10 --map-by ppr:10:node --bind-to none -np 10 \
  work/modern-build/benchmark_legacy_double_core_map_2d \
  examples/2d_double_core_one_map.nml
```

The executable accepts optional second and third arguments overriding the
metrics and sampled-residual output paths.

## Constrained iteration smoke test

`examples/2d_double_core_iteration_smoke.nml` exercises the complete repeated
map and legacy-Anderson path for two accepted updates. It evaluates the same
27 diagnostic points but admits only the 18 movable points in the two
half-core stencils to the Anderson vector. The center point and the complete
imported exterior remain fixed. This is an explicit temporary constraint:
the center suppresses translational drift, while the frozen exterior retains
the imported double-core orientation during this integration test.

```text
mkdir -p work/double-core-iteration-smoke
mpirun --host localhost:10 --map-by ppr:10:node --bind-to none -np 10 \
  work/strict-point-map-build/benchmark_legacy_double_core_map_2d \
  examples/2d_double_core_iteration_smoke.nml
```

Anderson history is now stored only for `update_point`. For this smoke test
that is 378 real values rather than the 7,585,221 values in the full imported
field. With history limit ten, the two main history matrices therefore occupy
about 65 KiB rather than about 1.24 GiB. Inactive values are copied exactly
from the current state after every update.

The first two updates gave:

```text
iteration  update maximum   mixing p   history
1          3.8069e-4        0.01       1
2          3.7932e-4        1.19       2
```

The sampled relative-L2 residual changed from `1.3820e-3` to `1.3763e-3`,
and maximum propagator-normalization error over both maps was `1.68e-15`.
The run deliberately reports `iteration_limit`, not convergence.

Set `final_field_file` to write the complete modern field state. A later run
can read it with `initialization_mode = 'restart'` and `restart_field_file`.
The restart is checked against the mesh inferred from the archived split
field before it is broadcast to all ranks.

The frozen exterior is not the final physical zero-mode treatment. Its role
is to verify repeated mapping, active-only Anderson storage, state broadcast,
restart output, and asymptotic-reference refresh before the active region is
expanded and translation/rotation tangent modes are projected explicitly.

The update region can now be expanded without changing the solver:

- `half_core_stencils` uses the two square diagnostic stencils;
- `half_core_disks` selects smooth disks of radius `update_radius` around the
  two half-core centers;
- `centered_disk` selects one disk about the vortex center for eventual
  inclusion of the complete hard and soft cores.

`examples/2d_double_core_patch_smoke.nml` is the first smooth expansion. Its
two disks have radius `0.4 xi0`, giving 26 movable points and 546 Anderson
values. A two-update ten-rank run completed with `p=0.01` and `p=1.23`. The
movable-region maximum residual changed from `4.330e-4` to `4.312e-4`, while
normalization stayed below `1.56e-15`. The sampled set contains 35 points and
the two maps took 31.9 and 34.3 seconds on the MacBook Pro.

This radius is intentionally small. Increasing `update_radius` provides a
measurable sequence from a machinery check to a hard-core calculation. A
centered disk reaching the asymptotic matching annulus is intended for a
cluster run, not as the next laptop test.

## Controlled active-radius ladder

`tools/run_double_core_radius_ladder.py` automates the next controlled step.
It runs strictly increasing paired-disk radii and carries the accepted values
from one radius into the next. The continuation file is sparse: it contains
only the coordinates and 21 real state values at the movable points. Each new
case first imports the unchanged archived 601 by 601 field and then overlays
the sparse state. The large reference field is therefore neither copied nor
silently modified between radii.

From the project directory, the short laptop check is:

```text
tools/run_double_core_radius_ladder.py \
  --no-build \
  --build-dir work/strict-point-map-build \
  --radii 0.4 0.6 \
  --iterations 2 \
  --ranks 10
```

Every radius receives its own input, log, metrics, iteration history, sampled
residuals, and sparse final state beneath one timestamped run directory. The
top level contains CSV and JSON summaries and, when Matplotlib is available,
a noninteractive PNG diagnostic. `--independent` disables continuation and is
the control for distinguishing a radius effect from a restart effect.

The verified continuation run gave:

```text
radius   movable points   Anderson values   imported sparse points   update max after map 2
0.4      26               546               0                        4.3117e-4
0.6      58               1218              26                       5.4919e-4
```

The 0.6 case therefore demonstrably consumed the complete 0.4 state before
adding the new annulus. Its core relative-L2 residual was `1.7714e-3`; the
maximum propagator-normalization error over the run was `1.56e-15`. Both cases
ended at the requested two-iteration limit, so these figures validate state
transfer and expansion, not physical convergence. The next solver change is
explicit projection of translational and orientation tangent modes before
larger hard- and soft-core radii are admitted.

## Template tangent-mode projection

The double-core update can now remove collective drift before the residual is
passed to Anderson mixing. The option is deliberately explicit and defaults
to `none`. The complete present constraint is selected with

```text
zero_mode_projection = 'gauge_translation_orientation'
```

For the initial state of each radius case, the solver constructs four tangent
vectors on the movable points:

```text
gauge:          (i A, 0)
translation x: (d_x A, d_x h)
translation y: (d_y A, d_y h)
orientation:   (i A + y d_x A - x d_y A,
                y d_x h - x d_y h)
```

Here `h` denotes the three current-related Fermi-liquid mean fields. The
orientation expression is the infinitesimal code-basis action
`exp(i theta) A(R(-theta) r)`. Its gauge factor leaves the winding-one bulk
field unchanged. It rotates the coordinate dependence but does not silently
apply an additional spin- or orbital-index rotation. If that physical action
is later required, it will be introduced as a separately named convention.

Spatial derivatives use three-point Lagrange differentiation on the general
rectilinear mesh, including nonuniform coordinates. The tangents are restricted
to the movable cells and twice orthonormalized in the same unweighted Euclidean
vector metric used by the current Anderson implementation. For residual `f`
and orthonormal template tangents `t_k`, the accelerator receives

```text
f_perp = f - sum_k t_k (t_k dot f).
```

The template remains fixed within a run. Thus the operation is a phase
condition, not a modification of the quasiclassical map. The history and
metrics retain both the unprojected sampled residual and the projected update
residual, together with all four removed coefficients.

The direct projected smoke input is
`examples/2d_double_core_projected_patch_smoke.nml`. The radius runner exposes
the same choice:

```text
tools/run_double_core_radius_ladder.py \
  --no-build \
  --build-dir work/strict-point-map-build \
  --radii 0.4 0.6 \
  --iterations 2 \
  --ranks 10 \
  --zero-mode-projection gauge_translation_orientation
```

In the verified `0.4 xi0` projected run all four modes were independent. Only
`1.46e-6` of the residual norm lay in their span, and the projected and
unprojected update maxima differed by `2.55e-10`. This is the expected control
result for the highly symmetric archived state: the projector is active but
does not manufacture an appreciable correction when drift is absent.

## First testable continuation version

The outer asymptotic halo can now be made a dependent part of the state rather
than a frozen copy of the imported solution. Select

```text
inactive_halo_policy = 'state_asymptotic'
```

together with `endpoint_policy = 'state_asymptotic'`. After every accepted
update, the root rank refits the bulk reference, matches the Eq. (20)-allowed
single power of every complex order-parameter component, and fits both radial
powers of all three Fermi-liquid mean fields. It then regenerates the dependent
halo and broadcasts it. Halo values do not enter the Anderson vector.

This deliberately leaves one constraint visible: points inside the matching
circle but outside the paired active disks are still the imported field. They
provide the source annulus for the present small-radius calculation. Thus this
is a runnable staged continuation and not yet a converged, unconstrained
double-core solution.

The shortest end-to-end test is:

```text
tools/run_first_double_core_2d.py
```

The launcher checks for a plotting-capable Python and, on a Homebrew or
MacPorts Mac, uses it automatically when the default Xcode Python has no
Matplotlib. On a cluster without Matplotlib the calculation and machine-readable
reports still complete; pass `--no-plot` to suppress the plot attempt.

If the strict build already exists, skip rebuilding it:

```text
tools/run_first_double_core_2d.py \
  --no-build \
  --build-dir work/strict-point-map-build
```

The default run uses ten MPI ranks, two updates at each of the paired-disk
radii `0.4` and `0.6 xi0`, all four tangent constraints, and the dependent
asymptotic halo. On the development MacBook Pro it takes approximately three
to four minutes. A successful run should show the movable set growing from 26
to 58 points, the second case applying 26 sparse restart points, 97,012
dependent halo points on the supplied 601 by 601 mesh, four retained tangent
modes, a propagator-normalization error near `1e-15`, and `"passed": true` in
the top-level JSON report.

The verified two-stage calculation gave:

```text
radius   update points   restart points   update max after map 2   tangent fraction
0.4      26              0                4.3117e-4                1.462e-6
0.6      58              26               5.4919e-4                7.269e-6
```

At radius `0.6`, the sample and core relative-L2 residuals were `1.6069e-3`
and `1.7714e-3`; the maximum normalization error was `1.67e-15`. The final
halo regeneration changed a halo value by at most `3.74e-7`. Both radii ended
at the requested two-iteration limit, as intended for this first test.

A stopped ladder can be resumed without repeating its smaller radii:

```text
tools/run_first_double_core_2d.py \
  --radii 0.8 \
  --initial-sparse-restart /absolute/path/to/r0p6/sparse_final_state.dat
```

For a more substantial calculation, increase both the radius sequence and the
iteration limit, for example `--radii 0.4 0.6 0.8 1.0 --iterations 20`. Such a
run still needs convergence and domain studies before it can be treated as a
new physical result.

## Dense-stencil result

With ten ranks, the historical 64 by 11 angular quadrature, eight Ozaki poles,
and a maximum trajectory interval of `0.25 xi0`, the asymptotic 27-point map
gave:

```text
sample relative L2 residual = 1.3820e-3
core relative L2 residual   = 1.7887e-3
outer relative L2 residual  = 3.2996e-6
maximum absolute residual   = 3.8069e-4
maximum normalization error = 1.5543e-15
maximum rank time           = 22.74 s
```

The largest point RMS occurs at `(x,y)=(0,-12.0)`, immediately beside one
half-core minimum. The fitted constant bulk reference has phase offset
`2.43e-8` and differs from the identity rotation by less than `4.6e-10`.
The normalized field on the `r=58 xi0` fitting circle nevertheless has RMS
deviation `8.90e-2` and maximum deviation `2.82e-1` from that constant bulk
matrix. This is the sizeable angle-dependent soft texture retained by the
ray-wise tail coefficients, not a constant far-field rotation.

Repeating the same dense stencil with `local_supplied_field` changes the
overall relative L2 residual by only `2.0e-14` and the outer relative L2
residual by `2.1e-11`. This verifies numerical continuity of the new endpoint
against the previous supplied-boundary calculation at the present probes.

The earlier single-point half-core check found a core relative L2 residual of
`1.4240e-3` at step `0.25 xi0`; refinement to `0.125 xi0` gave `1.4759e-3`.
A formal dense-stencil step-size convergence table remains to be completed.

The results are measurements, not an acceptance tolerance. The remaining
physics checks include the dense-stencil trajectory-step convergence table
and converged active-radius continuation after the tangent-mode constraints
are available.
