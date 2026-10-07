# Symmetry-constrained radial computation

`spatial_mode = 'radial_symmetry'` in the axisymmetric benchmark namelist
selects independent nodes on the non-negative x ray, including the origin.
The default `full_2d` is unchanged. The same spinful point map, quadratures,
Riccati propagation, MPI distribution and BB/Anderson engines are used.

## First test

### Temperature and runtime Ozaki generation

The radial example now generates its poles at solver startup:

```fortran
ozaki_mode = 'generate'
temperature = 0.30       ! T/Tc
ozaki_cutoff = 50.0
```

The same controls are accepted by the full-2D axisymmetric benchmark driver
and the double-core driver. Each rank generates the small table locally;
rank zero saves the actual table as `metrics.txt.ozaki.dat` beside the run
metrics. In generate mode `ozaki_file` is ignored and need not exist.
The radial runner also accepts `--temperature 0.30 --ozaki-cutoff 50`;
specifying --temperature enables generate mode even for an archived input.
Omitting these options retains the input's settings.

The generator is an explicit-real64 adaptation of the user-supplied colleague's
`/Users/mikael/Projects/Ozaki/ozaki.f90` and the conversion in `setpoles.f90`.
It retains the 200-entry cutoff table, selects the first last-pole value
strictly greater than cutoff/temperature, solves the same generalized
eigenproblem using LAPACK DSYGV, and scales poles by T/(2*pi). It does not
copy the interactive program's diagnostic Matsubara sums: Delta was used only
in those diagnostics, not in pole generation. The equality-at-maximum cutoff
edge case is now rejected rather than allowing count 201. Invalid or unsupported
temperatures/cutoffs fail explicitly. T=0 is not supported.

CMake now discovers LAPACK (Apple Accelerate on the tested Mac; site-provided
LAPACK/OpenBLAS/MKL on clusters). No Intel-specific compiler flags or paths are
required. Upstream authorship/license was not supplied; establish attribution
and redistribution permission before an external release of the adapted code.
The original files outside the project are unchanged.

Existing namelists default to `ozaki_mode='file'` for reproducibility.
To explicitly use a precomputed table, choose file mode and omit temperature
(its sentinel default is -1). A temperature supplied in file mode is rejected
to prevent silently using the wrong table.

### Automatic bulk gap

`bulk_gap_mode='auto'` (default) computes the bulk gap from the current Ozaki
quadrature when `ozaki_mode='generate'`. With old file-based inputs it preserves
the explicitly supplied `asymptotic_bulk_gap`, so historical runs are unchanged.
The radial example now computes both poles and gap from temperature at startup.
The resolved gap is printed and recorded as `bulk_gap` in metrics, and is set
before seed construction or asymptotic endpoint configuration in both drivers.

Explicit modes, also available through the radial runner's `--bulk-gap-mode`:

- `ozaki`: solve the isotropic bulk equation with exactly the vortex solver's
  current poles, weights and gap prefactor (also works with a saved table).
- `legacy`: reproduce the zero-flow `irep=1` equation in
  `new_src/bulkgap.f90`, with indices 0 through int(15/T+0.00001).
- `manual`: use `asymptotic_bulk_gap` unchanged.

The modern scalar solver brackets the positive root and bisects, avoiding
finite-difference Newton derivatives, mutable module state and unnecessary
angular integration. For irep=1 the normalized angular gap magnitude is unity.
This integration is for the isotropic B-phase bulk surrounding the vortex;
it does not import the separate irep=0 d-wave model, strong-coupling corrections,
or a general material model. Units remain Delta/(2*pi*kB*Tc).

At T/Tc=0.30, cutoff=50:
legacy gap = 0.279946155719; Ozaki-consistent gap = 0.279894824068.
The small difference reflects the energy-summation prescriptions and is not
silently treated as equality. The equation solved is
log(T) + T sum_j w_j [1/e_j - 1/sqrt(e_j^2+Delta^2)] = 0.
Vortex runs require 0 < T/Tc < 1 and a positive gap.

Reference profiles and restart fields are not regenerated or rescaled;
comparison to a reference at another temperature is not a validation.
Specify the tabulated Landau parameter directly as `fs1 = 5.4`, or pass
`--fs1 5.4` to the radial runner. Both modern drivers internally compute
`feedback_scale=Fs1/(1+Fs1/3)`. Explicit Fs1 takes precedence if a legacy
feedback_scale is also present; omitting Fs1 preserves old inputs. The radial
example now uses Fs1 rather than feedback_scale. Current drivers support
nonnegative Fs1. Both supplied Fs1 and resolved feedback are recorded in metrics.

### Matched BB scratch runs: 0plus and plus0

Launch the pair sequentially, using ten ranks per run:

```sh
bash /Users/mikael/Documents/Codex/3he-vortex-modernization/tools/run_radial_bb_seeds.sh
```

Pass `0plus` or `plus0` as the final argument to run just one seed.
Inputs are `examples/radial_bb_localized_0plus.nml` and
`examples/radial_bb_localized_plus0.nml`. Both use T=.30, Fs1=5.4,
60 independent radial points to 60 xi0, fitting at 45/55 xi0 and endpoints
at 120 xi0, matching the preceding radial benchmark.

The seed is a phase-wound normal-core background plus just one selected
harmonic with a compact smooth envelope of radius 5 xi0 and centre amplitude
equal to the computed bulk gap. No converged reference field or checkpoint is
used for initialization; reference data remain only for comparison plots.
BB uses initial/min/max step 1/.001/100, absolute curvature, no growth limit,
maximum-residual tolerance 2e-6 and a ceiling of 600 iterations. Every iteration
is checkpointed. Reaching the ceiling is not convergence.

The launcher uses a separate `work/radial-bb-seeds-build`, new dated run folders
and stops on a command failure. It waits for one complete run before starting
the other. It does not detect unrelated jobs: avoid starting it alongside
another ten-rank run. Inspect progress with the existing running-plot tool
using the exact run directory (or --latest radial-bb-scratch-0plus-T030-Fs54).

### Signed BB iteration diagnostics

Opt in with `save_iteration_diagnostics = .true.` or the runner flag
`--save-iteration-diagnostics`. This currently requires radial symmetry and
BB; other combinations fail explicitly. The default is off, and recording
does not select a different step, reject updates or add map evaluations.

For a fresh matched diagnostic run (wait until existing laptop jobs finish):

```sh
bash /Users/mikael/Documents/Codex/3he-vortex-modernization/tools/run_radial_bb_diagnostics.sh 0plus
# Or substitute plus0.
```

This uses the same T=.30, Fs1=5.4, 60-point scratch inputs and BB controls,
but a separate diagnostic build and a new dated run directory. Existing running
executables cannot acquire the new recording mid-run. Old scalar snapshots
cannot reconstruct the missing signed components retrospectively.

New files beside metrics:

- `metrics.txt.signed.layout.dat`: radial point ordering, mesh indices, radii
  and regions. There are 21 component-major blocks: real then imaginary part
  for each Cartesian A(spin,orbital), spin outer/orbital inner, then three
  mean-field components. Each block lists all N independent radial points.
- `metrics.txt.signed.mapNNNNNN.dat`: packed index, X_k, signed residual
  R_k=F(X_k)-X_k, and actual applied update X_next-X_k. These are precisely
  the independent active entries after the radial map projection; dependent
  halo/reconstructed plane entries are not duplicated. Files include a
  converged final evaluation, where the applied update is zero.
- `metrics.txt.bb_diagnostics.dat`: one row per update/evaluation, selected
  alpha, formula/status/clipping, validity flags, signed secant cosine,
  residual alignment, norms of s=X_k-X_previous and y=R_previous-R_k,
  raw signed BB1=(s.s)/(s.y), BB2=(s.y)/(y.y), and three regional RMS/max norms.

Candidates are before clipping and before the absolute-curvature convention.
In absolute mode, their absolute values are the corresponding unclipped
proposals when valid; status still determines whether a safeguard overrides
them. Invalid/missing quantities are zero placeholders with validity flags;
never interpret first-row zeros as measured curvature. Secant diagnostics are
also evaluated on the converged row before the usual early return.

Regions use the existing fine_region_half_width and medium_region_half_width:
core r<=fine, transition fine<r<=medium, outer r>medium. Norms are unweighted
over all 21 real entries per point, not spatial area integrals. An empty zone
has zero norm and can be identified from the layout. Each file is closed after
writing for live access. MPI rank zero alone writes diagnostics.

At 60 points the extra signed text data are approximately 64 MB for 600
iterations. The data enable signed mode correlations and slow/fast subspace
analysis; no adaptive or learned accelerator is introduced by this change.

### Launch

From any directory:

```sh
bash /Users/mikael/Documents/Codex/3he-vortex-modernization/tools/run_radial_symmetry.sh
```

This uses `examples/radial_symmetry_a_phase.nml`: T=0.30, F1s=5.4,
60 independent radial nodes from 0 to 60 xi0, a smooth sinh mesh,
quadratic radial interpolation, evolving fits at 45 and 55 xi0 and
trajectory endpoint radius 120 xi0. It starts from the preserved A-phase
reference, currently allows up to 400 BB updates on 10 ranks, and saves every update.
Only 60 maps are evaluated per iteration (1260 real mixing entries), not
the 14161 nodes in the reconstructed plotting grid. Runtime still depends
on trajectory lengths and quadratures; no fixed speedup is guaranteed.

To start from a seed and test an accelerator:

```sh
bash tools/run_radial_symmetry.sh --initialization-mode historical_aop --max-iterations 100
bash tools/run_radial_symmetry.sh --initialization-mode historical_aop --iteration-method anderson --max-iterations 100
```

Commands above with relative paths run from the project root. For normal-core
tests, use `--initialization-mode historical_nop`, but change the two reference
field paths in a copied input to the normal-core reference if you want meaningful
radial-reference comparison plots. The runner accepts that input via `--input`.
All runs have new dated output directories. No production run is started by
the implementation/tests. CongAcc is not yet an available engine.

Live plotting from the project root:

```sh
python3 tools/plot_running_axisymmetric.py --latest a-phase-radial-symmetry
```

The usual field and harmonic axis profiles are reconstructed over the plane.
Residual snapshots contain only independently evaluated ray points; their
spatial plots become radial profiles, not fictitious independently measured
2D residual maps. Scalar norms count only ray values. Area-weighted diagnostics
use annular areas rather than the Cartesian strip occupied by the ray.

## Representation and constraints

On a ray, each point carries nine complex tensor components and three real
mean-field components. At angle phi the harmonic amplitudes are multiplied by
exp[i(n-k-s)phi], then transformed back using the corrected physical tensor
convention. The planar mean field rotates as a vector and its z component
is unchanged. This does not impose the obsolete Axy=-Ayx transport projection.
Nonzero-winding harmonics are projected to zero at the origin; planar mean
fields are also zero there. Radial and axial vector components away from the
origin are retained, not artificially suppressed.

Trajectory sampling interpolates directly along the ray. Linear or continuous
blended-quadratic interpolation follows trajectory_interpolation_order. It
does not interpolate the reconstructed Cartesian cache. Outside the last
independent radius it uses the existing evolving state_asymptotic policy;
coefficients are sampled from the current ray, not frozen reference data.

The ray is taken from positive x mesh nodes up to active_radius. For an even
N-cell mesh with active_radius=half_width, the count is N/2+1; hence N=118
gives 60 points. With a smaller active_radius the last included node may be
slightly inside it. The matching radius must lie within that ray.

First version restrictions: centered unit-winding free vortex, positive
active radius, state_asymptotic endpoint policy, no 2D perturbation. Explicit
double-core, quadrupole and localized negative-index seeds are rejected.
The model cannot assess nonaxisymmetric stability. A restart deliberately
takes the +x ray and discards off-ray information; do not interpret this as
continuing an unconstrained double-core solution.

Full 2D storage/I/O and MPI array communication are retained for compatibility;
the expensive map count and accelerator vector, not all memory/communication,
are reduced. Coarse radial tests still need mesh/integration convergence and
comparison with new_src. Matching the radial symmetry does not make the two
discretizations or endpoint prescriptions identical.

## Verification

`test_radial_symmetry` checks polynomial interpolation, harmonic rotation,
agreement with the existing radial embedding, outside detection and origin
regularity. `tests/check_radial_symmetry.py BUILD_DIR` exercises two BB/AA
updates and compares one- and two-rank final fields. Existing full-2D tests
remain applicable. Use a separate build while an overnight calculation runs.
