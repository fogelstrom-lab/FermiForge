# Rotationally symmetric specular confinement in new_src

## First no-vortex test: radius 10

From the project root, run:

```sh
bash tools/run_cylinder_R10_no_vortex.sh smoke 10
caffeinate -i bash tools/run_cylinder_R10_no_vortex.sh relax 10
```

On Linux omit `caffeinate`. A separate angular-resolution check is:

```sh
bash tools/run_cylinder_R10_no_vortex.sh angular48 10
```

All three start from a uniform bulk B state (`istart=-1`), with zero winding,
T/Tc=0.30, F1s=0, and radius 10 xi_0, where xi_0 is the zero-temperature
coherence length. No vortex-core seed or restart is
used. Each gets a separate timestamped directory in runs; existing reference
data and new_src/qcv.inp are untouched. The matching T=0.30 Ozaki table is copied
by the runner. This radial path still uses tables, not the modern 2D startup
generator. The grid has 100 nodes, spacing 10/99, with reflected half-path
length 8000/99 (approximately 80.81).

The smoke run uses two azimuths and at most five NN startup updates plus one
Anderson update; it is only an execution/plotting check. The relaxation uses
32 azimuths, 11 polar nodes, all table poles, p_max=3, tolerance 2e-6, and up to
200 Anderson updates after startup. The angular48 case changes only the
azimuthal resolution. Reaching the iteration limit is not proof of convergence;
inspect the final residuals. The runner's completed status means process success.

The runner writes stdout.log (iteration history), stderr.log (MPI messages),
run_metadata.json, op_xyz, op_harm, curr, and Matplotlib plots under plots/.
Terminal output is not streamed each iteration: follow stdout.log with `tail -f`.
To regenerate plots during or after a run:

```sh
python3 tools/run_and_plot.py --solver new-src --plot-only runs/YOUR_RUN_DIRECTORY
```

For live plotting, wait until an update has finished: the legacy output files
are rewritten in place, so a read during a write may fail; retry afterwards.
Compare the wall-normal and tangential order-parameter components, regularity
at the centre, and the current diagnostic. At F1s=0 the stored current-related
mean field alone is not evidence of zero physical current. Do not interpret
the radial solver's density-like output as a computed superfluid-density tensor.
This is a symmetry-constrained B-like branch test, not a comparison of competing
confined phases or a validation of the still-unimplemented 2D cylinder boundary.
Angular agreement must be followed by radial-step and reflected-path checks.

## Modern framework: radial cylinder benchmark (2026-10-06)

The modern spinful point-map now accepts a specular cylinder geometry through
`mesh%cylinder`. The first supported field representation is the independent
radial ray; unconstrained Cartesian cylinder sampling is explicitly rejected
until an interior-only wall stencil is implemented. Existing free-vortex
drivers keep their previous default geometry and asymptotic behaviour.

### Bulk B or chiral bulk A, with no vortex

Select `--initialization -1` (default) for B, or `--initialization -2` for A:

```sh
caffeinate -i python3 tools/run_modern_cylinder.py --initialization -2 --ranks 10 --iterations 20
```

The Fortran namelist uses `initialization=-1` or `-2`, and `winding=0` in both
bulk cases. The additional `initialization=1` selects an A-phase-core vortex
and requires `winding=1`. Other combinations are rejected.
With imported profiles, choose the initialization matching the phase: it also
selects the appropriate radial symmetry class, not just a fresh seed.

The B seed is Delta times the identity. The A seed preserves the `new_src/bulkA`
choice A_zx=sqrt(2) Delta, A_zy=i sqrt(2) Delta, with all other entries zero:
d along z, orbital chirality p_x+i p_y, and no spatial phase winding. Delta is
the existing legacy B-gap scale used only to set the starting amplitude; this
is not a claim that the seed amplitude is the self-consistent bulk A gap.
Initial current fields are zero and the transport map can generate an
azimuthal current-related field during relaxation, including at F1s=0.
`amplitude_and_current_field.png` displays this field; its normalization is
not a calibrated physical mass-current density.

The radial reconstruction uses exp[i(m-s-k)phi]. Uniform B has symmetry label
m=0, whereas this uniform chiral A seed occupies C_(0,+) and needs m=1 to have
**zero Cartesian phase winding**. The code therefore separates the physical
winding (zero) from this intrinsic-chirality label, prints both, and uses the
same label for origin regularity throughout iteration. Reusing m=0 for A
would instead impose a spurious e^(-i phi) phase and suppress its centre.
The automated seed test checks spatial uniformity around a full circle and
preservation under origin projection. Do not equate the modern A physical
winding flag with an unchanged legacy harmonic reconstruction label when
setting up a new_src comparison. These radial classes constrain relaxation;
they do not explore arbitrary symmetry-breaking conversion between A and B.

Start with a single map from uniform bulk B, using all eight T=0.30 table poles:

```sh
python3 tools/run_modern_cylinder.py --ranks 10
```

Then try a short relaxation with the same physical parameters as the legacy
case (R=10 xi_0, no winding, T/Tc=0.30, F1s=0):

```sh
caffeinate -i python3 tools/run_modern_cylinder.py --ranks 10 --iterations 20
```

Each invocation builds the modern target, archives the executable and inputs
in a unique runs directory, streams iteration lines, captures MPI messages in
mpi.log, and makes Cartesian/harmonic profile and convergence PNGs. Use
`--initial-run runs/LEGACY_CYLINDER_RUN` to start from its op_xyz and curr instead
of bulk B, or pass a modern run directory to use its cylinder_final_op_xyz and
cylinder_final_curr. Anderson history starts fresh; original outputs are not
modified. The input must have the same uniform radial nodes and radius; the
driver checks coordinates and row counts. It does not infer temperature or
F1s from imported data: the convenience runner fixes T=0.30 and accepts
`--fs1` (default zero).
The lower-level namelist supports both parameters; its Ozaki table must match T.

For inspection during a run, after a completed update:

```sh
python3 tools/run_modern_cylinder.py --plot-only runs/YOUR_MODERN_CYLINDER_RUN
```

Output `cylinder_final_op_xyz` / `cylinder_final_curr` is the last **evaluated**
state, and `cylinder_mapped_*` its self-consistency map. The profile plot also
shows map(final), so a zero-update test is not mistaken for a relaxation.
`cylinder_history.dat` records absolute RMS, absolute maximum, relative L2,
normalization error and map time. These are not the legacy scaled maximum
residual. `--iterations N` means at most N AA updates followed by evaluation
of the final state (N+1 maps). There is no legacy five-update NN startup here.
The stopping criterion is absolute maximum residual <=2e-6 by default, with
p_max=3. Set `--tolerance` to change this absolute threshold.
An iteration limit is reported separately from convergence.

For a tighter continuation on the unchanged default B-cylinder configuration:

```sh
python3 tools/run_modern_cylinder.py --initial-run runs/YOUR_MODERN_CYLINDER_RUN \
  --tolerance 4e-9 --iterations 50 --ranks 10
```

Reapply any nondefault radius, resolution, trajectory or initialization options:
the restart imports fields, not the parent's run configuration. For A, retain
`--initialization -2`. The modern absolute tolerance must not be compared directly
with the legacy scaled maximum, which is N*max(abs(residual))/norm(state).

### Chiral-current reference: Sauls (2011)

J. A. Sauls, Physical Review B 84, 214509, provides surface-state and edge-current
benchmarks for chiral p-wave superfluids. Surface-normal pairing suppression,
specular-reflection continuity and spontaneous tangential current are relevant
to the A-cylinder branch. Its main analytic disk results assume a thin film,
a two-dimensional cylindrical Fermi surface and radius large compared with the
pair coherence length. Our spatially z-independent cylinder still integrates
over a three-dimensional momentum sphere; these are not quantitatively identical
models. The paper's temperature-dependent coherence length also differs from our
fixed zero-temperature length unit xi_0.

In particular, compare total current, not the bound-state contribution alone:
the continuum changes the zero-temperature angular momentum from the bound-state
estimate N*hbar to N*hbar/2 in the paper's specular large-disk model. Our stored
current-related field requires a normalization audit before interpreting it as
mass current or integrating angular momentum. The current B-cylinder comparison
is not a test of this chiral-current prediction.

### Reflections and controlled comparison

### A-phase-core vortex at R=20, F1s=5.4

From the project root on macOS:

```sh
caffeinate -i bash tools/run_cylinder_R20_a_core.sh
```

This starts from scratch at T/Tc=0.30, with unit circulation, the historical
`aop` seed in a B-phase background, 201 uniform radial nodes, step 0.1 xi_0,
reflected half-path 160 xi_0, 32 azimuths, 11 polar nodes, eight energy poles,
10 MPI ranks, AA p_max=3, and up to 300 updates. The absolute maximum residual
target is 2e-7. F1s=5.4 is converted to A1s=5.4/(1+5.4/3) internally.
The wall profile is relaxed, not clamped to bulk. This radial symmetry class
cannot develop a double core. The path length and discretization still require
refinement checks; this is an exploratory confined-vortex run.

Arguments appended to the script override defaults, e.g. `--iterations 50`.
For continuation, append `--initial-run runs/YOUR_MATCHING_RUN`; keep its grid
and physical settings unchanged. For live plots, use the existing `--plot-only`
command above after a completed iteration; retry if reading overlaps a write.
MPI messages go to mpi.log and profiles are saved after every map.

### Reflection settings

### Matched normal-core vortex

```sh
caffeinate -i bash tools/run_cylinder_R20_normal_core.sh
```

Uses the same settings as the R20 A-core script, changing only the seed to
historical `nop` (`initialization=0`, `winding=1`). It starts with a vanishing
order parameter at the centre and a phase-wound diagonal B background, without
the A-core filling components. Radial symmetry remains imposed, but no extra
projection is added to force the relaxed solution to remain normal-core.
Inspect the final harmonics to establish the branch reached. Outputs go to a
new dated `normal-core` directory; existing runs are untouched. The disk-plot
command below works for this case too. Omit caffeinate on Linux.

To reconstruct the saved radial solution over the disk, using the same angular
symmetry and blended-quadratic radial interpolation as the solver:

```sh
python3 tools/plot_cylinder_disk.py runs/YOUR_MODERN_CYLINDER_RUN
```

Seven PNGs are written to plots/disk: Cartesian and harmonic amplitude and
phase-cosine maps, density/current-related-field panels, and Cartesian and
harmonic profiles along both axes. The outside of the physical disk is masked.
The default 401-by-401 plotting grid adds no independent computed points.
The current-related field is not relabelled as calibrated physical mass current.

The default convenience-runner mode samples at uniform arc-length intervals
to compare with new_src. `--collision-aligned` adds both incoming and outgoing
wall limits at the same path coordinate; Riccati coherence is carried across
the zero-length interval without a reset. No finite integration interval then
straddles a momentum discontinuity. The two modes need not give identical maps
at a finite step. Compare them under step refinement:

```sh
python3 tools/run_modern_cylinder.py --collision-aligned --step 0.05 --iterations 0
```

`--half-length` controls the reflected propagation length on each side, not
the cylinder radius. Default step=10/99 and half-length=8000/99 match new_src.
The finite-path initialization remains bulk-coherence plus constant-endpoint
relaxation; there is **no asymptotic extrapolation or exterior reservoir**.
Length/step and initialization independence still need physical validation.
An exact wall target is evaluated at the same small interior displacement as
new_src; this convention also needs a refinement check.

The standalone check `python3 tests/check_modern_cylinder.py --build-dir
work/modern-cylinder-build` compiles an independent new_src oracle. It uses
two azimuths, 11 polar nodes and only one pole: a code check, not a physical
benchmark. At R=10 the bulk-seed and nonuniform quadratic finite-F1s maps agree
within about 3e-12 in the gap in the release build (below 7e-11 with the bounds-
checked debug build). At R=4 the release-build centre differs by about 1.6e-6 while
other nodes agree within 3e-12. Exactly coincident wall samples have a side
ambiguity in uniform sampling; matching the first incoming-wall convention
reduced the centre discrepancy from 5.8e-5. The centre difference also changes
with compiler optimization (about 4e-8 in the debug build). Do not call this bitwise agreement
or relax the physical convergence criterion to hide it. Collision-aligned
geometry is tested separately for specularity, position continuity, conserved
axial momentum, axial rays and the wall limit.

Next gate: compare **converged** radial cylinder profiles, then develop and
validate a wall-safe sampler before releasing the radial constraint in 2D.

The active radial solver now selects free-vortex sampling for `icyl=0` and
reflected cylinder sampling for `icyl=1`. This restores the legacy contained
transport path as a bounded first step toward cylindrical confinement in the
modern framework. Unconstrained 2D wall sampling is not yet supported.

## Input

The legacy ten-record `new_src` layout is still accepted. For `icyl=1`,
append an eleventh record containing the positive finite cylinder radius after
AA `p_max`. For example:

```text
0.30                  : T/Tc
1.0                   : vorticity
1                     : icyl
5.40                  : F1s
2                     : azimuthal directions (smoke only)
1.e-2                 : residual tolerance
1                     : A-core seed
3                     : Anderson mixing
0                     : main iteration limit
1.0                   : AA p_max
4.0                   : cylinder radius
```

All lengths are expressed in units of the zero-temperature coherence length
`xi_0`, as confirmed by Mikael on 2026-10-06. Thus radius `10.0` means
`10 xi_0`; no temperature-dependent rescaling is applied. This syntax differs
from the older fixed-form solver input.
Free-vortex inputs need no added record.

The Mac's explicit-grid layout is also accepted: insert `gridM` (free radial
extent) and `Rx` (free tangent-grid scale) immediately after `icyl`. This gives
12 records for a free vortex and 13 for a cylinder, with the cylinder radius
still last. For a cylinder the explicit wall radius controls the uniform grid;
the two free-grid settings do not change the contained geometry. With the
legacy layout the free-grid defaults are `gridM=70`, `Rx=10`. Free inputs may
omit the last AA parameter (9 or 11 records), defaulting to `p_max=1`.
Blank lines and full-line `!` or `#` comments are ignored. The Python runner
and Fortran reader use the same record-count convention; do not mix layouts
or append additional numeric records. Input is read through EOF, so use file
redirection for direct executable runs.

The supplied file is `examples/radial_specular_cylinder_smoke.inp`. It is not a
converged benchmark. In particular, `qcv` still runs up to five startup NN
updates even when its main iteration limit is zero. To exercise the complete
runner when that short calculation is desired:

```sh
python3 tools/run_and_plot.py --solver new-src \
  --input examples/radial_specular_cylinder_smoke.inp \
  --case-name radial-specular-cylinder --ranks 1
```

## Transport and interpolation

For confined calculations, the 100 radial nodes uniformly cover `0 <= r <= R`.
The integration step is `R/99`; the existing 800 samples on each half-ray give
a half-path length of approximately `8.08 R`. Free calculations retain their
original tangent-mapped grid, integration step, and asymptotic tail.

`intord_c` calls the module in `3DFS_MPICodes/trajectories.f90` to construct both
half-rays, including specular side-wall reflections. `poretraj.f90` is an
interactive geometry demonstration, not a library routine used by the solver.
The sampled momentum changes at reflections; its axial component is conserved
apart from the legacy epsilon handling of exactly transverse rays. Backward
sample momentum follows the legacy orientation toward the target.

At each sample the existing axial harmonic reconstruction and quadratic
interpolation supply all nine order-parameter components and the azimuthal
current-related mean field. Their contractions use the local reflected
momentum. No free-vortex extrapolation is used. The final three radial nodes
provide the wall interpolation stencil, avoiding an out-of-range node.

The Riccati calculation follows the reflected path without restarting the
coherence amplitudes at a wall. This is the source-compatible representation
of spin-independent specular continuity. The order parameter itself is not
forced to retain its bulk value at the wall; its surface response must emerge
from the self-consistency calculation.

## Checks and current limits

Run the bounded standalone checks:

```sh
python3 tests/check_cylindrical_trajectories.py
python3 tests/check_radial_cylinder.py
python3 tests/check_radial_input_formats.py
```

The geometry check compares the free-form routines with both original sources.
The radial check builds with bounds checks and floating-point exception traps,
compares reflected samples with an analytic quadratic gap and linear azimuthal
mean field, checks exact wall interpolation, and evaluates finite maps at two
radii using two azimuths, eleven polar nodes, and one energy pole. It also
rejects missing/invalid radii and a restart on an unrelated radial grid. These
standalone checks are separate from the existing CMake/CTest suite.

During implementation, free-vortex sampled self energies and one-pole
propagator contributions at four radial nodes were compared byte for byte with
an isolated pre-change build and agreed exactly. The complete radial executable
was built in an isolated directory. No production confined solution has yet
been accepted.

Before treating confined output as a physics benchmark:

- Check finite reflected-path length and endpoint-initialization independence.
  The original bulk-coherence start and constant-endpoint relaxation remain;
  they are initialization devices, not reservoirs attached to the cylinder.
- Check integration-step convergence at reflection discontinuities. The current
  routine samples at fixed arc-length intervals, and Riccati interpolation can
  straddle a reflection; collision-aligned intervals are a future refinement.
- Check the wall-node interior displacement (`1e-7` integration steps), grazing
  rays, angular/energy quadrature, radial resolution, and restart sensitivity.
- Compare a converged axisymmetric cylinder against the legacy contained
  reference, including surface suppression and current-related fields.

Cylinder restarts must already use the same 100-node uniform radial grid and
radius. Coordinate checks allow for the existing three-decimal text output;
changing the radius requires explicit field transfer, not relabelling old data.

After those gates, the 2D extension can reuse the reflected geometry with the
modern Cartesian field sampler, a circular physical-domain mask, explicit wall
stencils, and a cylinder boundary option replacing the free-vortex endpoint.
The full three-dimensional momentum quadrature remains necessary even when the
fields are translationally invariant along the cylinder axis.
