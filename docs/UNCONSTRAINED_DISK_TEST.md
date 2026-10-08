# First unconstrained 2D specular disk test (7 October 2026)

Experimental validation path, not yet grid-converged confined-vortex physics.
Fields depend on x,y, are invariant along z, and retain 3D momentum. No end caps.

## Launch

From the project root on macOS:

```sh
caffeinate -i bash tools/run_disk_2d_B_reference.sh
```

On Linux omit caffeinate and supply a copied radial-reference run directory
as the first argument (local runs are not in GitHub):

```sh
bash tools/run_disk_2d_B_reference.sh /path/to/radial-reference
```

Default reference: runs/261006-192650-696944-modern-cylinder-R10-n0-bulk-b.
It needs cylinder_final_op_xyz and cylinder_final_curr, with 100 uniform
radial nodes to radius10. Original reference files are untouched.

Settings: R=10 xi0, T/Tc=.30, Fs1=0, no vortex, 32 cells across each diameter,
797 independent disk points (1089 storage nodes), 32 azimuths, 11 polar nodes,
eight poles, step10/99, half-path8000/99, uniform collision sampling, 10 ranks,
AA pmax3/history10, at most five updates, absolute maximum tolerance2e-7.
The final saved state is evaluated. This coarse spatial grid makes the first
comparison economical; agreement to the radial stopping tolerance is not expected.

## What is unconstrained

The radial profile is reconstructed only once for initialization. Thereafter
all disk nodes enter the independent mixing vector, without angular or origin
projection. Exterior storage is passive and never sampled. Reflection occurs
at the exact circle, not a staircase boundary. Append --perturbation 0.001 for
a separate symmetry-breaking experiment; it adds gap*amplitude*(x/R)*(y/R)*
exp(-4*r²/R²) to Axx initially only. First use the unperturbed benchmark.

## Wall sampling

The first backend requires uniform square Cartesian cells. Interior nodes
within 3.5 spacings supply a six-term, total-degree-two moving least-squares
fit. Compact inverse-distance weights, capped near nodes for conditioning,
and reorthogonalized QR determine the sample. Exact node queries return the
stored value. No exterior ghost values, radial continuation or bulk wall clamp
are used. This differs from the free-vortex tensor-product quadratic sampler.
One-sided fits can overshoot; polynomial reproduction is not a physical error
estimate. Grid/support/trajectory refinement remains necessary. Deficient
stencils fail explicitly instead of falling back to symmetry.

## Inspect and restart

```sh
python3 tools/run_modern_cylinder.py --plot-only runs/YOUR_FULL2D_RUN
```

Plots include Cartesian/harmonic amplitude and phase-cosine maps, density and
current-related field, convergence, initial/final axis profiles and RMS/max
residual maps. disk_comparison.json records gap change from the initial state
and the maximum-residual location. Initial/final difference is a reference
comparison only if that initial state was the intended reference. It is not
the self-consistency residual. Exterior points are masked in plots.

Every map saves cylinder_final_fields_2d.dat and cylinder_mapped_fields_2d.dat;
cylinder_initial_fields_2d.dat remains intact. Files retain full rectangular
storage. Use --restart-2d PATH/TO/cylinder_final_fields_2d.dat with the same
grid/physics instead of --initial-run. AA history starts fresh; no remeshing.
Do not use plot_cylinder_disk.py here: that tool reconstructs radial fields.

## Verification and next gates

Manufactured tests cover quadratic reproduction, near-node and exact-wall
queries, exterior NaN poisoning, and reflected point-map agreement with the
radial polynomial solution at centre/off-axis/wall targets under both collision
policies. A perturbed one-update smoke is identical on one and two MPI ranks.
A 197-point/four-azimuth radial-reference smoke completed on ten ranks; this
is execution evidence, not angular or spatial convergence. Next: the user-run
B-cylinder benchmark above, refinement, then normal/A-core and physical
symmetry-breaking tests. The existing radial cylinder path remains available.
# Graded annular option (8 October 2026)

## Exploratory A texture seeds

Run one of these, not all concurrently on the laptop:

```bash
caffeinate -i bash tools/run_annular_A_texture.sh a_mermin_ho
caffeinate -i bash tools/run_annular_A_texture.sh a_panam
caffeinate -i bash tools/run_annular_A_texture.sh a_planar
```

Defaults: R10, T0.30, Fs1=0, 24 core/wall rings, stretch1.5, arc target1 xi0,
10 ranks,32 azimuths,20 AA updates,p_max3. Extra arguments override runner
options, e.g. `--iterations 5`. The seeds are evaluated independently at nodes;
no texture or symmetry restriction is imposed in relaxation.

`a_mermin_ho`: pure A, d=z, l tilting from axial at centre to radial at wall.
It carries one wall circulation quantum, despite zero *additional* imposed
winding in the namelist/run-name n0. `a_panam`: same orbital texture with
d=(cos eta,sin eta,0), eta=-pi*x*y/(2R²), a schematic hyperbolic-like spin trial,
not the minimized Takagi Fig.2 or an orbital Pan-Am texture. `a_planar`:
literal short-note formula with fixed d=z; mixed A/polar trial, not the
conventional planar phase. All use the old A seed's global amplitude scale
sqrt(2)*legacy_B_gap(T), not an independently solved A bulk gap.

`orbital_texture.png` shows normalized orbital chirality (unit l only for pure
A). `spin_texture.png` shows a principal spin director and spin rank-one
fraction; arrows are shown only above0.99. Usual matrix/harmonic plots remain.
The code does not yet include the paper's dipole/Zeeman/rotating-frame terms.
No GL equilibrium/stability comparison is implied. Restarts require
`--texture none --restart-2d FILE`, with otherwise identical grid settings.

## Normal-core vortex: resolve centre and surface together

```bash
caffeinate -i bash tools/run_annular_disk_normal_core.sh
```

This starts from the converged radial normal-core run
`runs/261007-103929-261905-modern-cylinder-R20-n1-normal-core-Fs5.4`.
Physics and transport settings match that reference: T/Tc=0.30, Fs1=5.4,
R=20 xi0, winding=1, 201 input radial points, trajectory step0.1 and
half-length160, 32 angular azimuths. Ten ranks, five updates, p_max3.

The new `--disk-radial-layout core_wall` uses a symmetric tanh radial map
with 32 rings, stretch2, target arc spacing1.25 xi0 and at least24 points
per ring. There are1877 independent nodes, seven rings inside radius2 xi0;
first/last radial steps0.104 xi0 and maximum step1.29 xi0 at half-radius.
The existing default `wall` layout is unchanged. These are trial spacings,
not a certified resolution. The solver does not constrain the core to stay
normal. For a standalone historical seed omit `--initial-run` when invoking
`run_modern_cylinder.py` directly and select `--initialization 0`.

The launcher accepts an alternative radial reference as its first positional
argument, with additional runner overrides afterwards. Output includes the
native grid and initial/final Cartesian and harmonic axis profiles.

## Vortex-free B-phase reference

The Cartesian mode described below remains available and is still the default.
For the first annular radial-reference comparison, from the repository root:

```bash
caffeinate -i bash tools/run_annular_disk_B_reference.sh
```

The first optional positional argument is an alternative radial reference run
directory. Further arguments override runner controls, for example:

```bash
bash tools/run_annular_disk_B_reference.sh runs/261006-192650-696944-modern-cylinder-R10-n0-bulk-b --disk-rings 24 --iterations 5
```

Defaults: T/Tc=0.30, Fs1=0, R=10 xi0, no vortex, 10 MPI ranks, p_max=3,
32 momentum azimuths, unchanged trajectory step and length, five AA updates.
The grid has 16 radial intervals, wall stretch 2 and target arc spacing 1.25 xi0:
585 independent nodes; radial steps range from 1.359 to 0.2084 xi0.
This is a trial B-phase wall grid, not a validated vortex grid.

Use `--disk-grid annular` with `--spatial-mode full_2d` for standalone runs.
`--disk-rings` controls radial intervals, `--disk-stretch` (0 to 4) clusters
radial points at the wall, and `--disk-tangent-spacing` sets target arc spacing.
`--radial-points` still controls the input reference ray, not the independent
annular grid. `--disk-cells` applies only to Cartesian mode.

Each point carries independent fields. There is one origin node and an exact
wall ring. Seven neighbouring rings and up to nine angular neighbours per ring
support the quadratic QR sampler, scaled separately in radial/tangential
directions. Native node checkpoints support same-grid restarts. The existing
runner's `--plot-only RUN_DIRECTORY` produces grid, Cartesian/harmonic maps,
axis profiles, current-related mean fields, density, residuals and convergence.
Display triangulation is not solver interpolation. Generic rectangular-map
plotters must not be used directly on these unstructured checkpoint files.

Control refinement by changing only radial intervals first, then only arc
spacing. Do not simultaneously change trajectory or momentum quadrature.
RMS and Anderson products are currently node-weighted, not area-weighted;
changing the node distribution changes that norm. Maximum errors and profiles
must also be checked. The hybrid Cartesian-interior/annular-wall option and
automatic refinement are not yet implemented.
