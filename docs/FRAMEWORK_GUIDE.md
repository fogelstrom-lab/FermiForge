# FermiForge: a hands-on map of the framework

An editable [LaTeX edition](technical_note/fermiforge_framework_guide.tex)
includes a vector flowchart and revision history. The two editions are
maintained together; they are not automatically synchronized.

Revision 0.5 (8 October 2026) adds the experimental **graded annular disk**:
independent 2D nodes, a wall-clustered radial layout, variable angular counts,
quadratic sampling and native-node plots. Start with
`tools/run_annular_disk_B_reference.sh`; see
[disk tests](UNCONSTRAINED_DISK_TEST.md) for resolution controls and limitations.
A Cartesian-interior/annular-wall hybrid and automatic refinement are future
steps, not enabled options.

Revision 0.6 adds `--disk-radial-layout core_wall` and
`tools/run_annular_disk_normal_core.sh`: a 1877-point unconstrained test from
the converged R20/T0.30/Fs1=5.4 radial normal-core cylinder, resolving both
centre and surface. It retains the default wall-only layout for vortex-free runs.

Revision0.7 adds `a_mermin_ho`, the schematic spin-texture `a_panam`, and
the author's original mixed A/polar seed `a_planar`. Use
`tools/run_annular_A_texture.sh NAME`. Texture diagnostics distinguish orbital
chirality from a unit l-vector and report principal spin direction/rank.
An orbital Pan-Am seed and dipole/Zeeman/rotation energetics remain future work.

LaTeX revision 0.2 (7 October 2026) adds Appendix A, **From COMMON blocks to
explicit simulation objects**, for readers coming from Fortran 77. It traces
the actual `3DFS_MPICodes/qcv.dat` COMMON declarations through
`new_src/global_dec.f90` to modern derived types, explaining allocation,
component access with `%`, precision, argument intent, type-bound procedures,
indexing, copying and MPI ownership, with a small worked example.

Code-oriented overview, 6 October 2026. Start here to navigate the package;
use the linked specialist notes for equations and detailed options. Historical
benchmark notes describe their own implementation date, not necessarily today's
capabilities. This guide does not change any physics or launch a calculation.

## 1. What we have built

The present solver is a spinful quasiclassical mean-field solver for 3He, with
fields depending on x and y and translation invariance along z. Momentum still
has three components: spatially 2D does **not** mean a 2D Fermi surface.
The modern solver is Fortran; Python prepares, archives, diagnoses and plots
runs. `new_src` is the separate radial reference, not the engine underneath
every modern run.

At an independent point the state contains nine complex values A(spin,orbital)
and three real current-related Fermi-liquid fields: 21 real mixing entries.
The propagator retains spin structure. General superconducting devices remain
the goal; the implemented self-consistency law here is the 3He model, not yet
a selectable library of arbitrary materials and interfaces.

| Calculation | Current route | What is constrained? |
|---|---|---|
| Free vortex, general 2D | Axisymmetric benchmark in `full_2d`, or double-core driver | No cylindrical field symmetry |
| Free vortex, radial | Axisymmetric benchmark in `radial_symmetry` | Independent +x ray; angular reconstruction imposed |
| Specular cylinder | Modern cylinder driver | Radial or experimental unrestricted uniform-grid disk; reflected trajectories |
| Legacy comparison | `new_src` | Separate radial solver and conventions |
| Manufactured demonstration | 2D demo | Prescribed fields, no self-consistency |

The historical name `benchmark_axisymmetric_core_2d` is misleading: in
`full_2d` it evolves an unconstrained planar field, even when initialized from
an axisymmetric reference. Conversely, plotting a radial solution over a plane
does not make its calculation unconstrained 2D.

## 2. The algorithm at a glance

```text
User: example namelist + explicit command-line overrides
                         |
Python runner: build -> unique run directory -> archive input/executable
                         |
Fortran driver: physics + quadratures + grid + initial/restarted fields
                         |
                  current state X_k
                         |
             refresh dependent exterior fields (free vortex)
                         |
     distribute independent target points over MPI ranks
                         |
     for each target: momentum directions x energy poles
          |
          +-- construct straight or reflected trajectory
          +-- sample A and mean fields along the trajectory
          +-- Riccati propagation from both ends
          +-- reconstruct spin-Nambu propagator
          +-- accumulate gap and Fermi-liquid self-consistency
                         |
              collect mapped state F(X_k)
                         |
       residual R_k = F(X_k)-X_k -> diagnostics/checkpoints
                         |
       converged? ---- yes ---> final output and plots
          |
          no -> chosen accelerator -> next state -> constraints -> repeat
```

This is a functional map, not a claim that every driver has identical control
flow. For example, zero-mode projection and external probes belong to specific
free-vortex workflows; the cylinder has its own short iteration driver.

MPI currently assigns whole target-point calculations to ranks. A rank loops
over that point's directions and poles locally. Fields are replicated; this
is not spatial domain decomposition with nearest-neighbour ghost exchange.
The accelerator update is coordinated and the resulting state communicated.
This avoids splitting one trajectory across ranks, but memory/communication
and uneven point cost can eventually limit scaling.

## 3. Where to find each layer

Paths below are relative to the repository root.

| Layer / question | Start reading here |
|---|---|
| Run preparation and archive | `tools/run_axisymmetric_core_benchmark.py` |
| Free-vortex solve orchestration | `app/benchmark_axisymmetric_core_2d.f90` |
| Double-core orchestration, probes and specialized diagnostics | `app/benchmark_legacy_double_core_map_2d.f90`, `tools/run_double_core_from_scratch.py` |
| Cylinder solve and launcher | `app/benchmark_radial_cylinder_2d.f90`, `tools/run_modern_cylinder.py` |
| Mesh coordinates and indexing | `src/cartesian_mesh_2d.f90` |
| Field storage and radial reconstruction | `src/spinful_state_2d.f90` |
| Cartesian interpolation | `src/bilinear_field_sampler_2d.f90` (name retained as options evolved) |
| Basis and harmonic conventions | `src/order_parameter_basis.f90` |
| Historical and localized core seeds | `src/historical_core_seed_2d.f90` |
| MPI / serial map | `src/he3_mpi_field_map_2d.f90`, `src/he3_serial_point_map.f90` |
| Straight / reflected rays | `src/straight_trajectory_2d.f90`, `src/specular_cylinder_2d.f90` |
| Free-vortex exterior | `src/free_vortex_asymptotic_2d.f90` |
| Riccati transport | `src/legacy_riccati_trajectory_2d.f90`, `src/riccati_segment.f90` |
| Propagator and spin algebra | `src/quasiclassical_propagator.f90`, `src/spin_matrix_2x2.f90` |
| Physical self-consistency integrand | `src/he3_self_consistency_integrand.f90` |
| Quadratures, poles and bulk gap | `src/he3_quadrature.f90`, `src/ozaki_generator.f90`, `src/he3_bulk_gap.f90` |
| Packing, residuals and updates | `src/new_src_iteration_layout_2d.f90`, `src/he3_nonlinear_iteration_2d.f90` |
| Accelerators | `src/legacy_anderson_mixing.f90`, `src/barzilai_borwein_mixing.f90`, `src/polyak_mixing.f90` |
| Regression coverage | `tests/`, registrations in `CMakeLists.txt` |

`examples/` contains inputs; `work/` contains disposable builds and working
products; `runs/` contains individually dated experiments; `benchmarks/` holds
reference cases. `2D_benchmarks/` contains separately maintained legacy data
and may be absent on another checkout. Do not assume it was uploaded to GitHub.

## 4. Grid, trajectories and the exterior are different objects

The full-2D field grid is a tensor-product Cartesian grid. Uniform, piecewise
multiscale and smoothly stretched choices exist in the appropriate drivers.
The smooth option grades node spacing through a sinh coordinate mapping;
it is not an unstructured finite-element mesh or a set of elliptical rings.
Active masks and asymptotic fit surfaces need not coincide with the outer box.

An **active point** receives a costly self-consistency map and contributes to
the iteration vector. A **dependent halo point** supplies interpolated fields
but may be filled from current-state asymptotics rather than independently
relaxed. External probes test this extrapolation; their role and whether they
are iterated depend on the selected driver/input. They must not be silently
counted as ordinary interior convergence points.

For a free vortex, the phase-wound bulk matrix A0 is the asymptotic reference.
Current-state fitting allows angle-dependent coefficients, including the
rotation-like order-parameter tail and inverse-radius mean-field terms.
Reference-based endpoint options also exist for controlled comparisons: they
are not equivalent to `state_asymptotic`. Read the input before interpreting a
run as an independent relaxation. See [evolving asymptotics](EVOLVING_ASYMPTOTIC_RESTART.md).

Trajectory sampling has its own step size, independent of field-node spacing.
Reducing the number of field points does not shorten a trajectory. In radial
mode, trajectory fields are interpolated directly from the independent ray;
they are not recovered by interpolation of the reconstructed plotting grid.

For a cylinder the wall is physical. Rays reflect specularly, stay inside it,
and use no free-vortex exterior fit. Reflected half-path length is propagation
distance, not cylinder radius. The cylinder supports uniform collision sampling
and collision-aligned sampling; finite-step results need not coincide.

## 5. Input: which knob changes what?

| Quantity | Meaning and main controls |
|---|---|
| Temperature | T/Tc; free-vortex startup can generate poles and bulk gap |
| `fs1` | Tabulated F1s, internally converted to F1s/(1+F1s/3) |
| `half_width` | Physical Cartesian extent, not resolution |
| `number_of_cells` | Cells across an axis; nodes = cells + 1 |
| `active_radius` | Extent of independently mapped region/ray |
| `mesh_kind`, `smooth_stretch` | Node placement; same count can give different core spacing |
| `trajectory_interpolation_order` | Field sampling order, not Riccati integration order |
| `trajectory_maximum_step` | Upper trajectory interval size |
| Asymptotic fit/matching/outer radii | Where exterior coefficients are extracted and rays end |
| `initialization_mode` | Reference, historical seed, localized harmonic or restart workflow |
| `iteration_method` | `anderson`, `bb`, or `polyak` in the free-vortex benchmark |
| `maximum_iterations`, `convergence_tolerance` | Work limit and stopping threshold |

In radial symmetry, an even N-cell grid with active radius equal to half-width
has N/2+1 independent nonnegative-ray nodes: 118 cells gives 60 points, not
118 radial points. The cylinder instead takes `radial_points` explicitly.

The cylinder runner is deliberately narrower: temperature is currently fixed
at .30 with its archived Ozaki table; `--fs1`, radius and discretization are
exposed. Its accelerator is Anderson, p_max=3, history length 10. Do not pass
free-vortex flags to it expecting them to work. The saved namelist is the
authoritative record of the run that was actually launched.

## 6. A short practical tour

The experimental unrestricted cylinder branch is now available through
`--spatial-mode full_2d` in the cylinder runner. See the
[first disk benchmark](UNCONSTRAINED_DISK_TEST.md) for the interior-only sampler,
limitations, radial-reference launch and independent 2D outputs.

Run commands from `/Users/mikael/Documents/Codex/3he-vortex-modernization`
(or the equivalent Linux checkout). Required tools are a Fortran compiler,
MPI, CMake, LAPACK, and Python with NumPy/Matplotlib. Avoid overlapping laptop
MPI runs. `--help` is safe and does not launch a solve.

### A. Explore the plotting pipeline without a physical solve

```sh
python3 tools/run_2d_demo.py --case-name orientation-demo --formats png
```

Inspect its input and field plots. These are manufactured fields: a test of
storage/interpolation/plotting, not an equilibrium vortex.

### B. Follow three radial updates through the real solver

```sh
python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/radial_bb_localized_0plus.nml \
  --case-name orientation-radial \
  --ranks 10 --max-iterations 3 --checkpoint-interval 1
```

This uses the existing compact 0plus scratch seed at T=.30, Fs1=5.4, 60 radial
points and BB. It is a short learning run, not a convergence test. The reference
file paths in that input support comparisons; they do not initialize this seed.
Use `--iteration-method anderson --anderson-pmax 3` for a separate AA experiment.
Keep all other settings fixed when comparing accelerators.

Inspect the run directory printed by the runner. To plot during a later run:

```sh
python3 tools/plot_running_axisymmetric.py --help
python3 tools/plot_running_axisymmetric.py --latest orientation-radial
```

Read after a completed checkpoint; retry if a file is being rewritten.
For an actual restart, use `--restart PATH/TO/final_fields_2d.dat` with a
matching input and an explicit iteration limit. A field restart is not proof
that accelerator history is preserved. Preserve parent directories.

### C. Understand the full-2D branch before launching it

Compare `examples/2d_double_core_Fs1_zero_smooth.nml` with the radial input.
Read [smooth-grid test](SMOOTH_GRID_TEST.md), then inspect:

```sh
python3 tools/run_double_core_from_scratch.py --help
```

Its `--iterations` differs from the axisymmetric runner's `--max-iterations`.
Do not simply change radial mode to full_2d at the same 118-cell resolution:
that changes the independent point count dramatically. Start from a known
small 2D example and use one map or very few updates first.

### D. The newly prepared confined-vortex experiment

```sh
caffeinate -i bash tools/run_cylinder_R20_a_core.sh
```

On Linux omit `caffeinate`. This is the longer user-run experiment, not a small
tutorial: R=20 xi0, T=.30, Fs1=5.4, 201 nodes, unit winding, A-phase-core seed
in B phase, 10 ranks and up to 300 AA updates. See the
[cylinder guide](RUNNING_RADIAL_SPECULAR_CYLINDER.md) for live plots and restart.
Initialization `-2` means bulk chiral A without a vortex; it is **not** the
A-phase-core vortex, which uses `1` in this cylinder driver.

## 7. Read a run like an experiment

1. Check saved input, initialization/restart source and executable provenance.
2. Confirm independent point count, temperature, Fs1 and quadrature.
3. Read terminal status: an iteration limit or single map is not convergence.
4. Inspect RMS and maximum fixed-point residuals, not only the applied step.
5. Locate the remaining residual: core motion, fit surfaces and mesh transitions
   suggest different next tests. A peak hopping between equivalent cores can
   be a tie-breaking effect, not motion.
6. Compare Cartesian and harmonic profiles, including imaginary components.
7. Refine grid, trajectory step, quadrature and outer treatment separately.

For packed independent state X, R=F(X)-X. Absolute RMS is sqrt(sum(R_i^2)/N),
maximum is max(abs(R_i)), and relative L2 is norm(R)/norm(X), subject to the
driver's zero-denominator safeguards. These mix gap and mean-field entries
without physical rescaling. Area-weighted and region-specific diagnostics
are separately labelled. Legacy printed scaled maxima are not directly the
modern absolute maximum. Accelerator update size need not equal residual size.

Propagator normalization checks the algebraic quasiclassical identity; values
near machine precision are useful but do not bound spatial, angular, energy
or endpoint errors. Likewise, `benchmark passed` only certifies the configured
checks; inspect convergence status and reference discrepancy independently.

The three stored current-related fields are not yet automatically calibrated
physical mass-current density. A gap-amplitude/pair-density plot is not a
computed superfluid-density tensor. Keep these distinctions in publications.

Cylinder files use `cylinder_*` names; free-vortex workflows generally use
`metrics.txt`, `iteration_history.dat`, `final_fields_2d.dat`, checkpoint files
and residual snapshots. The cylinder evaluates its final saved state; for
other drivers, check map/update indexing before attributing a residual to a
post-update checkpoint. A separate zero-update map is the unambiguous check.

## 8. Reading order and development boundaries

For hands-on use: this guide -> [radial mode](RADIAL_SYMMETRY.md) ->
[live plotting](PLOTTING_RUNNING_BENCHMARK.md) ->
[cylinder mode](RUNNING_RADIAL_SPECULAR_CYLINDER.md).
For theory: [implementation note sources](technical_note/README.md),
[mathematical conventions](architecture/MATHEMATICAL_CONVENTIONS.md), then
the actual map and integrand modules above. Older PDFs can lag new features.
For accelerator work: [BB](BB_ITERATION.md), [Polyak](POLYAK_ITERATION.md),
[AA cycling](CYCLED_ITERATION.md), and signed diagnostics in the radial guide.

Implemented does not mean validated for every parameter regime. We have
radial/full-2D cross-checks, converged vortex experiments and cylinder regression
checks. We do not yet have general-geometry 2D container-wall sampling, leads and
spin-active interfaces, a DG backend, GPU execution, or a learned accelerator.
Automatic mesh adaptation is not a universal solver feature: mesh utilities
and staged domain studies should not be confused with fully automatic error
control during every run. The GPU/HPC roadmap is in
[project plan](PROJECT_PLAN.md) and [scaling guide](hpc/SCALING.md).

When adding physics, keep interpolation, propagator algebra, self-consistency,
iteration and visualization separately testable. Begin with a small deterministic
test, compare one-rank and multi-rank results, then perform a controlled physics
benchmark. Record changes in `MODERNIZATION_LOG.md`; do not replace reference
data with a new solver's output merely because a residual is small.
