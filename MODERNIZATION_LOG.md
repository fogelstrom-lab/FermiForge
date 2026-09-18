# FermiForge Modernization Log

This file records analysis, decisions, changes, verification results, and open
questions for the modernization of the legacy Fortran/MPI vortex solver.

The working project name is **FermiForge**: a scalable quasiclassical
simulation framework for spinful superfluids and superconductors.

## Working rules

- Legacy file contents are immutable reference material. Organizational moves
  inside incoming/legacy-f77/ are allowed when recorded here.
- Modernized source, tests, build files, and documentation will be created
  outside incoming/.
- Physics-changing work will be separated from behavior-preserving
  modernization.
- Every accepted numerical change will be checked against a named legacy
  baseline.
- Unresolved conventions and suspected defects will be logged rather than
  silently corrected.

## Status

Current phase: modern radial-mesh foundation and legacy reproducibility audit.

No legacy source files have been modified.

Authoritative source baseline: the original `3DFS_MPICodes` fixed-form code,
preserved under `incoming/legacy-f77/`. The matching copy under
`new_src/3DFS_MPICodes/` confirms its provenance but is not a separate code
line to merge.

## 2026-09-14 - Initial intake

### Material received

- 74 files, approximately 848 KiB.
- Approximately 6,364 lines across the fixed-form Fortran sources, include
  files, Makefile, and primary input file.
- Main program: incoming/legacy-f77/qcv.f.
- Current Makefile build set: qcv.f, iter.f, gap.f, riccati.f, linpacks.f,
  trajectories.f, init_calc.f, mpicalls.f, and interpol.f.
- Additional source variants or utilities not in the current link:
  greenR.f, intord_v.f, intpol_new.f, iter_org.f, poretraj.f, and rescale.f.
- Shared state is held primarily in qcv.dat, iter.dat, and gap.dat through
  COMMON blocks.
- MPI uses the legacy mpif.h interface. Work is distributed over radial grid
  indices as iy = myid, nx, nproc, followed by reductions of nine
  order-parameter components and three vector components.
- The main driver uses five pre-iterations followed by runciterat calls with
  a sequence of mixing/subspace parameters.
- The present input selects t=0.30, unit vorticity, cylindrical mode,
  F_1^s=5.40 before the source's amplitude conversion, radius 20, 48
  azimuthal trajectories, tolerance 2e-6, an A-core initial state, and
  40 maximum iterations.
- Ozaki frequency tables and Gaussian quadrature tables are present.
- Reference output families are present for combinations labelled cyl, vort,
  a, and n. The current unsuffixed files contain 50 radial rows; several
  labelled baselines contain 100 rows, implying a different compile-time nx.

### Initial architectural observations

1. Cylindrical symmetry is embedded in array dimensions, interpolation,
   trajectory construction, input selection, and output layout; it is not
   confined to one routine.
2. The physical state contains a complex 3-by-3 order parameter represented
   as nine separate radial arrays, plus three complex vector arrays.
3. Fixed compile-time dimensions (nx, mx, frequency and quadrature limits)
   are propagated through include files and COMMON blocks.
4. -fdefault-real-8 determines the effective precision. A modern port must
   replace this compiler-wide behavior with explicit kinds before changing
   interfaces.
5. -fallow-argument-mismatch indicates legacy implicit interfaces or actual
   type/rank mismatches. These must be catalogued before explicit interfaces.
6. The bundled linpacks.f should initially be retained for reproducibility,
   then compared with a modern LAPACK replacement in a separate change.
7. Specular scattering appears in the trajectory layer. The exact
   reflection/interpolation path still needs a full trace.
8. The incoming directory includes local MPI compatibility headers. These
   should not be carried into the modern source tree.

### Reproducibility status

- The Makefile and primary input are present.
- GNU Fortran 15.2.0 and Open MPI 5.0.9 are installed under
  /opt/homebrew/bin. They are not in the task's default PATH, so build and run
  commands explicitly prepend that directory.
- Incoming reference outputs can still define file-level and field-level
  regression tests.
- The code expects ozaki.dat and gauss11.dat in the run directory and reads
  the main parameters from standard input.
- ozaki.dat is byte-for-byte identical to ozaki_T=0.3.dat and its header
  records T=0.30, so it is consistent with the supplied qcv.inp.
- The unsuffixed output files have exactly nx+1 = 50 rows and span radius
  0 to 20 in steps of about 0.408, consistent with the current compile-time
  grid and input radius.
- fort.30 contains an iteration diagnostic headed by average and maximum
  errors plus a field snapshot; fort.99 is empty. The surviving diagnostic
  begins at iteration 1 and is not by itself evidence of final convergence.

### Immediate next audit

- Build a complete routine call graph and COMMON-block ownership map.
- Trace the order-parameter phase winding and every cylindrical assumption.
- Trace specular reflection and interpolation of incoming/outgoing paths.
- Identify compile-time and runtime dimensions.
- Characterize each reference-output family and select a canonical baseline.
- Run static compiler checks once a Fortran/MPI toolchain is available.
- Propose module boundaries and a staged migration plan without editing the
  reference tree.

## 2026-09-14 - Toolchain verification and historical archive

### Toolchain

- Verified /opt/homebrew/bin/gfortran: GNU Fortran 15.2.0.
- Verified /opt/homebrew/bin/mpifort and mpirun: Open MPI 5.0.9.
- No broader filesystem access or Full Access mode is required. Commands can
  use the installed toolchain by setting PATH for the individual build/run.

### Baseline build

- Copied the uploaded tree to work/legacy-baseline-2026-09-14.
- Built the unchanged current Makefile target successfully.
- All nine objects compiled and the qcv executable linked.
- The compiler reported one legacy MPI diagnostic in mpicalls.f: MPI_BCAST is
  called with both LOGICAL and COMPLEX actual arguments through the implicit
  mpif.h interface. This is tolerated by -fallow-argument-mismatch and is a
  priority for the later explicit-interface conversion.
- No legacy source was edited to obtain the successful build.

### Historical-file archive

The following files were absent from the current Makefile and duplicate or
supersede routines in the active build. They were moved without content
changes to incoming/legacy-f77/archive/historical-source/:

- greenR.f
- intord_v.f
- intpol_new.f
- iter_org.f
- poretraj.f
- rescale.f

Unused alternate quadrature tables were moved to
incoming/legacy-f77/archive/historical-data/:

- Gaussian.dat
- gauss6.dat
- gauss10.dat

All supplied test outputs, current inputs, Ozaki tables, MPI compatibility
headers, and active build sources remain in their original locations.

## 2026-09-14 - GPU and multiscale target added

The modernization target now explicitly includes GPU portability and a
multiscale, self-consistent treatment of the full 3He order parameter in
general containers. The design direction is recorded in
docs/architecture/GPU_AND_MULTISCALE_DESIGN.md.

### Accelerator observations and decisions

- The legacy loop exposes independent work over evaluation point, momentum
  direction, and Ozaki energy. A Riccati integration remains sequential along
  its own ray, but many such rays can run concurrently.
- With the supplied input, that outer phase space contains 50 x 48 x 11 x 8 =
  211,200 independent `makeprops` calls per nonlinear iteration.
- The first modern implementation will remain a CPU reference, but its field
  layout and kernel interfaces will be GPU-compatible from the start.
- Standard modern Fortran will contain the physics. OpenMP target offload is
  the provisional portable accelerator path, with OpenACC or CUDA Fortran
  retained as optional measured backends rather than embedded requirements.
- MPI remains the inter-node layer. Initial GPU work will be batched within
  each MPI rank.
- Geometry and interpolation metadata will be separated from propagation and
  cached across nonlinear iterations where memory permits.

### Grid and self-consistency direction

- The leading mesh candidate is block-structured Cartesian AMR with an
  embedded-boundary or signed-distance geometry description.
- Fine blocks will resolve coherence-length core and surface structure; coarse
  blocks will cover the long-wavelength spin-orbit texture.
- The full complex 3-by-3 order parameter will initially remain an unconstrained
  field on the hierarchy. A phase/rotation/amplitude decomposition may later
  serve as a coarse correction or preconditioner, but will not be imposed near
  defects.
- The spatial field mesh will be separated from adaptive integration points
  along a trajectory. This avoids forcing every long ray to use the finest
  spatial step.
- Refinement will be driven by distance to defects and boundaries, field
  variation, texture variation, and the self-consistency residual, with
  resolution selected through convergence tests.
- Specular reflection will be expressed through a general geometry query and a
  separate boundary-scattering law, allowing later diffuse or mixed models.

### New open information

- Accelerator vendors, compiler toolchains, and per-GPU memory on the intended
  production systems.
- Whether the first symmetry-relaxed milestone should be a general 2D cross
  section or a fully 3D container.
- Which spin-orbit, dipole, magnetic, flow, and surface contributions must be
  active in the first self-consistent texture milestone.

## 2026-09-14 - Portability and first general geometry clarified

- The required implementation must not depend on a particular GPU family.
  Standard OpenMP target offload is the baseline accelerator model, paired with
  a CPU implementation of the same kernels.
- GPU-vendor APIs, types, and libraries will not enter physics-module
  interfaces. Vendor-specific build flags and optional tuned kernels may be
  added behind common interfaces and must have a portable reference path.
- The first symmetry-relaxed geometry will be a two-dimensional cross section
  with translational symmetry along the cylinder axis.
- Real-space fields will depend on `(x,y)`, but the full three-dimensional
  Fermi-surface momentum quadrature, including `p_z`, will be retained.
- Specular reflection from the extruded side wall changes the in-plane momentum
  and conserves `p_z`.
- The initial adaptive mesh will consequently be a block-structured quadtree
  with an embedded boundary. Interfaces will remain suitable for a later 3D
  backend.

## Open questions

1. Which labelled output family is the authoritative physics baseline?
2. Which compiler and MPI implementation most recently produced it?
3. What launch command and process count were used?
4. Are unsuffixed op_xyz, op_harm, curr, and xaxis completed, converged
   outputs for the supplied qcv.inp?
5. Is ozaki.dat intentionally the table for T=0.30, or should it be copied
   from ozaki_T=0.3.dat before a run?
6. Are poretraj.f and the alternate interpolation files part of intended
   work toward general geometry, or historical experiments?
7. Are there notes defining the units and meanings of vort, subvort, aa0,
   vx:vy:vz, and the a/n/d-core initial states?
8. Which GPU systems and Fortran compilers are available for production runs?
9. Which noncircular two-dimensional cross section should be the first general
   geometry test after the circular-cylinder regression?

## 2026-09-14 - Modern radial-mesh implementation initiated

The first modern source and test tree has been created alongside the immutable
legacy reference. No file below incoming/legacy-f77 was edited.

### Implemented

- Added a CMake Fortran 2018 build.
- Added an explicit double-precision kind module.
- Added a dynamically sized radial mesh with zero-based point indices to ease
  comparison with the legacy `0:nx` arrays.
- Added exact construction of the legacy 49-cell/50-point uniform grid.
- Added bisection-based refinement of independent core and surface regions.
- Added a two-to-one adjacent-cell balance rule, mirroring the intended later
  block-AMR constraint.
- Added binary cell lookup and three-point Lagrange interpolation for complex
  fields on a nonuniform radial grid.
- Added tests for legacy uniform coordinates, adaptive-mesh validity and
  balance, reduced point count relative to a globally fine grid, endpoint
  preservation, and exact interpolation of a complex quadratic field.

### Verification

- Configured and built with GNU Fortran 15.2.0.
- The normal debug build passes the `radial_mesh` test.
- A second build using `-Wall -Wextra -Wconversion-extra -fcheck=all
  -fbacktrace` is warning-clean and passes the same test.
- The uploaded legacy source and reference outputs remain unchanged.

### Low-temperature A-core benchmark definition

Added benchmarks/a_core_low_temperature_large/README.md. The provisional
starting family retains cylindrical symmetry and uses:

- A-core initialization;
- `T/Tc = 0.10`, for which an uploaded Ozaki table is available;
- primary radius `R = 80` in legacy length units;
- isolation comparisons at `R = 60` and `R = 100`;
- fine spacing target `0.125` in initially refined core and surface widths of
  `8.0`, with a coarse bulk limit of `4.0`.

These are deliberately labelled starting values. Core-surface overlap will be
declared negligible only after radius independence is demonstrated for both
localized profiles. The legacy length convention and the historically
problematic low temperature remain to be confirmed.

### Next implementation step

- Route uniform-grid interpolation through the modern lookup/interpolation
  interface and compare it pointwise with the legacy `intpol` implementation.
- Only then enable nonuniform nodes in the physics calculation and decouple the
  trajectory integration step from radial cell spacing.

### Dynamic state storage completed

- Added dynamically sized storage for all nine complex order-parameter
  components and the three complex exchange-field components. The spatial
  point is the contiguous array dimension for later batched CPU/GPU kernels.
- Added pack/unpack routines that exactly reproduce the component ordering in
  legacy `makearray` and `makeself`.
- Preserved and documented the legacy exchange-field convention: only its real
  part is packed, it is multiplied by ten in the iteration vector, and it is
  reconstructed as real after division by ten. This convention has not been
  silently reinterpreted.
- Added a round-trip test with and without the exchange field, including direct
  checks of the first `dxx` and `vx` positions in the packed vector.
- The normal and strict builds now contain two tests, `radial_mesh` and
  `radial_state`; both pass.

## 2026-09-14 - Reproducible MPI runner and plotting added

Added `tools/run_and_plot.py` for the current cylindrical reference solver.

### Runner behavior

- Builds a fresh copy of the active legacy source in
  `work/legacy-runner-build`; it does not add objects or executables to the
  immutable incoming tree.
- Uses ten MPI processes by default and accepts `--ranks` for other counts.
- Uses `slot:OVERSUBSCRIBE` placement with binding disabled because the desktop
  Open MPI environment advertises fewer scheduler slots than the requested
  local process count.
- Selects the temperature-labelled Ozaki table by reading and matching its
  header, then copies it into the run as `ozaki.dat`.
- Creates a unique timestamped directory below `runs/` and archives the input,
  executable, selected tables, stdout, stderr, elapsed time, command, host,
  and parsed input parameters.
- Explicitly records the legacy 49-cell radial spacing and warns that these
  runs do not yet exercise the adaptive mesh.
- Preserves MPI wrapper symlinks and supplies the Homebrew compiler directory
  in the build environment; resolving `mpifort` to its generic wrapper target
  or omitting `gfortran` from `PATH` prevents Open MPI from linking correctly.
- Interrupted runs created by the current script are marked `interrupted` in
  their metadata with return status 130.

### Plotting

- Plot-only mode accepts either an archived run or the uploaded reference-data
  directory.
- The Apple system Python lacks NumPy and Matplotlib, so the initial plotting
  implementation used the installed gnuplot executable. This was later
  superseded by the Homebrew Python/Matplotlib path documented below.
- Produces 3-by-3 Cartesian and harmonic order-parameter component figures,
  amplitude/exchange-field plots, and a logarithmic convergence plot.
- Parses both the five initial `it nr` records and the subsequent accelerated
  nonlinear-iteration records from captured stdout.
- Writes a normalized convergence table and JSON figure manifest beside the
  PNG files; the optional gnuplot backend also preserves its command files.

### End-to-end verification

- Plotting the uploaded reference output produced four PNG files successfully.
- A ten-process run of the supplied production input was healthy through 15
  recorded nonlinear evaluations before being intentionally interrupted; it
  was not used as a numerical regression result.
- Added `tests/data/qcv_workflow_smoke.inp`, which has deliberately loose and
  small parameters and is only an execution-path test.
- The ten-process workflow smoke run completed successfully in approximately
  2.1 seconds and generated all four plot families. It is not a physics
  benchmark.

## 2026-09-14 - Plotting migrated to Matplotlib

- A completed ten-process run in `runs/20260914-150620-cylindrical` exposed a
  gnuplot startup failure caused by the user's `GNUTERM=qt` environment setting
  and a gnuplot installation without Qt terminal support.
- The physics run itself completed normally in approximately 3,772 seconds.
  Its `op_xyz`, `op_harm`, `curr`, and `xaxis` files each contain all 50 radial
  points; the plotting error did not damage the numerical output.
- The runner now uses NumPy and Matplotlib with the noninteractive Agg backend
  by default. The available Homebrew Python 3.14 installation provides NumPy
  2.4.3 and Matplotlib 3.10.8.
- Added a project-local Matplotlib cache location so plotting does not depend
  on write access to the user's home-directory cache.
- Gnuplot is retained only behind `--plot-backend gnuplot`; that fallback now
  removes `GNUTERM` from its environment before startup.
- Replotted the completed run successfully with Matplotlib and visually checked
  all four figures.
- The convergence history contains 205 evaluations. The final reported average
  and maximum errors are approximately `9.18e-4` and `1.93e-3`, respectively,
  well above the effective requested tolerance printed by the solver
  (`5.6e-7`). The calculation completed by exhausting its scheduled iteration
  stages, not by satisfying the convergence criterion.

## 2026-09-14 - Review of the earlier `new_src` modernization

The newly supplied `new_src` directory was reviewed without modifying its
contents. It contains two different bodies of material:

- `new_src/3DFS_MPICodes/` is another legacy snapshot. Its active fixed-form
  source, Makefile, tables, MPI headers, and labelled reference-output families
  are byte-for-byte identical to the corresponding files already held under
  `incoming/legacy-f77/`. Its unsuffixed outputs and input differ and may be
  useful as an additional historical run.
- The top-level free-form `.f90` files are the substantive modernization
  attempt, dated in the main-program comments from October 2024 through January
  2025. Alternate and backup versions of several iteration and interpolation
  routines are present alongside the nominal build set.

### Useful design work to retain

- The radial field nodes use the analytic mapping
  `r(i) = Rx*tan(i*atan(gridM/Rx)/nx)`, while the trajectory propagation keeps
  a separate uniform step `dx = 2*gridM/mx`. Decoupling the field mesh from the
  Riccati propagation step is the correct design direction for adaptive work.
- Quadratic interpolation was generalized to nonuniform radial nodes and is
  performed in the symmetry-adapted harmonic representation. The origin parity
  and phase-winding transformations contain valuable physics knowledge that
  should be ported with focused tests.
- The top-level calculation is split into modules and uses array syntax and
  allocatable storage inside the nonlinear solvers. This is substantially more
  readable than the fixed-form program.
- Simple relaxation, Broyden, Barzilai--Borwein, and Anderson variants were
  separated conceptually, and average/global and largest-component convergence
  measures were recorded independently.
- Packing by component makes each radial component contiguous. That layout is
  potentially useful for vector and accelerator kernels, although it is not the
  same ordering as the legacy iteration vector.

### Build and correctness audit

- The supplied Makefile builds and links serially with GNU Fortran 15.2.0 and
  Open MPI 5.0.9. This success depends on implicit typing.
- Recompilation with implicit declarations disabled stops on undeclared names
  in `init_calc.f90`, `mpicalls.f90`, `iter.f90`, and `qcv.f90`. One important
  example is the declaration `NcCutOff` versus use of the legacy spelling
  `NcCutof`; these are different Fortran identifiers.
- `input_bcast` contains undeclared scalar buffers `inmax` and `pmax` passed to
  broadcasts with a count of ten. With implicit typing these calls compile but
  can access memory beyond the scalar. The representation flag `irep` is owned
  by a different module and is not actually broadcast through that module
  variable. Vorticity is broadcast twice.
- Several array broadcasts use literal counts smaller than the declared arrays.
  These happen not to affect the presently used low indices, but are unsafe as
  general interfaces. Modern broadcasts must derive counts from `size()` or
  transmit a configuration object using the type-safe `mpi_f08` interface.
- A parallel `make -j10` fails because the Makefile does not encode Fortran
  module dependencies. Serial build order only works because `SRC` happens to
  be listed in dependency order.
- The build enables OpenMP and links Apple's Metal framework, but there are no
  OpenMP, OpenMP-target, OpenACC, or other accelerator directives in the source.
  Linking Metal does not make the Fortran kernels GPU-capable.
- Precision still depends on `-fdefault-real-8`; primary dimensions remain
  compile-time constants; and most physical state remains mutable module-global
  storage. `use mpi` is an improvement over bundled `mpif.h`, but is not the
  type-safe `mpi_f08` target.
- The standalone `scsolver.f90` and `broyden.f90` files are not part of the
  Makefile and do not compile as Fortran 2018 without repairs. They should be
  treated as algorithm references, not production dependencies.

### Numerical and physical differences from the legacy cylinder solver

- The supplied mapping has 100 points over radius 70. Its spacing grows from
  approximately 0.144 at the core to 6.55 at the outer edge. It resolves the
  core preferentially but strongly under-resolves a surface layer; it is a
  static graded mesh rather than residual-driven adaptivity.
- The top-level source always calls the free-vortex interpolation path even
  though it still reads `icyl`. The cylindrical trajectory and specular
  reflection routine are absent from the active free-form implementation. The
  included input also selects `icyl=0`. Consequently its supplied outputs are
  not a cylinder or surface-state benchmark.
- The Riccati integration subdivision parameter was changed from `ndx=4` in the
  legacy source to `ndx=1`. This is a physics-numerics change and requires a
  step-convergence study before acceptance.
- The free-vortex far-field matching point is `xgrid(nx-5)`. On a nonuniform or
  adaptive mesh, selecting the fifth node from the end changes its physical
  location and is not a mesh-independent boundary rule.
- The interpolation routine diagnoses, but does not safely handle, a lookup at
  or beyond the final node before forming an `iz+1` index. Its current caller
  generally avoids the case through the far-field cutoff; the interpolation
  API itself is not safe at the boundary.
- The new packing scheme always includes 21 real component fields, omits the
  legacy factor-of-ten weighting of the exchange field, and changes the packed
  ordering from point-major to component-major. These may be sensible solver
  choices but are not behavior-preserving transformations.
- The dense Broyden implementation stores an `N` by `N` Jacobian. It is viable
  for the current 2,100-variable radial problem but cannot scale to a general
  two-dimensional adaptive mesh. Limited-memory Anderson mixing is the more
  promising algorithmic basis for the symmetry-relaxed problem.
- The included convergence histories are useful experimental evidence but use
  convergence definitions different from the legacy iterator and are not yet
  canonical regression data. For the recorded tolerance `2e-6`, the histories
  end after 20 simple, 83 Broyden, 250 Barzilai--Borwein, and 51 Anderson
  iterations; the Barzilai--Borwein history did not meet its recorded largest-
  component criterion.

### Integration decision

The top-level `new_src` implementation should be preserved as a scientific
design reference, not adopted wholesale as the new production base. The
recommended migration sequence is:

1. Add its tangent mapping as a named, tested static mesh option alongside the
   existing uniform and core-plus-surface-refined meshes.
2. Port the harmonic interpolation and origin transformations onto the generic
   modern mesh lookup, with pointwise uniform-grid comparison against the
   legacy routine and explicit endpoint tests.
3. Port the Riccati propagator behind a small kernel interface while initially
   retaining `ndx=4`; assess `ndx=1` only as a separate convergence experiment.
4. Restore the cylindrical/specular trajectory path from the legacy source
   before running the large-radius A-core surface-separation benchmark.
5. Introduce a common nonlinear-solver interface with simple relaxation and
   limited-memory Anderson first. Keep dense Broyden only as a small radial
   diagnostic if it proves useful.
6. Add type-safe `mpi_f08` distribution after the serial physics kernels and
   mesh-independent interpolation pass regression tests. Accelerator directives
   should then target batches of independent trajectories, not the globally
   coupled iteration driver.

## 2026-09-14 - Original Anderson-like accelerator ported

The fixed-form `citerat` algorithm in the original `3DFS_MPICodes/iter.f` was
confirmed as the authoritative nonlinear accelerator. The free-form Anderson
routine in the later `new_src` experiment is not being used as the baseline.

### Why the later experimental routine can be unstable

Comparison with the original exposed several numerical changes:

- Original `citerat` starts accelerated stages with back-coupling `p=0.01` and
  adapts it. The later routine performs two full-strength fixed-point updates
  before constructing its accelerated subspace.
- The original drops half the history when the new residual is worse than the
  previous *worst* retained residual, or when the residual-minimizing subspace
  makes insufficient progress. The later routine's test can drop history when
  a new residual is merely worse than the current *best* residual.
- The original mixing-ratio update uses `sqrt(usub/u)`; the later routine uses
  `usub/u` directly.
- The later least-squares construction doubles only diagonal matrix entries
  while doubling the right-hand side, changing the multi-vector minimization
  problem.
- Original `citerat` takes a deterministic stochastic step when its subspace
  system is singular or effectively singular. The later routine prints a
  LAPACK error and continues with the returned coefficient vector.
- The later routine also changes convergence normalization and the relative
  weighting of exchange-field components in the packed vector.

These differences explain why the later routine cannot be treated as a
free-form transcription of the stable original algorithm.

### Modern implementation

Added `src/legacy_anderson_mixing.f90` with:

- explicit `real64`-based kinds and `implicit none`;
- a dynamically allocated `legacy_anderson_t` state object instead of COMMON,
  EQUIVALENCE, and compile-time work arrays;
- an `update(current, mapped, tolerance, next, report)` interface that keeps
  MPI and the quasiclassical map evaluation outside the solver;
- the original residual `G(x)-x`, initial `p=0.01`, residual ordering, history
  rejection, residual-minimizing subspace, adaptive back-coupling, and
  deterministic singular-subspace fallback;
- explicit diagnostics for residual norms, mixing, history size, discarded
  vectors, convergence, and fallback use;
- a small partial-pivoting real solver in place of the bundled complex LINPACK
  storage overlay. This is mathematically equivalent for the real Gram system
  but is not expected to be bitwise identical near roundoff.

The state vector and map result are passed separately. This permits the future
driver to pack a `radial_state_t`, compute the expensive map collectively with
MPI or an accelerator, and perform the small history solve on the controlling
rank without embedding communications in the nonlinear algorithm.

### Regression verification

- Added `tests/test_legacy_anderson_mixing.f90`.
- Compiled the unchanged fixed-form `citerat`, its original history routines,
  and original LINPACK code in a one-rank deterministic harness.
- The first modern relaxation step is exact, the two-point update is checked
  directly against the original equations, and a four-update vector trace is
  compared with values generated by the unchanged fixed-form routine. The
  traces agree within `2e-9`; differences begin at roundoff in the small linear
  solve.
- Tests cover convergence to a known coupled fixed point, explicit reset,
  maximum-history rejection, and the singular-subspace fallback.
- Both normal and strict builds pass all three project tests. The strict build
  uses warnings, conversion diagnostics, interface diagnostics, bounds/runtime
  checking, and backtraces.

### Not yet connected

The new accelerator is not yet driving the legacy vortex solver. Integration
requires a modern self-consistency-map interface around `getnewop`. When that
is added, the original outer schedule must also be preserved initially: five
plain fixed-point pre-iterations followed by up to five fresh `citerat` stages
using the corresponding `pmax` and `inmax` input pairs. Any change to that
schedule will be a separately benchmarked convergence experiment.

## 2026-09-14 - Original accelerator connected to `new_src`

The preceding "not yet connected" status is now superseded for the experimental
`new_src` driver. `new_src/iter_AA.f90` has been replaced by a free-form module
that adapts the existing `getnewop`/`pacvector` map to the tested
`legacy_anderson_t` implementation in `src/legacy_anderson_mixing.f90`.
There is deliberately only one implementation of the algorithm: the adapter
owns the MPI-facing packing and unpacking, while the canonical module owns all
history, rejection, minimisation, adaptive-mixing, and fallback logic.

The adapter preserves the original `citerat` residual definition and reports
the unnormalised maximum-component residual to `qcv`. It initializes the
accelerator with the original first-step mixing `p=0.01`, ten-vector history,
and progress threshold `prog=0.1`. The current experimental driver supplies no
historical `pmax`/`inmax` stage arrays, so this first integration uses a
documented maximum mixing of `1.0` and one continuous solver stage. Restoring
the five original stage pairs from the legacy input is still required before
calling `new_src` a driver-level reproduction of the production calculation.
The newer driver's input tolerance also remains in its existing units; the
legacy driver multiplied its input tolerance by the bulk gap. These two
driver-level differences are intentionally recorded rather than hidden in the
solver adapter.

`new_src/iter.f90` now imports the adapter as a normal module instead of
textually including its source. Its undeclared scratch variables and the
obsolete `norm=0` assignment were also removed. `new_src/Makefile` now compiles
the shared kind and accelerator modules from `src` and declares the Fortran
module dependencies, permitting reproducible parallel builds.

### MPI defects exposed during integration

The first two-rank smoke run completed the Anderson step on rank zero and then
failed on another rank. This was traced to pre-existing broadcasts in
`new_src/mpicalls.f90`, not to the accelerator: the routine broadcast ten
elements from undeclared scalar `inmax` and `pmax` variables. Those overruns,
along with broadcasts of other unused or uninitialised names, were removed.
The remaining array counts now use either `size(array)` or the initialized
slice. `NcCutOff` was reconciled with the declared `NcCutof`, and
`implicit none` was added to the touched MPI units. An obsolete broadcast of
an undeclared `eps` at program startup was removed from `qcv.f90`.

### Verification

- A clean disposable build of `new_src` succeeds both serially and with
  `make -j10`.
- The changed adapter, dispatcher, MPI broadcast module, declarations, and
  main program compile under `-fimplicit-none`, interface/conversion warnings,
  and runtime checking. Only pre-existing real-to-complex and integer-to-real
  conversion warnings remain in the dispatcher and included BB routine.
- A deliberately tiny low-temperature calculation (`T=0.30`, A-core initial
  state, one trajectory, eight Ozaki frequencies) completes five preliminary
  fixed-point iterations and one Anderson iteration with both 2 and 10 MPI
  ranks. Both runs report the same first Anderson result:
  `max residual = 0.6175`, `p = 0.01`, and history size `1`.
- The normal and strict CMake configurations still pass all three tests,
  including the direct regression of the modern accelerator against the
  unchanged fixed-form `citerat` trace.

This verifies software integration and MPI consistency, not physical
convergence of the low-temperature A-core solution. The next convergence
benchmark should restore the legacy outer `pmax`/`inmax` schedule and run the
large-radius core/surface-separated case against archived production data.

## 2026-09-14 - Runtime Anderson `p_max` and distinct error diagnostics

A longer user run showed the adaptive mixing pinned at the adapter's compiled
limit `p=1.0`. The limit is now a runtime input named `aa_pmax` rather than a
parameter in `iter_AA.f90`:

- `new_src/qcv.inp` has an optional tenth line containing `AA p_max` after
  `itmax` (older nine-line inputs retain the default `1.0`);
- `init_calc.f90` reads and validates `aa_pmax >= 0.01` and prints the selected
  value at startup;
- `mpicalls.f90` broadcasts it to all ranks; and
- `iter_AA.f90` passes it to `legacy_anderson_t%initialize`.

This is one configurable cap for the current continuous AA stage. It does not
yet reintroduce the legacy driver's five separate `pmax`/`inmax` stages.

The first adapter revision returned the canonical accelerator's absolute
maximum residual in both `tolA` and `tolL`, causing the two printed errors to
be identical. The adapter now reproduces the diagnostics used by the other
`new_src` solvers for the same packed vectors:

```
tolA = ||G(x)-x||_2 / ||x||_2
tolL = N max_i |G_i(x)-x_i| / ||x||_2
```

These are the newer driver's global/average and scaled largest-component
measures; they are not the same as the original fixed-form driver's absolute
stopping residual. `qcv` owns the convergence decision using these values, so
the adapter disables the canonical iterator's internal absolute-tolerance
early return. Output now labels the pair as `error(avg,max)`.

### Verification

- A clean parallel `new_src` build succeeds.
- The touched units compile with implicit declarations disabled and runtime
  checking enabled. This also exposed and removed three remaining implicit
  local declarations in `init_calc.f90`; pre-existing conversion warnings
  remain.
- A 10-rank smoke calculation with `AA p_max=0.05` reaches and remains at
  `p=0.05`; the same calculation with `AA p_max=0.50` reaches `p=0.50` and then
  adapts down to `p=0.43`. This verifies that the selected input reaches the
  solver and is an active cap.
- The two reported errors are distinct in the smoke run (for example,
  `0.2114` and `49.68` on the first AA evaluation), and the normal and strict
  CMake test suites continue to pass all three tests.

## 2026-09-14 - `new_src` compiler and linker flags corrected

The active `new_src/Makefile` was reduced to flags and dependencies used by the
current solver. The previous line mixed compilation, checking, CPU threading,
and macOS framework options in every build even though the active source has
no OpenMP directives and makes no calls to Accelerate or Metal.

### Release build

The default `make -j10` and named `make -j10 release` builds now use:

```
-std=f2018 -fdefault-real-8 -fdefault-double-8 -fimplicit-none
-Wall -Wextra -Wimplicit-interface -O2
```

`-fdefault-real-8` remains a temporary compatibility requirement because most
of `new_src` still declares bare `real` and `complex` variables while its MPI
datatypes are double precision. `-fdefault-double-8` prevents GNU Fortran from
promoting explicit `double precision` declarations to 16 bytes as a side
effect. `-fimplicit-none` is now enabled project-wide and the complete active
source builds with it.

The release build intentionally omits fast-math, CPU-family-specific flags,
bounds checking, and accelerator-family flags. This keeps numerical behavior
conservative and the executable portable. OpenMP or GPU offload flags will be
introduced only with the corresponding source directives and regression
tests.

### Checked build

`make -j10 debug` cleans first and adds:

```
-O0 -g -fcheck=all -fbacktrace
-ffpe-trap=invalid,zero,overflow -finit-real=snan
```

Both named modes clean before rebuilding, preventing release and checked
objects from being mixed. `EXTRA_FFLAGS`, `EXTRA_LDFLAGS`, and `EXTRA_LDLIBS`
remain available for controlled experiments, and `make show-config` displays
the effective configuration.

Compilation and linking now consistently use the MPI Fortran wrapper. The
irrelevant `-fopenmp`, `-framework Accelerate`, `-framework Metal`, obsolete
vectorizer-verbosity option, and unused MKL fragments were removed. The
case-sensitive include spelling for `Iter_BR.f90` was also corrected for
non-macOS filesystems.

### Verification and warning retained

Clean release and checked builds both succeed, and a one-rank checked smoke
calculation completed five preliminary fixed-point iterations plus one AA
iteration without a bounds or floating-point trap. GNU Fortran still warns
that the local `sem(4,-mx:mx)` array in `getnewses.f90` is moved to static
storage. That is safe for the current one-thread-per-MPI-rank execution but
must be changed to per-thread/allocatable storage before adding threaded
trajectory evaluation.

The checked smoke calculation was run in `new_src` and therefore refreshed the
unsuffixed working outputs `op_xyz`, `op_harm`, `curr`, `xgrid.dat`, and
`error_log_AA.dat`. All suffixed historical/reference datasets were left
untouched; the previously unsuffixed working outputs were not retained.

## 2026-09-14 - Interpretation of the `p_max=5`, `md=10` run

The user reported an ongoing AA calculation in which the mixing parameter
mostly reaches `p_max=5` and the reported history size remains at `md=10`.
Inspection confirms that the subspace minimization is active:

- every AA evaluation inserts the new residual into residual-norm order;
- once the ten-vector history is full, the largest-residual retained vector is
  normally removed and the new vector is inserted, so the post-update size
  remains ten even though the contents turn over; and
- the Gram system for the residual-minimizing coefficients is rebuilt and
  solved on every iteration.

`md` is therefore the post-update history size, not a count of replacements or
a measure of subspace quality. It falls below ten only when a progress test or
a sufficiently bad new residual triggers a half-history rejection. The
reported trace itself shows active mixing adaptation before the cap dominates
(`5.00 -> 4.18 -> 2.90`, and later `5.00 -> 4.27`). From iteration 16 onward,
the uncapped formula evidently proposes values at least as large as five.

The global error falls from approximately `1.75e-2` to `8.64e-5` by iteration
44, while the scaled largest-component error falls from about `8.13` to
`5.94e-2` with substantial local fluctuations. This is convergence rather
than a frozen subspace, although the maximum-component trace is noisy.

One important limitation remains: `new_src` currently runs one continuous AA
stage. The original production driver ran five fresh `citerat` stages with
separate `pmax` and `inmax` pairs, resetting the history at each stage. That
outer restart/search schedule has not yet been restored. Future diagnostics
should expose the uncapped trial mixing, number of discarded vectors,
minimized-to-current residual ratio, and singular-system fallback flag; the
present `p` and `md` columns cannot show those mechanisms directly.

## 2026-09-15 - Project focus and staged roadmap

The project scope and development order are now defined in
`docs/PROJECT_PLAN.md`. The north-star target is a spinful quasiclassical
solver for general two-dimensional superfluid and superconducting devices,
with MPI-cluster and GPU execution. The immediate scientific acceptance target
is the free double-core vortex in 3He-B, not a general container.

The plan uses the fixed-form code as behavioral reference, `new_src` as the
transitional radial implementation, the attached `dcvlong.pdf` as the primary
double-core physics benchmark, SuperConga as an architectural comparison, and
the Seja-Lofwander DG method as the preferred but not yet accepted long-term
transport discretization.

The main implementation lane first builds a uniform-grid, full-spin 2D
CPU/MPI reference using the existing trajectory/Riccati physics. A bounded DG
lane begins after the field conventions are stable and must pass scalar,
two-dimensional, and matrix-spin equivalence tests before it can replace the
main transport method. Full adaptive infrastructure and GPU optimization will
not be implemented for both methods.

The project gates are now:

1. radial free/cylinder regression;
2. axisymmetric 2D embedding;
3. converged free double-core vortex with `a`, `C1`, and `C2` diagnostics;
4. evidence-based DG decision and adaptive-grid convergence;
5. CPU/GPU physical equivalence and measured acceleration; and
6. staged walls, interfaces, leads, weak links, and magnetic bilayers.

The previous cylindrical-trajectory request remains in G0 as a contained
regression task. It must not expand into general geometry work before the free
double-core milestone. A partially drafted trajectory module from the
interrupted implementation turn was removed rather than left unintegrated.

## 2026-09-15 - Unified run archive and panel summary

Extended `tools/run_and_plot.py` from a legacy-only wrapper to an isolated
runner for both `incoming/legacy-f77` and `new_src`. The `--solver new-src`
path copies the active free-form sources and shared modules into `work/` before
building, so test runs do not alter objects or unsuffixed results in `new_src`.
One-rank runs use Open MPI singleton execution; multi-rank runs retain the
ten-rank local default.

Each new archive now records the selected solver, input/Ozaki/executable SHA-256
hashes, a copied build log, mesh description, MPI command, timing, and terminal
process status. The legacy and new input formats are parsed separately. The
runner accurately warns that the legacy grid is fixed uniform and the
`new_src` tangent-mapped grid is graded but not error-adaptive.

Matplotlib output now includes `run_summary.png` and `run_summary.pdf`. The six
panels contain the magnitudes of all nine Cartesian order-parameter components,
the normalized quantity `Tr(A A^dagger)/3`, the `vx,vy,vz` solver fields, and a
linear-iteration/log-error convergence plot. The convergence parser recognizes
legacy NN/AA stages and modern NN, AA, BR, and BB output, with the active engine
and restarts marked on the plot. It falls back to the selected
`error_log_ENGINE.dat` when no captured standard output is available.

The density panel is intentionally labelled a pair-density proxy rather than
superfluid density. The latter is a response coefficient/tensor and cannot be
reconstructed from the current files. Likewise, `vx,vy,vz` are labelled a
current-related mean field because the code uses them to construct the
diagonal/Fermi-liquid self-energy; the physical current prefactor still needs
to be made explicit.

Verification performed:

- Python syntax and input/convergence parser checks pass.
- The `new_src` source snapshot builds successfully in an isolated directory
  with the current Fortran 2018 flags. The known `sem` static-storage warning
  remains.
- Completed legacy output and current `new_src` output both produce all detail
  plots plus the PNG/PDF panel summary without Qt or gnuplot.

## 2026-09-15 - G0 evidence, conventions draft, and milestone backlog

Added `benchmarks/radial_g0` with matched free/specular-cylinder smoke inputs
and matched full A-core reference candidates. Added
`tools/compare_radial_runs.py`, which compares `op_xyz`, `op_harm`, pair
amplitude, and `vx,vy,vz`, interpolating a candidate onto the reference radial
grid when necessary and emitting a JSON tolerance report.

Clean one-rank and ten-rank legacy smoke runs completed for both boundary
choices. Within the precision written to the result files, every compared value
is exactly equal between rank counts. The free and contained branches are not
accidentally identical: their relative L2 differences are approximately 0.1284
for the order parameter and 0.0733 for `curr` after normalizing its legacy
layout.

This does not close G0. The existing full cylinder run
`runs/20260914-150620-cylindrical` ended at average/maximum errors
`9.1833e-4/1.9316e-3`, well above its printed requested tolerance, despite a
zero process return code. It is now explicitly rejected as a canonical
baseline. Full converged free and contained references are still required, and
`new_src` currently reads `icyl` but always uses its free-vortex trajectory
routine.

Added `docs/architecture/MATHEMATICAL_CONVENTIONS.md`, cross-checked against
`dcvlong.pdf` and current source formulas. It separates transcribed source
content from proposed project conventions and leaves the length/energy units,
general tilde operation, propagator normalization, code-basis rotation, and
physical current prefactors open for explicit resolution.

Added `docs/MILESTONE_BACKLOG.md` with gate-linked task identifiers,
dependencies, acceptance evidence, and two focused review batches. The current
critical path is terminal convergence status, convention freeze, accepted
radial oracle runs, and only then the free/cylinder boundary strategy and the
uniform 2D kernel.

## 2026-09-15 - Shared order-parameter scale and collision-free plot names

The three spin-row panels in the run summary now use one common vertical scale,
computed from the largest magnitude among all nine Cartesian components. The
detailed Cartesian and harmonic component grids likewise use a common symmetric
scale across their nine panels, so amplitudes can be compared visually without
reading separate axes.

All derived plot, table, script, and manifest names now begin with
`yymmdd_case_initial-state`. Initialization targets are named `bulk-a`,
`bulk-b`, `n-core`, `a-core`, `d-core`, or `restart` from the archived input or
metadata. A serial suffix `_01`, `_02`, and so on is selected when a matching
stem already exists. Solver-native files (`op_xyz`, `op_harm`, and `curr`) keep
their historical names inside the already unique timestamped run directory so
restart and comparison compatibility is not broken.

The updated summary was regenerated for the current `new_src` A-core data and
visually checked in both PNG form and as a rendered one-page PDF. Python syntax,
date/case/initialization inference, and collision handling checks pass.

## 2026-09-15 - Colleague `scsolver.f90` Anderson audit

The requested `new_ses/sesolver.f90` path was not present. The only matching
colleague solver is `new_src/scsolver.f90`, authored in 2002--2005, and it was
reviewed as the likely intended file. It is not referenced by the current
Makefile or by any active source.

The central `MOBESC` update is a legitimate Anderson/Pulay-style affine
residual minimization. It ranks stored points by residual norm, discards stale
history, solves a small residual Gram system, and adapts the mixing parameter.
Its mixing update differs from the trace-validated `iter.f` translation: it
uses the new displacement estimate times a residual-ratio factor rather than
averaging that estimate with the previous mixing value. This makes it a useful
algorithmic benchmark, especially for the difficult low-temperature A-core
case, but not a drop-in replacement.

Direct integration is blocked by both source-modernization and correctness
issues. The file begins with a corrupted comment marker, is ISO-8859 text,
depends on default-real promotion and legacy syntax, hard-codes kind 8, and
stores all solver state in saved module variables. A strict Fortran 2018 syntax
check fails; after repairing the comment marker it compiles only with legacy
and default-real-8 compatibility flags.

More importantly, the documented residual is `g-G(g)` while every update adds
that residual; the actual caller convention must be established before any
comparison. The wrapper tests convergence using the proposed update size
rather than the true fixed-point residual, so a small damping value or stalled
step can give false convergence. The two-vector step has unguarded zero
denominators and no `pm_max` cap. The singular-system fallback repeatedly uses
only the newest residual instead of the indexed history vector and introduces
unseeded randomness. A negative LAPACK status can leave the next iterate
undefined. Several other norm divisions lack zero guards.

Decision: retain the file as provenance and an algorithm reference. For a fair
benchmark, extract its distinctive damping/history policy into the existing
typed accelerator interface, keep the nonlinear map and true residual
diagnostics identical, and compare map-evaluation counts, residual histories,
restarts, final fields, and MPI reproducibility against the current legacy-AA
engine. Do not link the historical module directly into `qcv`.

## 2026-09-15 - Working project name adopted

The project working name is **FermiForge**, with the descriptor "scalable
quasiclassical simulation of spinful superfluids and superconductors." Project
documentation now uses this name. Existing directory, executable, source, and
archived-run names remain unchanged to preserve scripts, paths, and numerical
provenance; those identifiers can be migrated deliberately when the package
layout is established.

## 2026-09-15 - Proposed first FermiForge 2D transition

The next proposed implementation milestone is a fixed, uniform Cartesian 2D
cross-section with the trusted radial nonlinear iterator retained initially.
The iterator will see the same physical components packed over
`npoint = nx*ny`, but the transport map requires new field-sampling and
trajectory-geometry layers; it is not only a vector-length change.

The first interpolation backend should be local bilinear (`Q1`) interpolation
of the complex Cartesian field components. It has constant-size four-node
stencils, exact reproduction of constant and affine fields, no high-order
overshoot, constant-time lookup on a uniform grid, and a layout suitable for
MPI and later GPU kernels. The sampling API must hide the backend so that an
element-local higher-order or adaptive finite-element evaluator can replace it
without changing the Riccati kernel. Complex components are interpolated
directly; phase and amplitude are not interpolated separately.

For a grid point `r0 = (x0,y0)` and Fermi-surface direction
`p = (px,py,pz)`, the sampled ray is `r_xy(s) = r0 + s*(px,py)`.
Translation invariance removes the spatial `z` coordinate, but `pz` remains in
the angular quadrature and in the contraction of the triplet order parameter
with momentum. Fixed-mesh ray locations and interpolation indices and weights
should be precomputed and reused for every Matsubara energy and nonlinear
iteration.

Acceptance is staged: manufactured interpolation tests; an embedded radial
field sampled along rays; checkpoint agreement for one Riccati trajectory; one
complete 2D map evaluation against the radial solver; and finally preservation
of an axial vortex under 2D relaxation. The double-core perturbation starts
only after those tests pass. Adaptive refinement will operate as an outer
solve-estimate-refine-transfer cycle and reset Anderson history whenever the
mesh changes.

## 2026-09-15 - First isolated 2D mesh and sampler implementation

Added strict-Fortran-2018 modules for a uniform Cartesian mesh, typed 2D
spinful fields, and local bilinear field sampling. The flat mesh index has `x`
vary fastest. The state stores the full complex `A(spin,orbital,point)` and
three real current-related mean fields. Interpolation stencils expose four
indices and weights for later trajectory precomputation. A separate sampler
contracts the interpolated order parameter with all three momentum components,
so translation invariance does not discard `pz`.

Added manufactured tests for point numbering, all four boundaries, explicit
outside detection, constant/affine/bilinear exactness, quadratic refinement,
cached-stencil equivalence, full-momentum contraction, and a smooth
singly-quantized axial-vortex field sampled along an oblique ray. The new
sources and tests pass strict syntax checking with Fortran 2018, implicit
typing disabled, conversion warnings enabled, and warnings promoted to errors.

Added an explicitly named `new_src` packing adapter. It preserves the current
component-major ordering and round-trips 21 real values per spatial point,
without importing the factor-of-ten convention from the older fixed-form
iterator.

The default developer-tool selection initially blocked object compilation
because the full Xcode licence has not been accepted. Selecting the already
installed Command Line Tools through `DEVELOPER_DIR` avoids that unrelated
toolchain path. A clean GNU Fortran 15.2 build then completed and all four test
executables passed, including the new Cartesian field suite and the existing
radial and Anderson suites. A second clean build with bounds and runtime checks,
floating-point traps, interface warnings, and conversion warnings also passed
all four tests. M1-01 is now `done`. No 2D code has yet been connected to the
trusted `new_src` Riccati or self-consistency calculation.

## 2026-09-16 - Axial reference embedding and cached 2D rays

Completed A-01 by transcribing all nine Cartesian-to-axial-harmonic formulas
from `new_src/init_calc.f90` into a pure explicit-kind module. Added known-field,
general complex round-trip, harmonic phase, and B-phase far-field tests. The
projection index order is explicit as `(+,0,-)`, and the source phase is
`exp(i*(M-sigma-ell)*phi)`.

The audit exposed a material distinction that must remain visible. The
harmonic-to-Cartesian formulas used by `new_src/interpol.f90` enforce
`Axy=-Ayx`; they are not the mathematical inverse of the general nine-component
forward transformation. FermiForge now exposes two explicitly named paths: a
true general inverse, and the source-faithful axial transport reconstruction.
The latter is used only by the radial-reference embedding. It has not been
silently generalized or used to constrain the unconstrained 2D state.

As a data check, the current `new_src/op_xyz` transformed back to harmonics
agrees with the separately printed `new_src/op_harm` to approximately
`1.01e-9`, consistent with output precision. The same data have a maximum
`|Axy+Ayx|` of approximately `5.04e-3` at radial point 33 (`r=5.159`), compared
with a maximum component magnitude of approximately `3.11e-1`. Thus the
restricted transport reconstruction can measurably project the printed
Cartesian data and cannot be treated as a harmless general inverse.

Added an axial radial-profile adapter that reproduces the source three-point
radial interpolation, including the negative-radius harmonic parity at the
origin, applies the vortex phase, rotates the azimuthal current-related field,
and embeds the result on a Cartesian mesh. A node beyond the available radial
support is reported explicitly rather than extrapolated. Manufactured parity,
embedding, oblique-ray refinement, and outside-support tests pass.

Added rectangular-domain straight-ray geometry for
`r_xy(s)=r0+s*(px,py)`. It retains the normalized full momentum including
`pz`, inserts the target exactly at `s=0`, bounds the path step independently
of the field mesh, and precomputes bilinear interpolation stencils. Rays
parallel to the invariant axis are represented by one spatial sample. Tests
cover entry/exit geometry, nonzero `pz`, the invariant-axis case, and cached
sampling of an affine spin-triplet field.

The strict runtime-checked GNU Fortran build now passes all seven tests with
bounds checks, floating-point traps, interface warnings, and conversion
warnings enabled. The next numerical gate remains M1-02: extract a
side-effect-free Riccati segment kernel and compare every propagation
checkpoint with the trusted radial source before attempting a 2D
self-consistency iteration.

## 2026-09-16 - Legacy Riccati kernel and 2D self-energy path

Completed M1-02 without importing `new_src` global state. The RK4 interval
routine in `new_src/riccati.f90:ricc` is now a pure explicit-kind procedure.
Its operation order and unusual half-step interpolation of linearly varying
coefficients are retained. A reference driver linked directly to the original
objects produced a nontrivial frozen interval result; the modern result agrees
within the rounding used by the checkpoint.

Added the source bulk-coherence initialization and a variable-interval
trajectory driver. A second original-source driver generated every forward and
reverse checkpoint for a five-point path. The modern entry-to-target and
exit-to-target branches reproduce all six stored source checkpoints. The
reference generators and exact compiler flags are retained under
`tests/reference` so that the origin of the frozen numbers is auditable.

Added the compatibility layer between cached 2D fields and propagation. It
constructs the triplet vector `d_alpha=sum_i A_alpha_i p_i`, projects the three
current-related fields onto the full momentum direction with an explicit
feedback scale, subtracts that value from the spectral energy, and reproduces
the distinct conjugation placement of source `setenergy1` and `setenergy2`.
The inactive spin-vector exchange self-energy remains explicit zero rather
than disappearing into an assumption. Affine manufactured tests cover both
the pair-field contraction and diagonal shift.

The strict clean build now compiles without warnings and passes all ten tests
with bounds checking, floating-point traps, interface warnings, and conversion
warnings enabled. M1-03 is the next gate: express and verify the 2 by 2
spin-matrix algebra, including normalization and tilde/particle-hole tests,
before using the path kernel in a complete 2D self-consistency map.

## 2026-09-16 - Git repository and guarded upload workflow

Prepared FermiForge for version control with a project-specific `.gitignore`.
The tracked set retains source, documentation, benchmark definitions, compact
comparison evidence, legacy reference inputs, and intentionally named
historical datasets. Local `work/`, `runs/`, and `tmp/` trees, compiler and
CMake products, active solver outputs, caches, executables, and credentials are
excluded without deleting them from the workstation.

Added `tools/upload_to_git.sh`. It supports either GitHub or GitLab SSH remotes
for `fogelstrom-lab/FermiForge`, runs the strict Fortran test gate by default,
checks staged whitespace and individual file size, verifies Git identity,
creates a commit, fetches an existing remote branch, and permits only a
fast-forward push. It contains no force-push path and stops on an unrelated or
diverged remote history. A dry-run mode performs verification and shows the
prospective import without staging, committing, or contacting a remote.

## 2026-09-16 - Runnable 2D demonstration and accessible field plots

Added a standalone Fortran driver and editable namelist for exercising the 2D
mesh, interpolation, field output, and visualization without representing the
result as a self-consistent solution. The driver constructs a smooth
manufactured singly quantized vortex, samples it at cell centers by default,
and writes a documented ASCII table containing `x,y`, the real and imaginary
parts of all nine Cartesian order-parameter components, pair density, and
three manufactured-current components.

Added a Python run wrapper that builds the driver, archives the input, logs,
hashes, data, and manifest in a unique timestamped run directory, then invokes
the plotter. Added a separate plot command for replotting any compatible 2D
map. The output contract is intended to remain valid when the manufactured
driver is replaced by the self-consistent map.

The amplitude figure is a 3 by 3 component grid with one shared white-to-indigo
scale. The phase figure displays `cos(arg(A_mn))` on a fixed blue-white-red
scale and marks phase-undefined low-amplitude regions grey. The observable
figure contains pair density, in-plane current magnitude with outlined vector
arrows, and signed axial current on its own symmetric scale. The palettes avoid
a red-green distinction and use strong lightness contrast and numerical
colorbars.

An end-to-end archived run produced the data and all six PNG/PDF files, which
were inspected visually. The strict runtime-checked build now passes all eleven
CTest cases, including a command-line smoke test of the new driver.

## 2026-09-16 - Normal Xcode toolchain restored

Confirmed that the Xcode licence is accepted and repeated a clean configuration
and strict GNU Fortran 15.2 build through the normally selected Apple developer
toolchain, without a `DEVELOPER_DIR` override. All eleven tests passed. Removed
the obsolete automatic Command Line Tools selection from the 2D runner and
renamed its default build directory to `work/fermiforge-2d-build`; the licence
workaround remains documented only as troubleshooting for another machine.
The revised one-command workflow then completed successfully and archived its
field map, logs, manifest, and six PNG/PDF figures under
`runs/20260916-162454-licensed-toolchain-check-2d`.

## 2026-09-16 - Explicit spin-matrix propagator started

Started M1-03 by adding explicit 2 by 2 Pauli and pair-potential algebra. The
forward conversion is exactly
`(gamma_0 I + gamma_alpha sigma_alpha) i sigma_y`; the reverse conversion
retains the distinct `Gt` signs in `new_src/riccati.f90:green` and is named as
a legacy reverse operation so that it cannot be mistaken for a fully resolved
general tilde convention.

Added reconstruction of the four 2 by 2 blocks of the dimensionless Nambu
propagator from independent forward and reverse coherence matrices. Tests now
cover coefficient round trips, the source triplet/singlet component layout,
matrix inversion, the spin-scalar limit, covariance of the triplet pair matrix
under an SU(2) spin rotation, the conjugate-coherence anomalous-block relation,
and `g_hat^2=-1`. The normalization test passes both manufactured nonsingular
coherences and coherence values frozen from the original radial trajectory.

A clean GNU Fortran 15.2 build with bounds checks, floating-point traps,
interface warnings, and conversion warnings completed without warnings. All
twelve tests pass. M1-03 remains active because the general momentum/energy
tilde and particle-hole convention still needs a source-backed acceptance
test. The next implementation slice is the pure one-target-point 3He
self-consistency integrand, using the explicit propagator without MPI or
nonlinear mixing.

## 2026-09-16 - One-trajectory 3He map contribution

Added a pure one-point accumulator for the nine complex gap-map components and
three real current-related mean fields. It receives the explicit spin-Nambu
propagator and preserves the exact component combinations and quadrature
factors in `new_src/riccati.f90:green`. The source `esum` factor is applied in a
separate finalization routine because the original applies it only to the gap
components.

Created a small driver that links the original source and freezes one
momentum/energy contribution. This exposed an existing language defect: the
legacy routine declares the accumulating `tem` array `intent(out)` while
reading its previous value. The modern code instead uses an explicitly
initialized `intent(inout)` accumulator and reproduces the zero-initialized
source result for all twelve outputs.

Connected the already tested 2D layers into a first complete vertical slice:
bilinear field sampling, triplet and current-feedback projection, two-sided
legacy Riccati propagation, explicit 2 by 2 coherence conversion, 4 by 4
propagator reconstruction, normalization diagnostics, and accumulation of one
quadrature contribution. An oblique-ray regression recreates the frozen source
coherences. During test construction it correctly detected an omitted diagonal
current-feedback field; including that field recovered the full reference
result.

A clean strict GNU Fortran 15.2 build completes without warnings and all
fourteen tests pass. No MPI, nonlinear update, or complete angular/energy loop
has been introduced. The next slice is a quadrature data model and a serial
one-target-point loop over every momentum direction and Ozaki pole, followed by
comparison with one radial `getnewop` point.

The user verification build inherited `-L/usr/local/opt/openblas/lib` in all
three CMake linker-flag cache entries. This was an obsolete Intel-Homebrew path,
not a FermiForge dependency or link failure; the installed Apple-silicon path
is `/opt/homebrew/opt/openblas`. Cleared the unused cached flags, rebuilt
without the linker warning, and repeated all fourteen tests successfully.

## 2026-09-16 - Complete serial one-point self-consistency map

Added typed angular and Ozaki quadrature objects and readers for the existing
`gauss11.dat` and `ozaki.dat` formats. The legacy conventions are explicit:
48 midpoint azimuths cover `[0,pi)`, the Gauss-Legendre polar weights are
halved, the direction ordering is azimuth outer/polar inner, and the source
gap prefactor is computed from the table poles and residues. Tests verify 528
directions, normalized weights, isotropic second moments, all eight poles, and
the table-derived prefactor.

Connected the complete serial one-target-point loop. A trajectory and its
interpolation cache are built once per direction; both Riccati branches are
then evaluated for every energy pole, reconstructed as an explicit spin-Nambu
propagator, accumulated, and finalized. Diagnostics report attempted,
accepted, and failed propagators, maximum normalization error, and separate
sampling, propagation, and accumulation times. The integration uses the
source energy-dependent substep floor `1+int(abs(pole))`, with a caller-set
minimum for refinement studies.

The real 48 by 11 angular table and all eight Ozaki poles were exercised on a
uniform B-phase field. All 4,224 contributions succeed, preserve rotational
component symmetry, produce zero equilibrium current, and satisfy the
propagator-normalization tolerance. The small manufactured case independently
expands the same loop through the lower-level trajectory connector.

## 2026-09-16 - Radial `getnewop` oracle comparison

Corrected the working `new_src/riccati.f90` accumulator declarations from the
historical invalid `intent(out)` to the intended `intent(inout)`; the numerical
statements are unchanged. Added an auditable legacy driver that evaluates the
complete original `getnewop` angular/energy loop at radial index 24, without
MPI reduction or nonlinear mixing. Its twelve reference values and exact
regeneration command are retained under `tests/reference`.

Added a source-format radial-state reader. It reconstructs the full-precision
tangent mesh hidden by the three-decimal `op_xyz` coordinate column, converts
the Cartesian axis values to the source axial harmonics, and retains the
positive-x-axis azimuthal mean field. The radial state is then embedded on a
Cartesian mesh and evaluated by the new point-map path.

At Cartesian spacing 0.25, the modern gap map differs from the independent
source result by `6.5113e-4` relative and the current-related mean field by
`2.5708e-5` absolute. Refining the spacing from 0.5 to 0.25 to 0.125 reduces
the relative gap error from `2.2945e-3` to `6.5113e-4` to `1.1983e-4`; the
mean-field error falls from `8.9095e-5` to `2.5708e-5` to `6.0659e-6`. This
identifies the remaining difference as Cartesian embedding/interpolation
error rather than a component, spin, or quadrature mismatch.

The default radial comparison is now an automated test with explicit
tolerances. A fresh GNU Fortran 15.2 build using warnings as errors, bounds
checks, floating-point traps, interface warnings, and conversion warnings
completed without warnings. All eighteen tests pass. The next implementation
slice is full-grid point batching with explicit ownership and reduction
interfaces, followed by packing the mapped field/residual for the existing
Anderson engine; the free-asymptotic endpoint policy remains a separate open
gate.

The transitional `new_src` executable also rebuilds successfully with the
Homebrew MPI wrapper after the accumulator-intent correction. Its previously
recorded `getnewses.f90` warning remains: the large local `sem` array is moved
to static storage. That warning is confined to the transitional radial code;
the modern library's fresh warnings-as-errors build is clean.

## 2026-09-16 - Full-grid map, MPI ownership, and Anderson connection

Extended the validated one-point map to an execution-neutral field-map layer.
It accepts an arbitrary, duplicate-checked list of flat Cartesian target
points, evaluates every direction and energy for each target, and returns a
sparse full-shaped mapped state plus a completion mask. The serial wrapper
evaluates all points. Aggregated diagnostics retain completed/failed point and
propagator counts, maximum normalization error, and sampling, propagation,
and accumulation times.

Added field residual diagnostics in the exact 21-value-per-point `new_src`
iteration layout. They report mean absolute, RMS, maximum absolute, relative
2-norm, legacy scaled maximum, and the spatial point with the largest local
RMS residual. Added an adapter that packs the current and mapped typed states,
applies the tested legacy Anderson implementation, and unpacks the next typed
state. Its first update reproduces the source 1% relaxation.

Added a separate MPI execution library using `mpi_f08`; MPI remains outside
the physics library. Points are assigned cyclically, while all directions and
energies for a target remain together. Each rank holds a replicated input
field, the complete mapped field is assembled by all-reduction, and rank zero
alone owns the Anderson history before broadcasting the next state and both
diagnostic reports.

The MPI regression uses all 528 directions and eight Ozaki poles at every
point and compares the gathered field with a serial evaluation. One, four,
and ten ranks give the identical checksum `6.26333449197233527`, maximum
normalization error `1.5543e-15`, and matching next Anderson state. Preliminary
strict-debug slowest-rank map times were 0.463, 0.133, and 0.064 seconds. This
is a functional scaling smoke test, not yet a production scaling report.

Added a runnable embedded-radial MPI driver and 17 by 17 example input. A full
one-iteration 10-rank run completed with slowest-rank map time 8.807 seconds,
RMS residual `1.4811e-3`, and maximum absolute residual `8.5121e-3`. It wrote
the common 2D field-map format, and the existing accessible Python plotter
successfully produced amplitude, phase, pair-density, in-plane mean-field,
and axial mean-field figures. The output is explicitly labelled an iteration
smoke test rather than a converged solution.

The strict warnings-as-errors build now passes all twenty-two tests, including
1-, 4-, and 10-rank MPI cases. The next physics-critical task remains the
controlled free-asymptotic trajectory endpoint, followed by an axially
symmetric full-grid convergence run before introducing a double-core seed.

## 2026-09-16 - Selectable free-vortex trajectory endpoint

Added an explicit free-vortex endpoint policy without changing the existing
local-endpoint regression. A ray can now be extended from the rectangular
field mesh to an enclosing outer circle. Exterior self energies are obtained
by radial matching to the numerical boundary field, followed by the Cartesian
tail powers already used in `new_src/interpol.f90`: mixed z/in-plane gap
departures decay as `1/r`, all other departures from phase-wound bulk B phase
decay as `1/r^2`, and the current-related mean field decays as `1/r`.

Separated the Riccati path propagator from the straight-trajectory stencil
object. The extended path therefore carries real exterior self energies and
coordinates rather than fake interpolation stencils. Its interior and tail
steps are independent controls. The bulk gap, winding, vortex center, outer
radius, and tail step are explicit input rather than hidden constants.

Tests now cover the source decay powers, exact continuity at the field
boundary, outer-circle geometry, preservation of the target and step bound,
uniform bulk-B invariance, and serial/MPI equality. The 1-, 4-, and 10-rank
MPI regression now exercises the free-vortex endpoint. All twenty-three tests
pass in the strict warnings-as-errors and runtime-checked build.

Added `examples/2d_mpi_free_vortex.nml` alongside the unchanged local endpoint
example. On the 17 by 17 embedded-radial field with ten ranks, one free-vortex
iteration to outer radius 24 gave RMS residual `1.1763e-3`, maximum residual
`6.3889e-3`, and slowest-rank map time 28.678 s. The matched local run gave
`1.4811e-3`, `8.5121e-3`, and 8.807 s. The maximum real/imaginary gap-column
difference between their output maps is `5.3006e-5`; this demonstrates that
the policy is active but is not yet an outer-cutoff convergence result.

The next gate is a controlled axial sequence in numerical box size, outer
radius, interior/tail trajectory steps, and Cartesian spacing. Only after that
sequence stabilizes the mapped field should multi-iteration axial convergence
and double-core seeding begin.

## 2026-09-17 - Normal-core 2D benchmark initiated

Started the first axial physics regression with the normal core, keeping the
driver core-independent so that an A-phase-core reference can use the same
path next. Added a dedicated MPI benchmark that reads a radial state, embeds it
on the full complex 3-by-3 Cartesian 2D field, evaluates the complete
spin-matrix quasiclassical map at every point, and records the residual against
the radial reference. It writes both input and mapped field maps plus
machine-readable metrics.

The radial reader now detects and reconstructs both source grid families. The
archived `op_xyz_n/curr_n` pair has 50 points on the older uniform radius-30
grid; the current `new_src` files have 100 points on a radius-70 tangent grid.
Their three-decimal coordinate columns are never used as the actual numerical
mesh. A regression covers both and rejects files matching neither family.

Factored the common 2D ASCII writer into `spinful_field_io_2d` and added an
organized Python runner. Its archive name contains `yymmdd`, case name, and
initialization target; the archive includes copied input, run/build logs,
hashes, metrics, input/mapped fields, and optional accessible Python plots.

The initial attempt to measure quarter-turn covariance by recovering axial
harmonics from an off-axis Cartesian field was rejected. The source transport
reconstruction enforces `A_xy=-A_yx`, whereas the stored radial values need not
satisfy that equality exactly. The projection is therefore lossy and cannot be
inverted off axis without changing the source convention. The benchmark now
compares every Cartesian point directly with the radial profile evaluated and
rotated by the authoritative source formulas. This decision is explicit in
the benchmark documentation.

The archived uniform normal-core run on a 17 by 17 grid and ten MPI ranks
passed its integration checks. The embedding errors were exactly zero, the
maximum propagator normalization error was `1.58e-15`, and the one-map
residuals were RMS `9.7592e-4`, maximum absolute `6.2002e-3`, and relative L2
`9.7204e-3`. The separately normalized maximum gap and current-field
differences were `2.1532e-2` and `2.7740e-2`. These are provisional historical
metrics, not final tolerances, because the data predate the current tangent
grid and circular asymptotic matching is still open.

A fresh isolated ten-rank `new_src` normal-core solve started from `istart=0`
on the current 100-point tangent grid with `T/Tc=0.30`, `F1s=5.4`, free-vortex
geometry, Anderson mixing, and `p_max=5`. It converged after five simple and 37
Anderson iterations. The final average and maximum residuals were `5.130e-9`
and `1.562e-6`, both below the declared `2e-6` tolerance. The input, complete
solver history, fields, executable hash, and file hashes are now the versioned
primary reference under
`benchmarks/normal_core_2d/reference/current_new_src_T0.30_F1s5.4/`.

The strict GNU Fortran build remains warning-clean and all 24 automated tests
pass, including one-, four-, and ten-rank MPI cases and the independent
normal-core point-map oracle. The next numerical gate is
the normal-core Cartesian-spacing, endpoint-radius, trajectory-step, and rank
matrix, followed by the identical benchmark for a converged A-phase core.

The 2D plotter now treats an axial current below `1e-12` of the in-plane scale
as numerical zero by default, labels its measured maximum, and displays it on
a nonsaturating scale. This prevents roundoff-level `j_z` in the normal core
from appearing as a physically structured red-blue field; the threshold is a
command-line control so genuinely weak axial currents in later cores remain
visible.

The primary 17 by 17 tangent-reference run gave RMS, maximum, and relative-L2
map residuals `9.8650e-4`, `6.2429e-3`, and `9.8111e-3`. Refining the same
half-width-eight box from spacing 1.0 to 0.5 only reduced relative L2 to
`9.5396e-3`; the maximum gap difference remained about `2.18e-2` relative,
although the current-field difference improved. Keeping 1,089 points but
doubling the half-width to 16 and outer radius to 48 reduced relative L2 to
`5.6912e-3`. This demonstrates a coupled core-resolution/domain requirement,
not a simple uniform-spacing error.

Added an axisymmetric radial-reference endpoint for benchmark use. Outside the
Cartesian box it samples the converged radial state directly to the source
matching radius `Rc=xgrid(94)=45.9993349989`, then applies the source `1/r` and
`1/r^2` tail. It is explicitly an oracle for axial cores, not a production
nonaxisymmetric boundary condition. The corresponding 17 by 17 run had
relative-L2 residual `9.7722e-3`, showing that the coarse full-grid discrepancy
is not primarily caused by the exterior tail implementation.

The decisive isolation test evaluates one point at radial index 24
(`r=3.60953947027`) in a large half-width-32 box with spacing 0.125. Against an
independent complete original `new_src:getnewop` map on the converged normal
core, the modern path differs by `2.0646e-5` maximum absolute gap,
`7.1891e-5` maximum relative gap, and `9.3457e-6` maximum absolute
current-field value; propagator normalization is `1.2213e-15`. The physics
kernel and conventions therefore reproduce `new_src`. The full-grid error is
the expected consequence of trying to combine a large incoming-trajectory
domain and coherence-scale core resolution on one uniform mesh. This makes the
next target a circular active region with fine core blocks and coarser outer
blocks; the same driver then extends directly to the A-phase core.

## 2026-09-17 - A-phase-core 2D benchmark

Accepted the user's clarification that the named 100-point
`op_xyz_a_converged` and `curr_a_converged` fields were first iterated to full
convergence and were subsequently reloaded with `istart=3` and zero requested
iterations as an immediate restart check. The restart input is preserved as
such and is not misidentified as the original generating input. Both field
files have identical tangent-grid coordinates from zero to radius 70.

Regenerated an independent original-`new_src:getnewop` result at radial index
24 directly from the named A-phase-core files. On a Cartesian embedding with
spacing 0.125, the modern 2D Fortran point path differs by `3.6068e-5` maximum
absolute gap, `1.1983e-4` maximum relative gap, and `6.0659e-6` maximum
absolute current-field value. Propagator normalization is `1.2212e-15`. The
A-phase-core source conventions and transport kernel therefore pass the same
independent oracle as the normal core.

The ten-rank 17 by 17 full-grid radial-reference run passed its integration
checks. Input embedding errors are exactly zero; the maximum normalization
error is `1.5543e-15`. The one-map residual is `1.1205e-3` RMS,
`6.3767e-3` maximum, and `1.1051e-2` relative L2. Separately normalized
maximum gap and current-field deviations are `2.3111e-2` and `4.6368e-2`.
These are uniform-grid convergence-study metrics, not final tolerances. The
result reinforces the need for separate resolution of the hard core, extended
soft-core texture, and outer trajectory domain.

The benchmark runner initially selected Apple Python for plotting even though
the invoking Homebrew interpreter had NumPy and Matplotlib. The physics run
was unaffected. Plotting was recovered from the archived fields, and the
runner now probes candidate interpreters for both dependencies before use.

## 2026-09-17 - First 2D multiscale mesh

Extended the Cartesian mesh to store explicit nonuniform axis coordinates
while retaining the uniform mesh as a special case. Added a symmetric
fine/medium/coarse constructor with exact origin and region-boundary nodes,
binary cell lookup, local cell widths in Q1 interpolation, and minimum/maximum
spacing diagnostics. Field storage and the 21-value-per-node nonlinear layout
remain contiguous and unchanged.

Added a circular active-point mask. Serial and MPI field maps evaluate only
active target points and copy the incoming state exactly through the inactive
outer halo. MPI work remains cyclic and balanced over the active list. The
halo is still sampled by trajectories, which makes the first axial benchmark
a circular relaxed region surrounded by a fixed radial-reference reservoir.
This is explicitly not a general nonaxisymmetric boundary law.

Residual reporting now accepts an active mask, preventing the frozen halo from
diluting convergence norms. The benchmark also divides the active region into
fine, medium, and outer zones. The Riccati step remains independent of local
field spacing. Mesh changes remain restricted to boundaries between nonlinear
stages and require trajectory-cache rebuild and Anderson-history reset.

The first normal-core multiscale run used half-width 12, region boundaries 3
and 7, target spacings 0.5, 1.0, and 2.5, and active radius 9. It stored 625
nodes and evaluated 429. Input embedding was exact, normalization error was
`1.78e-15`, and the active one-map RMS/maximum/relative-L2 residuals were
`1.0050e-3`, `6.2054e-3`, and `1.0575e-2`. The maximum current-field error
was reduced to `9.15e-3`. Fine/medium/outer relative-L2 residuals were
`9.78e-3`, `1.15e-2`, and `8.53e-3`.

The matched A-phase-core run also stored 625 nodes and evaluated 429. Its
normalization error was `1.33e-15`, with active RMS/maximum/relative-L2
residuals `1.0777e-3`, `6.3723e-3`, and `1.0848e-2`. Fine/medium/outer
relative-L2 residuals were `7.62e-3`, `1.21e-2`, and `1.09e-2`. The error has
therefore moved away from the fine hard-core region into the intermediate
soft-core band; the next convergence step should refine and extend that band
before further reducing the innermost spacing.

All 24 strict Fortran/MPI tests pass. New checks cover multiscale coordinate
construction, nonuniform bilinear exactness, circular masks, frozen-halo MPI
behavior, active-only residuals, and the existing one-, four-, and ten-rank
physics regressions. The implementation and its limitations are recorded in
`docs/architecture/MULTISCALE_RECTILINEAR_MESH.md`.

Added nodal control-volume weighting for comparisons between different point
distributions. The unweighted 21-value vector norm is retained because it is
the norm seen by Anderson; the area-weighted RMS and relative L2 norms are the
appropriate physical mesh-comparison diagnostics. A 25 by 25 uniform control
uses the same half-width, active radius, trajectory controls, and total stored
nodes as the multiscale case.

For the normal core, multiscale placement reduces the area-weighted fine-zone
relative L2 residual from `1.415e-2` to `9.844e-3` and the active-region value
from `1.030e-2` to `1.006e-2`. For the A-phase core, the fine and medium values
improve from `8.147e-3` and `1.290e-2` on the uniform control to `7.620e-3`
and `1.252e-2` on the multiscale mesh. Its coarse outer value worsens from
`1.004e-2` to `1.106e-2`, making the next adjustment unambiguous: retain the
fine core resolution but extend/refine the intermediate and outer active
regions rather than concentrating additional nodes at the origin.

## 2026-09-17 - Guarded GitLab publication helper

Configured `tools/upload_to_git.sh` to target
`git@gitlab.com:fogelstrom-lab/FermiForge.git` by default. The executable
helper initializes `main` when needed, runs the strict build and test gate,
stages only non-ignored project content, rejects unexpectedly large files,
shows the staged summary and exact destination, and requires confirmation
before committing and pushing. It checks remote ancestry and never
force-pushes. SSH and HTTPS GitLab authentication are selectable.

Extended the ignore rules to exclude extensionless executables and transient
products in nested legacy-source directories. Whitespace checks apply to the
maintained code and documentation while preserving imported legacy sources
and numerical reference data byte for byte.

## 2026-09-18 - GitHub established as the canonical repository

Confirmed that `fogelstrom-lab/FermiForge` is the existing GitHub project and
published the initial FermiForge import there. Changed the guarded publishing
helper and current documentation from the provisional GitLab default to
`git@github.com:fogelstrom-lab/FermiForge.git`; its HTTPS option now also
targets GitHub. The local `main` branch tracks `origin/main` without rewritten
history.

Recorded March 2027 as the start of collaborative FermiForge development. The
repository will be prepared before then with explicit licensing and provenance,
cross-platform onboarding, automated regression tests, contributor guidance,
and accepted two-dimensional normal- and A-phase-core workflows.
