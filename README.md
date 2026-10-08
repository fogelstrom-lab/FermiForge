# FermiForge

FermiForge is a scalable quasiclassical simulation framework for spinful
superfluids and superconductors. This workspace contains an immutable legacy
Fortran 77/MPI reference and a modern Fortran implementation being developed
alongside it.

The canonical repository is
<https://github.com/fogelstrom-lab/FermiForge>. It is being prepared as the
shared technical foundation for collaborative FermiForge development beginning
in March 2027.

The staged scientific and software roadmap is maintained in
`docs/PROJECT_PLAN.md`. Its immediate target is a validated two-dimensional,
spinful calculation of the free double-core vortex in 3He-B; general device
boundaries and GPU execution follow explicit validation gates.

The current equations, grid and core-seed figures, and complete standalone
run instructions are collected in the
[FermiForge implementation note](output/pdf/fermiforge_algorithm_and_architecture.pdf).
Its editable LaTeX and rebuild instructions are in
[docs/technical_note](docs/technical_note/README.md).

For alternating Anderson and simple updates, including a restart launcher for
the latest localized 0+ run, see [Cycled iteration](docs/CYCLED_ITERATION.md).
An independent safeguarded [BB engine](docs/BB_ITERATION.md) is also available.
The [Polyak momentum engine](docs/POLYAK_ITERATION.md) provides a matched radial
comparison using the method in SuperConga Appendix E.
For fast axisymmetric accelerator tests with only a radial ray of independent
points, see [Radial symmetry mode](docs/RADIAL_SYMMETRY.md).
NAISS preparation: [resource estimate](docs/hpc/NAISS_ALLOCATION_ESTIMATE.md),
[technical appendix](docs/hpc/NAISS_TECHNICAL_APPENDIX.md), and
[portable MPI scaling benchmark](docs/hpc/SCALING.md).

## Current status

Start with the [hands-on framework map](docs/FRAMEWORK_GUIDE.md), updated
6 October 2026: driver selection, algorithm flow, source navigation, input
controls, short exercises, diagnostics and present limitations.

The modern Fortran framework supports iterated full-2D free vortices, a cheaper
radial-symmetry mode using the same point-map physics, and a radial specular
cylinder driver. Python manages runs and plots; it is not the transport solver.
Free-vortex accelerator choices include Anderson, BB and Polyak. The cylinder
currently uses Anderson. General 2D walls, device leads, DG and GPU execution
remain development targets.

An experimental [unconstrained Cartesian/annular disk test](docs/UNCONSTRAINED_DISK_TEST.md)
now connects specular trajectories to interior-only quadratic sampling and
independent 2D iteration. General-geometry boundaries remain future work;
the disk backend still requires physical refinement benchmarks.

Standalone A-texture trials (`a_mermin_ho`, spin-texture `a_panam`, and the
original mixed A/polar `a_planar`) are available via
`bash tools/run_annular_A_texture.sh NAME`. On Linux omit the macOS-only
`caffeinate` prefix. The runner builds locally; no Mac build or local reference
run is needed for these texture trials. Radial-reference comparison launchers,
in contrast, need their referenced data directories transferred separately.

The following is an inventory of foundations and historical validation steps,
not a list of restrictions to the earliest single-map implementation:

- an exact uniform representation of the legacy 49-cell/50-point radial grid;
- core- and surface-localized adaptive refinement;
- two-to-one balancing between neighboring radial cells;
- binary cell lookup;
- quadratic interpolation on nonuniform nodes;
- dynamically sized storage for the complex 3-by-3 order parameter and
  three-component exchange field;
- tested packing and unpacking in the exact order expected by the legacy
  nonlinear iterator;
- a uniform Cartesian 2D mesh with deterministic flat point numbering;
- typed storage for `A(spin,orbital,point)` and the three real current-related
  mean fields;
- boundary-safe bilinear interpolation with reusable four-node stencils; and
- contraction of the interpolated order parameter with the full
  `(p_x,p_y,p_z)` momentum direction;
- exact packing and unpacking of the `new_src` component-major iteration layout
  with 21 real values per 2D node;
- tested Cartesian/axial-harmonic conversion, with the restricted legacy
  transport projection kept separate from the general inverse;
- source-faithful embedding of an axial radial profile on Cartesian nodes,
  including origin parity and vortex phase factors; and
- rectangular-domain straight-ray construction with bounded path steps and
  cached interpolation stencils;
- explicit trajectory self energies matching the active `new_src` triplet and
  current-feedback conventions; and
- a side-effect-free legacy Riccati interval/trajectory kernel checked against
  original forward and reverse propagation values at every checkpoint;
- explicit 2 by 2 spin-matrix and pair-potential algebra, including the exact
  legacy forward/reverse coherence conventions; and
- reconstruction of the full 4 by 4 spin-Nambu propagator with tested scalar,
  spin-rotation, conjugation, and normalization identities;
- a source-frozen one-sample 3He gap and current-feedback integrand; and
- an end-to-end one-trajectory path from interpolated 2D fields through
  two-sided Riccati propagation to a mapped mean-field contribution;
- typed legacy-compatible angular and Ozaki quadratures, including the
  source energy-dependent integration rule; and
- a complete serial one-point map over 528 directions and eight energy poles,
  checked against an independently generated radial `getnewop` value under
  Cartesian-grid refinement;
- a field-wide 2D map with RMS, mean-absolute, maximum-component, and
  worst-spatial-point residual diagnostics;
- cyclic MPI distribution over complete target points with replicated mapped
  fields and verified 1-, 4-, and 10-rank reproducibility;
- a root-owned Anderson update followed by broadcast of the next 2D state;
- a selectable free-vortex endpoint that continues rays outside the Cartesian
  box, matches the numerical boundary field continuously, applies the source
  `1/r` and `1/r^2` tails, and approaches an explicit bulk B-phase value;
- a source-preserving reader for both historical uniform and current
  tangent-mapped `new_src` radial profiles;
- a nonuniform fine/medium/coarse Cartesian mesh with local Q1 interpolation,
  binary cell lookup, a circular active region, and a fixed sampling halo;
- active-region MPI distribution and fine/medium/outer residual diagnostics;
- an axisymmetric radial-reference endpoint used as a benchmark oracle;
- converged normal- and A-phase-core references with automated radial-point
  oracles; and
- an organized full-grid MPI benchmark runner with accessible plots.

The trusted `new_src` radial calculation remains a numerical reference while
layers are validated independently. Subsequent work extended the original
single-map benchmarks to multi-iteration calculations, smooth graded grids,
quadratic sampling, evolving asymptotics, restarts and spatial residual
diagnostics. The grid remains structured Cartesian in full-2D mode; universal
automatic error-controlled refinement is not yet implemented. See the framework
guide for the distinction between tested cases and planned capabilities.

## Run the isolated 2D demonstration

From the project directory, run:

```text
/opt/homebrew/bin/python3 tools/run_2d_demo.py
```

This builds the manufactured-vortex driver, archives a collision-free run, and
generates accessible amplitude, phase, density, and current figures. See
`docs/RUNNING_THE_2D_DEMO.md` for input parameters, direct build commands, and
the output format. This demonstration tests infrastructure; it is not yet a
self-consistent quasiclassical result.

## Run the normal-core benchmark

From the project directory, run the organized ten-rank benchmark with:

```text
/opt/homebrew/bin/python3 tools/run_axisymmetric_core_benchmark.py \
  --case-name normal-core \
  --initialization current-new-src-tangent \
  --ranks 10
```

The runner builds the benchmark in Release mode, performs the MPI-distributed
2D map, generates accessible 3-by-3 Cartesian and spherical/harmonic
order-parameter figures together with density, current, and convergence
figures, and archives all inputs, logs, metrics, hashes, maps, and plots under
a collision-free timestamped directory in `runs/`. The default
uniform input retains the historical one-map comparison. The multiscale input
now runs a complete repeated self-consistency solve with the translated legacy
Anderson engine and restart checkpoints. See
`docs/RUNNING_NORMAL_CORE_2D_SOLVER.md` for the quick test, production solve,
restart, and two-cell qcv comparison.

Run the corresponding A-phase-core comparison with:

```text
/opt/homebrew/bin/python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_a_phase_core_benchmark.nml \
  --case-name a-phase-core \
  --initialization converged-radial-reference \
  --ranks 10
```

Its reference provenance and interpretation are documented in
`benchmarks/a_phase_core_2d/README.md`.

The first static multiscale comparisons use a fine core, an intermediate
soft-core region, a coarse outer halo, and a circular active update region:

```text
/opt/homebrew/bin/python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_a_phase_core_multiscale.nml \
  --case-name a-phase-core-multiscale \
  --initialization converged-radial-reference \
  --ranks 10
```

The design, boundary policy, and path to block-structured AMR are documented
in `docs/architecture/MULTISCALE_RECTILINEAR_MESH.md`.

Run the complete two-cell normal-core comparison with:

```text
/opt/homebrew/bin/python3 tools/run_normal_core_cell_study.py
```

It converges both cells, compares each with the radial `new_src/qcv` reference,
directly compares their common inner-grid values, and writes CSV and JSON
acceptance reports. The supplied calculation is intentionally a long
production run, not a smoke test.

## Rotationally symmetric cylindrical confinement

The active `new_src` radial solver supports `icyl=1` with a cylinder radius
appended after AA `p_max`. It uses the modernized specular trajectory routines
and interpolates the mean fields along the reflected rays. This is an
experimental axisymmetric path; the full 2D confinement backend and converged
surface benchmark remain to be validated. Input, short checks, and numerical
limits are in `docs/RUNNING_RADIAL_SPECULAR_CYLINDER.md`.

## Import a legacy full-2D state

Large legacy `op_x`, `op_y`, `op_z`, and `curr` datasets are kept outside Git.
When a local `2D_benchmarks` archive is available, inspect and convert a state
with `import_legacy_split_field_map_2d`. The reader infers the rectilinear
mesh, preserves all complex order-parameter and Fermi-liquid mean-field
components, and checks the redundant norm columns. See
`docs/IMPORTING_LEGACY_2D_FIELDS.md` for the command and the explicitly
different legacy harmonic convention. The same note documents
`tools/plot_legacy_double_core_reference.py`, which creates full-cell and
core-window versions of the modern reference plot set directly from the
ignored archive.

The first nonaxisymmetric transport check uses the complete imported field for
trajectory interpolation but maps only a configurable set of core and outer
probe points. Run `benchmark_legacy_double_core_map_2d` with
`examples/2d_double_core_one_map.nml`; interpretation and the initial ten-rank
result are documented in `docs/RUNNING_DOUBLE_CORE_2D_BENCHMARK.md`.

## Run the first double-core 2D continuation

With the local `2D_benchmarks` archive available, the first testable staged
continuation is one command:

```text
tools/run_first_double_core_2d.py
```

It builds the code, runs a ten-rank `0.4 -> 0.6 xi0` active-radius ladder,
transfers sparse state between the two radii, projects gauge/translation/
orientation tangent modes, regenerates the dependent asymptotic halo after
each accepted update, and writes a timestamped report and plot beneath
`runs/`. It is a two-update integration test, not a converged double-core
solution. Detailed controls, expected values, and restart instructions are in
`docs/RUNNING_DOUBLE_CORE_2D_BENCHMARK.md`.

## Start a double-core candidate from scratch

To start instead from the archived converged double-core solution, use
`python3 tools/run_converged_reference.py --iterations 20`. Omitting
`--iterations` measures the initial map defect without an update;
`--probe-only` provides a quick double-core check. The same wrapper accepts
`--core normal` and `--core a-phase` for the preserved converged `new_src`
radial references. See `docs/RUNNING_CONVERGED_REFERENCES.md` for frozen and
evolving exterior comparisons, mesh controls, and interpretation.

The modern solver can also generate a regularized London/two-half-core seed on
a newly constructed multiscale mesh, with zero initial Fermi-liquid mean field:

```text
tools/run_double_core_from_scratch.py
```

This path reads no archived 2D field. It begins MPI/Anderson iteration, writes
restartable state and machine-readable metrics, and produces the Cartesian,
harmonic, density/current, convergence, and symmetry-axis order-parameter
plots. The ordinary command now requests up to 20 updates and writes a sparse
restart checkpoint after every completed update; use `--iterations 2` for the
short integration test. Reaching the iteration limit is not the same as a
converged vortex. The seed equations, verified result, restart command, and
escalation path are documented in
`docs/RUNNING_DOUBLE_CORE_FROM_SCRATCH.md`.

For the first larger branch-stability calculation, use
`caffeinate -i tools/run_double_core_overnight.py`. This preset expands the
cell to `+/-40 xi0`, resolves the central region at `0.4 xi0`, relaxes a disk
of radius `22 xi0`, and starts the hard cores at `y=+/-10 xi0`. It records the
measured half-core separation after every update so collapse toward a
single-core state is visible during the run.

After that run, `caffeinate -i tools/run_double_core_adaptive_domain.py`
measures directional interface and boundary-collar diagnostics, enlarges only
the failing axes, and stops when the calculated state is contained. The
selected extent therefore responds to temperature and Fermi-liquid parameters
instead of using prescribed physical radii. Restart and diagnostic controls
are in
`docs/RUNNING_DOUBLE_CORE_FROM_SCRATCH.md`.

The original `nop`, `aop`, and `dop` initialization formulas are also
available as `historical_nop`, `historical_aop`, and `historical_dop`. The
large controlled comparison with the historical double-core seed is
`caffeinate -i tools/run_historical_dop_overnight.py`; see
`docs/HISTORICAL_CORE_SEEDS.md` for the exact conventions.

## Build and test

Configure with a Fortran 2018 compiler, then build and run the tests:

```text
cmake -S . -B work/modern-build -G Ninja \
  -DCMAKE_Fortran_COMPILER=/path/to/gfortran
cmake --build work/modern-build --parallel
ctest --test-dir work/modern-build --output-on-failure
```

The regression suite includes the embedded radial one-point comparison. Its
reference provenance and spatial-refinement evidence are documented in
`benchmarks/radial_point_map/README.md`.

The first complete MPI map/Anderson smoke driver is documented in
`docs/RUNNING_THE_2D_MPI_SMOKE.md`. It produces files accepted by the existing
accessible 2D plotter, but its supplied one-iteration run is not a converged
physical solution.

Two endpoint inputs are supplied. `examples/2d_mpi_smoke.nml` retains the
earlier local-endpoint regression, while `examples/2d_mpi_free_vortex.nml`
enables the controlled outer continuation. The latter is now the relevant
starting point for domain-, cutoff-, and trajectory-step convergence studies;
one iteration must still not be described as a converged vortex.

If a different Mac reports that its selected full Xcode installation has an
unaccepted licence, either accept that licence or temporarily select the
standalone Command Line Tools before configuring:

```text
export DEVELOPER_DIR=/Library/Developer/CommandLineTools
```

## Version control and upload

Generated builds, local run archives, plots, executables, compiler products,
and LaTeX intermediates are excluded by `.gitignore`. Source documentation is
tracked together with the canonical rendered PDF under `output/pdf`; transient
PDFs beside the LaTeX source are not tracked. Benchmark definitions, compact
comparison reports, legacy reference inputs, and intentionally named
historical datasets remain trackable.

The guarded upload helper requires the MPI execution layer, runs the strict
CMake/CTest suite, validates the maintained Python tools, previews and stages
changes, rejects individual files larger than 25 MiB, asks for final
confirmation, commits, checks the remote branch for fast-forward ancestry,
and pushes without ever using a force option. Its default destination is the
canonical FermiForge GitHub repository. The portable top-level shortcut is:

```text
./push_to_github.sh -m "Describe this reviewed change"
```

Use HTTPS instead of SSH if that is how GitHub authentication is configured:

```text
./push_to_github.sh --https -m "Describe this reviewed change"
```

The maintained implementation can also be called directly:

```text
tools/upload_to_git.sh -m "Describe this reviewed change"
```

Once a local Git repository exists, use `--dry-run` to run verification and
preview changes without staging, committing, or contacting the remote. The
script stops rather than guessing if the remote already contains unrelated
commits, and reports explicitly whether a local commit was created before any
failure. `--allow-no-mpi` permits an intentionally reduced test build, but it
should not be used for the normal FermiForge publication gate.

## Run and plot the cylindrical reference

The Python runner builds either the legacy solver or the transitional
`new_src` solver in an isolated work directory, launches ten MPI processes by
default, and archives the executable, input, build/run logs, hashes, and
metadata. Plotting uses NumPy and Matplotlib through a noninteractive backend,
so it does not require Qt or a display server. These packages are already
available in the Homebrew Python installation on the current Mac.

```text
python3 tools/run_and_plot.py
```

Run the current free-form source instead with:

```text
python3 tools/run_and_plot.py --solver new-src --case-name radial-aa
```

A fast end-to-end workflow check, not a physics benchmark, is available as:

```text
python3 tools/run_and_plot.py \
  --input tests/data/qcv_workflow_smoke.inp \
  --case-name workflow-smoke
```

Use a different input or MPI process count with `--input` and `--ranks`. Plot
existing output without running the solver using:

```text
python3 tools/run_and_plot.py --plot-only incoming/legacy-f77
```

The final summary PNG and PDF contain six panels: the magnitudes of all nine
Cartesian order-parameter components on a common vertical scale, a normalized
pair-density proxy, the current-related mean-field components, and a semilog
convergence history shaded and labelled by iteration engine. Derived filenames
use `yymmdd_case_initial-state_kind`, for example
`260915_radial-aa_a-core_run_summary.png`. If that stem already exists, a
serial suffix is added rather than overwriting it. The pair-density panel is
deliberately not labelled superfluid density: a true superfluid-density
response tensor is not present in the current output.

Gnuplot remains available only as an explicit fallback through
`--plot-backend gnuplot`. If another Python installation is used, install the
declared packages with `python3 -m pip install -r requirements.txt`.

Every new run is written below `runs/` with a timestamp. Metadata identifies the
legacy fixed uniform mesh or the `new_src` static tangent-mapped mesh; neither
must be mistaken for an error-driven adaptive-mesh result. Compare two radial
runs with `tools/compare_radial_runs.py`; the G0 inputs and current evidence are
documented in `benchmarks/radial_g0/README.md`.

## Layout

- `incoming/legacy-f77/`: immutable reference input, source, and data.
- `src/`: modern Fortran modules.
- `tests/`: unit and regression tests for modern components.
- `benchmarks/`: documented physics benchmark definitions.
- `docs/architecture/`: design decisions and migration constraints.
- `docs/architecture/MATHEMATICAL_CONVENTIONS.md`: source-preserving spin,
  Nambu, trajectory, phase, and normalization conventions.
- `docs/MILESTONE_BACKLOG.md`: gate-oriented executable backlog.
- `work/`: local build and reproducibility work areas.
- `MODERNIZATION_LOG.md`: chronological record of work and open questions.
