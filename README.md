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

## Current status

The modern support library now contains both the dynamically sized radial
reference infrastructure and the first isolated two-dimensional field-sampling
slice. It provides:

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

The complete map is now distributed over a full Cartesian grid with MPI and is
connected to the legacy-compatible Anderson update. The trusted `new_src`
radial calculation remains the numerical reference while each layer is
validated independently. The normal- and A-phase-core benchmarks each perform
one full 2D self-consistency map; neither is yet a converged multi-iteration 2D
vortex. Their comparisons demonstrate that a single uniform grid cannot
efficiently resolve both the coherence-length core and the large outer domain.
The first static multiscale/circular-active backend is now working; automatic
block refinement, field transfer, and a converged multi-iteration solution are
the next production gates.

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
2D map, generates the accessible 3-by-3 order-parameter, density, and current
figures, and archives all inputs, logs, metrics, hashes, maps, and plots under a
collision-free timestamped directory in `runs/`. The canonical converged
radial input and the quantitative acceptance evidence are described in
`benchmarks/normal_core_2d/README.md`. This benchmark deliberately separates
verification of the transport/map kernel from convergence of a future
adaptive 2D relaxation.

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

Generated builds, local run archives, plots, executables, and compiler products
are excluded by `.gitignore`. Source, documentation, benchmark definitions,
compact comparison reports, legacy reference inputs, and intentionally named
historical datasets remain trackable.

The guarded upload helper runs the strict CMake/CTest suite, previews and
stages changes, rejects individual files larger than 25 MiB, asks for final
confirmation, commits, checks the remote branch for fast-forward ancestry,
and pushes without ever using a force option. Its default destination is the
canonical FermiForge GitHub repository:

```text
tools/upload_to_git.sh -m "Initial FermiForge import"
```

Use HTTPS instead of SSH if that is how GitHub authentication is configured:

```text
tools/upload_to_git.sh --https -m "Initial FermiForge import"
```

After `origin` has been configured, subsequent uploads need only a message:

```text
tools/upload_to_git.sh -m "Describe this reviewed change"
```

Once a local Git repository exists, use `--dry-run` to run verification and
preview changes without staging, committing, or contacting the remote. The
script stops rather than guessing if the remote already contains unrelated
commits.

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
