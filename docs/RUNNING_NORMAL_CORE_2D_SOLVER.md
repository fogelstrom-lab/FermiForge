# Running the self-consistent normal-core 2D solver

## What is implemented

`benchmark_axisymmetric_core_2d` can now perform a complete repeated
self-consistency solve on a uniform or static multiscale Cartesian mesh. Each
nonlinear evaluation:

1. distributes the active Cartesian points cyclically over MPI ranks;
2. samples the full complex 3 by 3 order parameter and three Fermi-liquid
   current-related mean fields along every trajectory;
3. propagates the spin-matrix Riccati amplitudes and accumulates the 3He map;
4. applies the translated legacy Anderson accelerator on rank zero; and
5. broadcasts the next field to all ranks.

The inactive halo is held at the embedded converged radial `new_src/qcv`
solution. This supplies a controlled axisymmetric reference boundary for this
benchmark. It is not the eventual boundary condition for a double-core vortex.
The mesh is fixed during a solve; automatic adaptation is a later stage.

## Quick verification

From the FermiForge project directory, run:

```text
/opt/homebrew/bin/python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_normal_core_convergence_smoke.nml \
  --case-name normal-core-solver-smoke \
  --initialization radial-qcv-reference \
  --ranks 10
```

This deliberately uses a very loose tolerance and is only an integration
test. It should report `terminal status: converged`, produce initial/final
field plots and a logarithmic convergence plot, and archive all inputs,
outputs, logs, and hashes below `runs/`.

The strict automated suite contains the same ten-rank solver path:

```text
ctest --test-dir work/strict-point-map-build --output-on-failure
```

## Production normal-core solve

The first physically meaningful multiscale run is:

```text
/opt/homebrew/bin/python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_normal_core_multiscale.nml \
  --case-name normal-core-multiscale-converged \
  --initialization radial-qcv-reference \
  --ranks 10
```

The supplied input allows 60 evaluations, uses the same `2e-6` unnormalised
maximum-component stopping tolerance as the radial calculation, and uses the
legacy Anderson controls `md=10`, progress threshold `0.1`, and `p_max=5`.
This is a substantially longer calculation than the quick verification.

Each run directory contains:

- `input_fields_2d.dat`: the exact initial 2D field;
- `mapped_fields_2d.dat`: the last quasiclassical map;
- `final_fields_2d.dat`: the final accepted/mixed state;
- `checkpoint_fields_2d.dat`: the latest restart state;
- `iteration_history.dat`: residuals, Anderson state, normalization, and time;
- `metrics.txt`: terminal status and comparisons with the radial qcv field;
- `run.log`, `input.nml`, `manifest.json`, and accessible plots.

For both the initial and final fields, the plotter writes common-scale 3 by 3
amplitude and phase-cosine figures in two representations:

- Cartesian components `Axx,...,Azz`;
- spherical/harmonic components `A++,...,A--`, ordered exactly as in
  `new_src/op_harm`.

The harmonic figures are obtained with the authoritative transformation in
`src/order_parameter_basis.f90`. The plotting implementation checks pointwise
that the transformation preserves the total order-parameter norm.

`current_mean_field` and the plotted `j_x,j_y,j_z` columns retain the present
code's current-related Fermi-liquid mean-field convention. They should not yet
be presented as an independently normalized physical mass-current observable.

## Restarting

To continue a stopped run, point the runner at its checkpoint:

```text
/opt/homebrew/bin/python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_normal_core_multiscale.nml \
  --case-name normal-core-restart \
  --initialization saved-2d-checkpoint \
  --restart runs/YYMMDD-HHMMSS-CASE/checkpoint_fields_2d.dat \
  --ranks 10
```

The restart reader verifies the complete point count and every mesh coordinate.
The Anderson history is intentionally reset, because only fields—not an
accelerator subspace—are valid restart state.

## Computational-cell comparison

The requested cell-size test is organized by
`benchmarks/normal_core_2d/cell_study.json`. Run both supplied cells with:

```text
/opt/homebrew/bin/python3 tools/run_normal_core_cell_study.py
```

The study first solves a half-width-12 cell with active radius 9 and then a
half-width-18 cell with active radius 15. Both retain identical fine and
medium meshes inside radius 7. The generated CSV and JSON reports contain:

- convergence and normalization checks for each cell;
- area-weighted final-field errors relative to the converged radial qcv state;
- direct final-field differences on common Cartesian nodes inside radius 7;
- iteration counts, mesh sizes, timings, and links to both archived runs.

The provisional tolerances are explicit in the JSON file and can be tightened
after the first full production pair. A result is not accepted merely because
the executable reached its iteration limit: the required terminal status is
`converged`.

To change the cell study, copy the JSON file and pass the copy with `--config`.
This keeps the comparison itself reproducible and prevents command-line
settings from being lost.
