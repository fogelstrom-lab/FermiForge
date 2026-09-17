# Running the first MPI 2D self-consistency iteration

## Scope

This program embeds the trusted axial `new_src` radial state on a uniform 2D
Cartesian grid, evaluates the complete 3He self-consistency map at every grid
point, and applies the modernized legacy Anderson mixer. Target points are
distributed cyclically over MPI ranks; every momentum direction and energy
for one target remains on the same rank.

The supplied input performs one iteration on a 17 by 17 nodal grid. It is an
integration and scaling smoke test, not a converged 2D solution and not yet a
double-core calculation. The stored three-component vector is the
current-related mean field used by self-consistency, not the physical mass
current.

There are now two endpoint choices. `examples/2d_mpi_smoke.nml` uses the
original local mesh-edge initialization. `examples/2d_mpi_free_vortex.nml`
continues every non-axial ray to a caller-selected outer circle. Outside the
box it radially matches the numerical boundary field, applies the same
Cartesian `1/r` and `1/r^2` departure-from-bulk powers as the free-vortex tail
in `new_src/interpol.f90`, and approaches an explicit phase-wound bulk B phase.

Run all commands from:

```text
/Users/mikael/Documents/Codex/3he-vortex-modernization
```

## Build

```text
cmake -S . -B work/mpi-2d-build \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_BUILD_TYPE=Release
cmake --build work/mpi-2d-build --parallel 10 \
  --target fermiforge_2d_mpi_smoke
```

## Run with ten MPI processes

```text
mkdir -p work/my-mpi-smoke
mpiexec -n 10 work/mpi-2d-build/fermiforge_2d_mpi_smoke \
  examples/2d_mpi_smoke.nml work/my-mpi-smoke/fields_2d.dat \
  | tee work/my-mpi-smoke/run.log
```

If Open MPI cannot infer local process slots on macOS, use:

```text
mpiexec --host localhost:10 --map-by ppr:10:node --bind-to none -n 10 \
  work/mpi-2d-build/fermiforge_2d_mpi_smoke \
  examples/2d_mpi_smoke.nml work/my-mpi-smoke/fields_2d.dat \
  | tee work/my-mpi-smoke/run.log
```

Plot the output with the established accessible palettes:

```text
/opt/homebrew/bin/python3 tools/plot_2d_fields.py \
  work/my-mpi-smoke/fields_2d.dat \
  --output-dir work/my-mpi-smoke/plots \
  --prefix mpi-smoke
```

To exercise the free-vortex continuation instead, substitute the second input:

```text
mpiexec --host localhost:10 --map-by ppr:10:node --bind-to none -n 10 \
  work/mpi-2d-build/fermiforge_2d_mpi_smoke \
  examples/2d_mpi_free_vortex.nml \
  work/my-mpi-smoke/free-vortex-fields-2d.dat
```

## Input controls

The namelist controls the domain, cell count, trajectory step, minimum Riccati
substeps, boundary-relaxation distance, number of nonlinear iterations,
Anderson history and maximum mixing, source data paths, and output path.
`endpoint_policy` is either `local` or `free_vortex`. The free-vortex choice
also requires the vortex center, bulk gap, winding, outer radius, and maximum
tail step. Its outer circle must strictly enclose every corner of the Cartesian
box. Keep `iteration_count=1` for the first check. Increasing it is a solver
experiment, not evidence of convergence, until outer-radius/domain/step
convergence and gauge/orientation constraints are completed.

Each iteration reports RMS and maximum absolute residuals, the worst spatial
point and its RMS residual, Anderson mixing/history, and the slowest-rank map
time. A converged production driver will additionally write a structured
iteration-history file and an explicit terminal status.

The first matched 17 by 17, ten-rank comparison gave the following one-step
diagnostics in the strict build:

| endpoint | RMS residual | maximum residual | slowest rank |
|---|---:|---:|---:|
| local | `1.4811e-3` | `8.5121e-3` | 8.807 s |
| free vortex, outer radius 24 | `1.1763e-3` | `6.3889e-3` | 28.678 s |

The maximum change in any real/imaginary gap output column was `5.30e-5`; the
maximum current-related mean-field change was `7.97e-5`. These numbers show
that the endpoint is active and quantify its initial cost. They are not a
cutoff-convergence result.
