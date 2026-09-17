# Free-vortex endpoint benchmark

This benchmark isolates the effect and cost of the first controlled outer
continuation for the two-dimensional trajectory solver.

## Inputs

- Local endpoint: `examples/2d_mpi_smoke.nml`
- Free-vortex endpoint: `examples/2d_mpi_free_vortex.nml`
- Initial field: `new_src/op_xyz` and `new_src/curr`
- Grid: 16 by 16 cells, 17 by 17 nodes, half-width 8
- Angular/energy quadrature: 528 directions and eight Ozaki poles
- MPI ranks: 10
- Free-vortex outer radius: 24
- Interior and exterior maximum trajectory step: 0.25

## First one-iteration comparison

| endpoint | RMS residual | maximum residual | worst point | slowest-rank map time |
|---|---:|---:|---:|---:|
| local | `1.4811e-3` | `8.5121e-3` | 289 | 8.807 s |
| free vortex | `1.1763e-3` | `6.3889e-3` | 226 | 28.678 s |

Across the written maps, the maximum change in an individual real or
imaginary order-parameter column is `5.3006e-5`, the maximum pair-density
change is `1.6804e-6`, and the maximum current-related mean-field change is
`7.9747e-5`.

This is activation evidence, not a converged physics result. Acceptance of the
endpoint requires a sequence in Cartesian box size, outer radius, interior
trajectory step, exterior trajectory step, spatial resolution, and angular
and energy quadrature. The local endpoint stays in the suite as a regression
control.
