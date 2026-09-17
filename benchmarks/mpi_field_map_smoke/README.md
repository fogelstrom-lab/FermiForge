# MPI field-map smoke benchmark

The deterministic test problem uses a 4 by 4 nodal field, 528 momentum
directions, eight Ozaki poles, and a nonuniform manufactured spinful state.
Every MPI result is compared directly with the serial full-field map. It also
checks a root-owned first Anderson update and broadcast of the next state.

Strict debug-build results on the development Mac were:

| MPI ranks | points per rank | field checksum | slowest-rank map time |
|---:|---:|---:|---:|
| 1 | 16 | 6.26333449197233527 | 0.463 s |
| 4 | 4 | 6.26333449197233527 | 0.133 s |
| 10 | 1-2 | 6.26333449197233527 | 0.064 s |

All rank counts had maximum propagator-normalization error `1.5543e-15` and
agreed with the serial state to the test tolerances. These timings establish
functional point-level parallelism only; they are too small for a production
scaling claim.

A separate 17 by 17 embedded-radial smoke run with 10 ranks, the full
quadrature, and one Anderson update took 8.807 s on its slowest rank. Its first
map residual was RMS `1.4811e-3`, maximum absolute `8.5121e-3`. The result is
not converged and is retained only as end-to-end evidence.
