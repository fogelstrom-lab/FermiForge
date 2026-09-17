# Current `new_src` A-phase-core reference

This directory preserves the 100-point radial A-phase-core state supplied by
the user for the 2D regression. The controls recorded by the restart check are
`T/Tc=0.30`, unit winding, free-vortex geometry, `F1s=5.4`, 48 azimuths,
Gauss-11 polar integration, eight Ozaki poles, and a declared tolerance of
`2e-6`.

The state was first iterated to full convergence. The supplied input then
loaded the resulting `op_xyz` and `curr` through `istart=3` with zero requested
iterations to verify immediate restart behavior. Consequently,
`restart_check.inp` is provenance for the reload check, not the generating
iteration history.

`point_24.expected` is an independent complete original-`new_src` map of this
state at radial index 24. `reference-point-24.txt` retains the complete output
of its generating driver. The modern 2D point path is compared against this
fixture on a fine Cartesian embedding.
