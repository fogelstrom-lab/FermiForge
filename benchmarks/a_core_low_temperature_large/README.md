# Low-temperature A-core benchmark

Status: benchmark design; the adaptive mesh is not yet connected to the full
legacy physics solver.

## Purpose

This benchmark will test the reported convergence difficulty of the A-phase
vortex core at low temperature while separating the core and surface boundary
layers sufficiently that their overlap is numerically negligible.

## Initial parameter family

- Cylindrical symmetry retained.
- A-core initial condition (`istart = 1` in the legacy input).
- Temperature `T/Tc = 0.10`, selected because an uploaded Ozaki table exists
  at this temperature. This is provisional and does not assert that 0.10 was
  the temperature of the historical convergence problem.
- Primary radius `R = 80` in the legacy solver's length units.
- Radius-isolation checks at `R = 60` and `R = 100`.
- Initial fine radial spacing target `0.125` near the core and surface.
- Initial coarse bulk spacing limit `4.0`.
- Initially refined widths of `8.0` from the origin and the surface.

These values are starting points for a convergence study, not accepted
production parameters. In particular, the correspondence between the legacy
length unit and the chosen coherence-length convention must be confirmed.

## Required comparisons

1. Repeat each radius at successively smaller fine spacing.
2. Compare the full complex order-parameter components in a fixed core window.
3. Compare the same components as functions of distance from the outer wall in
   a fixed surface window.
4. Monitor the order-parameter and Fermi-liquid self-consistency residuals
   separately in the core, bulk, and surface regions.
5. Declare core-surface overlap negligible only when increasing the radius no
   longer changes both localized solutions within a recorded tolerance.
6. Record iteration histories so that apparent convergence caused by global
   averaging cannot hide a poorly converged core or surface cell.

## Implementation sequence

1. Verify the modern mesh in uniform mode against the legacy 50-point grid.
2. Replace uniform radial lookup with nonuniform lookup while holding the
   physical point set uniform.
3. Make field arrays and nonlinear iteration vectors dynamically sized.
4. Enable the core/surface refined mesh and transfer all fields consistently.
5. Decouple spatial spacing from the Riccati integration step.
6. Run the temperature, radius, resolution, and nonlinear-mixing study.
