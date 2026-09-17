# Normal-core point-map comparison

The complete original `new_src` angular/Ozaki map was evaluated independently
at radial index 24 of the converged normal-core reference
(`r=3.60953947027193633`). The modern 2D point map embedded the same radial
state in a square of half-width 32 with Cartesian spacing 0.125.

| metric | difference |
|---|---:|
| maximum absolute gap-component error | `2.0646e-5` |
| maximum relative gap-component error | `7.1891e-5` |
| maximum absolute current-field error | `9.3457e-6` |
| maximum propagator normalization error | `1.2213e-15` |

This is the decisive normal-core kernel comparison: component conventions,
spin-matrix propagation, quadrature, gap prefactor, and current feedback agree
with the source at the tested point. The much larger full-grid residual on the
17-by-17 and 33-by-33 uniform meshes is therefore a domain/resolution problem,
not a normal-core physics mismatch.

The full-grid experiments expose the reason an adaptive mesh is required. A
box large enough to reproduce the incoming trajectories must extend through
the long radial tail, while the vortex core still needs coherence-scale
spacing. Refining only the small box from spacing 1.0 to 0.5 barely changed the
relative L2 residual (`9.811e-3` to `9.540e-3`). Doubling the box half-width
from 8 to 16 at 1.0 spacing reduced it to `5.691e-3`, but did not remove the
largest local component difference. A globally uniform mesh that combines
half-width 32 with spacing 0.125 would require 513 by 513 points and is not an
appropriate production route.

The next implementation target is consequently a circular active domain with
fine core blocks and coarser outer blocks, while retaining the radial-reference
endpoint as an axisymmetric oracle. The same point and field comparisons then
transfer unchanged to the A-phase core.
