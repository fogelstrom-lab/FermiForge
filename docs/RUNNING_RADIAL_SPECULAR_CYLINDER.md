# Rotationally symmetric specular confinement in new_src

The active radial solver now selects free-vortex sampling for `icyl=0` and
reflected cylinder sampling for `icyl=1`. This restores the legacy contained
transport path as a bounded first step toward cylindrical confinement in the
modern 2D solver. The full 2D solver has not been changed.

## Input

The first ten `new_src` records retain their existing meanings. For `icyl=1`,
append an eleventh record containing the positive finite cylinder radius after
AA `p_max`. For example:

```text
0.30                  : T/Tc
1.0                   : vorticity
1                     : icyl
5.40                  : F1s
2                     : azimuthal directions (smoke only)
1.e-2                 : residual tolerance
1                     : A-core seed
3                     : Anderson mixing
0                     : main iteration limit
1.0                   : AA p_max
4.0                   : cylinder radius
```

The radius uses the existing solver's dimensionless length units. Their physical
conversion remains an open convention question; this input does not label the
radius as `xi_0`. This syntax differs from the older fixed-form solver input.
Free-vortex inputs need no added record.

The supplied file is `examples/radial_specular_cylinder_smoke.inp`. It is not a
converged benchmark. In particular, `qcv` still runs up to five startup NN
updates even when its main iteration limit is zero. To exercise the complete
runner when that short calculation is desired:

```sh
python3 tools/run_and_plot.py --solver new-src \
  --input examples/radial_specular_cylinder_smoke.inp \
  --case-name radial-specular-cylinder --ranks 1
```

## Transport and interpolation

For confined calculations, the 100 radial nodes uniformly cover `0 <= r <= R`.
The integration step is `R/99`; the existing 800 samples on each half-ray give
a half-path length of approximately `8.08 R`. Free calculations retain their
original tangent-mapped grid, integration step, and asymptotic tail.

`intord_c` calls the module in `3DFS_MPICodes/trajectories.f90` to construct both
half-rays, including specular side-wall reflections. `poretraj.f90` is an
interactive geometry demonstration, not a library routine used by the solver.
The sampled momentum changes at reflections; its axial component is conserved
apart from the legacy epsilon handling of exactly transverse rays. Backward
sample momentum follows the legacy orientation toward the target.

At each sample the existing axial harmonic reconstruction and quadratic
interpolation supply all nine order-parameter components and the azimuthal
current-related mean field. Their contractions use the local reflected
momentum. No free-vortex extrapolation is used. The final three radial nodes
provide the wall interpolation stencil, avoiding an out-of-range node.

The Riccati calculation follows the reflected path without restarting the
coherence amplitudes at a wall. This is the source-compatible representation
of spin-independent specular continuity. The order parameter itself is not
forced to retain its bulk value at the wall; its surface response must emerge
from the self-consistency calculation.

## Checks and current limits

Run the bounded standalone checks:

```sh
python3 tests/check_cylindrical_trajectories.py
python3 tests/check_radial_cylinder.py
```

The geometry check compares the free-form routines with both original sources.
The radial check builds with bounds checks and floating-point exception traps,
compares reflected samples with an analytic quadratic gap and linear azimuthal
mean field, checks exact wall interpolation, and evaluates finite maps at two
radii using two azimuths, eleven polar nodes, and one energy pole. It also
rejects missing/invalid radii and a restart on an unrelated radial grid. These
standalone checks are separate from the existing CMake/CTest suite.

During implementation, free-vortex sampled self energies and one-pole
propagator contributions at four radial nodes were compared byte for byte with
an isolated pre-change build and agreed exactly. The complete radial executable
was built in an isolated directory. No production confined solution has yet
been accepted.

Before treating confined output as a physics benchmark:

- Check finite reflected-path length and endpoint-initialization independence.
  The original bulk-coherence start and constant-endpoint relaxation remain;
  they are initialization devices, not reservoirs attached to the cylinder.
- Check integration-step convergence at reflection discontinuities. The current
  routine samples at fixed arc-length intervals, and Riccati interpolation can
  straddle a reflection; collision-aligned intervals are a future refinement.
- Check the wall-node interior displacement (`1e-7` integration steps), grazing
  rays, angular/energy quadrature, radial resolution, and restart sensitivity.
- Compare a converged axisymmetric cylinder against the legacy contained
  reference, including surface suppression and current-related fields.

Cylinder restarts must already use the same 100-node uniform radial grid and
radius. Coordinate checks allow for the existing three-decimal text output;
changing the radius requires explicit field transfer, not relabelling old data.

After those gates, the 2D extension can reuse the reflected geometry with the
modern Cartesian field sampler, a circular physical-domain mask, explicit wall
stencils, and a cylinder boundary option replacing the free-vortex endpoint.
The full three-dimensional momentum quadrature remains necessary even when the
fields are translationally invariant along the cylinder axis.
