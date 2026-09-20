# First 2D multiscale mesh

Status: implemented and benchmarked, 2026-09-17.

## Scope

The first multiscale backend is a symmetric, nonuniform tensor-product
Cartesian mesh. It introduces three independently controlled spatial zones:

- a fine central region for coherence-length core structure;
- an intermediate region for the A-phase/soft-core structure; and
- a coarse outer region used as a trajectory-sampling halo.

This is an incremental backend, not the final quadtree/AMR implementation. It
was selected because it exercises nonuniform lookup, interpolation, MPI work
distribution, fixed-boundary handling, residual localization, output, and
plotting through the existing full spin-matrix physics path without creating a
second transport solver.

## Coordinate construction

The origin and all requested region boundaries are exact mesh nodes. Within
each positive-axis segment the number of cells is the ceiling of the segment
width divided by its requested maximum spacing. The segment is then divided
uniformly and reflected about the origin. Consequently, actual cell spacings
never exceed their requested fine, medium, or coarse targets.

The mesh stores explicit contiguous coordinate arrays. Cell lookup uses a
binary search, and Q1 interpolation uses the local cell widths. Uniform meshes
remain a special case of the same type and retain their existing numbering and
iteration-vector layout.

## Circular active region and fixed halo

The physical update can be restricted to nodes inside a circle. MPI distributes
only those active target points. Nodes outside the circle remain available to
trajectory interpolation but are copied unchanged from the incoming state to
the mapped state. For the axial benchmarks this creates a circular relaxed
region surrounded by a fixed radial-reference halo inside the rectangular
storage mesh.

This is a benchmark boundary policy, not a general container boundary
condition. The nonaxisymmetric free-vortex calculation will replace the fixed
radial halo with the angle-dependent asymptotic continuation described below.

Inactive nodes have exactly zero nonlinear residual. Residual reports can take
an active mask, and the benchmark additionally reports fine-, medium-, and
outer-zone norms so a large number of quiet halo nodes cannot hide a poorly
resolved core or texture.

Two norm families are reported. The ordinary 21-value-per-node norm is the
Euclidean vector norm used by Anderson. Mesh-to-mesh physics comparisons use
nodal control-volume weights formed from half the adjacent cell widths in each
axis. This prevents a fine region from receiving extra importance merely
because it contains more nodes.

## Separation from trajectory integration

The Riccati integration step remains an independent input. A coarse outer
field cell therefore does not force a large propagation step, and a fine core
cell does not force that step along the complete ray. Cached trajectory
stencils carry local nonuniform interpolation weights while the propagation
coordinate remains independently bounded.

## Nonlinear iteration rule

The mesh and active mask are fixed during one Anderson solve. Frozen-halo
values pass through the nonlinear map unchanged, so their residual and update
are zero. Any later mesh change invalidates both cached stencils and Anderson
history; fields must be transferred and the accelerator reset before
restarting the solve.

## Angle-dependent asymptotic exterior

For a free nonaxisymmetric vortex, points outside the self-consistent matching
surface do not remain frozen. Following the double-core calculation, the
order parameter will be continued as

```text
A(r,phi) = A0(phi)
         + A1(phi) (Rc/r)
         + A2(phi) (Rc/r)^2 + O(r^-3),
```

The coefficient functions are evaluated independently at each spatial azimuth.
Equation (20) of `dcvlong` fixes their Cartesian block structure: `A1` has
only the mixed z/in-plane entries, while `A2` has the in-plane 2 by 2 block
and the zz entry. Thus each order-parameter component has one allowed leading
power, not a fitted mixture of both. The functions are not restricted to the
two-constant hydrostatic `C1,C2` form. The Fermi-liquid current-related mean
field has the corresponding, unrestricted large-distance form

```text
nu(r,phi) = nu1(phi) (Rc/r) + nu2(phi) (Rc/r)^2 + O(r^-3).
```

Equivalently, the factors of `Rc` may be absorbed into coefficients multiplying
`1/r` and `1/r^2`. The order parameter is matched continuously at `Rc` using
the single power permitted for each block. The mean field is sampled at two
radii, `Rin` and `Rc`, and both coefficients are solved independently for all
three real components. This retains the slow A-core channels without imposing
axial symmetry on a future double-core state.

The endpoint has two source modes. The radial-reference mode extracts the
coefficients from a converged `new_src` profile and is the controlled normal-
and A-core validation path. The state-fit mode extracts them from the current
2D field at the same polar angle; its regression test uses different
coefficients at `phi=0` and `phi=pi`. The latter is the mechanism intended for
nonaxisymmetric iteration.

For an anisotropic active domain, `Rin(phi)` and `Rc(phi)` are radial
intersections with two nested ellipses. Both surfaces lie a configurable
number of local mesh intervals inside the active ellipse. This keeps every
interpolation stencil in the self-consistently updated region. A circular
special case remains available for axisymmetric calculations.

The outer nodes are dependent asymptotic values rather than independent
Anderson unknowns. After each accepted update and before the next trajectory
map, the benchmark regenerates every node outside the active ellipse from the
fitted expansion. The current axisymmetric
validation deliberately refits the preserved radial solution; switching the
driver to the state-fit source will regenerate the coefficients from each new
interior iterate for the double core.

Two-radius interpolation is the minimum identifiable fit for the Fermi-liquid
field and makes continuity easy to test. It should eventually be upgraded to
an overdetermined annular least-squares fit with radial-window and angular
resolution diagnostics; the Eq. (20) order-parameter block constraints must
remain explicit in that upgrade.

## Initial benchmark

The first normal- and A-phase-core runs use a square halo of half-width 12,
zone boundaries at radii 3 and 7, target spacings 0.5, 1.0, and 2.5, and a
circular active radius of 9. There are 625 stored nodes, of which 429 are
active and 196 form the dependent asymptotic halo in the A-core test.

For the A-phase core the one-map relative L2 residuals are approximately
`7.62e-3`, `1.21e-2`, and `1.09e-2` in the fine, medium, and outer active
zones. The largest relative error has therefore moved to the intermediate
soft-core region. The next static refinement study should reduce the medium
spacing and extend its radius before reducing the innermost spacing.

Against a uniform 25 by 25 control with identical storage extent and active
radius, multiscale placement improves the area-weighted A-core fine and medium
relative L2 residuals from `8.15e-3` and `1.29e-2` to `7.62e-3` and
`1.25e-2`. The coarse outer zone degrades from `1.00e-2` to `1.11e-2`; the
next parameter study must therefore extend or refine that zone.

## Remaining path to adaptive blocks

The tensor-product mesh still refines complete coordinate bands and therefore
becomes inefficient for displaced half cores, walls, and devices. The next
backend should reuse the validated contracts while adding:

1. dense rectangular blocks refined by factors of two;
2. finest-valid-block interpolation and cached block/stencil identifiers;
3. conservative field restriction and tested prolongation;
4. gradient, curvature, and local-residual indicators;
5. two-to-one balancing and refinement hysteresis; and
6. solve-estimate-transfer-restart cycles with Anderson history reset.

The present implementation supplies reference behavior for each of those
steps and remains useful for axial regression and controlled convergence
studies.
