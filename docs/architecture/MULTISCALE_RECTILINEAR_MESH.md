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
condition. A future nonaxisymmetric calculation must replace the radial halo
with a controlled asymptotic, material, or domain-decomposition policy.

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

## Initial benchmark

The first normal- and A-phase-core runs use a square halo of half-width 12,
zone boundaries at radii 3 and 7, target spacings 0.5, 1.0, and 2.5, and a
circular active radius of 9. There are 625 stored nodes, of which 429 are
active and 196 form the frozen halo.

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
