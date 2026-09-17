# GPU and Multiscale Design Direction

Status: provisional architecture decision, 2026-09-14.

This note records requirements and design constraints for modernizing the
superfluid 3He vortex solver. It is not a change to the physical model or to
the legacy reference implementation.

## Requirements

1. Preserve the present cylindrical calculation as a numerical reference.
2. Make GPU execution possible without maintaining a separate physics code.
3. Retain MPI for multi-node calculations.
4. First support general two-dimensional container cross sections with
   translational symmetry along the cylinder axis and a replaceable boundary
   scattering model, beginning with specular reflection. Preserve a route to
   later fully three-dimensional containers.
5. Resolve coherence-length structure at vortex cores and surfaces while also
   treating the much longer spin-orbit texture self-consistently.
6. Do not impose a bulk-B-phase rotation parameterization in regions where the
   full complex order parameter may leave that manifold.

## What the legacy loop tells us

`getnewop` currently distributes radial evaluation points over MPI ranks and,
for every point, loops over momentum directions and Ozaki energies. Each call
to `makeprops` performs two causal Riccati integrations along a trajectory and
then contributes to a small local self-consistency sum.

The useful accelerator decomposition is therefore:

- one logical work item per evaluation point, momentum direction, and energy;
- sequential propagation along the individual ray inside that work item;
- hierarchical reduction over energy and direction into the order-parameter
  and Fermi-liquid fields.

For the supplied input this already amounts to 50 evaluation points times 48
azimuthal directions times 11 polar quadrature directions times 8 energies,
or 211,200 independent `makeprops` calls per nonlinear iteration.

The sequential integration along one ray is not an obstacle to GPU use because
the number of independent rays is large. The legacy code should first be
re-expressed in this form on the CPU and checked against the reference output.

## Programming model and data layout

The physics layer should use modern standard Fortran with explicit kinds,
modules, explicit interfaces, allocatable contiguous arrays, and small
side-effect-free numerical kernels. GPU directives belong in a thin execution
layer rather than in field definitions or physics interfaces.

The baseline accelerator route will be standard OpenMP target offload, with a
CPU implementation of the same kernels. No GPU vendor API, device data type,
or library will appear in a physics-module interface. Compiler- and
vendor-specific settings belong in build presets and the thin execution layer.
An OpenACC or CUDA Fortran implementation may be used later as an optional,
measured optimization, but it must not become the only implementation of the
Riccati or self-consistency mathematics.

OpenMP source portability does not imply identical compiler support or
performance on every device. Continuous verification should therefore include at least
one CPU compiler and, when access is available, one compiler for each intended
GPU family. Device-specific tuning is permitted behind common interfaces;
device-specific physics code is not.

Fields and ray batches should be stored in flat, contiguous arrays whose
leading dimension corresponds to the accelerator lane for the hot loops.
Pointer-linked cells, polymorphism inside kernels, per-ray allocation, and
derived types containing many separately allocated components should be kept
out of the device path. The final ordering will be selected using measured
CPU and GPU access patterns, not assumed in advance.

Geometry-dependent information is constant during a nonlinear field update.
Boundary intersections, surface normals, reflection events, cell traversal,
and interpolation metadata should therefore be built or cached separately
from the Riccati propagator and reused across iterations. Storage must be
bounded: cache compact ray topology and weights in batches rather than an
unbounded dense table for every possible ray sample.

## First symmetry-relaxed problem

The first general-container problem has translational symmetry along the
cylinder axis. The fields are functions of `(x,y)`, but the quasiparticle
momentum remains three-dimensional. In particular, the quadrature over `p_z`
must be retained. A trajectory with momentum `(p_x,p_y,p_z)` follows its
projected path in the cross section, while `p_z` changes the in-plane speed per
unit trajectory length.

For a wall normal lying in the cross section, specular reflection changes only
the in-plane momentum components; `p_z` is conserved. Directions nearly
parallel to the cylinder axis require an explicit numerical policy rather than
the small-number substitution used in the legacy cylindrical trajectory
routine.

This first problem therefore uses:

- a two-dimensional block hierarchy rather than a three-dimensional octree;
- vortex lines parallel to the invariant axis;
- two-dimensional boundary intersections and normals;
- the full three-dimensional momentum quadrature and full complex 3-by-3 order
  parameter.

The interfaces will nevertheless use dimension-independent concepts such as
points, directions, intersections, and field components so that a later 3D
geometry backend does not require a new transport solver.

## Spatial discretization

The recommended first mesh is a block-structured Cartesian quadtree/AMR
hierarchy, normally refined by factors of two, with the physical cross section
represented by an embedded boundary or signed-distance description.

An implemented intermediate backend now uses explicit nonuniform rectilinear
coordinates with fine, medium, and outer zones plus a circular active-region
mask. It exercises nonuniform interpolation, active-only MPI work distribution,
frozen-halo sampling, and area-weighted cross-mesh diagnostics without changing
the contiguous field layout. This is a validation bridge, not a replacement
for the block hierarchy described below; details and benchmark results are in
`MULTISCALE_RECTILINEAR_MESH.md`.

This choice offers:

- regular, dense blocks suitable for GPU kernels;
- local coherence-length resolution at vortex cores and surfaces;
- coarse coverage of slowly varying bulk texture;
- straightforward prolongation, restriction, and multilevel corrections;
- a geometry interface that does not assume a cylinder.

Mapped multiblock grids may be valuable for selected smooth container families,
but should be an optional geometry backend. A fully unstructured mesh is not
the first choice because irregular gathers and ray-element traversal make both
the quasiclassical sweep and GPU execution substantially harder.

The field mesh and the one-dimensional Riccati integration mesh must be
distinct concepts. A ray should choose integration steps from local field
variation and boundary proximity while sampling the spatial field through a
documented interpolation operator. It must not inherit a single global step
size from the finest spatial cell.

Initial refinement indicators should include:

- distance to a vortex line and to a material boundary;
- gradients or curvature of the complex order-parameter tensor;
- the local self-consistency residual;
- gradients of a diagnosed spin-orbit rotation or texture field.

The number of points per coherence length and per texture length will be set by
convergence studies. They should be input policies rather than compile-time
constants.

## Full order parameter and slow texture

The safe initial representation is the full complex tensor
`A(spin, orbital, point)` on the adaptive hierarchy. Far from cores and walls,
the solution may be diagnosed in terms of amplitude, phase, and a spin-orbit
rotation. That decomposition is useful as a coarse correction or
preconditioner, but it should not constrain the fine-grid tensor near defects.

The weakly varying spin-orbit sector can otherwise dominate nonlinear
convergence. The solver should therefore support a multilevel update:

1. evaluate quasiclassical and self-consistency residuals on the active mesh;
2. relax short-scale amplitude and core/surface components on refined blocks;
3. restrict the long-wavelength rotational residual to coarser levels;
4. solve or precondition the texture correction there;
5. prolongate the correction and repeat until both fine- and coarse-scale
   residuals converge.

Anderson or Broyden mixing can provide the first nonlinear accelerator. A
nonlinear multigrid or physics-based coarse-space correction can be introduced
after the modern CPU reference is stable. Global gauge and spin-rotation zero
or near-zero modes need explicit conventions so that they do not appear as
false non-convergence.

## Geometry and boundary scattering

The container interface should provide ray intersection points, outward unit
normals, region membership, and material or boundary identifiers. Specular
reflection then uses the local normal through

`p_out = p_in - 2 * dot_product(p_in, normal) * normal`.

The reflection law should be separate from the geometry query so that diffuse,
partially specular, and spin-active surface models can be introduced without
rewriting the trajectory integrator.

## MPI and GPU scaling stages

The current MPI decomposition over spatial evaluation points should be retained
for the first behavior-preserving port. Within each rank, directions and
energies form GPU batches. Replicating the complete field on every rank is a
reasonable first implementation and closely matches the legacy communication
pattern.

For meshes too large to replicate, two scalable alternatives must be measured:

- block-distributed directional sweeps that exchange Riccati boundary data;
- direction/energy decomposition with distributed reductions of spatial
  self-consistency fields.

Specularly reflected rays can create long paths and cyclic dependencies, so a
distributed sweep should not be selected before its communication and
convergence behavior has been prototyped.

## Verification gates

1. Modern CPU cylinder reproduces the selected legacy baseline.
2. Batched CPU implementation reproduces the scalar modern implementation.
3. GPU and CPU agree in physical norms, symmetries, conserved quantities, and
   converged observables; bitwise identity is not required.
4. Uniform-grid general geometry reproduces the cylinder.
5. Adaptive-grid results converge to the uniform-grid result.
6. Texture and full-tensor residuals satisfy separate reported tolerances.
7. Specular reflection tests cover planar, cylindrical, and curved surfaces,
   including grazing incidence.

## Open decisions

- Primary accelerator vendors and the compilers available on the production
  cluster.
- The first non-cylindrical target is fixed as a two-dimensional cross section
  with translational symmetry along the cylinder axis; the first concrete
  noncircular geometry remains to be selected.
- Authoritative legacy baseline, launch process count, and convergence status.
- Required dipole, magnetic, flow, and surface terms for the first texture
  calculation.
- Maximum practical replicated-field problem size on the target GPUs.
