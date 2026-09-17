# FermiForge Project Plan

Scalable quasiclassical simulation of spinful superfluids and superconductors.

Status: living project plan, version 1, 2026-09-15

Execution is tracked in `MILESTONE_BACKLOG.md`; mathematical definitions are
kept in `architecture/MATHEMATICAL_CONVENTIONS.md`. This document remains the
scope and gate authority.

## 1. Project focus

### North-star outcome

Develop a validated, open, modern quasiclassical solver that can compute
self-consistent spinful superfluid and superconducting states in general
two-dimensional device cross sections, run efficiently on MPI clusters, and
use GPUs without making the physics implementation depend on one accelerator
family.

The framework should ultimately cover:

- the full complex spin-triplet order parameter of superfluid 3He;
- singlet, triplet, and mixed spin structures in superconductors;
- Fermi-liquid, exchange, impurity, electromagnetic, and spin-dependent
  self-energies as selected by the physical model;
- vortices, material regions, surfaces, weak links, leads, and spin-active
  interfaces;
- equilibrium self-consistency first, with spectroscopy and nonequilibrium
  extensions added only after the equilibrium foundation is reliable.

### Immediate scientific target

The first symmetry-relaxed target is the free double-core vortex in 3He-B,
with translation symmetry along the vortex axis and no container wall. This is
the project's primary acceptance problem because it simultaneously requires:

- a fully two-dimensional solution that can break axial symmetry;
- the full complex 3 by 3 order parameter and the Fermi-liquid feedback;
- accurate coherence-length resolution around two half cores;
- a large domain or controlled asymptotic continuation for the long 1/r
  spin-rotation tail;
- robust nonlinear iteration in the presence of gauge, orientation, and
  metastability issues.

The existing axially symmetric A-phase-core vortex remains a control problem.
The cylindrical wall implementation remains a valuable regression test, but
general walls and leads are not on the critical path until the free double-core
calculation has passed its validation gates.

### Scope guardrail

For now, do not combine all of the following in one implementation step:
two-dimensional fields, adaptive meshes, discontinuous Galerkin (DG), full
spin structure, general interfaces, MPI domain decomposition, and GPU
offload. Each feature must first agree with a simpler reference at a named
validation gate.

## 2. Technical basis and implications

### Existing 3He codes

The fixed-form `incoming/legacy-f77` code is the behavioral reference. The
free-form `new_src` code is the starting point for extracting modern physics
kernels, but it is not yet a two-dimensional solver and its adaptive radial
mesh foundation is not yet connected to propagation.

The useful inheritance from these codes is:

- mature 3He order-parameter and Fermi-liquid physics;
- angular and Ozaki quadratures;
- Riccati or equivalent causal transport along a trajectory;
- MPI distribution over independent self-consistency evaluations;
- accumulated experience with Anderson-like iteration accelerators.

The parts that must not define the new architecture are compile-time array
limits, mutable global state, implicit precision, radial-only indexing, and a
geometry model embedded directly in the trajectory routine.

### Double-core reference calculation

The attached `dcvlong.pdf` supplies the scientific benchmark and several
important numerical constraints. In particular, it describes self-consistent
Eilenberger calculations on a square lattice inside a radius `Rc`, followed by
a symmetry-allowed outer expansion containing the 1/r spin-rotation and 1/r^2
amplitude tails. The reported calculation uses `Rc = 38 xi0`, spacing
`d = 0.2 xi0`, an 11-point polar quadrature, and an azimuthal sum. It identifies
the following comparison observables:

- the two pair-density minima and mass-current stagnation points;
- the half-core separation `a`;
- the asymptotic coefficients `C1` and `C2`;
- the 2m'm' symmetry of the converged vortex;
- the strong dependence on temperature and `F1s`;
- the A-phase-core vortex as an axially symmetric comparison solution.

These are more useful acceptance tests than component-by-component equality
with one historical output file, because the paper's transport method is not
identical to the present Riccati implementation.

### SuperConga

SuperConga demonstrates several architectural ideas worth adopting:

- configuration-driven preprocessing, a compute backend, and separate
  postprocessing;
- distinct geometry, boundary-condition, order-parameter, Riccati,
  quadrature, observable, and nonlinear-accelerator components;
- batching over grid points, momentum directions, and energies;
- extensive numerical and physics tests;
- a thin portability boundary between CUDA and experimental HIP backends.

It is not a drop-in base for this project. The inspected development version
is designed around spin-degenerate/singlet fields, a regular dense grid, and
CUDA/HIP device containers; it has no MPI execution layer. We should reuse
ideas and validation cases before considering source-level sharing. Any reused
LGPL source must retain its licensing and attribution.

### Discontinuous Galerkin method

The Seja-Lofwander method is a discontinuous Galerkin finite-element method,
not a finite-difference method. Its key attraction is that the coherence
amplitudes and self-energies live on the same mesh, removing the repeated
interpolation between a field grid and separate trajectory stepping points.
The upwind numerical flux preserves the causal transport direction, and each
energy and Fermi-momentum direction remains an independent transport problem.

The published implementation treats a scalar, spin-degenerate Riccati
amplitude. The paper states that full spin gives four coupled unknowns, but a
robust matrix-valued implementation, its nonlinear solver, its memory cost,
and its accelerator behavior still have to be demonstrated for this project.
DG is therefore the preferred long-term hypothesis, subject to an early
prototype and a formal go/no-go gate.

References:

- K. M. Seja and T. Lofwander, "A finite element method for the
  quasiclassical theory of superconductivity," Phys. Rev. B 106, 144511
  (2022), <https://doi.org/10.1103/PhysRevB.106.144511>.
- SuperConga source, <https://gitlab.com/superconga/superconga>, inspected at
  development commit `81299835461788335572db7e26972070cf4e9e15`.
- SuperConga documentation, <https://superconga.gitlab.io/superconga-doc/>.

## 3. Architectural rules from the beginning

### Physics is independent of execution

The same mathematical kernels must be callable from a scalar CPU reference,
an MPI batch executor, a threaded CPU executor, and a GPU executor. MPI calls,
GPU directives, and vendor-specific data types must not appear in the physical
model interfaces.

### Full spin is a data-model requirement, not a later patch

The core representation should support matrix-valued quasiclassical
propagators and self-energies from the first 2D implementation. A spin-scalar
problem may select a smaller optimized path, but scalar assumptions must not
be built into the interfaces.

For 3He-B the principal field is naturally represented as
`A(spin, orbital, dof)`, with three spin and three orbital components. The
transport layer should consume the resulting 2 by 2 spin-space gap matrix,
not know how `A` was generated. This separation is what later permits singlet
and triplet superconducting models to use the same transport code.

### Mesh, geometry, and boundary scattering are separate concepts

The mesh owns degrees of freedom, element connectivity, refinement, and
field transfer. Geometry owns material-region and boundary labels. A boundary
model maps incoming to outgoing coherence data. The transport solver uses all
three without embedding a circle, a specular law, or a material model in its
integration kernel.

### Adapt only between nonlinear stages initially

The mesh should be held fixed during one self-consistency solve. Refine or
coarsen after a converged or deliberately stopped stage, transfer the fields,
reset the Anderson history, and reconverge. Changing the vector space inside
an Anderson history makes residual comparisons and subspace vectors
ill-defined.

### Configuration and output are versioned

Every run records the code revision, input, mesh, quadratures, physical model,
solver tolerances, rank/thread/device layout, residual history, and restart
provenance. Output must contain enough metadata to reject an incompatible
restart rather than silently reinterpret it.

## 4. Target component model

The intended dependency direction is:

```text
run configuration
       |
       v
physical model ---> field/state representation <--- mesh + geometry
       |                      |
       v                      v
self-energy map <--- transport operator ---> boundary/interface model
       |
       v
nonlinear solver ---> diagnostics/observables ---> restart + postprocessing
       ^
       |
execution backend: scalar CPU | MPI/CPU | MPI/GPU
```

Recommended module responsibilities:

1. `models`: 3He weak-coupling/Fermi-liquid model first; later singlet,
   triplet, magnetic, spin-orbit, and impurity terms.
2. `fields`: explicit-kind, contiguous storage and pack/unpack views for all
   self-consistent fields.
3. `mesh`: uniform Cartesian reference, then one selected adaptive element
   hierarchy.
4. `geometry`: region membership, boundary facets, normals, and stable IDs.
5. `transport`: scalar reference and matrix Riccati kernels; later DG
   assembly/solve behind the same high-level interface.
6. `boundaries`: bulk/asymptotic, periodic, specular, transmissive, reservoir,
   and eventually spin-active scattering laws.
7. `quadrature`: Fermi-surface and energy rules with convergence controls.
8. `nonlinear`: simple mixing, the legacy Anderson implementation, and later
   multilevel/preconditioned variants.
9. `execution`: work batching, MPI collectives, CPU threading, and GPU data
   residence.
10. `io`: versioned configuration, checkpoints, observables, and Python
    analysis.

## 5. Roadmap and gates

The effort estimates below are order-of-magnitude estimates for one
full-time researcher. They are sequencing aids, not deadlines.

### M0 - Lock references and finish the radial bridge (2-4 weeks)

Goal: make the existing calculation a dependable oracle before expanding its
dimensionality.

Work:

- select canonical free A-core and cylindrical-wall input/output baselines;
- finish the modern free/cylinder switch as a regression facility, while
  keeping general wall work out of scope;
- document every order-parameter, self-energy, normalization, coordinate,
  phase-winding, and momentum convention;
- convert the active interfaces to explicit kinds and shapes far enough to
  remove compiler-dependent ABI assumptions;
- report average and maximum residuals, quadrature settings, trajectory step,
  and nonlinear-accelerator state in machine-readable logs;
- add deterministic unit tests for interpolation, reflection, Riccati
  propagation, packing, and the Anderson history.

Exit gate G0:

- free and cylindrical radial cases run from clean directories;
- results are rank-count reproducible within a declared tolerance;
- the modernized radial calculation agrees with the selected baseline in
  fields and integrated observables;
- all deliberate numerical differences are listed rather than hidden.

### M1 - Uniform-grid 2D spinful CPU reference (6-10 weeks)

Goal: remove cylindrical symmetry without adding adaptivity or device walls.

Work:

- introduce a uniform Cartesian 2D mesh and full-field storage;
- implement a matrix-valued spinful transport kernel with a scalar reference
  mode;
- evaluate fields along straight trajectories in an effectively unbounded
  system using a controlled outer/asymptotic condition;
- implement the 3He gap and Fermi-liquid self-consistency map on all 2D
  points;
- define gauge fixing and an orientation convention for comparison without
  artificially imposing axial symmetry;
- retain simple mixing as the debugging reference and add Anderson only after
  the map is independently tested.

Exit gate G1:

- an axially symmetric field sampled on the 2D grid reproduces the radial
  solver after spatial, angular, energy, and trajectory-step convergence;
- normalization, particle-hole relations, axial symmetry, and current
  diagnostics pass pointwise or norm-based tests;
- a complete self-consistency iteration is deterministic and restartable.

### M2 - Free double-core vortex and first MPI scaling (8-16 weeks)

Goal: converge the primary broken-symmetry physics benchmark.

Work:

- implement at least two independent seeds: a London/two-half-core seed and
  a symmetry-breaking perturbation of an axial vortex;
- use continuation in temperature, `F1s`, domain size, and resolution rather
  than attempting the hardest low-temperature point directly;
- distribute independent point/direction/energy batches over MPI while
  initially replicating the mean fields on each rank;
- diagnose and pin only the global gauge and continuous orientation zero
  modes needed for stable nonlinear iteration;
- compute pair density, current/stagnation points, half-core separation, and
  asymptotic `C1`, `C2` coefficients automatically.

Exit gate G2:

- a stable two-half-core solution is obtained without imposing axial
  symmetry;
- the solution satisfies 2m'm' symmetry within tolerance and is distinct in
  free energy or stability from the A-core control;
- `a`, `C1`, and `C2` are converged with respect to domain, mesh, trajectory
  step, energy cutoff, and angular quadrature, and are consistent with the
  trends in `dcvlong.pdf`;
- the result is independent of reasonable seed variations and MPI rank count;
- a scaling report identifies compute, interpolation, reduction, and memory
  bottlenecks.

### Parallel method track D - Decide the spatial transport method

This bounded prototype starts once the M1 equations and field conventions are
stable. It must not become a second production solver indefinitely.

D0. Reproduce the one-dimensional scalar benchmark from Seja-Lofwander with
piecewise constant, linear, and quadratic DG bases, including the expected
error-order study.

D1. Solve a clean two-dimensional scalar test on the same mesh with the
trajectory/Riccati and DG methods. Compare fields, observables, cost, memory,
and convergence.

D2. Extend the DG residual and Jacobian action to the 2 by 2 spin coherence
matrix and compare it against the spinful trajectory solution on a uniform
mesh.

Decision gate GD:

- select DG as the production transport method only if it matches the
  reference physics, remains robust for matrix Riccati amplitudes, supports
  the required boundary maps, and has a credible matrix-free or scalable
  sparse-solve path;
- otherwise retain adaptive trajectory integration as the production method
  and keep DG as a research branch;
- do not build full adaptivity and GPU optimization for both approaches.

### M3 - Production adaptive 2D solver (8-16 weeks after GD)

Goal: resolve the half cores and long tail at substantially lower cost than a
globally fine grid.

Work:

- implement h-refinement on the mesh selected at GD;
- refine using distance to half cores, field gradients/curvature, and the
  self-consistency residual; retain enough coarse coverage to represent the
  long spin-rotation tail;
- implement conservative field transfer and symmetry-aware error norms;
- freeze the mesh during nonlinear stages and reset accelerator history after
  transfer;
- add multilevel correction or preconditioning for the long-wavelength
  spin-rotation sector if Anderson mixing alone stalls;
- compare an explicit asymptotic outer region with simply enlarging the
  numerical domain.

Exit gate G3:

- adaptive and uniform sequences converge to the same `a`, `C1`, `C2`, free
  energy, and current observables;
- the low-temperature/high-`F1s` double core converges reliably from a
  documented continuation sequence;
- refinement reduces time or memory by a measured amount without shifting the
  physical solution beyond the error budget;
- no mesh-dependent drift of gauge or vortex orientation remains.

### M4 - Portable accelerator execution (8-16 weeks)

Goal: accelerate the validated production method without forking the physics.

Work:

- first add threaded CPU batching and profile the actual hot kernels;
- keep fields resident for a complete transport/self-consistency batch;
- introduce standard OpenMP target offload as the first Fortran portability
  path, with compiler-specific build presets outside the physics modules;
- distribute energy/direction batches across MPI ranks and use one or more
  devices per rank only through the execution layer;
- compare the practical portability tradeoffs with SuperConga's CUDA/HIP
  backend separation;
- add continuous tests for at least one CPU compiler and each accelerator
  family actually used by the collaboration.

Exit gate G4:

- CPU and GPU converge to the same physical solution and diagnostics within a
  stated floating-point tolerance;
- GPU execution gives a measured benefit on a production-size double-core
  case, including data-transfer costs;
- multi-node runs have no device-specific code in the physics interface;
- scaling limits and memory per field, angle, energy, and mesh degree of
  freedom are documented.

OpenMP source portability does not guarantee equal compiler maturity or
performance on all GPUs. If a vendor-specific optimized kernel is eventually
needed, it must remain optional behind the same tested interface.

### M5 - Materials, walls, interfaces, and leads (12+ weeks per major layer)

Goal: turn the vortex solver into a general spinful device framework.

Implement in this order:

1. material-region labels with the same material everywhere as a null test;
2. a circular specular wall that reproduces the radial cylinder calculation;
3. general closed 2D boundaries with specular reflection;
4. transparent and partially transmitting nonmagnetic interfaces;
5. bulk/reservoir leads with phase or current constraints;
6. singlet/triplet material models and superconducting-magnetic bilayers;
7. spin-active scattering matrices and spin-orbit interface terms;
8. diffuse or mixed-specularity surfaces;
9. nonequilibrium distribution functions only after equilibrium current
   conservation and interface tests pass.

Exit gate G5 for the first device release:

- the circular wall agrees with the radial reference;
- current is conserved through a weak link to the discretization tolerance;
- normal-state, bulk-superfluid, transparent-interface, opaque-interface, and
  spin-rotation limiting cases are recovered;
- at least one weak-link problem and one magnetic-bilayer problem have mesh,
  quadrature, and solver convergence studies;
- checkpoints and visualization identify material regions and boundary laws.

## 6. Validation matrix

Every physics result must pass the applicable rows before being treated as a
benchmark.

| Test class | Required evidence |
|---|---|
| Algebra | Green-function normalization, conjugation, particle-hole, and spin-rotation identities |
| Radial regression | Free A-core and specular cylinder against selected legacy outputs |
| 2D embedding | Axisymmetric 2D solution against the radial solver |
| Double core | 2m'm' symmetry, two half cores, `a`, `C1`, `C2`, currents, and free-energy/stability comparison |
| Discretization | Mesh, polynomial order or trajectory step, domain size, angular quadrature, and energy cutoff |
| Nonlinear solve | Both L2/average and maximum residuals, restart consistency, seed/continuation sensitivity |
| Parallelism | Rank/thread/device-count reproducibility and reduction-order sensitivity |
| Boundaries | Reflection/transmission conservation, limiting transparencies, grazing incidence, and corner policy |
| Devices | Local current conservation and compatible reservoir/interface conditions |
| Performance | Time and memory split by transport, interpolation/assembly, reduction, nonlinear update, and I/O |

Bitwise equality is not a general requirement. Each benchmark must define
absolute and relative tolerances for fields and, separately, for physically
important observables.

## 7. Project organization and decision discipline

Use the repository as follows:

- `incoming/legacy-f77/`: immutable historical reference;
- `new_src/`: transitional radial application used for regression;
- `src/`: reusable modern modules, with no application-global assumptions;
- `tests/`: unit and small deterministic integration tests;
- `benchmarks/`: versioned scientific cases, expected observables, and
  convergence studies;
- `docs/architecture/`: architecture-decision records and mathematical
  conventions;
- `runs/`: generated runs, never canonical source;
- `MODERNIZATION_LOG.md`: chronological record of work and findings.

Every substantial change should be labeled as exactly one primary category:

1. behavior-preserving refactor;
2. numerical-method change;
3. physical-model change;
4. performance-only change.

Changes from categories 2-4 require a before/after benchmark. Avoid mixing
them in one review because a faster result is not useful if the source of a
physical shift cannot be identified.

At the end of each milestone, record:

- which exit criteria passed or failed;
- the exact benchmark revisions and run metadata;
- known numerical or physical limitations;
- the next milestone selected at the gate;
- features explicitly deferred.

## 8. Main risks and mitigations

### Nonlinear convergence and metastability

The double core may drift in orientation, collapse to another vortex, or
converge only along particular temperature and `F1s` paths. Use explicit gauge
and orientation diagnostics, continuation, independent seeds, controlled
noise, and free-energy/stability checks. Treat the nonlinear accelerator as a
replaceable module rather than part of the physics map.

### Long- and short-length scales in one field

Core and surface structure vary on coherence lengths while spin rotation has
a long tail. Use adaptive resolution plus either an asymptotic outer solution
or a multilevel coarse correction. Do not infer convergence from a small local
residual alone.

### Premature commitment to DG

DG is promising, but the scalar published example does not settle robustness
or performance for matrix spin structure. Enforce gate GD before replacing
the trajectory mainline or building a large DG-specific infrastructure.

### MPI scaling and memory

Replicated mean fields simplify the first MPI solver but eventually limit
problem size. First distribute independent angle/energy/point batches; measure
memory; introduce mesh domain decomposition only when the measured limit
requires it. A DG implementation should prefer matrix-free Jacobian actions
or carefully evaluated sparse solvers rather than a dense global Jacobian.

### GPU portability

Portable directives and portable performance are different goals. Keep a
canonical CPU path, isolate the execution layer, test multiple compilers, and
allow optional optimized kernels without duplicating the physical equations.

### Scope expansion into devices

Walls, weak links, bilayers, spin-active interfaces, electromagnetic response,
and nonequilibrium transport can each become a project. Admit them one at a
time only after their simpler limiting cases pass.

## 9. Immediate work package

The next focused work package is complete when G1 is in sight, not when every
future abstraction exists.

1. Finish G0's free/cylinder radial regression and select canonical baseline
   runs.
2. Write the mathematical conventions note for spin, Nambu, orbital,
   trajectory direction, tilde operation, phase winding, and normalization.
3. Define the 2D field/state API with full spin capacity and explicit kinds.
4. Add a uniform Cartesian mesh and tested interpolation independent of the
   radial globals.
5. Extract a side-effect-free matrix Riccati propagation kernel and compare it
   with `new_src` on identical sampled trajectories.
6. Implement one 2D self-consistency map and reproduce an embedded axial
   vortex before enabling a symmetry-breaking seed.
7. Add MPI batching only after the scalar CPU iteration is numerically stable.
8. Start DG track D0 as a bounded prototype once the field conventions are
   frozen.

### Explicitly deferred from this work package

- general walls and leads;
- magnetic bilayers;
- diffuse or spin-active boundary scattering;
- production DG integration;
- GPU offload;
- three-dimensional fields or loss of translation symmetry along the vortex;
- microwave response and bound-state spectroscopy beyond diagnostics needed
  to validate the equilibrium vortex.

This sequencing keeps the project aimed at one decisive result: a reproducible,
self-consistent free double-core vortex from a modern 2D spinful solver. Once
that result and its convergence evidence exist, the choices of adaptive method,
DG promotion, accelerator backend, and device boundary framework can be made
from measured evidence rather than architectural guesswork.
