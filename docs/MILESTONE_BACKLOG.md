# Milestone backlog

This is the executable backlog for `PROJECT_PLAN.md`. The plan defines scope;
this file defines the next reviewable pieces of work. Status values are
`done`, `ready`, `active`, `blocked`, and `deferred`.

## Current critical path

| ID | Status | Deliverable | Depends on | Acceptance evidence |
|---|---|---|---|---|
| G0-01 | done | Isolated legacy and `new_src` run/archive wrapper | - | Builds outside source trees; timestamped metadata, logs, inputs, hashes, plots |
| G0-02 | done | Free/cylinder smoke and rank-count check | G0-01 | 1- and 10-rank files agree exactly; boundary switch produces a nonzero difference |
| G0-03 | ready | Converged legacy free A-core reference | G0-01 | Both residuals meet tolerance; resolution and restart recorded |
| G0-04 | ready | Converged legacy specular-cylinder A-core reference | G0-03 | Same controls as free case; core and surface convergence recorded |
| G0-05 | ready | Boundary-strategy interface in radial code | C-01 | Explicit free-asymptotic and specular-cylinder implementations; no hidden `cyl` global branch |
| G0-06 | blocked | `new_src` free/cylinder field regression | G0-03, G0-04, G0-05 | Comparison report for all complex fields and current-related fields |
| G0-07 | ready | Solver terminal status and machine-readable residual log | G0-01 | Distinguishes converged, iteration limit, invalid state, and interruption |
| G0-08 | done | Converged current-`new_src` free normal-core reference | G0-01 | 100-point tangent grid; final average/max residuals `5.130e-9`/`1.562e-6`; exact input, history, fields, hashes, and status archived |
| G0-09 | done | Converged current-`new_src` free A-phase-core reference | G0-01 | Consistent 100-point tangent fields archived; zero-iteration restart provenance recorded; independent original/modern point maps agree to `1.20e-4` relative gap error |
| C-01 | active | Freeze mathematical conventions | draft note | Resolve unit, tilde, propagator-normalization, and code-basis questions; approve tests |
| A-01 | done | Pure Cartesian/harmonic conversion module | C-01 | Exact source forward transform, general inverse, axial transport projection, phase, and far-field tests |
| A-02 | active | Explicit-kind 2D field/state API | C-01 | `A(3,3,npoint)` plus diagonal/spin self-energies; accelerator-independent packing tests |
| M1-01 | done | Uniform Cartesian mesh and interpolation | A-02 | Constant/affine/bilinear exactness, quadratic and axial-vortex refinement, boundary-safe lookup |
| M1-02 | done | Side-effect-free legacy Riccati kernel | C-01, M1-01 | Frozen forward/reverse trajectory agrees with `new_src` at every integration checkpoint |
| M1-03 | active | 2 by 2 spin-matrix Riccati kernel | M1-02 | Scalar limit, source normalization, coefficient conversion, conjugate coherence, and spin rotation pass; general tilde/particle-hole convention still open |
| M1-04 | active | Unbounded 2D trajectory sampler | M1-01, M1-03 | Selectable free-vortex continuation, source tail powers, boundary continuity, uniform-bulk invariance, and serial/MPI tests pass; outer-radius/domain/step convergence remains |
| M1-05 | active | 2D 3He self-consistency map | A-02, M1-04 | Full-grid serial/MPI map, spatial residual diagnostics, packing, Anderson adapter, and free-vortex endpoint pass; converged axial 2D regression remains |
| M1-06 | active | Embedded normal-core and A-core 2D regression | G0-08, G0-09, M1-05 | Both independent point oracles and uniform-grid map residuals recorded; field tolerances remain to be met under multiscale mesh/endpoint/step/rank refinement |
| M1-07 | done | Static multiscale rectilinear mesh and circular active region | M1-06 | Fine/medium/coarse coordinates, nonuniform Q1 lookup, active-only MPI mapping, frozen halo, masked/zonal residuals, normal/A-core maps, and strict tests pass |
| M2-01 | deferred | Double-core seeds and gauge/orientation constraints | M1-06 | London and perturbed-axial seeds; only zero modes pinned |
| M2-02 | deferred | Temperature/domain/resolution continuation | M2-01 | Stable two-half-core solution from at least two seed paths |
| M2-03 | deferred | `a`, `C1`, `C2`, pair density, current diagnostics | M2-02 | Automated asymptotic windows and uncertainty from convergence studies |
| M2-04 | active | First MPI batching/scaling report | M1-05 | Point batching is identical for 1/4/10 ranks with preliminary timing; production-size compute/reduction/interpolation/I/O report remains |
| D0-01 | deferred | Scalar 1D DG reproduction | C-01, M1-02 | Published convergence order for basis degree 0, 1, and 2 |
| D1-01 | deferred | 2D scalar trajectory-vs-DG comparison | D0-01, M1-06 | Same mesh/physics, field error, cost, memory, robustness |
| D2-01 | deferred | Spin-matrix DG prototype and method gate | D1-01, M1-03 | Matrix equivalence and documented keep/reject decision |
| M3-01 | active | Adaptive 2D production mesh | M1-07, method gate | Dense block hierarchy, core/tail indicators, conservative transfer, two-to-one balance, hysteresis, convergence study |
| M4-01 | deferred | Portable GPU execution layer | M3-01 | CPU/GPU physics equivalence before performance claims |
| M5-01 | deferred | Walls, interfaces, leads, and devices | G0-05, M3-01 | Boundary conservation and limiting-case suite before weak links/bilayers |

## Next two focused review batches

### Batch A - close the radial oracle

1. G0-07: make termination status reliable before spending hours on reference
   runs.
2. C-01: review and freeze the mathematical convention questions.
3. G0-03: run the free A-core reference with a continuation/restart schedule;
   do not accept `itmax` as convergence.
4. G0-04: repeat the matched contained case and add surface-window diagnostics.
5. G0-05 and G0-06: port the boundary law only after the two oracle outputs are
   accepted.

### Batch B - make the first 2D kernel testable

1. A-02: finish the field-state API; A-01 basis conversion is complete.
2. M1-02 is complete: the source transport arithmetic and forward/reverse
   checkpoints are frozen independently of global state.
3. M1-03: add spin matrices with algebraic invariants.
4. M1-04 and M1-05: evaluate an axial field in 2D before self-consistent
   symmetry breaking.

## Rules for changing priority

- A task enters `active` only with a named acceptance test.
- A numerical-method, physical-model, or performance change needs a before and
  after benchmark; do not combine these categories in one change.
- General device geometry and GPU kernels remain deferred until their incoming
  gate is passed.
- A failed or nonconverged run is retained as evidence but cannot become a
  canonical baseline.
