# Sparse outer extension test

Run from the project root:

```bash
caffeinate -i env OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm \
  python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_a_phase_extended_outer.nml \
  --build-dir work/quadrupole-seed-build \
  --case-name a-phase-extended-outer --initialization core-transfer \
  --ranks 10 --formats png
```

The example performs up to five updates, checkpointing every update. No
--restart option is necessary: the input selects restart_core and names the
50-update overnight final field. The archived, corrected radial data initialize
the newly active outer points and remain available for independent comparison.
They do not fix those points during iteration.

## Grid and transfer

`mesh_kind='extended_smooth'` preserves the entire old 65-node coordinate
axis on [-20,20]. Twelve new coordinates on each side extend it to [-60,60].
Outer spacings grow with the configurable ratio 1.15, with a small common
rescaling to hit the outer endpoint without a short last cell. The grid is
89x89 (7921 points), with 7453 active inside radius 58. Minimum spacing is
unchanged at 0.29142; maximum outer spacing is 6.41672.

This is a tensor-product extension with a circular active region, not a polar
annular discretization. Cross-products of old and new coordinates create more
points than a true annular grid. A polar or multi-patch method needs separate
interpolation/trajectory work. The common mesh constructor can be reused for
double-core development, but that driver's mesh input is not extended here.

`restart_core` copies every old active point at r<=16 exactly, verifies every
required point exists and rejects duplicates/nonmatching grids. New active
points start from the corrected radial reference without smoothing across the
transfer circle. All active points then relax with the same equations.

Matching moves from r=14 to r=52; the Fermi-liquid two-power fit uses r=40
and 52. The bilinear fit stencil at r=52 reaches at most r=57.1942, inside
the active radius. Trajectory quadratic interpolation and step=0.25 are
unchanged. The outer trajectory termination radius increases to 100, enclosing
the enlarged square. Very remote reference-map corners beyond the radial
data's r=70 support are filled by a radial-reference tail for display only;
the full active comparison region is inside the supplied data support.

## What this tests

Compare core and outer profiles, both separately and against the corrected
radial solution. Global errors now cover a larger area, so they are not
directly comparable with the old radius-16 global norm. Core transfer is exact,
but the initialized outer solution and changed boundary closure can give a
transient residual. Coarse outer interpolation and matching-radius adequacy
remain numerical questions to check, not assumed successes.

Cost: active points grow by 2.43 times and trajectories become longer. Use the
first measured map time to budget subsequent runs; this is not a negligible
cost increase. Preparation and regression tests were run, not a full-domain
trajectory map or production relaxation.
