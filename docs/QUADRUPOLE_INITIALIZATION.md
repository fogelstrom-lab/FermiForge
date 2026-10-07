# Fourth core seed: quadrupole

**September 24 update:** the extended non-winding tail in this original seed
is not suitable for the new controlled initialization tests. Prefer the four
[localized single-harmonic options](LOCALIZED_HARMONIC_SEEDS.md). The original
mode below is retained to reproduce earlier runs, not silently redefined.

Select `initialization_mode = 'historical_qop'` (shared seed kind `qop`).
Both 2D drivers and Python runners accept it. Despite the compatibility name
`historical`, this is a new initializer derived from the user's specification
and Fogelstrom and Kurkijarvi, JLTP 116 (1999), Fig. 1 and pp. 3-5, not a
transcription from the old nop/aop/dop Fortran routines.

In the code's (+,0,-) harmonic order, the added core field is

    C_0- = +sqrt(2) Delta f(r),  C_-0 = -sqrt(2) Delta f(r)
    C_0+ = C_+0 = 0

where f(r)=tanh(s)/s, s=r/(3 xi_seed), xi_seed=1/sqrt(1-t^2), f(0)=1.
This uses the existing A-phase seed envelope and its core amplitude scale.
Both occupied core components have NO azimuthal phase winding. Applying
the axisymmetric harmonic phase rule here would instead give exp(2 i phi)
and would describe the wrong core. Full harmonic-to-Cartesian conversion gives
A_zx=+Delta f, A_zy=-i Delta f, A_xz=-Delta f, A_yz=+i Delta f.
The diagonal phase-wound B background is unchanged. All initial Fermi-liquid
mean fields are zero. The seed is finite and direction-independent at r=0.

These are initialization choices only: no harmonic is held zero during
iteration, no axial symmetry is imposed, and the common state_asymptotic
boundary scheme remains unchanged. Relaxation must establish whether the
seed reaches a metastable quadrupole state; this is not yet demonstrated.
The A-phase radial reference can remain as a comparison, but is not a
quadrupole convergence target. Existing double-core separation diagnostics
are not a quadrupole shape diagnostic.

For later exploration, from the project root (do not run alongside a busy
overnight calculation unless resources permit):

```bash
caffeinate -i env OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm \
  python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_quadrupole_fresh.nml \
  --initialization historical-qop \
  --case-name quadrupole-fresh \
  --build-dir work/quadrupole-seed-build --ranks 10 --formats png
```

This retains the current exploration parameters T/Tc=0.30, F1s=5.4 and up to
50 iterations; it does not reproduce Fig. 1's T/Tc=0.60, F1s=6.0 parameters.

## Compact F1s=0 exploration (2026-09-23)

Use `examples/2d_quadrupole_Fs1_0.nml` for T/Tc=0.30, F1s=0,
up to 50 updates, with a checkpoint every update:

```bash
caffeinate -i env OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm \
  python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_quadrupole_Fs1_0.nml \
  --initialization historical-qop --case-name quadrupole-T030-Fs0 \
  --build-dir work/quadrupole-seed-build --ranks 10 --formats png
```

The smooth rectilinear grid has 65 x 65 nodes and half-width 22 xi0;
3481 nodes inside r=20 xi0 are active. The fit radii are 14 and 17 xi0.
Trajectories continue through the evolving asymptotic field to radius 70 xi0,
with maximum path step 0.25 xi0. The small computational cell is an exploratory
choice, not a demonstrated domain-convergence result for this core.

`feedback_scale=F1s/(1+F1s/3)=0` turns off Fermi-liquid feedback.
Temperature is read from the Ozaki file (0.30), not inferred from filenames.
The driver still loads an A-phase radial reference with F1s=5.4 for its
benchmark comparison output. That reference does NOT initialize this seed
or set its evolving asymptotics. Its radial-difference numbers and benchmark
pass flag must not be interpreted as quadrupole accuracy/physical stability
tests. Inspect residuals and the evolving harmonics instead.

The existing running-checkpoint plotter can inspect this run too; radial
overlays are only a comparison to a different core at a different F1s.
Avoid running two ten-rank calculations simultaneously on the laptop.

Preparation-only validation succeeded with one MPI rank (no nonlinear maps).
This confirms input/mesh/seed construction, not quadrupole convergence.

### Inspect while it runs

After the first checkpoint is saved, use a second terminal:

```bash
python3 tools/plot_running_axisymmetric.py --latest quadrupole \
  --no-radial-comparison
```

This selects the newest dated quadrupole run and creates a separate dated
inspection folder containing field/harmonic maps, axis profiles, convergence,
and spatial residual dynamics. Repeat manually after later updates. No solver
is launched or stopped and previous inspections remain available. The command
reports clearly if no quadrupole run or checkpoint exists yet. This mode omits
the unrelated A-phase radial comparisons. An explicit run directory can be
used instead of `--latest quadrupole` to inspect an older calculation.
