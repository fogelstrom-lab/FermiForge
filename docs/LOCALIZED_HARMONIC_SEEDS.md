# Four independent compact harmonic seeds

These are the preferred controlled tests following the September 23 quadrupole
run. The old `historical_qop` initializer remains available for reproducibility,
but its non-winding 1/r core filling extends into the fitting region. It is not
used by these tests. The September 23 PDF describes that older preset; this
note documents the new options added September 24.

The code's harmonic index order is **spin, orbital**, with projections (+,0,-).

| Test name | `initialization_mode` | Only initially occupied mixed harmonic | Axisymmetric outer phase for n=1 |
|---|---|---|---|
| `0plus` | `localized_0plus` | C(0,+1) | 1 |
| `plus0` | `localized_plus0` | C(+1,0) | 1 |
| `0minus` | `localized_0minus` | C(0,-1) | exp(2 i phi) |
| `minus0` | `localized_minus0` | C(-1,0) | exp(2 i phi) |

The usual diagonal unit-winding B background is retained. Only the selected
mixed harmonic receives an addition, initially real and positive:

```
C_selected(r) = seed_amplitude * Delta_B * (1 - r^2/R_seed^2)^2, r < R_seed
              = 0,                                         r >= R_seed
```

Defaults: `seed_core_radius=5.0` (in xi0), `seed_amplitude=1.0` (in bulk-gap
units). The envelope and its first derivative vanish at its edge. All four
tests use the same amplitude; no opposite-sign partner is added automatically.
The other three mixed harmonics and all initial current-related mean fields
are zero. This does not mean the entire B-phase background is zero.

No angular winding is attached to the localized core addition. Outside the
support all four mixed harmonics initially vanish, including at the asymptotic
fit surfaces. They can develop during iteration. The existing state-dependent
tail retains the 1/r amplitude rule for mixed z/in-plane components and their
angle-dependent matching data. It does **not** project the phase onto a single
angular harmonic: exp[i(n-s-k)phi] is the expected rotationally symmetric
far-field law to check, not a new hard constraint on the broken-symmetry 2D
solution. Localizing the seed removes the imposed wrong tail; it does not by
itself prove the eventual solution has the desired asymptotic phase.

## Launch one test

From the repository root, choose `0plus`, `plus0`, `0minus`, or `minus0`:

```bash
seed=0plus
caffeinate -i env OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm \
  python3 tools/run_axisymmetric_core_benchmark.py \
  --input "examples/2d_localized_${seed}_Fs1_0.nml" \
  --initialization-mode "localized_${seed}" \
  --initialization "localized-${seed}" \
  --case-name "single-harmonic-${seed}-T030-Fs0" \
  --build-dir work/localized-seed-build --ranks 10 --formats png
```

Repeat for the other three choices, preferably sequentially on the laptop.
All inputs retain T/Tc=0.30, Fs1=0, 3481 active points in r<=20 xi0, the
65x65 smooth storage grid of half-width22, quadratic ray sampling, fit radii
14/17, endpoint radius70 and up to50 updates with a checkpoint every update.
The dedicated build directory avoids replacing the prior run's executable.

Inspect without stopping a run:

```bash
python3 tools/plot_running_axisymmetric.py \
  --latest single-harmonic-0plus --no-radial-comparison
```

Replace `0plus` with the selected test. The unrelated radial A-phase reference
is still a driver dependency, not a physical convergence target for these
calculations. No expensive four-run batch is launched automatically.

## Restrictions and checks

These modes are exposed through `benchmark_axisymmetric_core_2d` and its Python
runner, not yet the separate double-core driver. They require a centred,
unit-winding vortex and `state_asymptotic` endpoints, with seed radius strictly
inside the active and inner fitting radii. No restart or extra perturbation may
be specified simultaneously. Keep every fit stencil outside the compact seed
when changing meshes or radii; the supplied inputs have ample separation.

Regression tests check the complete harmonic transform at every grid point,
exact zero additions at and beyond the cutoff, unchanged B background and zero
initial mean field for all four modes. Historical seed tests remain unchanged.
