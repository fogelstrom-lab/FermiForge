# Starting from converged reference fields

Run these commands from `/Users/mikael/Documents/Codex/3he-vortex-modernization`.
The wrapper prepares an input and source SHA-256 manifest under a unique
`work/converged-reference-inputs/` directory, then uses the existing timestamped
run/plot/checkpoint wrappers. Source data are read without modification.

## Double core

First evaluate the map at 11 probe points without taking an Anderson step:

```sh
python3 tools/run_converged_reference.py --probe-only
```

To retain the original 601 by 601 grid for this check:

```sh
python3 tools/run_converged_reference.py --probe-only --native-grid
```

The normal full-domain diagnostic evaluates the initial map on the active
ellipse (28 by 32 xi0 radii). It can take substantially longer than the probe
check:

```sh
python3 tools/run_converged_reference.py
```

To start a 20-update relaxation from the imported solution:

```sh
caffeinate -i python3 tools/run_converged_reference.py --iterations 20
```

This uses twice the azimuthal resolution of the preceding 32-direction
overnight run, plus a larger mesh/update region. Allow more time per update;
`--iterations 2` is a useful initial timing and drift check before a long run.

Defaults use the archive at T/Tc=0.30, input F1s=5.4, feedback
5.4/(1+5.4/3)=1.9285714285714286, bulk gap 0.2799462, and the archived
quadrature files with 64 azimuthal directions and 11 polar nodes. The modern
multiscale mesh spans +/-60 xi0, with 0.4 spacing to 16, 0.8 to 36, and
2.0 to 60 xi0 (155 by 155 nodes). The Fortran `legacy_resampled` initializer
uses the tested bilinear sampler to transfer the full complex matrix and all
three mean fields. Requests outside the source grid fail explicitly. Use
`--fine-spacing 0.2` to test central resolution, or `--native-grid` to retain
the original full grid; a full active-domain iteration on that grid is costly.

The first comparison freezes the imported exterior and uses the archived
field at radii 50--58 xi0 to supply the asymptotic endpoint. No zero-mode
projection is applied, so residuals expose the complete map defect. It is a
controlled test of an imported reference with fixed exterior data, not a
fully relaxed exterior texture. Afterward, repeat with:

```sh
caffeinate -i python3 tools/run_converged_reference.py \
  --iterations 20 --halo state_asymptotic
```

This second case regenerates the exterior from evolving field values at
radii 24--27 xi0 inside the active ellipse. Compare its initial halo change,
residual, and core-size evolution with the frozen-exterior case. The fitting
radii are test controls, not a claim that this close-in fit is asymptotically
accurate. Existing probe/shadow checks can then be enabled in the prepared
input. They are disabled in these initial reference tests to isolate the map.

Every iteration writes a full-state sparse checkpoint. Resume with the
existing `run_double_core_from_scratch.py --input PATH/TO/input.nml
--sparse-restart PATH/TO/sparse_checkpoint_state.dat --iterations 20` route.
The legacy initializer recreates the base mesh before overlaying the complete
checkpoint. `--prepare-only` writes a reviewable input without running MPI;
`--no-build` reuses the current executable.

## Normal and A-phase cores

The preserved converged `new_src` radial datasets are already under
`benchmarks/normal_core_2d/reference/current_new_src_T0.30_F1s5.4/` and
`benchmarks/a_phase_core_2d/reference/current_new_src_T0.30_F1s5.4/`.
The A-phase reference is the user's 100-point converged calculation with an
istart=3 reload check. No fresh radial calculation is needed for these tests.

```sh
python3 tools/run_converged_reference.py --core normal
python3 tools/run_converged_reference.py --core a-phase
```

Both rotate/embed the radial state on a Cartesian mesh, retain the radial
reference at trajectory endpoints, and measure the initial map defect. The
default mesh spans +/-32 xi0 with 0.25 central spacing and active radius 9.
It tests the inner region; it does not establish convergence of the entire
A-phase soft core. Use `--iterations 20` to relax the interior from the same
reference and `--fine-spacing` for a resolution check. The prepared namelist
allows larger active radii and different grids for a subsequent study.

The radial benchmark currently uses its existing fixed 48-by-11 angular
quadrature; this differs from the double-core archive's 64-by-11 quadrature.
Its radial matching radius is selected from the outer reference grid, rather
than inherited from the short-radius asymptotic smoke template.

## Interpretation and first checks

An old converged solution need not be a fixed point of a newly discretized
map. Compare unprojected residuals before mixing, final-minus-initial fields,
core separation, and resolution dependence. A successful executable exit
and a roundoff-level propagator normalization error do not establish physical
fixed-point agreement.

The first ten-rank double-core probe checks completed with all 704 directions
and eight energy poles. On the original grid the maximum residual was
3.8069e-4 and the relative L2 residual 6.4715e-4, with the largest point RMS
at a half core. The normalization error was 1.1145e-15. On the resampled mesh
the corresponding maxima were 8.8104e-5 and 1.2925e-4 (normalization
1.2212e-15). These are not a pointwise grid-convergence pair: the coarse mesh
samples the half-core probes at +/-11.6 while the native grid samples +/-11.8.
The smaller coarse-grid residual cannot be interpreted as greater accuracy.
The full transferred state agrees with an independent Python interpolation
of all matrix and mean-field components to machine precision.

A subsequent `--probe-only --fine-spacing 0.2` check uses a 235 by 235
multiscale mesh and recovers a=23.6 xi0. At the same half-core coordinates
as the native grid, the maximum residual is 3.8096e-4 and core relative L2
residual 1.4238e-3, essentially unchanged from the native-grid core result.
This points beyond transfer error alone; it is not yet a diagnosis of the
remaining discrepancy. The outer relative residual is more sensitive to the
coarser exterior sampling (3.1238e-5 versus 3.2996e-6 on the native grid).

The default normal-core map also completed on ten ranks (2037 active points,
295 seconds). Input embedding errors were zero; its map relative L2 residual
was 1.0544e-2 and maximum residual 6.2595e-3. A compact A-phase smoke test
with only 109 active points gave relative L2 residual 7.979e-3 and maximum
residual 3.9477e-3, again with exact embedding. These finite residuals need
resolution and map-convention checks before claiming agreement with the
radial fixed point. The A-phase smoke test does not validate the default
larger mesh. The two-update double-core integration check updated only three
probe points; a full-domain reference relaxation is still to be run.
