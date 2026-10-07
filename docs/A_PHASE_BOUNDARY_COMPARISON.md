# Preserved A-phase reference: boundary diagnostic

Both runs use T/Tc=0.30 and F1s=5.4 from the preserved 100-point radial
reference, a 65x65 smooth Cartesian mesh over [-20,20]^2, active radius 16,
quadratic trajectory interpolation and the same Anderson settings. This
is an inner-region diagnostic, not a claim that radius 20 resolves the
entire A-phase soft texture. The driver uses its existing 48x11 angular
quadrature and archived energy poles, unlike the double-core 64x11 route.

Run each separately so they do not compete for the ten MPI ranks:

```sh
cd /Users/mikael/Documents/Codex/3he-vortex-modernization
caffeinate -i python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_a_phase_smooth_radial_reference.nml \
  --build-dir work/double-core-scratch-build \
  --case-name a-phase-smooth-reference --ranks 10 --max-iterations 10
```

Then:

```sh
caffeinate -i python3 tools/run_axisymmetric_core_benchmark.py \
  --input examples/2d_a_phase_smooth_radial_asymptotic.nml \
  --build-dir work/double-core-scratch-build \
  --case-name a-phase-smooth-asymptotic --ranks 10 --max-iterations 10
```

Use `--max-iterations 0` first if only comparing the initial maps. The
runner writes independent timestamped directories and overrides output
paths. Ten updates are a maximum; early convergence is allowed. These
inputs explicitly allow iteration-limit completion without asserting
physical convergence. Existing benchmark defaults still require convergence.

In the first run, exterior mesh values remain the embedded radial reference;
trajectory continuation uses that radial reference with its automatic outer
matching radius. In the second, continuation is fitted at radii 12 and 15
and the exterior halo is reconstructed from that fit. These fitted
coefficients come from the preserved radial profile, not the evolving 2D
solution. This deliberately tests the known-reference continuation before
testing a fully self-consistent 2D fit. It changes the exterior policy, not
just a flag with identical halo fields. All other controls are matched.

Compare residual maps near the axes and active/halo interface and changes
relative to the radial reference. Per-map snapshots are pre-update defects
and produce residual_maps.png and residual_peak_track.png in plots. For
live inspection between completed map writes:

```sh
python3 tools/plot_residual_locations.py runs/YOUR-RUN --output runs/YOUR-RUN/plots
```

The driver now records per-rank MPI shutdown in
metrics.txt.shutdown.rankNNNNNN.txt. Each should contain
MPI_Finalize_return_code=0 and MPI_Finalized=T. See MPI_SHUTDOWN_INVESTIGATION.md.
