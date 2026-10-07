# Compact F1s=0 legacy-start baseline

From the project directory:

```sh
cd /Users/mikael/Documents/Codex/3he-vortex-modernization
caffeinate -i python3 tools/run_double_core_from_scratch.py \
  --input examples/2d_double_core_Fs1_zero_reference.nml \
  --run-name double-core-Fs1-zero-reference --ranks 10 --iterations 20
```

The runner builds the executable, makes a timestamped run directory, imports
the original archived state without modifying it, and checkpoints every
update. Twenty is the maximum update count; it may stop earlier on its
convergence criterion. No run has been pre-launched for you.

The mesh spans [-20,20]^2 xi0 with spacings 0.2 inside +/-6, 0.4 through
+/-12, and 0.8 outside. The iterated disk has radius 16; the remaining
points form an asymptotically reconstructed halo, not additional independently
iterated points. The fitting radii are 12 and 15, within the iterated disk.
These choices require validation, particularly since the fit intersects
coarse cells. Interpolation is unchanged (bilinear). Angular quadrature is
64 by 11, with the archived energy poles. The imported half cores are at
y=+/-3.2. feedback_scale=0 eliminates the transport self-energy coupling;
the stored current-related vector is retained, not incorrectly set to zero.

Every evaluated state produces a spatial snapshot and peak-location record:
`sampled_residuals.dat.mapNNNNNN.dat` and `sampled_residuals.dat.peaks.dat`.
They record the raw residual before the corresponding Anderson update.
For live inspection between completed map writes:

```sh
python3 tools/plot_residual_locations.py runs/YOUR-RUN --output runs/YOUR-RUN/plots
```

At completion all standard plots plus residual maps/peak tracks are generated.
Passive exterior ray probes also assess the continuation, but do not establish
independent convergence of the exterior. There is no shadow relaxation here.

To measure the true residual after the final update, run a separate map:

```sh
python3 tools/run_double_core_from_scratch.py \
  --input runs/YOUR-RUN/input.nml \
  --sparse-restart runs/YOUR-RUN/sparse_final_state.dat \
  --iterations 0 --ranks 10 --run-name Fs1-zero-final-check
```

To resume updates instead, use the same command with `--iterations 20`.
This restores fields but starts a fresh Anderson history. Later comparisons
should vary one of interpolation, mesh grading, fit radii, or halo policy
at a time; this baseline alone cannot exclude discretization error.
