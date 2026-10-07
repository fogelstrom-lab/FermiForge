# Smaller smooth-grid test

Run from the project root:

```sh
cd /Users/mikael/Documents/Codex/3he-vortex-modernization
caffeinate -i python3 tools/run_double_core_from_scratch.py \
  --input examples/2d_double_core_Fs1_zero_smooth.nml \
  --run-name double-core-Fs1-zero-smooth --ranks 10 --iterations 10
```

The production run is left for the user. It starts again from the F1s=0
legacy data, with quadratic trajectory interpolation, 704 directions,
eight archived poles, unchanged trajectory integration settings and
unchanged asymptotic/Anderson controls. Output goes to a new timestamped
directory, with per-update checkpoints and spatial residual snapshots.

The new mesh option is `scratch_mesh_kind='smooth'`. In each Cartesian
direction it uses x(u)=L sinh(s u)/sinh(s), -1<=u<=1, with uniform u.
Here L=20, s=2.3 and 64 cells (65 nodes) per axis. The s=0 limit is uniform.
The existing piecewise multiscale mesh remains unchanged and available.
The old fine/medium spacing parameters in the copied input are ignored
when smooth mode is selected.

The full grid has 4225 points and the radius-16 active disk has 3065,
compared with 12321 and 9309 respectively in the prior quadratic run.
The central spacing is 0.2914 xi0, the largest outer spacing 1.4163 xi0,
and adjacent outward cell widths grow by no more than about 7.3%.
Point-count scaling suggests roughly ten rather than thirty minutes per
map, but actual timing must be measured; interpolation, paths, and load
balance also affect it.

This economical test changes both grading and resolution. It is not a
matched-resolution convergence study. Keep an eye on the residual between
the half cores and near the fit surfaces at radii 12/15. Core separation
is currently measured at grid nodes: the nearest node to y=3.2 is 3.16896,
so an initial separation near 6.338 instead of 6.4 is a sampling effect,
not evidence of physical shrinkage. Compare field profiles as well.

To refine the smooth mesh in a later test, increase scratch_number_of_cells
to an even number (e.g. 80 or 96), leaving the stretch and physics fixed.
This is still a Cartesian tensor-product grid, not yet elliptical rings.
The smooth grid does not remove the separate asymptotic matching surface.

For live residual plotting between completed map writes:

```sh
python3 tools/plot_residual_locations.py runs/YOUR-RUN --output runs/YOUR-RUN/plots
```

Snapshots are pre-update. A zero-update restart from the final checkpoint
is needed to measure the final state's full residual, as in
RUNNING_FS1_ZERO_REFERENCE.md.
