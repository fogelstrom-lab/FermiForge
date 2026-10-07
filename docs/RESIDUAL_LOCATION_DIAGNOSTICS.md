# Spatial residual diagnostics

The double-core solver now saves `sampled_residuals.dat.map000001.dat`,
etc., after every successful map evaluation, before mixing or zero-mode
projection. Each is the raw defect F(X)-X of the evaluated state, not the
Anderson step. Map 1 evaluates the initial state; map k evaluates the state
after k-1 updates. These snapshots do not entail extra trajectory solves.
The final `sampled_residuals.dat` remains the last evaluated map; after an
update a separate zero-update restart is necessary to measure the final state.

`sampled_residuals.dat.peaks.dat` is flushed by closing it after every map.
It records both maximum point-RMS and maximum absolute real-component
residual, with their coordinates, restricted to the update mask. The RMS
includes 18 real order-parameter components and three real Fermi-liquid
components, divided by 21. No component rescaling is applied. Ties select
the first point; near symmetry can make the selected peak switch between
cores without indicating actual motion of the error distribution.

Snapshot columns preserve the existing ten columns and append `update_mask`.
The plotting tool accepts old ten-column files using their region flag as
the interior mask. It plots all snapshots when there are at most three;
otherwise maps show first and last, while the peak track includes every map.
Exterior passive probes are excluded from these interior maps.

During or after a run, execute:

```sh
python3 tools/plot_residual_locations.py PATH/TO/RUN --output PATH/TO/RUN/plots
```

Run between map writes to avoid reading a partially written snapshot.
Completed runs generate these two PNGs automatically. Indigo intensity uses
a shared power-law scale within each residual type, with zero white and
unmeasured cells grey. The peak JSON gives exact locations and source files.
Full snapshot storage grows linearly with iteration count.

For the completed 260920-183515 validation batch, the three available maps
show the largest component defect at (0,-12.6), (-1.2,-11.8), then (0,11.6),
in xi0 units. Point RMS instead ends at (-1.2,-11.8). The residual remains
concentrated around the half cores; bands near the fine-grid transition at
|y|=16 are also visible. Their numerical origin needs a mesh-sensitivity
test, not an assumption that they represent physical structure. The old
run did not save every iteration's spatial map, so missing frames cannot
be reconstructed from its scalar history.
