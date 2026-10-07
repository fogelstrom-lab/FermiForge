# Inspecting a running radial-reference benchmark

From the repository directory, run:

```sh
python3 tools/plot_running_axisymmetric.py runs/260923-135645-a-phase-extended-outer-core-transfer
```

Requires the same NumPy, Matplotlib and Pillow environment as the existing
plotting tools. No MPI launch, rebuild or solver evaluation is performed.
The command creates a unique `inspection-<timestamp>-updateNNNN/plots`
directory inside the run, keeping all earlier inspections. Repeat whenever
a new checkpoint is available.

The saved checkpoint, radial reference, input, history and completed residual
maps are copied before plotting. Files changing during a read are rejected;
if checkpoint writing is in progress, retry shortly. The checkpoint grid is
validated against its declared dimensions. Residual/history entries beyond
the checkpoint update are excluded.

Outputs include Cartesian and harmonic amplitude/phase maps, pair-density
and current-related mean-field plots, axis profiles against the preserved
radial reference, component difference maps, convergence, and spatial RMS/max
residual dynamics (including animated PNGs).

Residual map k describes the field **before** update k; the checkpoint is
**after** update k. Small iteration residuals are not proof of agreement with
the radial solution or spatial convergence. Density/current panels retain
their existing proxy interpretation, not a superfluid-density response.

The snapshot's `final_fields_2d.dat` is just a compatibility filename for the
checkpoint copy, not a claim that the running calculation is finished.
`snapshot_manifest.json` records its original source, update and SHA-256 hash.
The synthetic `metrics.txt` only supplies the active plotting radius; it is
not a solver report. Original run outputs are never overwritten.

Every inspection also includes `harmonic_axis_profiles.png`: real and imaginary
local C(s,k) components on the positive x and y axes, using the same shared
normalization and reversed-x/forward-y layout as the Cartesian axis profiles.
The spin/orbital order is (+,0,-); no vortex or angular phase is removed.
Imaginary channels use open-circle markers, with row colours and column line
styles distinguishing components. As with the Cartesian plot, channels below
0.002 of the bulk normalization are hidden by default. For a different cutoff,
run `tools/plot_double_core_axis_profiles.py` with `--basis harmonic` and
`--visibility-threshold` explicitly. Both inspection modes (with and without
radial comparison) produce this checkpoint harmonic plot.
