# Continue the completed 100-update BB run

Run from any directory:

```sh
bash /Users/mikael/Documents/Codex/3he-vortex-modernization/tools/continue_0plus_bb_100.sh
```

The launcher starts from final_fields_2d.dat in
runs/260928-232830-single-harmonic-0plus-bb-bb-radial-style, using that run's
input file. It requests 100 additional updates on 10 ranks with a checkpoint
every update. BB remains absolute-curvature, initial alpha=1, maximum=100,
minimum=0.001 and growth damping off. Temperature, feedback, mesh, trajectory
sampling and asymptotics are unchanged. Convergence can stop the run early.
Allow approximately 15–16 hours at the previous measured rate.

The executable is rebuilt in work/bb-continuation-build; results go into a
new dated single-harmonic-0plus-bb-continued directory. No prior field,
history or plot is overwritten. The macOS PMIx shutdown mitigation and
diagnostics are applied by the normal runner.

BB history is reset, not restored. This is a physical-field continuation,
not a bitwise uninterrupted 200-step BB trajectory.

## Concatenated output

The new directory contains:

- parent_iteration_history.dat: archived copy of the original 100-row history.
- iteration_history.dat: new local history (1–100).
- combined_iteration_history.dat: generated after the solver returns, containing
  old 1–100 and new 101–200 (or fewer if convergence stops early).
- plots/combined_convergence.png: the combined RMS/max convergence and mixing
  plot with a dashed restart marker before iteration 101.

The combined file preserves all diagnostic columns and adds local_iteration,
segment and history_restart. The new manifest identifies the parent and its
archived history checksum. Launcher failure remains a failure even when
saved products can be postprocessed; combining histories does not certify
MPI shutdown or physical convergence.

The ordinary Cartesian and harmonic axis profiles are produced too. Each
segment's fields and spatial-residual snapshots remain separate: unlike a
scalar convergence history, those files must not simply be appended.

While the calculation is running, inspect its checkpoint as before:

```sh
python3 tools/plot_running_axisymmetric.py --latest single-harmonic-0plus-bb-continued --no-radial-comparison
```

These live plots show the new segment's local iteration numbers. The combined
plot is created by the completion postprocessing.

For later continuations, the runner's --history-parent-run option can use
the parent's combined history when present, extending the chain without
duplicating earlier segments. It requires a matching saved final/checkpoint
field and a completed parent segment.
