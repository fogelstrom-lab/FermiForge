# Alternating Anderson and simple mixing

## Restart-only comparison (2026-09-27)

Set `--anderson-cycle 20 --simple-cycle 0` to reset Anderson at updates
21, 41, 61, etc., without any intervening simple-update block. The reset
retains the physical field but discards the entire accelerator history;
the first update of each block still uses the legacy AA startup factor 0.01.

Run `bash tools/compare_0plus_aa_restart.sh` from the project directory.
It uses the **same initial field** as the previous AA20/simple3 experiment,
not its final field, to isolate the change in schedule. It requests up to
60 updates on 10 ranks with a checkpoint every update in a new run directory.
At the latest measured 554 seconds/map, allow roughly 9.2 hours.
No production calculation is launched by preparing this script.

To continue the improved field instead (not a same-start comparison), append:

```sh
--restart runs/260926-205451-single-harmonic-0plus-cycled-aa20-simple3/final_fields_2d.dat
```

The regression command `python3 tests/check_iteration_cycle.py
work/cycled-iteration-build --simple-steps 0` tests the restart-only schedule.
See also [BB code review](BARZILAI_BORWEIN_REVIEW.md).

The full 2D solver driven by `run_axisymmetric_core_benchmark.py` now supports
a repeating schedule, independent of the chosen vortex seed:

- Updates 1–20: legacy Anderson mixing.
- Updates 21–23: simple damped updates with no Anderson extrapolation.
- Updates 24–43: Anderson with completely reset history.
- Updates 44–46: simple damped updates, and so on.

The simple step is `X_next = X + p (F(X) - X)`, with a conservative initial
choice `p=0.1`. This is a configurable choice, not an experimentally optimized
value. Use `--simple-mixing 1` for undamped direct substitution if desired.
All active order-parameter and Fermi-liquid components are treated together.
Inactive points retain the existing halo policy. State-dependent asymptotics
are refreshed after either type of update as before.

Simple steps do not enter the Anderson history. At each new AA block, all
history, adaptive mixing and internal restart state are reset; the first
AA step uses the legacy initial factor 0.01. The maximum-component convergence
test still applies to the raw map residual in both engines, and can stop the
run before a block finishes.

## Continue the latest 0+ test

From the project directory:

```sh
bash tools/continue_0plus_cycled.sh
```

This rebuilds in Release mode, preserves the existing run, and starts a new
timestamped run from its final field, with 10 ranks, 46 updates and a checkpoint
every update. At roughly 20 minutes per map, budget about 15 hours.
To request another duration, append e.g. `--max-iterations 69`.
No production run has been launched automatically.

For other inputs, use the existing runner with:

```text
--anderson-cycle 20 --simple-cycle 3 --simple-mixing 0.1
```

The corresponding namelist keys are `anderson_cycle_iterations`,
`simple_cycle_iterations` and `simple_mixing`. Cycling is disabled by default
(`anderson_cycle_iterations=0`); old inputs remain unchanged.
The runner can add these optional keys to archived inputs.
Restart files contain fields, not iteration-engine state: every new invocation
starts a fresh AA block at iteration 1.

The terminal labels each update AA or simple. History adds `engine`
(1=AA, 0=simple, -1=single-map) and `cycle_position` columns. Convergence plots
shade and hatch simple intervals and mark fresh-history transitions.
Existing histories without these columns remain readable.
Checkpoint plotting works as before; for this non-axisymmetric state use
`plot_running_axisymmetric.py --latest single-harmonic-0plus-cycled --no-radial-comparison`.

This schedule is implemented in the general 2D axisymmetric-benchmark driver
used by the recent localized-seed runs. It does not change the legacy radial
iterator or the separate historical double-core driver.

## Verification

The field-adapter regression checks damping of both fields, frozen points,
unchanged AA history during simple steps, convergence detection, and fresh AA
initialization. A bounded integration check is available:

```sh
python3 tests/check_iteration_cycle.py work/cycled-iteration-build
```

It runs 24 tiny-grid updates on two MPI ranks and verifies the actual
20 AA / 3 simple / 1 fresh AA sequence and its convergence plot.
It tests control flow, not vortex accuracy or improved convergence.
