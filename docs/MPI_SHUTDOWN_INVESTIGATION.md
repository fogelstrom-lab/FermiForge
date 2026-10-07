# MPI abnormal-exit investigation (2026-09-22)

## 2026-09-29: acknowledgement-timeout reproducer and scoped mitigation

All ten audit files from the completed 100-update run
260928-232830-single-harmonic-0plus-bb-bb-radial-style report successful
MPI_Finalize return and MPI_Finalized=T, despite the launcher reporting
rank 9 as an improper exit. That rules out a missing application finalize
call for this occurrence.

Installed versions: Open MPI 5.0.9, PMIx 5.0.9, bundled PRRTE
3.0.12a12025-10-30. Executable linkage resolves to Homebrew Open MPI.
PMIx reports a default finalize acknowledgement timeout of 10 seconds.
In its source, the timeout callback wakes the waiter and finalization
continues to return PMIX_SUCCESS:
https://github.com/openpmix/openpmix/blob/v5.0.9/src/client/pmix_client.c

The isolated tools/check_mpi_shutdown.py experiment produced:

| Probe | Launcher result | PMIx timeout messages | Rank audits |
|---|---|---|---|
| No delay, 10-second timeout | 0 | 0 | finalized successfully |
| Launcher paused 15 seconds, timeout 10 | 1, PMIX_ERR_UNREACH | 2 | finalized successfully |
| Same pause, timeout 60 | 0 | 0 | finalized successfully |

Evidence is retained in work/mpi-shutdown-dvdfnf65/results.json and logs.
**Fault injection explicitly sets OMPI_MCA_async_mpi_finalize=1 only in the
test**, bypassing Open MPI's initial finalize fence so a launcher pause can
exercise the later ACK stage. Without this test setting, the earlier probe
waited at the fence and exited cleanly. The injected error is NOT the exact
production improper-exit message. This demonstrates a launcher failure despite
successful MPI_Finalize audits, not a reproduction of the exact production
failure or proof of its trigger. Production does not set that flag.
Do not conclude that App Nap, sleep, or network changes caused the production
delay; those have not been established.

The current 2D Python runner now sets, on macOS only and only if not already
specified by the user:

    PMIX_MCA_pmix_finalize_timeout=60
    PMIX_MCA_pmix_client_base_verbose=2

This is a bounded mitigation and diagnostic trial. It does not suppress
warnings, skip finalization, or turn launcher failure into success. Healthy
finalization returns as soon as the acknowledgement arrives; it does not
always wait 60 seconds. Client startup/shutdown messages go into run.log and
are also displayed by the streaming runner. They distinguish
"finwait_cbfunc received" from "finwait timeout fired".

The audit now records the configured PMIx controls and time inside
MPI_Finalize. The run manifest records the selected settings too. Existing
audit files and scientific data are unchanged.

The next real long run must validate this mitigation. If the warning recurs
without a timeout, or despite the extended allowance, inspect the new logs
before claiming this is fixed or changing the installed MPI stack.
The generic double-core runner and direct radial mpiexec invocations have
not been automatically changed. To test the same controls there, prefix the
existing command with:

    env PMIX_MCA_pmix_finalize_timeout=60 PMIX_MCA_pmix_client_base_verbose=2

This does not require changing .bashrc or reinstalling MPI. Do not use timeout
0 as "unlimited": this PMIx implementation schedules it as an immediate timer.

To repeat the bounded two-rank fault-injection experiment:

    python3 tools/check_mpi_shutdown.py

It suspends/resumes only the launcher it creates, never an existing solver.
It takes roughly 40 seconds and retains a unique report folder.

## Earlier investigation

The recent production runs wrote outputs and subsequently received an
Open MPI/PRRTE improper-exit warning. The warning is not a numerical
convergence diagnostic. However, saved rank-zero outputs do not prove
that all workers completed MPI_Finalize; the earlier wrapper's recovery
comment attributing this to sleep/resume was unsupported and removed.

Source inspection found explicit MPI_Finalize calls on the normal paths
of both drivers, after a collective result broadcast. The field-map MPI
communication uses blocking collectives, with no outstanding nonblocking
requests found in that module. The linked libraries are Homebrew Open MPI,
and mpifort reports Open MPI 5.0.9. These checks do not identify the root
cause of the long-run warning. We have not changed or upgraded the toolchain,
or hidden the warning using -quiet.

Both drivers now write a separate audit file for every rank, alongside the
metrics file. It records the runtime library version, entry to MPI_Finalize,
its return code, and MPI_Finalized after return. The audit changes local
logging, not the collective algorithm. MPI_Finalized is legal after
MPI_Finalize; see the official MPI_Finalize documentation:
https://docs.open-mpi.org/en/main/man-openmpi/man3/MPI_Finalize.3.html

Interpretation on the next occurrence:

- Missing file: that rank did not reach the audited shutdown point, or
  could not write it; inspect preceding errors.
- Entry only: that rank entered but did not record return from finalization.
- Every rank reports successful return and finalized=true, but launcher
  reports failure: focus investigation on the runtime/launcher or process
  exit after finalization, rather than a missing application finalize call.

The double-core runner still preserves and plots completed products on a
late launcher error, but now reports mpi_launcher_clean=false and overall
passed=false with exit status 1. calculation_passed is separately recorded.
This is deliberately stricter than the previous recovery behavior. An old
report saying passed=true despite launcher_warning is not proof of clean
MPI termination. The axisymmetric runner already treats launcher failure
as failure; raw files and audit records remain available.

An instrumented short test can establish that the audited path works; it
cannot rule out the failure seen only after multi-hour production runs.
