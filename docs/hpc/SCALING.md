# Portable CPU/MPI scaling benchmark

The harness is `tools/benchmark_mpi_scaling.py`. It never compiles, submits a
scheduler job or performs a self-consistency relaxation. Default invocation
only PREPARES a uniquely named campaign. `--execute` is required to launch.
No laptop MPI transport restrictions are injected; remove any inherited
`OMPI_MCA_btl=self,sm` before multi-node use. Existing production inputs,
builds and runs are not changed.

## 1. Build with the centre's toolchain

Obtain a compute allocation and load mutually compatible compiler, MPI,
CMake and Python modules following centre documentation. In the project root:

```sh
cmake -S . -B work/hpc-build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=mpifort -DFERMIFORGE_REQUIRE_MPI=ON
cmake --build work/hpc-build --parallel 8
ctest --test-dir work/hpc-build --output-on-failure
```

The centre may require its compiler wrapper instead of mpifort. CTest's MPI
launcher must likewise match the installation. Do not reuse a build tree made
with a different compiler/MPI or one serving an active run.

## 2. Prepare without launching

```sh
python3 tools/benchmark_mpi_scaling.py \
  --executable work/hpc-build/benchmark_axisymmetric_core_2d \
  --ranks 1 2 4 8 \
  --launcher 'srun --exclusive --nodes=1 --ntasks={ranks} --cpus-per-task=1 --cpu-bind=cores'
```

This archives the input/dependencies and CMake cache, hashes the executable,
records Git state and selected MPI/scheduler environment, and writes a plan
under work/hpc-scaling. It executes no solver. Paths inside the archived
namelist point to the new archive; prepare it on the target filesystem.
The plan is for inspection: a later --execute invocation creates a fresh
campaign rather than reusing possibly stale prepared files.

## 3. Run a bounded smoke campaign

Use the same command with `--execute --timeout 300 --budget-seconds 3000`
INSIDE a compute allocation with at least eight physical cores. The default
example has 113 active points, short trajectories, full angular/frequency
tables and a single map. It is an infrastructure test, not a scientifically
converged vortex or a basis for extrapolating large-node scaling.

`tools/hpc_scaling.slurm` is a submission template for this one-node test.
Supply the real account, CPU partition and resources on the sbatch command
line; adjust MPI plugin/binding options with centre support. It has not been
submitted. The Slurm allocation is held for the entire serial campaign,
including repetitions with fewer active ranks, so allocation charges can
exceed sum(ranks * step elapsed).

Each rank count has one excluded warm-up launch and three measured launches
by default. These are fresh processes mapping the SAME initial field; warm-up
primarily removes first-launch/filesystem effects, not an in-process cache
warm-up. A single launch is bounded by --timeout and the campaign by
--budget-seconds (plus a short termination allowance). Failure or disagreement
stops the campaign and preserves logs; it is not silently treated as timing data.

## 4. Representative scaling, then multiple nodes

After smoke validation, supply --input with a representative full-2D namelist
and increase the budget only deliberately. All physical/numerical settings
are retained except maximum_iterations=0, checkpoint_interval=0,
prepare_only=false, require_convergence=false. No mixing update is performed.
For a restart, transfer the referenced field and repair its path before
preparing the campaign. Input data paths are resolved relative to the project
directory, consistent with the solver. The four reference/quadrature paths
must be explicitly quoted in the input.

Suggested stages: 1/2/4/8 ranks on the small case; then 8/16/32/64 on a
representative case; then suitable full-node and multi-node sizes based on
measured utilization. These are suggestions, not pre-approved allocations.
Use a site-appropriate launcher template with a deliberate node/task layout.
The harness does NOT infer node counts or oversubscribe a scheduler allocation.
For a fixed two-node layout, for example, replace --nodes=1 with --nodes=2
and choose rank counts compatible with the requested allocation. Compare
campaigns with the same input/executable hashes and documented placement.
Do not spend many nodes on the 60-point radial problem.

## Products and interpretation

- manifest.json: plan, hashes, build cache provenance, environment, final status.
- Per launch: run.log, mpi.log, inputs/outputs, metrics and solver shutdown audit.
- measurements.json: incremental raw measurements, including warm-up launches.
- summary.csv: medians/ranges, speedups and relative efficiencies against the
  smallest tested rank count; also rank-hours per complete launch.

`map_seconds` is the maximum rank-local point-map time and excludes subsequent
MPI reductions. Complete elapsed time includes startup, input, collectives,
mandatory field/diagnostic writing and finalization. Neither isolates the
nonlinear accelerator cost, which this fixed-map test deliberately omits.
Use later short iteration tests to measure full production-loop performance.
Field agreement defaults to atol=1e-12 plus rtol=1e-10, checked on all numeric
columns of mapped fields (not just scalar residual norms).

Mapped-field outputs are checked for finite values; a
nonzero MPI exit invalidates a measurement even if metrics say passed.
Scheduler accounting is authoritative for billable resources. On Slurm:

```sh
sacct -j JOB_ID --units=G \
  --format=JobID,State,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ExitCode
```

Record node type, number of nodes, compiler/MPI versions, binding, tasks/node,
physical-core versus SMT usage, and filesystem. MaxRSS is generally a
per-task maximum, not aggregate replicated memory; interpret it with centre
guidance and inspect the solver steps, not only the batch-shell record.
Output-file size is not RAM usage. Replace laptop-based allocation estimates
only after representative cluster timings and charging rules are available.

## Local harness validation (2026-10-01)

Five unit tests passed, including rejection of a failed launcher. A bounded
113-point single-map check passed on one and two local MPI ranks with identical
written mapped fields. Map times were 10.61 and 5.57 seconds; complete-launch
times were 10.72 and 5.72 seconds. This one-repetition smoke test validates the
measurement path, not cluster scaling or production performance. Its records
are under work/hpc-scaling/261001-083646-781960. No scheduler job was submitted.
