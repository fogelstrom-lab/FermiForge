#!/usr/bin/env bash
# Template: submit with site-specific account, partition and node resources.
# No allocation/account is guessed and this file does not submit itself.
# Example (edit values first):
# sbatch --account=YOUR_ACCOUNT --partition=YOUR_CPU_PARTITION --nodes=1 \
#   --ntasks=8 --cpus-per-task=1 --time=01:00:00 tools/hpc_scaling.slurm
set -euo pipefail
: "${SLURM_JOB_ID:?Run through sbatch in a compute allocation}"
project_dir="${FERMIFORGE_ROOT:-${SLURM_SUBMIT_DIR}}"
cd "$project_dir"
# Load the site's mutually compatible compiler, MPI, CMake and Python modules
# before submission, or add the module commands supplied by centre support here.
# Build beforehand in an interactive compute allocation; never modify a build
# tree used by an active job.
exe="${FERMIFORGE_EXE:-work/hpc-build/benchmark_axisymmetric_core_2d}"
if [[ ! -x "$exe" ]]; then
  echo "Missing executable: $exe. See docs/hpc/SCALING.md" >&2
  exit 1
fi
# Do not inherit the laptop's self,sm-only transport restriction on a cluster.
if [[ "${OMPI_MCA_btl:-}" == "self,sm" ]]; then
  echo 'Unset the laptop OMPI_MCA_btl=self,sm setting before cluster submission.' >&2
  exit 1
fi
python3 tools/benchmark_mpi_scaling.py --executable "$exe" \
  --ranks 1 2 4 8 --warmups 1 --repeats 3 \
  --launcher 'srun --exclusive --nodes=1 --ntasks={ranks} --cpus-per-task=1 --cpu-bind=cores' \
  --timeout 300 --budget-seconds 3000 --execute
