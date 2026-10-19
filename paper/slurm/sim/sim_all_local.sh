#!/bin/bash
# Run simulation stages 0-3 without Slurm by executing the same launchers locally.
# Usage (from paper/):  bash slurm/sim/sim_all_local.sh [--parallel N]
# N array tasks run at once (default 1); each uses `threads` cores (sim/settings.R).
# Order: data -> sim_ind -> sim_main -> sim_ll -> sim_BIC + sim_dense. Each
# stage's task range is read from its "#SBATCH --array=A-B" line.
set -euo pipefail
cd "$(cd "$(dirname "$0")/../.." && pwd)"   # paper/
mkdir -p sim/output

PARALLEL=1
case "${1:-}" in
  --parallel|-j) PARALLEL="${2:?usage: bash slurm/sim/sim_all_local.sh [--parallel N]}";;
  "")            : ;;
  *)             echo "usage: bash slurm/sim/sim_all_local.sh [--parallel N]   (N = array elements at once)"; exit 2;;
esac

run_stage() {  # run a launcher over its array range, PARALLEL tasks at a time
  local script=$1 range A B
  range=$(grep -oE '#SBATCH --array=[0-9]+-[0-9]+' "$script" | grep -oE '[0-9]+-[0-9]+')
  A=${range%-*}; B=${range#*-}
  echo ">>> $script  (array $A-$B, $PARALLEL at a time)"
  seq "$A" "$B" | xargs -P "$PARALLEL" -I{} bash -c \
    'SLURM_ARRAY_TASK_ID={} bash "$1" || { echo "FAILED: $1 [element {}]"; exit 255; }' _ "$script"
}

# stage 0: data
bash slurm/sim/data_generation.sh

# stage 1: independent fits, then the joint fits that start from them
run_stage slurm/sim/main/sim_ind.sh
run_stage slurm/sim/main/sim_main.sh

# stage 2: likelihoods
run_stage slurm/sim/main/sim_ll.sh

# stage 3: BIC and dense predictions, concurrently
run_stage slurm/sim/main/sim_BIC.sh   &  pid_bic=$!
run_stage slurm/sim/main/sim_dense.sh &  pid_dense=$!
wait "$pid_bic"; wait "$pid_dense"   # separate waits so either failure stops the run

echo "All stages finished. Logs and results: sim/output/"
