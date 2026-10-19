#!/bin/bash
# LMC stage 2 (sim/lmc/sim_lmc_pred.R), one task per replicate. Submit from paper/.
#SBATCH --job-name=sim_lmc_pred
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=3:59:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=sim/output/lmc-pred-%A_%a.out
#SBATCH --array=0-14

set -euo pipefail
export THISJOBVALUE=$(( ${SLURM_ARRAY_TASK_ID:-0} + 1 ))

module load r/4.4.1 2>/dev/null || module load r 2>/dev/null
module load gsl/2.8 2>/dev/null || module load gsl 2>/dev/null   # GSL for the CH family

# Run from paper/ (the submit directory) with the repository's renv library.
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../../.." && pwd)}"
export RENV_PROJECT="$(cd "$PROJ_ROOT/.." && pwd)"   # repo root holds renv/
export R_PROFILE_USER="$RENV_PROJECT/.Rprofile"
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
cd "$PROJ_ROOT"
test -f sim/lmc/sim_lmc_pred.R
dataset_tag=${SSGGP_LMC_DATASET:-default}
[[ "$dataset_tag" =~ ^[a-z_]+$ ]] || exit 1
mkdir -p sim/output
out_file_lmc="sim/output/lmc_${dataset_tag}_pred_${THISJOBVALUE}_${SLURM_ARRAY_JOB_ID:-local}.Rout"

# BLAS threads match --cpus-per-task.
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

# source() parses the whole script before running it.
Rscript --no-save --no-restore -e 'source("sim/lmc/sim_lmc_pred.R", echo = TRUE)' \
  > "$out_file_lmc" 2>&1
