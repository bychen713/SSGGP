#!/bin/bash
# Dense-fit prediction (sim/main/sim_dense.R), one task per replicate. Submit from paper/.
#SBATCH --job-name=sim_pred_dense
#SBATCH --cpus-per-task=11
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=3:59:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=sim/output/slurm-%j.out
#SBATCH --array=0-14
VALUES=({1..1000})
export THISJOBVALUE=${VALUES[$SLURM_ARRAY_TASK_ID]}
export out_file_main=sim/output/sim_pred_dense$((THISJOBVALUE)).out

module load r/4.4.1 2>/dev/null || module load r 2>/dev/null
module load gsl/2.8 2>/dev/null || module load gsl 2>/dev/null   # GSL for the CH family

# Run from paper/ (the submit directory) with the repository's renv library.
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../../.." && pwd)}"
export RENV_PROJECT="$(cd "$PROJ_ROOT/.." && pwd)"   # repo root holds renv/
export R_PROFILE_USER="$RENV_PROJECT/.Rprofile"
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
cd "$PROJ_ROOT"

R CMD BATCH --no-save --no-restore sim/main/sim_dense.R $out_file_main
