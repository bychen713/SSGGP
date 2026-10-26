#!/bin/bash
# Figure S2, p = 9: likelihood error against L, one replicate per task.
# Submit from paper/: SSGGP_FAMILY=mixed sbatch slurm/sim/sensitivity_L/p9.sh
#SBATCH --job-name=sensL-p9
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --time=1:00:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --array=1-15
#SBATCH --output=sim/output/sensL-p9-%A_%a.out

set -euo pipefail
module load r/4.4.1 2>/dev/null || module load r 2>/dev/null
module load gsl/2.8 2>/dev/null || module load gsl 2>/dev/null

# Run from paper/ (the submit directory) with the repository's renv library.
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../../.." && pwd)}"
export RENV_PROJECT="$(cd "$PROJ_ROOT/.." && pwd)"
export R_PROFILE_USER="$RENV_PROJECT/.Rprofile" RENV_CONFIG_SANDBOX_ENABLED=FALSE
cd "$PROJ_ROOT"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}" OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}" OMP_PROC_BIND=close OMP_PLACES=cores
export SSGGP_REPLICATE="${SLURM_ARRAY_TASK_ID:?submit as an array}"
export SSGGP_FAMILY="${SSGGP_FAMILY:-mixed}"
Rscript --no-save --no-restore sim/sensitivity_L/p9.R
