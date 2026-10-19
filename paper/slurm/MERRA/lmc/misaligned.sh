#!/bin/bash
# LMC for the misaligned design; mode pilot, fit, or continue.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=3:59:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=MERRA/output/%x-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from paper/}"
module load r/4.4.1 2>/dev/null || module load r
module load gsl/2.8 2>/dev/null || true
# The repository renv library (the repository root is the parent of paper/).
export RENV_PROJECT="$(cd .. && pwd)"
export R_PROFILE_USER="$RENV_PROJECT/.Rprofile" RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=1
Rscript --no-save --no-restore MERRA/lmc/misaligned.R "${1:?mode required}"
