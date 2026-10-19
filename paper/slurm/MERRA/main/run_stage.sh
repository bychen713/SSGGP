#!/bin/bash
# One SSGGP stage (MERRA/main/): independent, estimate, likelihood, or predict.
# Submitted by submit_all.sh from paper/.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=3:59:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=MERRA/output/%x-%A_%a.out
set -euo pipefail
stage=${1:?stage required}
cd "${SLURM_SUBMIT_DIR:?submit from paper/}"
module load r/4.4.1 2>/dev/null || module load r
module load gsl/2.8 2>/dev/null || true
# The repository renv library (the repository root is the parent of paper/).
export RENV_PROJECT="$(cd .. && pwd)"
export R_PROFILE_USER="$RENV_PROJECT/.Rprofile" RENV_CONFIG_SANDBOX_ENABLED=FALSE
export MERRA_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
case "$stage" in
  independent) Rscript --no-save --no-restore MERRA/main/independent.R ;;
  estimate|likelihood|predict) Rscript --no-save --no-restore MERRA/main/joint.R "$stage" ;;
  *) echo "unknown stage: $stage" >&2; exit 2 ;;
esac
