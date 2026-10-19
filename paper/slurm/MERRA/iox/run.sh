#!/bin/bash
# IOX baseline (MERRA/iox/run.R), one task per fold: estimation, then
# prediction. Submit from paper/:
#   MERRA_DESIGN=colocated sbatch --array=0-2 slurm/MERRA/iox/run.sh
#SBATCH --job-name=me-iox
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=3:59:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=MERRA/output/%x-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from paper/}"
: "${MERRA_DESIGN:?MERRA_DESIGN must be colocated or misaligned}"
export IOX_LAYOUT=$MERRA_DESIGN IOX_REP=$((SLURM_ARRAY_TASK_ID + 1))
export IOX_MCMC=${IOX_MCMC:-10000} IOX_NEIGHBORS=${IOX_NEIGHBORS:-15}
export IOX_PRED_DRAWS=${IOX_PRED_DRAWS:-$((IOX_MCMC / 4))}
module load gcc/8.5.0 r/4.4.1 2>/dev/null || module load r
export RENV_PROJECT="$(cd .. && pwd)"
export R_PROFILE_USER="$RENV_PROJECT/.Rprofile" RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_PROC_BIND=close OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
echo "IOX MERRA $IOX_LAYOUT fold=$IOX_REP iterations=$IOX_MCMC"
IOX_STAGE=estimate Rscript --no-save --no-restore MERRA/iox/run.R
IOX_STAGE=predict  Rscript --no-save --no-restore MERRA/iox/run.R
