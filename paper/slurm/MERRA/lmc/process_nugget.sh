#!/bin/bash
# Process-specific-nugget LMC (MERRA/lmc/process_nugget.R); mode pilot, fit,
# or continue. Needs the shared-nugget LMC fit of the same design first.
#   MERRA_DESIGN=colocated sbatch --array=0-2 slurm/MERRA/lmc/process_nugget.sh pilot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=3:59:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=MERRA/output/%x-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from paper/}"
: "${MERRA_DESIGN:?MERRA_DESIGN must be colocated or misaligned}"
module load r/4.4.1 2>/dev/null || module load r
module load gsl/2.8 2>/dev/null || true
export RENV_PROJECT="$(cd .. && pwd)"
export R_PROFILE_USER="$RENV_PROJECT/.Rprofile" RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=1
Rscript --no-save --no-restore MERRA/lmc/process_nugget.R "${1:?mode required}"
