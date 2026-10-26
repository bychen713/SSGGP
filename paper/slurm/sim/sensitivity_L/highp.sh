#!/bin/bash
# Figure S2, p = 25: stage integrated, then simulated (which needs the first).
# Submit from paper/: SSGGP_P=25 SSGGP_FAMILY=mixed SSGGP_STAGE=integrated sbatch slurm/sim/sensitivity_L/highp.sh
#SBATCH --job-name=sensL-hp
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64G
#SBATCH --time=3:59:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --array=1-15
#SBATCH --output=sim/output/sensL-hp-%A_%a.out

set -euo pipefail
module load r/4.4.1 2>/dev/null || module load r 2>/dev/null
module load gsl/2.8 2>/dev/null || module load gsl 2>/dev/null

# Run from paper/ (the submit directory) with the repository's renv library.
PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../../.." && pwd)}"
export RENV_PROJECT="$(cd "$PROJ_ROOT/.." && pwd)"
export R_PROFILE_USER="$RENV_PROJECT/.Rprofile" RENV_CONFIG_SANDBOX_ENABLED=FALSE
cd "$PROJ_ROOT"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}" OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}" OMP_PROC_BIND=close OMP_PLACES=cores
export SSGGP_REPLICATE="${SLURM_ARRAY_TASK_ID:?submit as an array}"
export SSGGP_P="${SSGGP_P:?set SSGGP_P to 25 or 50}"
export SSGGP_FAMILY="${SSGGP_FAMILY:?set SSGGP_FAMILY to matern or mixed}"
export SSGGP_STAGE="${SSGGP_STAGE:?set SSGGP_STAGE to integrated or simulated}"
Rscript --no-save --no-restore sim/sensitivity_L/highp.R
