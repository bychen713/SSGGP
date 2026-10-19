#!/bin/bash
# IOX baseline (sim/iox/sim_iox.R) for one setting, one task per replicate:
# estimation, then prediction. Submit from paper/:
#   IOX_SETTING=random_neg_coloc sbatch slurm/sim/iox/sim_iox.sh
#SBATCH --job-name=sim_iox
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=3:59:00
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --array=0-14
#SBATCH --output=sim/output/%x_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from paper/}"
: "${IOX_SETTING:?IOX_SETTING must be set}"
if [[ $IOX_SETTING == *_miss* ]]; then export IOX_LAYOUT=misaligned; else export IOX_LAYOUT=colocated; fi
export IOX_REP=$((SLURM_ARRAY_TASK_ID + 1))
export IOX_MCMC=${IOX_MCMC:-10000} IOX_NEIGHBORS=${IOX_NEIGHBORS:-15}
export IOX_PRED_DRAWS=${IOX_PRED_DRAWS:-$((IOX_MCMC / 4))}
module load gcc/8.5.0 r/4.4.1 2>/dev/null || module load r 2>/dev/null
export RENV_PROJECT="$(cd .. && pwd)"
export R_PROFILE_USER="$RENV_PROJECT/.Rprofile" RENV_CONFIG_SANDBOX_ENABLED=FALSE
export OMP_PROC_BIND=close OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
echo "IOX $IOX_LAYOUT setting=$IOX_SETTING replicate=$IOX_REP iterations=$IOX_MCMC"
IOX_STAGE=estimate Rscript --no-save --no-restore sim/iox/sim_iox.R
IOX_STAGE=predict  Rscript --no-save --no-restore sim/iox/sim_iox.R
