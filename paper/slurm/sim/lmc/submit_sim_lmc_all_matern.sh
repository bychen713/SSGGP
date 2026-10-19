#!/bin/bash
# Submit both LMC stages for the four tridiagonal all-Matern settings (15
# replicates each; sim_lmc_structure.R uses six latent factors for these).
# Results go to sim/output/simulation/<setting>/.
# Usage (from paper/):
#   bash slurm/sim/lmc/submit_sim_lmc_all_matern.sh [--dry-run|--continue]
set -euo pipefail
if (( $# > 1 )) || [[ ${1:-} != "" && ${1:-} != --dry-run && ${1:-} != --continue ]]; then
  echo "Usage: bash slurm/sim/lmc/submit_sim_lmc_all_matern.sh [--dry-run|--continue]" >&2
  exit 1
fi
continue_lmc=0
[[ ${1:-} == --continue ]] && continue_lmc=1
project_root=$(cd "$(dirname "$0")/../../.." && pwd -P)
cd "$project_root"
datasets=(
  tridiagonal_neg_coloc tridiagonal_neg_miss
  tridiagonal_pos_coloc tridiagonal_pos_miss
)

# Check every input and output before submitting anything.
for dataset in "${datasets[@]}"; do
  for replicate in {1..15}; do
    [[ -f "sim/data/$dataset/$replicate.RData" ]] || {
      echo "Missing input: sim/data/$dataset/$replicate.RData" >&2; exit 1; }
  done
  if (( continue_lmc )); then
    for replicate in {1..15}; do
      test -f "sim/output/simulation/$dataset/lmc_structure_$replicate.RData"
      if [[ ! -f "sim/output/simulation/$dataset/lmcscov_$replicate.RData" &&
            ! -f "sim/output/simulation/$dataset/lmc_checkpoint_$replicate.rds" ]]; then
        echo "No LMC state to continue for $dataset replicate $replicate" >&2
        exit 1
      fi
    done
  elif compgen -G "sim/output/simulation/$dataset/lmc*_*.RData" > /dev/null ||
       compgen -G "sim/output/simulation/$dataset/lmc_checkpoint_*.rds" > /dev/null; then
    echo "LMC results/checkpoints already exist in sim/output/simulation/$dataset" >&2
    exit 1
  fi
done
if [[ ${1:-} == --dry-run ]]; then
  printf '%s: replicates 1-15\n' "${datasets[@]}"
  exit 0
fi

mkdir -p sim/output
for dataset in "${datasets[@]}"; do
  if (( continue_lmc )); then
    prediction=$(sbatch --parsable --array=0-14 --job-name="lmc-mat-cont-$dataset" \
      --export="ALL,SSGGP_LMC_DATASET=$dataset,SSGGP_LMC_CONTINUE=true" \
      slurm/sim/lmc/sim_lmc_pred.sh)
    prediction=${prediction%%;*}
    echo "$dataset: continuation $prediction"
  else
    structure=$(sbatch --parsable --array=0-14 --job-name="lmc-mat-structure-$dataset" \
      --export="ALL,SSGGP_LMC_DATASET=$dataset" slurm/sim/lmc/sim_lmc_structure.sh)
    structure=${structure%%;*}
    prediction=$(sbatch --parsable --array=0-14 --job-name="lmc-mat-pred-$dataset" \
      --dependency="aftercorr:$structure" --kill-on-invalid-dep=yes \
      --export="ALL,SSGGP_LMC_DATASET=$dataset" slurm/sim/lmc/sim_lmc_pred.sh)
    prediction=${prediction%%;*}
    echo "$dataset: structure $structure -> prediction $prediction"
  fi
done
