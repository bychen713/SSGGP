#!/bin/bash
# Submit both LMC stages for every simulated setting in sim/data/ (16
# settings, 15 replicates each). Each prediction task waits for its own
# structure task. Results go to sim/output/simulation/<setting>/.
# Usage (from paper/):
#   bash slurm/sim/lmc/submit_sim_lmc_all.sh [--dry-run|--continue]
set -euo pipefail
if (( $# > 1 )) || [[ ${1:-} != "" && ${1:-} != --dry-run && ${1:-} != --continue ]]; then
  echo "Usage: bash slurm/sim/lmc/submit_sim_lmc_all.sh [--dry-run|--continue]" >&2
  exit 1
fi
continue_lmc=0
[[ ${1:-} == --continue ]] && continue_lmc=1
project_root=$(cd "$(dirname "$0")/../../.." && pwd -P)
cd "$project_root"

# Check every input and output before submitting anything.
datasets=()
for dataset_path in sim/data/*; do
  [[ -d "$dataset_path" ]] || continue
  dataset=${dataset_path##*/}
  if [[ ! "$dataset" =~ ^(random|tridiagonal)_(pos|neg)_(coloc|miss)(_mix)?$ ]]; then
    echo "Unexpected dataset directory: $dataset" >&2; exit 1
  fi
  for replicate in {1..15}; do
    [[ -f "$dataset_path/$replicate.RData" ]] || {
      echo "Missing input: $dataset_path/$replicate.RData" >&2; exit 1; }
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
  datasets+=("$dataset")
done
[[ ${#datasets[@]} == 16 ]] || { echo "Expected 16 settings, found ${#datasets[@]}" >&2; exit 1; }
if [[ ${1:-} == --dry-run ]]; then
  printf '%s: replicates 1-15\n' "${datasets[@]}"
  exit 0
fi

mkdir -p sim/output
for dataset in "${datasets[@]}"; do
  if (( continue_lmc )); then
    prediction=$(sbatch --parsable --array=0-14 --job-name="lmc-cont-$dataset" \
      --export="ALL,SSGGP_LMC_DATASET=$dataset,SSGGP_LMC_CONTINUE=true" \
      slurm/sim/lmc/sim_lmc_pred.sh)
    prediction=${prediction%%;*}
    echo "$dataset: continuation $prediction"
  else
    structure=$(sbatch --parsable --array=0-14 --job-name="lmc-structure-$dataset" \
      --export="ALL,SSGGP_LMC_DATASET=$dataset" slurm/sim/lmc/sim_lmc_structure.sh)
    structure=${structure%%;*}
    prediction=$(sbatch --parsable --array=0-14 --job-name="lmc-pred-$dataset" \
      --dependency="aftercorr:$structure" --kill-on-invalid-dep=yes \
      --export="ALL,SSGGP_LMC_DATASET=$dataset" slurm/sim/lmc/sim_lmc_pred.sh)
    prediction=${prediction%%;*}
    echo "$dataset: structure $structure -> prediction $prediction"
  fi
done
