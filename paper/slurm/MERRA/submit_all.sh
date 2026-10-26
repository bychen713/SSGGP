#!/bin/bash
# Submit every MERRA-2 fit. Run from paper/ after MERRA/prepare_folds.R:
#   bash slurm/MERRA/submit_all.sh
# Continue only nonconverged/interrupted LMC fits with:
#   bash slurm/MERRA/submit_all.sh --continue-lmc
# When all jobs have finished, collect the results with
#   Rscript MERRA/result/summarize.R
set -euo pipefail
if (( $# > 1 )) || [[ ${1:-} != "" && ${1:-} != --continue-lmc ]]; then
  echo "Usage: bash slurm/MERRA/submit_all.sh [--continue-lmc]" >&2
  exit 1
fi
root=$(cd "$(dirname "$0")/../.." && pwd -P)
cd "$root"

data_root=${MERRA_DATA_ROOT:-MERRA/data}
out_root=${MERRA_OUTPUT_ROOT:-MERRA/output}
test -f "$data_root/preprocessed_P${MERRA_P:-40}_uncentered.rds"

if [[ ${1:-} == --continue-lmc ]]; then
  printf 'design\tstage\tjob_id\tdependency\tarray\n'
  for design in misaligned colocated; do
    script=slurm/MERRA/lmc/$design.sh
    lmc_dir=$out_root/$design/lmc
    for fold in 1 2 3; do
      test -f "$data_root/$design/fold${fold}.rds"
      test -f "$lmc_dir/pilot_fold${fold}.rds"
      if [[ ! -f "$lmc_dir/fit_fold${fold}.rds" &&
            ! -f "$lmc_dir/checkpoint_fold${fold}.rds" ]]; then
        echo "No LMC state to continue for $design fold $fold" >&2
        exit 1
      fi
    done
    job=$(sbatch --parsable --export="ALL,MERRA_DESIGN=$design" --array=0-2 \
      --job-name="mfs-${design:0:1}lmc-cont" \
      --output="$out_root/${design}-lmc-continue-%A_%a.out" "$script" continue)
    job=${job%%;*}
    printf '%s\tlmc_continue\t%s\t-\t0-2\n' "$design" "$job"
  done
  exit 0
fi

for design in misaligned colocated; do
  for fold in 1 2 3; do
    test -f "$data_root/$design/fold${fold}.rds"
  done
  test ! -e "$out_root/$design/ssggp"
  test ! -e "$out_root/$design/lmc"
  test ! -e "$out_root/$design/iox"
done
mkdir -p "$out_root"
printf 'design\tstage\tjob_id\tdependency\tarray\n'
record() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5"; }

# SSGGP: independent fits -> 45 estimates (3 folds x 15 penalties) ->
# likelihoods -> dense/BIC prediction.
submit_ssggp() {
  local design=$1 label=$2
  local log=$out_root/${design}-ssggp
  local common="ALL,MERRA_DESIGN=$design"
  local ind est ll pred
  ind=$(sbatch --parsable --export="$common" --array=0-2 \
    --job-name="${label}-ind" --output="$log-ind-%A_%a.out" \
    slurm/MERRA/main/run_stage.sh independent)
  ind=${ind%%;*}; record "$design" independent "$ind" - 0-2
  est=$(sbatch --parsable --export="$common" --array=0-44 \
    --dependency="afterok:$ind" --job-name="${label}-est" \
    --output="$log-est-%A_%a.out" slurm/MERRA/main/run_stage.sh estimate)
  est=${est%%;*}; record "$design" estimate "$est" "afterok:$ind" 0-44
  ll=$(sbatch --parsable --export="$common" --array=0-44 \
    --dependency="afterok:$est" --job-name="${label}-ll" \
    --output="$log-ll-%A_%a.out" slurm/MERRA/main/run_stage.sh likelihood)
  ll=${ll%%;*}; record "$design" likelihood "$ll" "afterok:$est" 0-44
  pred=$(sbatch --parsable --export="$common" --array=0-2 \
    --dependency="afterok:$ll" --job-name="${label}-pred" \
    --output="$log-pred-%A_%a.out" slurm/MERRA/main/run_stage.sh predict)
  pred=${pred%%;*}; record "$design" predict "$pred" "afterok:$ll" 0-2
}

submit_ssggp misaligned mfs-m
submit_ssggp colocated mfs-c

# LMC: pilot -> fit and predict.
submit_lmc() {
  local design=$1 script=$2 label=$3
  local log=$out_root/${design}-lmc
  local common="ALL,MERRA_DESIGN=$design"
  local pilot fit
  pilot=$(sbatch --parsable --export="$common" --array=0-2 \
    --job-name="${label}-pilot" --output="$log-pilot-%A_%a.out" "$script" pilot)
  pilot=${pilot%%;*}; record "$design" lmc_pilot "$pilot" - 0-2
  fit=$(sbatch --parsable --export="$common" --array=0-2 \
    --dependency="afterok:$pilot" --job-name="${label}-fit" \
    --output="$log-fit-%A_%a.out" "$script" fit)
  fit=${fit%%;*}; record "$design" lmc_fit "$fit" "afterok:$pilot" 0-2
}

submit_lmc misaligned slurm/MERRA/lmc/misaligned.sh mfs-mlmc
submit_lmc colocated slurm/MERRA/lmc/colocated.sh mfs-clmc

# IOX: fit and predict, both designs.
for design in misaligned colocated; do
  iox=$(sbatch --parsable --export="ALL,MERRA_DESIGN=$design" --array=0-2 \
    --job-name="mfs-iox-$design" --output="$out_root/${design}-iox-%A_%a.out" \
    slurm/MERRA/iox/run.sh)
  iox=${iox%%;*}; record "$design" iox_fit "$iox" - 0-2
done

