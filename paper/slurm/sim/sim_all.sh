#!/bin/bash
# Submit simulation stages 1-3 with Slurm dependencies. Run from paper/ after
# data generation:  bash slurm/sim/sim_all.sh
set -euo pipefail
cd "$(cd "$(dirname "$0")/../.." && pwd)"
mkdir -p sim/output

JOB1_MAIN="slurm/sim/main/sim_main.sh"
JOB1_IND="slurm/sim/main/sim_ind.sh"
JOB2="slurm/sim/main/sim_ll.sh"
JOB3_BIC="slurm/sim/main/sim_BIC.sh"
JOB3_DENSE="slurm/sim/main/sim_dense.sh"

# The joint fits start from the independent fits.
jid1_ind=$(sbatch --parsable "$JOB1_IND" | cut -d';' -f1)
echo "Submitted $JOB1_IND as job $jid1_ind"

jid1_main=$(sbatch --parsable --dependency=afterok:"$jid1_ind" "$JOB1_MAIN" | cut -d';' -f1)
echo "Submitted $JOB1_MAIN as job $jid1_main (depends on $jid1_ind)"

jid2=$(sbatch --parsable --dependency=afterok:"$jid1_main" "$JOB2" | cut -d';' -f1)
echo "Submitted $JOB2 as job $jid2 (depends on $jid1_main)"

# Dense prediction needs only the estimates; BIC prediction also needs the likelihoods.
jid3_dense=$(sbatch --parsable --dependency=afterok:"$jid1_main" "$JOB3_DENSE" | cut -d';' -f1)
echo "Submitted $JOB3_DENSE as job $jid3_dense (depends on $jid1_main)"

jid3_bic=$(sbatch --parsable --dependency=afterok:"$jid2" "$JOB3_BIC" | cut -d';' -f1)
echo "Submitted $JOB3_BIC as job $jid3_bic (depends on $jid2)"

echo "Submitted chain:"
echo "  JOB1 main: $jid1_main"
echo "  JOB1 ind : $jid1_ind (independent)"
echo "  JOB2     : $jid2"
echo "  JOB3 BIC : $jid3_bic"
echo "  JOB3 dense: $jid3_dense"
