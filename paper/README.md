# Reproducing the paper

This directory holds the scripts behind the simulation study and the MERRA-2
application. To use SSGGP on your own data, see the top-level README and
`../demo.R`.

**Run everything from this directory** (`cd paper`). Generated files go to
`sim/data/`, `sim/output/`, and generated subdirectories of `MERRA/`; none are
versioned.

## Setup

From the repository root, to reproduce everything:

```bash
bash setup_env.sh --with-merra --with-iox
cd paper
```

`bash setup_env.sh` alone installs SSGGP and its dependencies, enough for the
SSGGP and LMC fits. The options add to it and can be combined: `--with-merra`
for MERRA-2 preprocessing and figures, and `--with-iox` for the bundled `spiox`
used by IOX (see [`iox/README.md`](iox/README.md)). `paper/.Rprofile` and the
Slurm launchers activate the same renv library as `../demo.R`.

## Layout

Both studies share one structure: settings and data preparation at the top,
then one folder per method and a `result/` folder for tables and figures.

```text
paper/
├── sim/                    simulation study
│   ├── settings.R          the only file to edit to change the simulation
│   ├── data_generation.R   simulated datasets
│   ├── data/               generated datasets, one folder per setting
│   ├── output/             results, tables, figures, and Slurm logs
│   ├── main/               SSGGP fits, penalty selection, prediction
│   ├── lmc/                LMC baseline
│   ├── iox/                IOX baseline
│   └── result/             tables and figures
├── MERRA/                  MERRA-2 application
│   ├── settings.R
│   ├── preprocess_and_eda.R, prepare_folds.R
│   ├── data/               preprocessed data and folds
│   ├── output/             results and Slurm logs
│   ├── main/  lmc/  iox/   SSGGP and the two baselines
│   └── result/             tables and figures
├── iox/                    IOX code shared by sim/ and MERRA/, and spiox/
└── slurm/                  Slurm launchers, mirroring the folders above
```

`sim/lmc/lmc_model.R` and `MERRA/lmc/lmc_model.R` hold the same LMC
implementation. Every R script can be run directly with `Rscript`; the `slurm/`
launchers only add scheduler directives, array indexing, and job chaining.
Adjust them for your cluster or pass options at submission
(`sbatch -A ACCOUNT ...`). Slurm logs go directly to `sim/output/` or
`MERRA/output/`; Slurm opens its log before the script starts, so that folder
must exist before calling `sbatch` directly.

## Tables and figures

Every table and figure is made from saved results, with no fitting. Run each
command from `paper/` with the setting variables of the listed settings (see
[Simulation](#simulation)). Simulation tables go to `sim/output/tables/<setting>/`
and figures to `sim/output/figure/<setting>/`; MERRA-2 outputs go to
`MERRA/output/` and `MERRA/output/figure/`. Table entries are mean (SD) over
replicates (simulation) or folds (MERRA-2) unless noted.

| Paper | Settings | Command | Output |
|---|---|---|---|
| Table 1, Fig. 1 | `random_{pos,neg}_coloc_mix` | `Rscript sim/result/sim_result.R` | `covariance_estimation_error.csv`; `correlation.png`, `precision.png`, `sigmainv_prec_{BIC,dense}.png` |
| Tables 2, 3 | `random_{pos,neg}_{coloc,miss}_mix` | `sim_result.R`; `sim_iox_result.R` | `prediction.csv` (Penalized = BIC, Dense, Ind, LMC); `prediction_iox.csv` |
| Fig. 2 | the four `random_*_mix` settings | `sim_result.R` (drawn once all four exist) | `sim/output/figure/pr_curves_random.png` |
| Fig. S1 | `random_pos_coloc_mix` | `sim_result.R` | `processes_1_21.png` |
| Fig. S2 | – | `Rscript sim/result/plot_sensitivity_L.R RESULTS_DIR` | `sim/output/figure/sensitivity_L_mixed_p{9,25,50}.png` |
| Fig. S3 | `random_{pos,neg}_coloc_mix` | `SSGGP_COMPUTE_DIAGNOSTIC_LOGLIK=true Rscript sim/result/sim_result.R` | `diagnostics/*.png` |
| Table S1 | all | `sim_result.R` | `computation_time.csv` |
| Table S2 | `random_{pos,neg}_{coloc,miss}_mix` | `sim_result.R` | `graph_recovery.csv` |
| Table S3 | `random_{pos,neg}_miss_mix` | `sim_result.R` | `covariance_estimation_error.csv` |
| Table S4, Fig. S4 | `tridiagonal_{pos,neg}_coloc_mix` | `sim_result.R` | `covariance_estimation_error.csv`; figures as Fig. 1 |
| Table S5 | `tridiagonal_{pos,neg}_{coloc,miss}_mix` | `sim_result.R`; `sim_iox_result.R` | `prediction.csv`; `prediction_iox.csv` |
| Tables S6, S7 | all-Matérn (`SSGGP_FAMILY=matern`): random (S6), tridiagonal (S7) | `sim_result.R`; `sim_iox_result.R` | `estimation_error.csv`; `prediction.csv`; `prediction_iox.csv` |
| Table 4 | – | `Rscript MERRA/result/summarize.R` | `prediction_summary.csv` |
| Fig. 3 | – | `summarize.R` | `data_maps.png` |
| Fig. 4 | – | `summarize.R` | `precision_colocated_{sample,bic,dense}.png` |
| Table S8 | – | `summarize.R` | `estimation_summary.csv` (medians over folds of within-fold means) |

Fig. S2 plots saved results of the likelihood-error experiment: `RESULTS_DIR`
holds `mixed_p9/`, `mixed_p25/` and `mixed_p50/` (see the script header).

## Simulation

All settings are in `sim/settings.R`. Every stage runs one setting, chosen with
`SSGGP_PREC` (`random`/`tridiagonal`), `SSGGP_SIGN` (`neg`/`pos`), `SSGGP_MISS`
(`miss`/`coloc`), and `SSGGP_FAMILY` (`mixed`/`matern`); the first value of each
is the default, and `SSGGP_FAMILY=matern` gives the all-Matérn design.

| Stage | Script | Slurm launcher | Output |
|---|---|---|---|
| 0 | `sim/data_generation.R` | `slurm/sim/data_generation.sh` | `sim/data/<setting>/` |
| 1a | `sim/main/sim_ind.R` | `slurm/sim/main/sim_ind.sh` | independent baseline |
| 1b | `sim/main/sim_main.R` | `slurm/sim/main/sim_main.sh` | penalty-path fits |
| 2 | `sim/main/sim_ll.R` | `slurm/sim/main/sim_ll.sh` | marginal likelihoods |
| 3a | `sim/main/sim_BIC.R` | `slurm/sim/main/sim_BIC.sh` | BIC-selected fit and prediction |
| 3b | `sim/main/sim_dense.R` | `slurm/sim/main/sim_dense.sh` | unpenalized fit and prediction |
| 4 | `sim/result/*.R` | run with `Rscript` | tables and figures |

```bash
# With Slurm
mkdir -p sim/output
sbatch slurm/sim/data_generation.sh          # stage 0; wait for it to finish
bash slurm/sim/sim_all.sh                    # stages 1–3, chained

# Without Slurm (runs the same launchers locally, one task at a time)
bash slurm/sim/sim_all_local.sh

Rscript sim/result/sim_result.R              # stage 4
```

To run another setting, export its variables once so that every stage sees
them (`sbatch` passes the environment on to each job). For example, the
all-Matérn, tridiagonal, colocated setting:

```bash
export SSGGP_FAMILY=matern SSGGP_PREC=tridiagonal SSGGP_MISS=coloc
sbatch slurm/sim/data_generation.sh          # wait for it to finish
bash slurm/sim/sim_all.sh
Rscript sim/result/sim_result.R
```

The penalized observed log-likelihood trace needs a likelihood evaluation at
iterations 50, 100, …, 500 of every selected fit. These are computed once with
`SSGGP_COMPUTE_DIAGNOSTIC_LOGLIK=true Rscript sim/result/sim_result.R` and
cached; without the cache, that one plot is skipped.

### Baselines

```bash
bash slurm/sim/lmc/submit_sim_lmc_all.sh          # LMC, all 16 settings
bash slurm/sim/lmc/submit_sim_lmc_all_matern.sh   # fixed-Matérn LMC, 4 tridiagonal all-Matérn settings
IOX_SETTING=random_neg_coloc sbatch slurm/sim/iox/sim_iox.sh   # IOX, one setting
SSGGP_FAMILY=matern SSGGP_MISS=coloc Rscript sim/result/sim_iox_result.R   # same setting
```

`submit_sim_lmc_all.sh` needs all 16 settings (both precision structures, signs,
and observation designs, each mixed-family and all-Matérn) with replicates
1–15. Run only one of the two LMC launchers on the same outputs.

## MERRA-2

See [`MERRA/README.md`](MERRA/README.md) for the raw input file and the command
sequence.
