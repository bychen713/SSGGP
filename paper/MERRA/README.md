# MERRA-2 analysis

## Setup

Download `MERRA2_400.inst3_3d_aer_Nv.20220715.nc4` from NASA GES DISC and place
it in `MERRA/raw/`. Then, from the repository root, run

```bash
bash setup_env.sh --with-merra --with-iox
```

`--with-merra` adds the preprocessing and figure packages (`ncdf4` and `sf` need
GDAL, GEOS, PROJ, udunits, and netCDF unless your R already provides them), and
`--with-iox` installs the bundled `spiox` for the IOX baseline.

Please cite the MERRA-2 data as: Global Modeling and Assimilation Office [GMAO]
(2015). MERRA-2 inst3_3d_aer_Nv: 3d, 3-Hourly, Instantaneous, Model-Level,
Assimilation, Aerosol Mixing Ratio V5.12.4. Goddard Earth Sciences Data and
Information Services Center (GES DISC), Greenbelt, MD, USA.
[doi:10.5067/LTVB4GPCOTK2](https://doi.org/10.5067/LTVB4GPCOTK2).


## Layout

```text
MERRA/
├── settings.R              all configuration
├── preprocess_and_eda.R    raw NetCDF -> data/preprocessed_P40.rds
├── prepare_folds.R         -> data/{misaligned,colocated}/fold[1-3].rds
├── main/                   SSGGP: independent.R, joint.R
├── lmc/                    LMC: misaligned.R, colocated.R
├── iox/                    IOX: run.R
└── result/                 summarize.R (tables and figures)
```

`raw/`, `data/`, and `output/` are generated and not versioned. `output/` holds
fits and predictions in `output/<design>/<method>/`, the tables, `figure/`, and
the Slurm logs. IOX fits can reach gigabytes per fold.

## Running

From `paper/`:

```bash
Rscript MERRA/preprocess_and_eda.R
Rscript MERRA/prepare_folds.R
bash slurm/MERRA/submit_all.sh          # all fits (Slurm)
Rscript MERRA/result/summarize.R        # after the jobs finish
```

- **Folds:** three train/test splits of 500 training and 200 test locations,
  with disjoint test sets. The misaligned design masks 40% of each process's
  training values.
- **`submit_all.sh`:** for both designs, submits SSGGP (independent fits, 3
  folds × 15 penalties, likelihood selection, dense and BIC prediction), LMC,
  and IOX. It refuses to overwrite existing results, so clear `output/` before
  a new run.

## Outputs

`summarize.R` writes to `output/`:

- `prediction_by_fold.csv` and `prediction_summary.csv` (Table 4)
- `estimation_by_fold.csv` and `estimation_median.csv`: fitted parameters by fold
  and medians across folds for IND, dense and BIC SSGGP, LMC, and IOX
- `estimation_summary.csv` (Table S8)
- `figure/correlation.png`: empirical residual correlation
- `figure/precision_colocated_{sample,bic,dense}.png` (Figure 4): fold-median
  colocated precision matrices
- `figure/data_maps.png` (Figure 3): raw and processed values for processes 1, 20, and 40
