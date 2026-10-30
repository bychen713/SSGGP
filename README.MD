# Sparse Spectral Graphical Gaussian Processes

This repository contains a research release of the **SSGGP** R
package and the scripts used for its simulation and MERRA-2 analyses. SSGGP
fits and predicts highly multivariate Gaussian processes using random Fourier
features, stochastic-approximation EM, and a graphical-lasso penalty for sparse
cross-process dependence.

## Repository layout

```text
.
├── pkg/             # SSGGP R package: R code, C++ kernels, docs
├── demo.R           # small end-to-end example (start here)
├── run_demo.sh      # runs demo.R inside the project environment
├── setup_env.sh     # dependency setup helper (renv)
├── renv.lock        # reproducible R dependency snapshot
└── paper/           # scripts reproducing the paper's simulations and MERRA-2 analysis
```

If you want to **use SSGGP**, you only need `pkg/` and `demo.R`. If you want to
**reproduce the paper**, see [`paper/README.md`](paper/README.md).

## Requirements and installation

- R 4.4.x
- a C/C++/Fortran compiler
- GNU GSL
- for the MERRA-2 preprocessing only: GDAL, GEOS, PROJ, udunits, and netCDF
  (or an R installation that already provides `sf` and `ncdf4`)

From the repository root:

```bash
bash setup_env.sh     # restores dependencies and installs SSGGP from pkg/
bash run_demo.sh      # runs demo.R (about a minute)
```

Add `--with-merra` to also set up the MERRA-2 preprocessing packages (`ncdf4`,
`sf`, `patchwork`, `maps`), and `--with-iox` to install the bundled `spiox` used
by the IOX baseline (`paper/iox/spiox`). After editing files under `pkg/`,
reinstall with `bash setup_env.sh` (or `R CMD INSTALL pkg` from the repository
root).

After `demo.R`, the vignette `pkg/vignettes/ssggp.Rmd` gives a longer
walkthrough: simulation, fitting, graph recovery, prediction, penalty
selection, and adding a spectral family.

## Reproducing the paper

The simulation study, the LMC and IOX baselines, and the MERRA-2 application
live under `paper/`. Run them from that directory. See
[`paper/README.md`](paper/README.md) for the pipeline and the scripts behind
each table and figure. Generated data, fitted objects, logs, and figures are
not versioned.

## Checking the package

```bash
R CMD build pkg
R CMD check SSGGP_*.tar.gz --no-manual
```

## License

The SSGGP package and the scripts in this repository are released under the
[MIT License](LICENSE). The bundled `paper/iox/spiox` is spiox (GPL >= 2) and
keeps its own license.
