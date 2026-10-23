# IOX baseline

`sim/iox/sim_iox.R` and `MERRA/iox/run.R` fit spiox's latent IOX model with
Matérn marginals and a fixed zero mean, started from `spiox::autostart()`
(`latent_autostart.R`), sampled with the single-site latent sampler
(`sampling = 2`), and estimating each process's Matérn range and smoothness.
Each runs in two stages, `IOX_STAGE=estimate` then `IOX_STAGE=predict`, and
writes `init_`, `fit_` and `pred_fold<k>.rds`. Prediction uses the final
quarter of the iterations (2,500 of 10,000). A fit file holds the full latent
chain, about 2 GB per simulation replicate.

## The bundled spiox

`spiox/` is based on upstream spiox 0.8.1 (M. Peruzzi, GPL >= 2) from
<https://github.com/mkln/spiox>, commit `5357bb7` (2026-05-27), with these
changes (`spiox.patch`):

1. The Matérn smoothness upper bound is raised from 2.1 to 2.9, so that the
   smoothness of 2.5 used in the simulations can be reached.
2. The intercept check in `spiox.h` loops over the p columns of `X` rather
   than the q outcomes, so that `X` can be a single column of zeros.
3. The response model's missing-data update (`SpIOX::sample_Y_misaligned`) no
   longer clips imputed values to ±3 residual SDs, and updates the residuals
   after each row. Only the response model calls this update; the latent model
   used here does not, so the results are unaffected.
4. The version is set to 0.8.1.9001, to tell the bundled copy apart from
   upstream 0.8.1.

Install it with `bash setup_env.sh --with-iox` from the repository root. Both
IOX runners stop unless exactly version 0.8.1.9001 is installed.
