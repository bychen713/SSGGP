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

`spiox/` is spiox 0.8.1 (M. Peruzzi, GPL >= 2) from
<https://github.com/mkln/spiox>, commit `5357bb7` (2026-05-27), with the Matérn smoothness upper bound raised
from 2.1 to 2.9 (`spiox.patch`), so that the smoothness of 2.5 used in the
simulations can be reached. A patch lets `spiox.h` accept `X` as a single column of zeros.
