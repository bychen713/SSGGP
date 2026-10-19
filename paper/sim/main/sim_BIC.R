# Select the penalty by BIC and predict at the test locations with the selected fit.
#   sbatch slurm/sim/main/sim_BIC.sh         (one task per replicate)
# The first penalty (the dense fit) is excluded; see sim_dense.R.
# Writes BIC_<setting>_<sim>.RData.

rm(list = ls())
if (!file.exists("settings.R")) {  # work from sim/, where settings.R lives
  if (file.exists(file.path("sim", "settings.R"))) setwd("sim") else
  if (file.exists(file.path("..", "settings.R"))) setwd("..") else
    stop("run this script from paper/")
}
source("settings.R")
setup_compute()

j <- job_index(1L)                    # replicate index, 1..n_sim
stopifnot(j >= 1, j <= n_sim)

# Rows of grid_params for this replicate, one per penalty, excluding the dense fit.
rows <- which(grid_params$sim == j)
rows_penalised <- rows[grid_params$rho[rows] != rho_vals[1]]   # drop the dense one

ll_k_all <- n_par_all <- rep(NA_real_, length(rows_penalised))
for (i in seq_along(rows_penalised)) {
  f <- file.path(out_dir, paste0("ll_", rows_penalised[i], ".RData"))
  if (file.exists(f)) { load(f); ll_k_all[i] <- ll_k; n_par_all[i] <- n_par }
}

# BIC = -2 logLik + log(n) k, with k the number of selected edges and n the
# number of training locations.
BICs <- -2 * ll_k_all + log(n_train) * n_par_all
if (all(!is.finite(BICs))) stop("no usable likelihood for replicate ", j)

best     <- which.min(BICs)
best_row <- rows_penalised[best]
rho_best <- grid_params$rho[best_row]

cat(sprintf("replicate %d: selected rho = %g (BIC %.3f, %d edges)\n",
            j, rho_best, BICs[best], n_par_all[best]))
print(data.frame(rho = grid_params$rho[rows_penalised],
                 logLik = ll_k_all, edges = n_par_all, BIC = BICs))

load(file.path(out_dir, paste0("est", best_row, ".RData")))   # result_est

time_a <- Sys.time()
pred <- predict(result_est$fit, result_est$grid_test,
                L = L_pred, MC = MC_pred, g = proposal_pred,
                summarise = "blockdiag", keep_chain = FALSE,
                threads = threads, seed = j)
time_b <- Sys.time()
time_pred <- difftime(time_b, time_a, units = "mins")

# Predictive means in process-fastest order and the per-location covariance blocks.
pred_var_gibbs <- Matrix::bdiag(lapply(seq_len(dim(pred$cov_block)[1]),
                                       function(i) pred$cov_block[i, , ]))
result_final <- list(
  Y          = result_est$Y,
  grid       = result_est$grid,
  Y_test     = result_est$Y_test,
  Y_mat_test = result_est$Y_mat_test,
  grid_test  = result_est$grid_test,
  a_est      = result_est$a_est,
  nu_est     = result_est$nu_est,
  sigma_inv_est      = result_est$sigma_inv_est,
  sigma_inv_est_vals = result_est$sigma_inv_est_vals,
  D_est_vals = result_est$D_est_vals,
  rho        = rho_best,
  prediction_gibbs = as.vector(t(pred$mean)),
  pred_var_gibbs   = pred_var_gibbs,
  cov_block        = pred$cov_block,   # per-location P x P, for the sum
  var_sum          = var_sum(pred),
  BICs   = BICs,
  llMCEMs = ll_k_all,
  n_par  = n_par_all,
  ll_MCEM = ll_k_all[best])

save(time_approx, time_pred, result_final,
     file = file.path(out_dir, paste0("BIC_", setting_name, "_", j, ".RData")))
cat("wrote BIC_", setting_name, "_", j, ".RData\n", sep = "")
