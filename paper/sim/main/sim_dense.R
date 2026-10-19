# Dense baseline: predict with the fit at the first penalty (rho = 0), without
# selection, so every edge is kept.
#   sbatch slurm/sim/main/sim_dense.sh       (one task per replicate)
# Writes dense_<setting>_<sim>.RData with the same fields as the BIC stage.

rm(list = ls())
if (!file.exists("settings.R")) {  # work from sim/, where settings.R lives
  if (file.exists(file.path("sim", "settings.R"))) setwd("sim") else
  if (file.exists(file.path("..", "settings.R"))) setwd("..") else
    stop("run this script from paper/")
}
source("settings.R")
setup_compute()

j <- job_index(1L)
stopifnot(j >= 1, j <= n_sim)

row <- which(grid_params$sim == j & grid_params$rho == rho_vals[1])[1]
est_file <- file.path(out_dir, paste0("est", row, ".RData"))
if (!file.exists(est_file)) stop("no dense fit at ", est_file)

load(est_file)                        # result_est, time_approx
cat(sprintf("replicate %d: dense baseline, rho = %g\n", j, result_est$rho))

# Likelihood of the dense fit, if sim_ll.R computed it.
ll_file <- file.path(out_dir, paste0("ll_", row, ".RData"))
ll_MCEM <- if (file.exists(ll_file)) { load(ll_file); ll_k } else NA_real_

time_a <- Sys.time()
pred <- predict(result_est$fit, result_est$grid_test,
                L = L_pred, MC = MC_pred, g = proposal_pred,
                summarise = "blockdiag", keep_chain = FALSE,
                threads = threads, seed = j)
time_b <- Sys.time()
time_pred <- difftime(time_b, time_a, units = "mins")

# Block-diagonal predictive covariance in process-fastest order, as in sim_BIC.R.
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
  rho        = result_est$rho,
  ll_MCEM    = ll_MCEM,
  prediction_gibbs = as.vector(t(pred$mean)),
  pred_var_gibbs   = pred_var_gibbs,
  cov_block        = pred$cov_block,
  var_sum          = var_sum(pred))

save(time_approx, time_pred, result_final,
     file = file.path(out_dir, paste0("dense_", setting_name, "_", j, ".RData")))
cat("wrote dense_", setting_name, "_", j, ".RData\n", sep = "")
