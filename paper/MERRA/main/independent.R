#!/usr/bin/env Rscript
fold <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")) + 1L
source("MERRA/settings.R")
suppressPackageStartupMessages(library(SSGGP))
dir.create(merra_output_dir, recursive = TRUE, showWarnings = FALSE)
x <- readRDS(merra_fold_file(fold))
stopifnot(identical(dim(x$Y_mat), c(merra_n_train, merra_p)),
          all(colSums(!is.na(x$Y_mat)) ==
                as.integer(round(merra_n_train * (1 - merra_missing_prop)))))
t0 <- Sys.time()
fit <- fit_independent(x$Y_mat, x$grid, par_init = c(NA, 3, 1, 0.01),
  lower = log(c(1 / 100, 1 / 10, 0.01, 0.001)),
  upper = log(c(1 / 0.001, 1 / 0.01, 3, 100)), threads = merra_threads)
pred_mean <- pred_var <- matrix(NA_real_, merra_n_test, merra_p)
for (j in seq_len(merra_p)) {
  keep <- !is.na(x$Y_mat[, j]); a <- fit$f_par[[j]][["a"]]; nu <- fit$f_par[[j]][["nu"]]
  Coo <- SSGGP:::matern_cov(fields::rdist(x$grid[keep, , drop = FALSE]), a = a,
    nu = nu, sigma = fit$variance[j], nugget = fit$nugget[j])
  Cno <- fit$variance[j] * fields::Matern(fields::rdist(x$grid_test,
    x$grid[keep, , drop = FALSE]), smoothness = nu, alpha = a)
  ch <- chol(Coo)
  pred_mean[, j] <- as.vector(Cno %*% backsolve(ch,
    backsolve(ch, x$Y_mat[keep, j], transpose = TRUE)))
  v <- backsolve(ch, t(Cno), transpose = TRUE)
  pred_var[, j] <- pmax(fit$variance[j] + fit$nugget[j] - colSums(v^2), 0)
}
err <- pred_mean - x$Y_mat_test
metrics <- data.frame(fold, method = "IND", rho = NA_real_, mse = mean(err^2),
  norm = mean(sqrt(rowSums(err^2))), coverage = mean(abs(err) <= 1.96 * sqrt(pred_var)),
  pi_length = mean(3.92 * sqrt(pred_var)), sum_mse = mean(rowSums(err)^2),
  sum_coverage = mean(abs(rowSums(err)) <= 1.96 * sqrt(rowSums(pred_var))),
  sum_pi_length = mean(3.92 * sqrt(rowSums(pred_var))),
  elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
saveRDS(list(fold = fold, data_md5 = unname(tools::md5sum(merra_fold_file(fold))),
  fit = fit, prediction = pred_mean, pred_var = pred_var, metrics = metrics),
  merra_ind_file(fold))
print(metrics)
