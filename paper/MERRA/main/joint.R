#!/usr/bin/env Rscript
# Joint SSGGP fit, run in three stages (first argument):
#   estimate    one fit per (fold, penalty); array index = (rho_index - 1) * 3 + fold - 1
#   likelihood  marginal log-likelihood of each fit, for BIC
#   predict     dense (first penalty) and BIC-selected prediction for each fold
args <- commandArgs(trailingOnly = TRUE)
stage <- match.arg(args[1L], c("estimate", "likelihood", "predict"))
index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
source("MERRA/settings.R")
suppressPackageStartupMessages(library(SSGGP))
dir.create(merra_output_dir, recursive = TRUE, showWarnings = FALSE)
save_atomic <- function(value, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid()); saveRDS(value, tmp)
  if (!file.rename(tmp, path)) stop("could not write ", path)
}
key <- function(i) list(fold = i %% merra_n_folds + 1L,
                        rho_index = i %/% merra_n_folds + 1L)
# Fold data, its independent fit, and the Matern proposal derived from that fit.
inputs <- function(fold) {
  x <- readRDS(merra_fold_file(fold)); ind <- readRDS(merra_ind_file(fold))
  stopifnot(identical(ind$data_md5, unname(tools::md5sum(merra_fold_file(fold)))))
  aa <- vapply(ind$fit$f_par, `[[`, numeric(1), "a")
  nn <- vapply(ind$fit$f_par, `[[`, numeric(1), "nu")
  pp <- c(a = median(aa), nu = min(nn))
  stopifnot(all(is.finite(pp)), pp[["a"]] > 0, pp[["nu"]] > 0)
  list(x = x, ind = ind, proposal_par = pp,
       proposal = matern_proposal(a = pp[["a"]], nu = pp[["nu"]]))
}

if (stage == "estimate") {
  z <- key(index); stopifnot(z$rho_index %in% seq_along(merra_rho))
  output <- merra_est_file(z$fold, z$rho_index)
  if (file.exists(output)) {
    old <- tryCatch(readRDS(output), error = function(e) NULL)
    if (!is.null(old) && identical(old$data_md5, unname(tools::md5sum(merra_fold_file(z$fold)))) &&
        (isTRUE(old$fit$converged) || old$fit$iterations == merra_iter)) {
      cat("valid completed estimate exists:", output, "\n"); quit(status = 0L)
    }
  }
  inp <- inputs(z$fold)
  checkpoint <- if (identical(Sys.getenv("MERRA_CHECKPOINT", "0"), "1"))
    file.path(merra_output_dir, "checkpoints",
      sprintf("est_fold%d_rho%d.rds", z$fold, z$rho_index)) else NULL
  t0 <- Sys.time()
  # Default starting precision, as in the simulation: the inverse training
  # sample covariance when colocated, the identity when misaligned.
  fit <- SSGGP::fit(Y = inp$x$Y_mat, locs = inp$x$grid, f_fun = "matern",
    f_par = inp$ind$fit$f_par, g = inp$proposal, rho = merra_rho[z$rho_index],
    nugget = "by_process", nugget_init = mean(inp$ind$fit$nugget),
    L = merra_L, MC = merra_MC, burnin = merra_burnin, iter = merra_iter,
    lag = merra_lag, tol = merra_tol, diag_pen = FALSE,
    lower = log(c(a = 0.1, nu = 0.01)), upper = log(c(a = 100, nu = 3)),
    threads = merra_threads, seed = z$fold, verbose = TRUE,
    checkpoint_file = checkpoint, checkpoint_every = 10L)
  save_atomic(list(fold = z$fold, rho_index = z$rho_index,
    rho = merra_rho[z$rho_index], data_md5 = inp$ind$data_md5,
    proposal_par = inp$proposal_par, fit = fit,
    elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins"))), output)
}

if (stage == "likelihood") {
  z <- key(index); inp <- inputs(z$fold); est <- readRDS(merra_est_file(z$fold, z$rho_index))
  ll <- logLik(est$fit, L = merra_L_ll, seed = 999L + z$fold,
               g = inp$proposal, threads = merra_threads)
  save_atomic(list(fold = z$fold, rho_index = z$rho_index, rho = est$rho,
    logLik = as.numeric(ll), edges = as.integer(attr(ll, "df")),
    # BIC = -2 logLik + log(n) k, n = training locations, k = selected edges
    bic = -2 * as.numeric(ll) + log(merra_n_train) * attr(ll, "df")),
    file.path(merra_output_dir, sprintf("ll_fold%d_rho%d.rds", z$fold, z$rho_index)))
}

if (stage == "predict") {
  if (merra_pred_single_method) {
    stopifnot(index >= 0L, index < 2L * merra_n_folds)
    fold <- index %% merra_n_folds + 1L
    selected_method <- c("dense", "BIC")[index %/% merra_n_folds + 1L]
  } else {
    stopifnot(index >= 0L, index < merra_n_folds)
    fold <- index + 1L
    selected_method <- c("dense", "BIC")
  }
  inp <- inputs(fold)
  ll_files <- file.path(merra_output_dir,
    sprintf("ll_fold%d_rho%d.rds", fold, seq_along(merra_rho)))
  stopifnot(all(file.exists(ll_files)))
  scores <- lapply(ll_files, readRDS); bic <- vapply(scores, `[[`, numeric(1), "bic")
  chosen <- c(dense = 1L, BIC = which.min(replace(bic, 1L, Inf)))
  for (method in selected_method) {
    ri <- chosen[[method]]; est <- readRDS(merra_est_file(fold, ri)); t0 <- Sys.time()
    pred <- predict(est$fit, inp$x$grid_test, L = merra_L_pred, MC = merra_MC_pred,
      burnin = merra_MC_pred %/% 2L, g = inp$proposal, summarise = "blockdiag",
      keep_chain = FALSE, threads = merra_threads, seed = fold)
    variance <- t(vapply(seq_len(merra_n_test),
      function(i) diag(pred$cov_block[i, , ]), numeric(merra_p)))
    sum_variance <- vapply(seq_len(merra_n_test),
      function(i) sum(pred$cov_block[i, , ]), numeric(1))
    err <- pred$mean - inp$x$Y_mat_test
    metrics <- data.frame(fold, method, rho = est$rho, logLik = scores[[ri]]$logLik,
      edges = scores[[ri]]$edges, mse = mean(err^2),
      norm = mean(sqrt(rowSums(err^2))),
      coverage = mean(abs(err) <= 1.96 * sqrt(variance)),
      pi_length = mean(3.92 * sqrt(variance)), sum_mse = mean(rowSums(err)^2),
      sum_coverage = mean(abs(rowSums(err)) <= 1.96 * sqrt(sum_variance)),
      sum_pi_length = mean(3.92 * sqrt(sum_variance)),
      elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    save_atomic(list(metrics = metrics, prediction = pred$mean,
      cov_block = pred$cov_block, bic = bic, proposal_par = inp$proposal_par,
      MC_pred = merra_MC_pred, L_pred = merra_L_pred,
      data_md5 = inp$ind$data_md5),
      file.path(merra_output_dir,
        sprintf("pred_%s_fold%d%s.rds", method, fold, merra_pred_suffix)))
    print(metrics)
  }
}
