# IOX results for one setting: prediction summaries, parameter
# estimation error (Matern processes), and graph recovery.
#   SSGGP_PREC=random SSGGP_SIGN=neg SSGGP_MISS=coloc Rscript sim/result/sim_iox_result.R
# The setting must match the IOX_SETTING used for sim/iox/sim_iox.R.

rm(list=ls())
library(magrittr)
setwd(this.path::this.dir())
source("../settings.R")

folder_results <- out_dir
n_folds <- n_sim
alpha   <- 0.05

avg_norm_resid_folds <- numeric(n_folds)
rmse_sum_folds       <- numeric(n_folds)
avg_coverage_ind_folds  <- numeric(n_folds)
coverage_sum_folds   <- numeric(n_folds)
avg_intlen_ind_folds <- numeric(n_folds)
intlen_sum_folds     <- numeric(n_folds)

theta_mean_folds <- list()   # 4 x q matrix of posterior means
theta_sd_folds   <- list()   # 4 x q matrix of posterior sds
Sigma_mean_folds <- list()   # q x q matrix of posterior means
Sigma_sd_folds   <- list()   # q x q matrix of posterior sds
var_id_folds     <- list()   # length-q identified marginal variance
nugget_mean_folds <- list()
coverage_ind_proc <- list()  # length-q per-process coverage, one per fold
intlen_ind_proc   <- list()  # length-q per-process interval length, one per fold

time <- c()
for (fold in 1:n_folds) {
  # Outputs of sim/iox/sim_iox.R: predictive draws, and the fit's chains with
  # theta in (phi, sigmasq, nu, nugget) rows; spiox() returns (phi, nu, alpha)
  # with the marginal variance in Sigma and the nugget in Ddiag.
  iox_dir <- file.path(folder_results, "iox")
  pred  <- readRDS(file.path(iox_dir, paste0("pred_fold", fold, ".rds")))
  saved <- readRDS(file.path(iox_dir, paste0("fit_fold", fold, ".rds")))
  data_env <- new.env()
  load(file.path(data_dir, paste0(fold, ".RData")), envir = data_env)
  Y_out <- as.matrix(data_env$Y_mat_test)
  estim_time <- saved$initialization_elapsed + saved$elapsed
  spiox_predicts <- list(Y = pred$draws)
  theta_chain <- array(1, c(4L, dim(saved$fit$Theta)[2:3]))
  theta_chain[1L, , ] <- saved$fit$Theta[1L, , ]
  theta_chain[3L, , ] <- saved$fit$Theta[2L, , ]
  theta_chain[4L, , ] <- saved$fit$Ddiag
  spiox_out <- list(theta = theta_chain, Sigma = saved$fit$Sigma)
  rm(saved, pred, data_env)

  time <- c(time,estim_time["elapsed"]/60)
  Y_samples <- spiox_predicts$Y
  n_samples <- dim(Y_samples)[3]
  
  Yhat_mean <- apply(Y_samples, 1:2, mean)
  
  # Mean residual norm over test locations
  resid_mat <- Y_out - Yhat_mean                        # n_test x q
  norms     <- apply(resid_mat, 1, function(r) sqrt(sum(r^2)))  # n_test
  avg_norm_resid_folds[fold] <- mean(norms)
  
  # RMSE of sum of processes
  Y_sum_true    <- rowSums(Y_out)      # n_test
  Yhat_sum_mean <- rowSums(Yhat_mean)  # n_test
  rmse_sum_folds[fold] <- sqrt(mean((Y_sum_true - Yhat_sum_mean)^2))
  
  # 95% intervals from the 2.5% and 97.5% posterior predictive quantiles
  # (type 7), for each process and for the sum across processes within each draw.
  lo_ind <- apply(Y_samples, 1:2, quantile, probs = 0.025, type = 7, names = FALSE)
  hi_ind <- apply(Y_samples, 1:2, quantile, probs = 0.975, type = 7, names = FALSE)
  covered_ind <- (Y_out >= lo_ind) & (Y_out <= hi_ind)  # n_test x q
  avg_coverage_ind_folds[fold] <- mean(covered_ind)
  avg_intlen_ind_folds[fold] <- mean(hi_ind - lo_ind)
  coverage_ind_proc[[fold]] <- colMeans(covered_ind)          # length q
  intlen_ind_proc[[fold]] <- colMeans(hi_ind - lo_ind)

  Y_sum_samples <- apply(Y_samples, c(1, 3), sum)             # n_test x draws
  lo_sum <- apply(Y_sum_samples, 1, quantile, probs = 0.025, type = 7, names = FALSE)
  hi_sum <- apply(Y_sum_samples, 1, quantile, probs = 0.975, type = 7, names = FALSE)
  covered_sum <- (Y_sum_true >= lo_sum) & (Y_sum_true <= hi_sum)
  coverage_sum_folds[fold] <- mean(covered_sum)
  intlen_sum_folds[fold] <- mean(hi_sum - lo_sum)

  # Posterior summaries of theta (rows phi, sigmasq, nu, nugget) and Sigma
  theta_mcmc <- spiox_out$theta

  # Last 5000 draws, the same for theta and Sigma (their product is used below).
  last5000 <- tail(seq_len(dim(theta_mcmc)[3]), 5000)
  theta_mcmc_last <- theta_mcmc[, , last5000, drop = FALSE]   # 4 x q x 5000
  theta_mean_folds[[fold]] <- apply(theta_mcmc_last, 1:2, mean)   # 4 x q
  theta_sd_folds[[fold]]   <- apply(theta_mcmc_last, 1:2, sd)     # 4 x q

  Sigma_mcmc <- spiox_out$Sigma[, , last5000, drop = FALSE]   # q x q x 5000
  Sigma_mean_folds[[fold]] <- apply(Sigma_mcmc, 1:2, mean)  # q x q
  Sigma_sd_folds[[fold]]   <- apply(Sigma_mcmc, 1:2, sd)    # q x q

  # IOX marginal variance is theta[2, j] * Sigma[j, j].
  Sigma_diag <- apply(Sigma_mcmc, 3, diag)                  # q x 5000
  var_id_folds[[fold]] <- rowMeans(theta_mcmc_last[2, , ] * Sigma_diag)
  nugget_mean_folds[[fold]] <- theta_mean_folds[[fold]][4, ]

  cat(sprintf("Fold %d done.\n", fold))
}
mean(time)

# Prediction summaries: mean and SD across replicates of each replicate's value
fold_metrics <- list(
  "Residual norm"                = avg_norm_resid_folds,
  "Individual 95% coverage"      = avg_coverage_ind_folds,
  "Individual 95% int. length"   = avg_intlen_ind_folds,
  "Sum RMSE"                     = rmse_sum_folds,
  "Sum 95% coverage"             = coverage_sum_folds,
  "Sum 95% int. length"          = intlen_sum_folds
)

results <- data.frame(
  metric = names(fold_metrics),
  mean   = vapply(fold_metrics, mean, numeric(1)),
  sd     = vapply(fold_metrics, stats::sd, numeric(1)),
  stringsAsFactors = FALSE
)
rownames(results) <- NULL

print(results, digits = 6)

cat("\n-- mean (sd) over", n_folds, "folds --\n")
iox_tab <- data.frame(metric = results$metric,
                      IOX = sprintf("%.3f (%.3f)", results$mean, results$sd))
print(iox_tab, row.names = FALSE)
table_dir <- here("sim", "output", "tables", setting_name)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(iox_tab, file.path(table_dir, "prediction_iox.csv"), row.names = FALSE)

# Parameter summaries averaged across folds
theta_param_names <- c("phi", "sigmasq", "nu", "nugget")

theta_mean_avg <- Reduce("+", theta_mean_folds) / n_folds  # 4 x q

rownames(theta_mean_avg) <- theta_param_names

cat("\n--- Theta: posterior mean (averaged over folds) ---\n")
print(round(rowMeans(theta_mean_avg), 2))

cat("\n--- Theta: posterior SD (averaged over folds) ---\n")
row_sds <- apply(combined_df <- do.call(cbind, theta_mean_folds), 1, sd, na.rm = TRUE)
print(round(row_sds,2))

# Estimation error for the Matern processes; IOX fits a Matern to every process,
# so the others have no comparable true parameters.

# Parameter bounds of the spiox sampler; `outside` counts true values beyond them.
spiox_bounds <- list(phi = c(0.3, 200), var = NULL,
                     nu = c(1e-5, 2.9), nugget = NULL)

family_true <- unlist(f_fun_true)
matern_idx  <- which(family_true == "matern")
Sigma_true  <- solve(true_precision())

cat("\n\n================ estimation error: Matern processes only ",
    "================\n", sep = "")
cat(length(matern_idx), " of ", P, " processes were generated Matern; the ",
    "other ", P - length(matern_idx), " (",
    paste(sprintf("%d %s", table(family_true[family_true != "matern"]),
                  names(table(family_true[family_true != "matern"]))),
          collapse = ", "),
    ") are fitted Matern too but have no true phi/nu to score against and are ",
    "left out.\n", sep = "")

if (length(matern_idx)) {
  truth_by_par <- list(
    phi    = vapply(f_par_true[matern_idx], function(p) p[["a"]],  numeric(1)),
    var    = diag(Sigma_true)[matern_idx],
    nu     = vapply(f_par_true[matern_idx], function(p) p[["nu"]], numeric(1)),
    nugget = rep(nugget_true, length(matern_idx))
  )

  est_long <- do.call(rbind, lapply(seq_len(n_folds), function(fold) {
    est_by_par <- list(phi    = theta_mean_folds[[fold]][1, matern_idx],
                       var    = var_id_folds[[fold]][matern_idx],
                       nu     = theta_mean_folds[[fold]][3, matern_idx],
                       nugget = nugget_mean_folds[[fold]][matern_idx])
    do.call(rbind, lapply(names(est_by_par), function(par) {
      data.frame(fold = fold, process = matern_idx, par = par,
                 truth = truth_by_par[[par]], estimate = est_by_par[[par]],
                 stringsAsFactors = FALSE)
    }))
  }))
  est_long$error <- est_long$estimate - est_long$truth

  # Pooled over Matern processes and replicates
  err_pooled <- do.call(rbind, lapply(names(truth_by_par), function(par) {
    x <- est_long[est_long$par == par, , drop = FALSE]
    b <- spiox_bounds[[par]]
    data.frame(par       = par,
               n         = nrow(x),
               truth_rng = sprintf("%.2f-%.2f", min(x$truth), max(x$truth)),
               mean_est  = mean(x$estimate),
               bias      = mean(x$error),
               rmse      = sqrt(mean(x$error^2)),
               outside   = if (is.null(b)) NA_integer_ else
                 sum(truth_by_par[[par]] < b[1] | truth_by_par[[par]] > b[2]),
               stringsAsFactors = FALSE)
  }))

  cat("\n-- pooled over ", length(matern_idx), " processes x ", n_folds,
      " folds --\n", sep = "")
  print(transform(err_pooled,
                  mean_est = round(mean_est, 4),
                  bias     = round(bias, 4),
                  rmse     = round(rmse, 4)),
        row.names = FALSE, na.print = "-")

  if (any(err_pooled$outside > 0, na.rm = TRUE)) {
    hit <- err_pooled[which(err_pooled$outside > 0), ]
    cat("\nNOTE: ", paste(sprintf("%d true %s value(s) lie outside spiox's [%s] ",
                                  hit$outside, hit$par,
                                  vapply(hit$par, function(p)
                                    paste(spiox_bounds[[p]], collapse = ", "),
                                    character(1))),
                          collapse = "and "),
        "sampling bound, so the posterior cannot reach them. The bias on that\n",
        "row is partly the bound, not the fit.\n", sep = "")
  }

  # One row per Matern process; `sd` is the spread of the posterior mean across replicates.
  err_by_process <- do.call(rbind, lapply(names(truth_by_par), function(par) {
    x <- est_long[est_long$par == par, , drop = FALSE]
    parts <- split(x, x$process)
    do.call(rbind, lapply(parts, function(p) {
      data.frame(process  = p$process[1], par = par, truth = p$truth[1],
                 mean_est = mean(p$estimate), sd = stats::sd(p$estimate),
                 bias     = mean(p$error),
                 rmse     = sqrt(mean(p$error^2)),
                 stringsAsFactors = FALSE)
    }))
  }))
  err_by_process <- err_by_process[order(err_by_process$par,
                                         err_by_process$process), ]
  rownames(err_by_process) <- NULL

  print_by_process <- TRUE
  if (print_by_process) {
    cat("\n-- per Matern process, mean over folds --\n")
    old_width <- options(width = 200)
    num <- vapply(err_by_process, is.numeric, logical(1))
    show <- err_by_process
    show[num] <- lapply(show[num], function(z) round(z, 4))
    print(show, row.names = FALSE)
    options(old_width)
  }

} else {
  message("No true Matern processes: parameter-estimation errors are not applicable.")
}

Sigma_mean_avg <- Reduce("+", Sigma_mean_folds) / n_folds  # q x q
Sigma_sd_avg   <- Reduce("+", Sigma_sd_folds)   / n_folds  # q x q
Prec_mean_avg <- Reduce("+", lapply(Sigma_mean_folds, solve)) / n_folds  # q x q


# Graph-recovery metrics relative to the generating precision matrix.
sparse <- c()
Sigma_inv <- true_precision()

off_diag_mask <- row(Sigma_inv) != col(Sigma_inv)
true_pos_mask <- (Sigma_inv != 0) & off_diag_mask
true_neg_mask <- (Sigma_inv == 0) & off_diag_mask

tol <- 0.1

n_true_pos    <- sum(true_pos_mask)
n_true_neg    <- sum(true_neg_mask)

for(fold in 1:n_folds){
  est_mat     <- solve(Sigma_mean_folds[[fold]])
  est_nonzero <- abs(est_mat) > tol
  
  TP <- sum(est_nonzero & true_pos_mask)
  FP <- sum(est_nonzero & true_neg_mask)
  FN <- sum(!est_nonzero & true_pos_mask)
  
  precision = TP/(TP+FP)
  recall = TP/(TP+FN)
  
  FPR <- if (n_true_neg > 0) FP / n_true_neg else NA_real_
  TPR <- if (n_true_pos > 0) 1 - (FN / n_true_pos) else NA_real_
  FDR <- if ((FP + TP) > 0) FP / (FP + TP) else 0
  F1 <- 2*TP/(2*TP+FP+FN)
  F2 <- 5*precision*recall / (4*precision+recall)
  
  sparse <- rbind(sparse, c(TPR,FDR,FPR))
}
colnames(sparse) <- c('TPR','FDR','FPR')
print(round(colMeans(sparse),2))

