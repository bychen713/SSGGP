#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
mode <- match.arg(args[1L], c("pilot", "fit", "continue"))
fold <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")) + 1L
root <- normalizePath(Sys.getenv("SLURM_SUBMIT_DIR", getwd()), mustWork = TRUE)
setwd(root)
# LMC for the colocated design: "pilot" creates starting values, "fit" runs
# maximum likelihood and prediction, and "continue" warm-restarts a saved fit.
if (!nzchar(Sys.getenv("MERRA_DESIGN"))) Sys.setenv(MERRA_DESIGN = "colocated")
source("MERRA/settings.R")
stopifnot(merra_design == "colocated")
source("MERRA/lmc/lmc_model.R")

out_dir <- Sys.getenv("MERRA_LMC_OUT_DIR", file.path(merra_output_root, "lmc"))
q_fixed <- as.integer(Sys.getenv("MERRA_LMC_Q", "10"))
stopifnot(length(q_fixed) == 1L, is.finite(q_fixed),
          q_fixed >= 1L, q_fixed <= merra_p)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pilot_file <- function(i) file.path(out_dir, sprintf("pilot_fold%d.rds", i))
fit_file <- function(i) file.path(out_dir, sprintf("fit_fold%d.rds", i))
checkpoint_file <- function(i) file.path(out_dir, sprintf("checkpoint_fold%d.rds", i))
save_atomic <- function(value, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(value, tmp)
  if (!file.rename(tmp, path)) stop("could not save ", path)
}

stopifnot(fold %in% seq_len(merra_n_folds))
data_path <- merra_fold_file(fold)
x <- readRDS(data_path)
stopifnot(identical(x$design, "colocated"),
          identical(dim(x$Y_mat), c(merra_n_train, merra_p)),
          identical(dim(x$Y_mat_test), c(merra_n_test, merra_p)),
          all(is.finite(x$Y_mat)), all(is.finite(x$Y_mat_test)))
Y <- as.matrix(x$Y_mat)
D <- fields::rdist(x$grid)
Dx <- fields::rdist(x$grid_test, x$grid)

if (mode == "pilot") {
  t0 <- Sys.time()
  bins <- lmc_bins(x$grid)
  vgm <- lmc_covariogram(Y, bins)
  nug_est <- lmc_shared_nugget(Y, vgm)
  nugget_hat <- nug_est$nugget
  S0 <- stats::cov(Y)
  ld <- lmc_loadings(S0, nugget_hat, q_max = q_fixed,
                     var_target = 0.80, q_min = q_fixed)
  q <- ld$q
  Zs <- Y %*% ld$V
  fits <- lapply(seq_len(q), function(k)
    lmc_score_fit("matern", Zs[, k], ld$d[k], nugget_hat, bins$D))
  if (any(!vapply(fits, `[[`, logical(1), "converged")))
    stop("a Matérn PCA-score initialization failed")
  par_raw <- lapply(fits, `[[`, "par")
  save_atomic(list(fold = fold, data_md5 = unname(tools::md5sum(data_path)),
    q = q, ld = ld, nugget_hat = nugget_hat,
    nugget_by_process = nug_est$per_process,
    nugget_method = "marginal_covariogram_wls_median", min_shared = nrow(Y),
    fam = rep("matern", q), par_raw = par_raw,
    start_method = "pca_score_gaussian_likelihood",
    selection_method = "Fixed Matern family; PCA-score Gaussian likelihood initialization",
    candidate_fits = fits,
    elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins"))),
    pilot_file(fold))
  cat("fold", fold, "q", q, "signal fraction", ld$frac,
      "pilot nugget", nugget_hat, "score-likelihood Matérn starts\n")
}

if (mode %in% c("fit", "continue")) {
  pilot <- readRDS(pilot_file(fold))
  stopifnot(pilot$fold == fold,
            pilot$data_md5 == unname(tools::md5sum(data_path)),
            pilot$q == q_fixed,
            identical(pilot$fam, rep("matern", q_fixed)),
            identical(pilot$nugget_method,
                      "marginal_covariogram_wls_median"),
            identical(pilot$start_method,
                      "pca_score_gaussian_likelihood"))
  bnd <- lmc_bounds(pilot$fam)
  nug0 <- max(pilot$nugget_hat, exp(tail(bnd$lower, 1L)))
  nll <- lmc_nll_coloc_factory(Y, D, pilot$ld$V, pilot$ld$d, pilot$fam)
  # Everything the fit depends on; a saved fit or checkpoint is reused only
  # when this matches exactly.
  identity <- list(design = "colocated", fold = fold, data_md5 = pilot$data_md5,
                   pilot_md5 = unname(tools::md5sum(pilot_file(fold))),
                   family = pilot$fam, lower = bnd$lower, upper = bnd$upper)
  previous <- NULL
  if (file.exists(fit_file(fold))) {
    previous <- readRDS(fit_file(fold))$optimization
    if (is.null(previous))
      stop("existing fit predates resumable fits; delete it to refit: ",
           fit_file(fold), call. = FALSE)
  }
  opt <- lmc_optimize(nll, lmc_pack(pilot$par_raw, nug0), bnd$lower, bnd$upper,
                      identity, checkpoint_file(fold),
                      resume = mode == "continue", previous = previous,
                      elapsed_before = pilot$elapsed_min,
                      save_every = as.numeric(Sys.getenv("SSGGP_LMC_CHECKPOINT_SECONDS", "300")))
  if (is.null(opt)) {
    cat("converged LMC fit already exists for fold", fold, "\n")
    quit(status = 0L)
  }
  fit <- opt$fit
  est <- lmc_unpack(fit$par, pilot$fam)
  pred <- lmc_predict_coloc(Y, D, Dx, pilot$ld$V, pilot$ld$d,
                            pilot$fam, est$par, est$nugget)
  variance <- t(vapply(pred$cov, diag, numeric(merra_p)))
  sum_variance <- vapply(pred$cov, sum, numeric(1))
  err <- pred$mean - x$Y_mat_test
  metrics <- data.frame(fold, method = "LMC_Matern", q = pilot$q,
    signal_fraction = pilot$ld$frac, mse = mean(err^2),
    norm = mean(sqrt(rowSums(err^2))),
    coverage = mean(abs(err) <= 1.96 * sqrt(variance)),
    pi_length = mean(3.92 * sqrt(variance)), sum_mse = mean(rowSums(err)^2),
    sum_coverage = mean(abs(rowSums(err)) <= 1.96 * sqrt(sum_variance)),
    sum_pi_length = mean(3.92 * sqrt(sum_variance)),
    pilot_nugget = pilot$nugget_hat, fitted_nugget = est$nugget,
    nll_start = opt$start_value, nll_final = fit$value,
    convergence = fit$convergence, optimization_passes = opt$passes,
    evaluations = opt$evaluations, elapsed_min = opt$elapsed_min)
  save_atomic(list(fold = fold, data_md5 = pilot$data_md5, pilot = pilot,
    fit = fit, optimization = opt[names(opt) != "fit"],
    est = est, metrics = metrics,
    prediction = pred$mean, prediction_cov = pred$cov), fit_file(fold))
  print(metrics)
}

