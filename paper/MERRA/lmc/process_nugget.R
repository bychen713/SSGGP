#!/usr/bin/env Rscript
# LMC with process-specific nuggets, the LMC reported for MERRA-2. Runs after
# the shared-nugget LMC (colocated.R or misaligned.R, by MERRA_DESIGN), whose
# pilot and fit give the starting values:
#   pilot      per-process variogram nuggets (clamped inside the nugget box),
#              loadings from the covariance minus the per-process nugget
#              diagonal, and the shared-nugget MLE spatial parameters
#   fit        joint maximum likelihood over the factor parameters and the P
#              nuggets, then prediction; "continue" resumes an unfinished fit
mode <- match.arg(commandArgs(trailingOnly = TRUE)[1L], c("pilot", "fit", "continue"))
fold <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")) + 1L
root <- normalizePath(Sys.getenv("SLURM_SUBMIT_DIR", getwd()), mustWork = TRUE)
setwd(root)
source("MERRA/settings.R")
source("MERRA/lmc/lmc_model.R")
source("MERRA/lmc/lmc_process_nugget.R")
design <- merra_design
stopifnot(fold %in% seq_len(merra_n_folds))

shared_dir <- Sys.getenv("MERRA_LMC_OUT_DIR", file.path(merra_output_root, "lmc"))
out_dir <- file.path(merra_output_root, "lmc_process_nugget")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pilot_file <- function(i) file.path(out_dir, sprintf("pilot_fold%d.rds", i))
fit_file <- function(i) file.path(out_dir, sprintf("fit_fold%d.rds", i))
checkpoint_file <- function(i) file.path(out_dir, sprintf("checkpoint_fold%d.rds", i))
save_atomic <- function(value, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(value, tmp)
  if (!file.rename(tmp, path)) stop("could not save ", path)
}
q_fixed <- as.integer(Sys.getenv("MERRA_LMC_Q", "10"))
families <- rep("matern", q_fixed)

path <- merra_fold_file(fold)
x <- readRDS(path)
Y <- as.matrix(x$Y_mat)
stopifnot(identical(dim(Y), c(merra_n_train, merra_p)),
          identical(dim(x$Y_mat_test), c(merra_n_test, merra_p)),
          if (design == "colocated") !anyNA(Y) else all(colSums(!is.na(Y)) > 0),
          all(is.finite(x$Y_mat_test)))
D <- fields::rdist(x$grid)
Dx <- fields::rdist(x$grid_test, x$grid)
data_md5 <- unname(tools::md5sum(path))

if (mode == "pilot") {
  shared_pilot_file <- file.path(shared_dir, sprintf("pilot_fold%d.rds", fold))
  shared_fit_file <- file.path(shared_dir, sprintf("fit_fold%d.rds", fold))
  sp <- readRDS(shared_pilot_file)
  sf <- readRDS(shared_fit_file)
  stopifnot(identical(sp$data_md5, data_md5),
            identical(sf$data_md5, data_md5), length(sf$est$par) == q_fixed,
            identical(sf$pilot$fam, families))
  nug_raw <- sp$nugget_by_process
  stopifnot(length(nug_raw) == ncol(Y), all(is.finite(nug_raw)))
  # The marginal variogram pilot can put a nugget at zero. Clamp it strictly
  # inside the same [1e-4, 10] box used for the final likelihood optimization.
  nug0 <- pmin(pmax(nug_raw, 1e-3), 9.9)
  if (design == "colocated") {
    S0 <- stats::cov(Y)
  } else {
    O <- !is.na(Y)
    Y0 <- sweep(Y, 2L, colMeans(Y, na.rm = TRUE))
    Y0[!O] <- 0
    Nsh <- crossprod(O * 1)
    S0 <- crossprod(Y0) / Nsh
    S0[!is.finite(S0)] <- 0
    S0 <- (S0 + t(S0)) / 2
  }
  ld <- lmc_loadings(S0, nug0, q_max = q_fixed,
                     var_target = if (design == "colocated") 0.80 else 1,
                     q_min = q_fixed)
  theta0 <- lmc_pn_pack(sf$est$par, nug0)
  bounds <- lmc_pn_bounds(families, ncol(Y))
  start_clamped <- which(theta0 <= bounds$lower | theta0 >= bounds$upper)
  theta0 <- pmin(pmax(theta0, bounds$lower + 1e-6), bounds$upper - 1e-6)
  stopifnot(all(theta0 > bounds$lower), all(theta0 < bounds$upper))
  save_atomic(list(design = design, fold = fold, data_md5 = data_md5,
    shared_pilot_md5 = unname(tools::md5sum(shared_pilot_file)),
    shared_fit_md5 = unname(tools::md5sum(shared_fit_file)),
    q = q_fixed, fam = families, ld = ld, nugget_raw = nug_raw,
    nugget_start = nug0, spatial_start = sf$est$par, theta0 = theta0,
    start_clamped = start_clamped,
    start_method = "per-process marginal variogram nugget; shared-nugget MLE spatial parameters",
    loading_method = "PCA of covariance minus process-specific nugget diagonal"),
    pilot_file(fold))
  cat(design, "fold", fold, "pilot nugget range", range(nug0),
      "signal fraction", ld$frac, "clamped starts", length(start_clamped), "\n")
}

if (mode %in% c("fit", "continue")) {
  pilot <- readRDS(pilot_file(fold))
  stopifnot(identical(pilot$data_md5, data_md5), pilot$q == q_fixed,
            identical(pilot$fam, families))
  bounds <- lmc_pn_bounds(families, ncol(Y))
  nll <- lmc_pn_nll_factory(Y, D, pilot$ld$B, families)
  # Everything the fit depends on; a saved fit or checkpoint is reused only
  # when this matches exactly.
  identity <- list(design = design, fold = fold, data_md5 = data_md5,
                   pilot_md5 = unname(tools::md5sum(pilot_file(fold))),
                   family = families, lower = bounds$lower, upper = bounds$upper)
  previous <- NULL
  if (file.exists(fit_file(fold))) {
    previous <- readRDS(fit_file(fold))$optimization
    if (is.null(previous))
      stop("existing fit predates resumable fits; delete it to refit: ",
           fit_file(fold), call. = FALSE)
  }
  opt <- lmc_optimize(nll, pilot$theta0, bounds$lower, bounds$upper,
                      identity, checkpoint_file(fold),
                      resume = mode == "continue", previous = previous,
                      save_every = as.numeric(Sys.getenv("SSGGP_LMC_CHECKPOINT_SECONDS", "300")))
  if (is.null(opt)) {
    cat("converged process-nugget LMC fit already exists for", design, "fold", fold, "\n")
    quit(status = 0L)
  }
  fit <- opt$fit
  est <- lmc_pn_unpack(fit$par, families, ncol(Y))
  pred <- lmc_pn_predict(Y, D, Dx, pilot$ld$B, families, est$par, est$nugget)
  variance <- t(vapply(pred$cov, diag, numeric(ncol(Y))))
  sum_variance <- vapply(pred$cov, sum, numeric(1))
  stopifnot(all(is.finite(pred$mean)), all(is.finite(variance)),
            all(variance > 0), all(sum_variance > 0))
  err <- pred$mean - x$Y_mat_test
  metrics <- data.frame(fold, method = "LMC_Matern", q = q_fixed,
    signal_fraction = pilot$ld$frac, mse = mean(err^2),
    norm = mean(sqrt(rowSums(err^2))),
    coverage = mean(abs(err) <= 1.96 * sqrt(variance)),
    pi_length = mean(3.92 * sqrt(variance)),
    sum_mse = mean(rowSums(err)^2),
    sum_coverage = mean(abs(rowSums(err)) <= 1.96 * sqrt(sum_variance)),
    sum_pi_length = mean(3.92 * sqrt(sum_variance)),
    pilot_nugget_median = median(pilot$nugget_start),
    fitted_nugget_min = min(est$nugget),
    fitted_nugget_median = median(est$nugget),
    fitted_nugget_max = max(est$nugget),
    nll_start = opt$start_value, nll_final = fit$value,
    convergence = fit$convergence, optimization_passes = opt$passes,
    evaluations = opt$evaluations, elapsed_min = opt$elapsed_min)
  save_atomic(list(design = design, fold = fold, data_md5 = data_md5,
    pilot = pilot, fit = fit, optimization = opt[names(opt) != "fit"],
    est = est, metrics = metrics,
    prediction = pred$mean, prediction_cov = pred$cov,
    nugget_model = "process_specific_diagonal; jointly estimated by Gaussian likelihood"),
    fit_file(fold))
  print(metrics)
}
