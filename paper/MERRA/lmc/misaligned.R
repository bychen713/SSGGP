#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
mode <- match.arg(args[1L], c("pilot", "fit", "continue"))
fold <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0")) + 1L
root <- normalizePath(Sys.getenv("SLURM_SUBMIT_DIR", getwd()), mustWork = TRUE)
setwd(root)
# LMC for the misaligned design: "pilot" creates starting values, "fit" runs
# maximum likelihood and prediction, and "continue" warm-restarts a saved fit.
if (!nzchar(Sys.getenv("MERRA_DESIGN"))) Sys.setenv(MERRA_DESIGN = "misaligned")
source("MERRA/settings.R")
stopifnot(merra_design == "misaligned")
source("MERRA/lmc/lmc_model.R")

out_dir <- Sys.getenv("MERRA_LMC_OUT_DIR", file.path(merra_output_root, "lmc"))
skip_wls <- identical(tolower(Sys.getenv("MERRA_LMC_SKIP_WLS", "true")), "true")
q_fixed <- as.integer(Sys.getenv("MERRA_LMC_Q", "10"))
nu_upper <- as.numeric(Sys.getenv("MERRA_LMC_NU_UPPER", "Inf"))
stopifnot(length(q_fixed) == 1L, is.finite(q_fixed), q_fixed >= 1L, q_fixed <= merra_p)
stopifnot(length(nu_upper) == 1L, !is.na(nu_upper), nu_upper > 0)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pilot_file <- function(i) file.path(out_dir, sprintf("pilot_fold%d.rds", i))
fit_file <- function(i) file.path(out_dir, sprintf("fit_fold%d.rds", i))
checkpoint_file <- function(i) file.path(out_dir, sprintf("checkpoint_fold%d.rds", i))
save_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid()); saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("could not save ", path)
}

stopifnot(fold %in% seq_len(merra_n_folds))
data_path <- merra_fold_file(fold)
x <- readRDS(data_path)
expected_observed <- as.integer(round(merra_n_train * (1 - merra_missing_prop)))
stopifnot(identical(dim(x$Y_mat), c(merra_n_train, merra_p)),
          all(colSums(!is.na(x$Y_mat)) == expected_observed),
          identical(dim(x$Y_mat_test), c(merra_n_test, merra_p)))
Y <- as.matrix(x$Y_mat)
D <- fields::rdist(x$grid)
Dx <- fields::rdist(x$grid_test, x$grid)

if (mode == "pilot") {
  t0 <- Sys.time()
  bins <- lmc_bins(x$grid)
  vgm <- lmc_covariogram(Y, bins)
  nug_est <- lmc_shared_nugget(Y, vgm)
  # lag-0 covariance from the locations where both processes are observed
  O <- !is.na(Y)
  Y0 <- sweep(Y, 2L, colMeans(Y, na.rm = TRUE)); Y0[!O] <- 0
  Nsh <- crossprod(O * 1)
  S0 <- crossprod(Y0) / Nsh; S0[!is.finite(S0)] <- 0
  S0 <- (S0 + t(S0)) / 2
  ld <- lmc_loadings(S0, nug_est$nugget, q_max = q_fixed,
                     var_target = 1, q_min = q_fixed)
  B <- ld$B; q <- ld$q; P <- ncol(Y)
  X <- vapply(seq_len(q), function(k) as.vector(tcrossprod(B[, k])), numeric(P * P))
  rho_hat <- matrix(NA_real_, vgm$n_bins, q)
  for (l in seq_len(vgm$n_bins)) {
    w <- as.vector(vgm$N[[l]]); y <- as.vector(vgm$C[[l]]); Xw <- X * w
    rho_hat[l, ] <- tryCatch(solve(crossprod(Xw, X), crossprod(Xw, y)),
                             error = function(e) rep(NA_real_, q))
  }
  rho_hat[] <- pmin(pmax(rho_hat, -1), 1)
  wl <- vapply(vgm$N, sum, numeric(1))
  if (skip_wls) {
    fits <- NULL
    par_raw <- replicate(q, LMC_STARTS$matern[[1L]], simplify = FALSE)
  } else {
    fits <- lapply(seq_len(q), function(k) {
      ok <- is.finite(rho_hat[, k]) & is.finite(wl) & wl > 0
      lmc_wls("matern", rho_hat[ok, k], wl[ok], vgm$h[ok],
              sill_max = 1, fix_sill = TRUE)
    })
    if (any(!vapply(fits, `[[`, logical(1), "converged")))
      stop("a Matérn latent-curve pilot fit failed")
    par_raw <- lapply(fits, `[[`, "par")
  }
  save_atomic(list(fold = fold, data_md5 = unname(tools::md5sum(data_path)),
    q = q, ld = ld, nugget_hat = nug_est$nugget,
    nugget_by_process = nug_est$per_process, min_shared = min(Nsh),
    fam = rep("matern", q), par_raw = par_raw, skip_wls = skip_wls,
    rho_hat = rho_hat, pilot_fits = fits,
    elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins"))),
    pilot_file(fold))
  cat("fold", fold, "q", q, "signal fraction", ld$frac, "skip WLS", skip_wls,
      "pilot nugget", nug_est$nugget, "minimum shared", min(Nsh), "\n")
}

if (mode %in% c("fit", "continue")) {
  pilot <- readRDS(pilot_file(fold)); data_path <- merra_fold_file(fold)
  stopifnot(pilot$fold == fold, pilot$data_md5 == unname(tools::md5sum(data_path)),
            pilot$q == q_fixed)
  bnd <- lmc_bounds(pilot$fam)
  if (is.finite(nu_upper)) {
    nu_position <- 2L * seq_len(pilot$q)
    bnd$upper[nu_position] <- pmin(bnd$upper[nu_position], log(nu_upper))
  }
  nug0 <- max(pilot$nugget_hat, exp(tail(bnd$lower, 1L)))
  nll <- lmc_nll_miss_factory(Y, D, pilot$ld$B, pilot$fam)
  # Everything the fit depends on; a saved fit or checkpoint is reused only
  # when this matches exactly.
  identity <- list(design = "misaligned", fold = fold, data_md5 = pilot$data_md5,
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
  pred <- lmc_predict_miss(Y, D, Dx, pilot$ld$B, pilot$fam, est$par, est$nugget)
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
    nu_upper = nu_upper, fit = fit, optimization = opt[names(opt) != "fit"],
    est = est, metrics = metrics,
    prediction = pred$mean, prediction_cov = pred$cov), fit_file(fold))
  print(metrics)
}

