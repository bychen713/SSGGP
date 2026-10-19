# LMC stage 1: choose the model structure for one replicate -- shared nugget,
# loadings B with a fixed number of factors, and a covariance family with
# starting values for each latent factor.
#   sbatch slurm/sim/lmc/sim_lmc_structure.sh   (one task per replicate)
# SSGGP_LMC_DATASET selects an existing dataset other than the active setting.
# Writes lmc_structure_<sim>.RData.

rm(list = ls())

if (file.exists("sim/settings.R")) {
  source("sim/settings.R")
} else if (file.exists("../settings.R")) {
  source("../settings.R")
} else {
  stop("run from paper/ or paper/sim/lmc/", call. = FALSE)
}
source(here("sim", "lmc", "lmc_model.R"))
invisible(list2env(lmc_dataset_config(project_root, setting_name), envir = environment()))

array_id <- job_index()
cat("setting:", setting_name, " replicate:", array_id, "\n")

load(file.path(data_dir, paste0(array_id, ".RData")))
P <- ncol(Y_mat)
n <- nrow(grid)
stopifnot(nrow(Y_mat) == n)

t_start <- Sys.time()

# 1. Empirical covariograms and the shared nugget
bins <- lmc_bins(grid)
vgm  <- lmc_covariogram(Y_mat, bins)
cat(sprintf("bins: %d, cutoff %.3f, width %.3f, lags %.3f..%.3f\n",
            vgm$n_bins, bins$cutoff, bins$width, min(vgm$h), max(vgm$h)))

nug_est <- lmc_shared_nugget(Y_mat, vgm)
nugget_hat <- nug_est$nugget
cat(sprintf("shared nugget (median over %d processes): %.4f  [IQR %.4f, %.4f]\n",
            sum(!is.na(nug_est$per_process)), nugget_hat,
            quantile(nug_est$per_process, .25, na.rm = TRUE),
            quantile(nug_est$per_process, .75, na.rm = TRUE)))

# 2. Lag-0 covariance
if (identical(missing_TF, "coloc")) {
  if (anyNA(Y_mat)) stop("missing_TF is 'coloc' but Y_mat has NAs", call. = FALSE)
  S0 <- stats::cov(Y_mat)
  min_shared <- n
} else {
  # Each pair uses the locations where both processes are observed.
  O  <- !is.na(Y_mat)
  Y0 <- sweep(Y_mat, 2, colMeans(Y_mat, na.rm = TRUE)); Y0[!O] <- 0
  Nsh <- crossprod(O * 1)
  min_shared <- min(Nsh)
  if (min_shared < 30) {
    warning("some process pair shares only ", min_shared,
            " locations; its lag-0 covariance is poorly determined")
  }
  S0 <- crossprod(Y0) / Nsh
  S0[!is.finite(S0)] <- 0
  S0 <- (S0 + t(S0)) / 2
}
cat(sprintf("lag-0 covariance assembled; min shared locations per pair: %d\n",
            min_shared))

# 3. Loadings with a fixed number of latent factors
q_default <- if (identical(prec_structure, "tridiagonal")) 6L else 10L
q_fixed <- as.integer(Sys.getenv("SSGGP_LMC_Q", as.character(q_default)))
stopifnot(length(q_fixed) == 1L, !is.na(q_fixed), q_fixed >= 1L, q_fixed <= P)
require_all_families <- !identical(Sys.getenv("SSGGP_LMC_REQUIRE_ALL_FAMILIES", "1"), "0")
subtract_loading_nugget <- Sys.getenv("SSGGP_LMC_SUBTRACT_LOADING_NUGGET", "1")
if (!subtract_loading_nugget %in% c("0", "1")) {
  stop("SSGGP_LMC_SUBTRACT_LOADING_NUGGET must be 0 or 1")
}
loading_nugget <- if (subtract_loading_nugget == "1") nugget_hat else 0
ld <- lmc_loadings(S0, loading_nugget, q_max = q_fixed, var_target = 0.80,
                   q_min = q_fixed)
B <- ld$B; V <- ld$V; d <- ld$d; q_fit <- ld$q
cat(sprintf("q = %d (fixed by %s structure; empirical 80%% rule wanted %d); retained variance %.3f\n",
            q_fit, prec_structure, ld$q_rule, ld$frac))
cat(sprintf("loading covariance: %s pilot nugget %.4f\n",
            if (loading_nugget == 0) "raw (not subtracting)" else "subtracting",
            nugget_hat))
cat("  d_k:", sprintf("%.2f", d), "\n")

# 4. Covariance family and starting values for each latent factor
fit_families <- if (family_mix) LMC_FAMILIES else "matern"
candidate_fits <- vector("list", q_fit)

if (identical(missing_TF, "coloc")) {
  # Colocated: ML fit to each principal-component score.
  selection_method <- if (family_mix) {
    "PCA-score Gaussian likelihood AIC (conditional on B and pilot nugget)"
  } else "Fixed Matern family; PCA-score Gaussian likelihood initialization"
  Zs <- Y_mat %*% V                        # n x q scores
  for (k in seq_len(q_fit)) {
    candidate_fits[[k]] <- setNames(lapply(fit_families, lmc_score_fit,
      z = Zs[, k], dk = d[k], nugget = nugget_hat, D = bins$D), fit_families)
  }
  rho_hat <- NULL
} else {
  # Misaligned: latent correlation curves from the cross-covariograms, by
  # weighted regression on the loadings at each lag.
  X <- vapply(seq_len(q_fit), function(k) as.vector(tcrossprod(B[, k])),
              numeric(P * P))
  rho_hat <- matrix(NA_real_, vgm$n_bins, q_fit)
  for (l in seq_len(vgm$n_bins)) {
    w <- as.vector(vgm$N[[l]]); y <- as.vector(vgm$C[[l]])
    Xw <- X * w
    sol <- tryCatch(solve(crossprod(Xw, X), crossprod(Xw, y)),
                    error = function(e) rep(NA_real_, q_fit))
    rho_hat[l, ] <- sol
  }
  # Moment estimates need not lie in [-1, 1].
  rho_hat[] <- pmin(pmax(rho_hat, -1), 1)

  selection_method <- if (family_mix) {
    "WLS Gaussian working-likelihood AIC (not spatial-data AIC)"
  } else "Fixed Matern family; WLS correlation-curve initialization"
  wl <- vapply(vgm$N, sum, numeric(1))     # location pairs in each bin
  for (k in seq_len(q_fit)) {
    y <- rho_hat[, k]
    ok <- is.finite(y)
    # Latent sill fixed at 1; the loadings carry the scale.
    candidate_fits[[k]] <- setNames(lapply(fit_families, lmc_wls,
      y = y[ok], w = wl[ok], h = vgm$h[ok], sill_max = 1, fix_sill = TRUE),
      fit_families)
  }
  Zs <- NULL
}

sel_detail <- selection_aic <- family_selection <- NULL
if (family_mix) {
  sel_detail <- lapply(candidate_fits, function(fits) {
    vapply(fits, `[[`, numeric(1), "aic")
  })
  selection_aic <- do.call(rbind, sel_detail)
  if (require_all_families) {
    family_selection <- lmc_select_families(selection_aic)
    fam <- family_selection$family
    stopifnot(all(LMC_FAMILIES %in% fam))
  } else {
    if (any(rowSums(is.finite(selection_aic)) == 0L))
      stop("a latent factor has no finite candidate AIC")
    scores <- selection_aic
    scores[!is.finite(scores)] <- Inf
    index <- apply(scores, 1L, which.min)
    fam <- LMC_FAMILIES[index]
    family_selection <- list(family = fam, index = index,
      total_aic = sum(scores[cbind(seq_len(q_fit), index)]),
      aic_increase = 0, changed = integer(0))
    selection_method <- paste(selection_method, "without all-family coverage")
  }
} else {
  fam <- rep("matern", q_fit)
}
par_raw <- lapply(seq_len(q_fit), function(k) candidate_fits[[k]][[fam[k]]]$par)
if (any(!vapply(seq_len(q_fit), function(k)
    isTRUE(candidate_fits[[k]][[fam[k]]]$converged), logical(1)))) {
  stop("a selected latent-family initialization failed", call. = FALSE)
}
cat("selection criterion:", selection_method, "\n")
if (family_mix) {
  cat(sprintf("family coverage: total AIC %.4f, increase over unconstrained %.4f; changed factors: %s\n",
              family_selection$total_aic, family_selection$aic_increase,
              paste(family_selection$changed, collapse = ", ")))
}
cat("latent families:", paste(sprintf("%d:%s", seq_len(q_fit), fam), collapse = " "), "\n")
cat("  selected:", paste(names(table(fam)), table(fam), collapse = ", "), "\n")

time_lmc_est <- difftime(Sys.time(), t_start, units = "mins")
cat(sprintf("stage 1 took %.2f min\n", as.numeric(time_lmc_est)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save(B, V, d, q_fit, fam, par_raw, sel_detail, nugget_hat, nug_est,
     loading_nugget,
     require_all_families,
     candidate_fits, selection_aic, selection_method, family_selection,
     rho_hat, vgm, bins, ld, min_shared, time_lmc_est,
     file = file.path(out_dir, paste0("lmc_structure_", array_id, ".RData")))
cat("wrote", file.path(out_dir, paste0("lmc_structure_", array_id, ".RData")), "\n")
