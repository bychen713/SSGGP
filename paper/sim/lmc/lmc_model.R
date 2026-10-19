# Linear model of coregionalization (LMC), the baseline method.
# Y(s) = B w(s) + eps(s): B is P x q, the latent processes w_k are independent
# with unit variance, and eps has a single nugget variance shared by all processes.

suppressMessages({
  library(fields)
  library(SSGGP)
})

# Data and output directories for one simulation setting. SSGGP_LMC_DATASET
# selects an existing dataset other than the one configured in settings.R.
lmc_dataset_config <- function(project_root, default_setting,
                               dataset = Sys.getenv("SSGGP_LMC_DATASET", default_setting)) {
  if (!grepl("^(random|tridiagonal)_(pos|neg)_(coloc|miss)(_mix)?$", dataset)) {
    stop("invalid SSGGP_LMC_DATASET: ", dataset, call. = FALSE)
  }
  parts <- strsplit(dataset, "_", fixed = TRUE)[[1]]
  data_dir <- file.path(project_root, "sim", "data", dataset)
  if (!dir.exists(data_dir)) stop("dataset directory not found: ", data_dir, call. = FALSE)
  out_dir <- file.path(project_root, "sim", "output", "simulation", dataset)
  list(setting_name = dataset, prec_structure = parts[1], offdiag_sign = parts[2],
       missing_TF = parts[3], family_mix = length(parts) == 4L,
       data_dir = data_dir, out_dir = out_dir,
       fig_dir = file.path(project_root, "sim", "output", "figure", dataset))
}

# Candidate latent covariance families.
LMC_FAMILIES <- c("matern", "ch", "se")

# Two starting values per family, one short-range and one long-range.
LMC_STARTS <- list(
  matern   = list(c(a = 2, nu = 1), c(a = 0.5, nu = 2.5)),
  ch       = list(c(beta = 0.5, nu = 1, alpha = 3), c(beta = 2, nu = 2, alpha = 6)),
  se       = list(c(a = 2), c(a = 0.5))
)

# Parameter names and log-scale bounds from the SSGGP family registry.
lmc_par_names <- function(family) SSGGP:::get_spectral(family)$par_names
lmc_lower     <- function(family) unname(SSGGP:::get_spectral(family)$par_lower)
lmc_upper     <- function(family) {
  out <- unname(SSGGP:::get_spectral(family)$par_upper)
  # The GSL Tricomi U evaluation is unstable at alpha = 20.
  if (identical(family, "ch")) out[3] <- log(19.9)
  out
}
lmc_npar      <- function(family) length(lmc_par_names(family))

# Assign a family to each latent factor, minimizing the total AIC subject to
# every family being used at least once. `aic` has one row per factor and one
# column per family; the constrained optimum is found by dynamic programming
# over the set of families used so far.
lmc_select_families <- function(aic, families = LMC_FAMILIES) {
  aic <- as.matrix(aic)
  nf <- length(families)
  if (!is.numeric(aic) || nf < 1L || nf > 10L || anyDuplicated(families) ||
      ncol(aic) != nf || nrow(aic) < nf) {
    stop("family coverage requires at least one factor per candidate family",
         call. = FALSE)
  }
  if (!is.null(colnames(aic))) {
    if (anyDuplicated(colnames(aic)) || !setequal(colnames(aic), families)) {
      stop("AIC columns must match the candidate families", call. = FALSE)
    }
    aic <- aic[, families, drop = FALSE]
  }
  aic[!is.finite(aic)] <- Inf
  q <- nrow(aic)
  if (any(rowSums(is.finite(aic)) == 0L)) {
    stop("at least one factor has no valid candidate fit", call. = FALSE)
  }
  independent <- apply(aic, 1L, which.min)
  index <- independent
  if (length(unique(index)) < nf) {
    ns <- bitwShiftL(1L, nf)
    bits <- bitwShiftL(1L, seq_len(nf) - 1L)
    cost <- rep(Inf, ns); cost[1L] <- 0
    parent <- choice <- matrix(NA_integer_, q, ns)
    for (k in seq_len(q)) {
      next_cost <- rep(Inf, ns)
      for (state in which(is.finite(cost))) {
        for (f in seq_len(nf)) {
          next_state <- bitwOr(state - 1L, bits[f]) + 1L
          value <- cost[state] + aic[k, f]
          if (value < next_cost[next_state]) {
            next_cost[next_state] <- value
            parent[k, next_state] <- state
            choice[k, next_state] <- f
          }
        }
      }
      cost <- next_cost
    }
    if (!is.finite(cost[ns])) {
      stop("no finite-AIC allocation can include every family; check candidate fits",
           call. = FALSE)
    }
    state <- ns
    for (k in rev(seq_len(q))) {
      index[k] <- choice[k, state]
      state <- parent[k, state]
    }
  }
  total <- sum(aic[cbind(seq_len(q), index)])
  unconstrained <- sum(aic[cbind(seq_len(q), independent)])
  list(family = families[index], index = index, total_aic = total,
       unconstrained_aic = unconstrained, aic_increase = total - unconstrained,
       changed = which(index != independent),
       counts = table(factor(families[index], levels = families)))
}

# CH correlation at arbitrary lags. ch_cov_cpp evaluates the lower triangle of a
# square matrix, so the lags are packed into the smallest such triangle.
ch_cor_at <- function(hv, par) {
  L <- length(hv)
  m <- ceiling((-1 + sqrt(1 + 8 * L)) / 2)      # smallest m with m(m+1)/2 >= L
  M <- matrix(0, m, m)
  idx <- which(lower.tri(M, diag = TRUE), arr.ind = TRUE)[seq_len(L), , drop = FALSE]
  M[idx] <- hv
  out <- SSGGP:::ch_cov_cpp(M, beta = par[1], nu = par[2], alpha = par[3],
                           variance = 1)[idx]
  lmc_ch_repair(hv, out, par)
}

# Replace CH correlations that GSL failed to evaluate by numerical integration of
#   rho(h) = int_0^1 exp(-z u / (1 - u)) dbeta(u; alpha, nu) du,  z = h^2 / (2 beta^2).
lmc_ch_repair <- function(h, out, par) {
  bad <- which(!is.finite(out) | out < 0 | out > 1)
  if (!length(bad)) return(out)
  bad_h <- unique(as.vector(h)[bad])
  values <- vapply(bad_h, function(lag) {
    if (lag == 0) return(1)
    z <- lag^2 / (2 * par[1]^2)
    integrate(function(u) {
      value <- exp(-z * u / (1 - u) + dbeta(u, par[3], par[2], log = TRUE))
      value[u == 0 | u == 1] <- 0
      value
    }, 0, 1, rel.tol = 1e-10, abs.tol = 1e-12, subdivisions = 300L)$value
  }, numeric(1))
  out[bad] <- values[match(as.vector(h)[bad], bad_h)]
  out
}

# Correlation of one family at lags `h`: a vector, a symmetric distance matrix, or
# a rectangular cross-distance matrix. Use symmetric = TRUE only for a symmetric
# matrix, because the CH branch then reads just its lower triangle.
lmc_cor <- function(family, par, h, symmetric = FALSE) {
  par <- as.numeric(par)
  if (family == "matern") {
    fields::Matern(h, alpha = par[1], smoothness = par[2])
  } else if (family == "se") {
    exp(-(par[1] * h)^2 / 2)
  } else if (family == "ch") {
    if (symmetric) {
      stopifnot(is.matrix(h), nrow(h) == ncol(h))
      out <- SSGGP:::ch_cov_cpp(h, beta = par[1], nu = par[2], alpha = par[3],
                               variance = 1)
      lmc_ch_repair(h, out, par)
    } else {
      out <- ch_cor_at(as.vector(h), par)
      if (is.matrix(h)) matrix(out, nrow(h), ncol(h)) else out
    }
  } else {
    stop("unknown LMC family '", family, "'", call. = FALSE)
  }
}

# NULL when the correlation cannot be evaluated, so optimizers treat the point as infeasible.
lmc_cor_safe <- function(family, par, h, symmetric = FALSE) {
  r <- tryCatch(lmc_cor(family, par, h, symmetric = symmetric),
                error = function(e) NULL)
  if (is.null(r) || any(!is.finite(r))) NULL else r
}

# A square root L with L L' = R. Falls back to an eigendecomposition when Cholesky
# fails on a numerically singular correlation (e.g. a smooth squared exponential kernel).
lmc_cor_root <- function(R) {
  cs <- tryCatch(chol(R), error = function(e) NULL)
  if (!is.null(cs)) return(t(cs))
  ee <- eigen(R, symmetric = TRUE)
  tolerance <- 100 * .Machine$double.eps * max(1, max(abs(ee$values)))
  if (any(!is.finite(ee$values)) || min(ee$values) < -tolerance) {
    stop("latent correlation is indefinite beyond numerical roundoff", call. = FALSE)
  }
  sweep(ee$vectors, 2, sqrt(pmax(ee$values, 0)), `*`)
}

# Parameters on a common scale for reporting: `a` is an inverse range (1 / beta
# for CH) and `nu` the smoothness (NA for the squared exponential family).
lmc_canonical <- function(family, par) {
  par <- as.numeric(par)
  switch(family,
    matern   = c(a = par[1],     nu = par[2]),
    se = c(a = par[1],     nu = NA_real_),
    ch       = c(a = 1 / par[1], nu = par[2]),
    stop("unknown LMC family '", family, "'", call. = FALSE))
}

# Empirical covariogram

# Distance bins shared by all process pairs: 15 equal-width bins up to one third
# of the bounding-box diagonal, as in gstat::variogram().
lmc_bins <- function(grid, n_bins = 15L) {
  D <- fields::rdist(grid)
  cutoff <- sqrt(sum(apply(grid, 2, function(x) diff(range(x)))^2)) / 3
  width  <- cutoff / n_bins
  bin <- ceiling(D / width)
  bin[D <= 0 | D > cutoff] <- 0L          # 0 = "not in any bin", incl. s == t
  list(D = D, bin = bin, n_bins = n_bins, cutoff = cutoff, width = width)
}

# Binned empirical cross-covariograms for all process pairs. With W_l the
# indicator matrix of lag bin l, O the observed mask, and Y0 the centered data
# with missing entries set to zero,
#   C_l = (Y0' W_l Y0) / (O' W_l O)   (elementwise),
# and N_l = O' W_l O is the number of location pairs behind each entry.
lmc_covariogram <- function(Y_mat, bins) {
  O  <- !is.na(Y_mat)
  mu <- colMeans(Y_mat, na.rm = TRUE)
  Y0 <- sweep(Y_mat, 2, mu)
  Y0[!O] <- 0
  Ostore <- O * 1

  nb <- bins$n_bins
  C <- N <- vector("list", nb)
  hbar <- rep(NA_real_, nb)
  for (l in seq_len(nb)) {
    W <- (bins$bin == l) * 1
    N[[l]] <- crossprod(Ostore, W %*% Ostore)
    Cl     <- crossprod(Y0, W %*% Y0) / N[[l]]
    # Empty cells are 0/0; set them to 0 (their pair count is also 0).
    Cl[!is.finite(Cl)] <- 0
    C[[l]] <- Cl
    hbar[l] <- if (sum(W) > 0) sum(bins$D * W) / sum(W) else NA_real_
  }
  keep <- which(!is.na(hbar) & vapply(N, function(x) max(x) > 0, logical(1)))
  list(C = C[keep], N = N[keep], h = hbar[keep], n_bins = length(keep))
}

# Latent starting values

# Gaussian ML fit of one family to a principal-component score z with eigenvalue
# dk and the given nugget. The AIC is used to compare families.
lmc_score_fit <- function(family, z, dk, nugget, D) {
  obj <- function(lp) {
    R <- lmc_cor_safe(family, exp(lp), D, symmetric = TRUE)
    if (is.null(R)) return(1e12)
    S <- dk * R; diag(S) <- diag(S) + nugget
    cs <- tryCatch(chol(S), error = function(e) NULL)
    if (is.null(cs) || any(!is.finite(diag(cs))) || min(diag(cs)) <= 0) return(1e12)
    u <- backsolve(cs, z, transpose = TRUE)
    v <- sum(log(diag(cs))) + 0.5 * sum(u^2)
    if (is.finite(v)) v else 1e12
  }
  best <- NULL
  for (st in LMC_STARTS[[family]]) {
    o <- tryCatch(optim(log(st), obj, method = "L-BFGS-B",
                        lower = lmc_lower(family), upper = lmc_upper(family)),
                  error = function(e) NULL)
    if (!is.null(o) && o$convergence == 0L && is.finite(o$value) &&
        o$value < 1e12 && (is.null(best) || o$value < best$value)) best <- o
  }
  np <- lmc_npar(family)
  if (is.null(best)) {
    return(list(family = family, par = NA, nll = Inf, npar = np,
                aic = Inf, converged = FALSE))
  }
  nll <- best$value + 0.5 * length(z) * log(2 * pi)
  list(family = family, par = setNames(exp(best$par), lmc_par_names(family)),
       nll = nll, npar = np, aic = 2 * nll + 2 * np, converged = TRUE)
}

# AIC of a weighted least-squares curve fit under a working model with
# Var(e_l) = sigma^2 / w_l. Covariogram bins are dependent, so this is a working
# criterion for comparing families, not the likelihood of the spatial data.
lmc_wls_aic <- function(wss, w, npar) {
  L <- length(w)
  if (!L || any(!is.finite(w) | w <= 0) || !is.finite(wss) || wss < 0) {
    stop("WLS AIC requires finite WSS and positive finite bin weights", call. = FALSE)
  }
  variance <- max(wss / L, .Machine$double.xmin)
  nll <- 0.5 * (L * (log(2 * pi) + 1 + log(variance)) - sum(log(w)))
  list(nll = nll, npar = npar + 1L, aic = 2 * nll + 2 * (npar + 1L),
       residual_variance = variance, nobs = L)
}

# Weighted least-squares fit of s * rho(h; par) to a covariogram curve. The sill s
# is profiled out and capped at `sill_max`; fix_sill = TRUE sets s = 1 for a
# latent correlation curve.
lmc_wls <- function(family, y, w, h, sill_max = Inf, fix_sill = FALSE) {
  stopifnot(length(y) == length(w), length(y) == length(h))
  ok <- is.finite(y) & is.finite(w) & w > 0 & is.finite(h)
  y <- y[ok]; w <- w[ok]; h <- h[ok]
  np <- lmc_npar(family) + as.integer(!fix_sill)
  failed <- list(family = family, par = NA, sill = NA_real_, wss = Inf,
                 nll = Inf, npar = np + 1L, aic = Inf, score = Inf,
                 nobs = length(y), converged = FALSE)
  if (length(y) <= np) return(failed)
  profile_sill <- function(r) {
    if (fix_sill) return(1)
    den <- sum(w * r * r)
    if (!is.finite(den) || den <= 0) return(NA_real_)
    min(max(sum(w * r * y) / den, 0), sill_max)
  }
  obj <- function(lp) {
    r <- lmc_cor_safe(family, exp(lp), h)
    if (is.null(r)) return(1e12)
    s <- profile_sill(r)
    if (!is.finite(s)) return(1e12)
    v <- sum(w * (y - s * r)^2)
    if (!is.finite(v)) 1e12 else v
  }
  best <- NULL
  for (st in LMC_STARTS[[family]]) {
    o <- tryCatch(
      optim(log(st), obj, method = "L-BFGS-B",
            lower = lmc_lower(family), upper = lmc_upper(family)),
      error = function(e) NULL)
    if (!is.null(o) && o$convergence == 0L && is.finite(o$value) &&
        o$value < 1e12 && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(failed)
  r <- lmc_cor_safe(family, exp(best$par), h)
  s <- if (is.null(r)) NA_real_ else profile_sill(r)
  if (!is.finite(s)) return(failed)
  ic <- lmc_wls_aic(best$value, w, np)
  c(list(family = family,
         par = setNames(exp(best$par), lmc_par_names(family)),
         sill = s, wss = best$value, score = ic$aic, converged = TRUE), ic)
}

# Family with the smallest WLS AIC for one curve.
lmc_wls_select <- function(y, w, h, sill_max = Inf, fix_sill = FALSE,
                           families = LMC_FAMILIES) {
  fits <- lapply(families, lmc_wls, y = y, w = w, h = h,
                 sill_max = sill_max, fix_sill = fix_sill)
  fits[[which.min(vapply(fits, `[[`, numeric(1), "aic"))]]
}

# Nugget and loadings

# Shared nugget: for each process, fit s_i rho(h) to its covariogram over lags
# h > 0 and take var(y_i) - s_i. The median over processes is returned because
# individual estimates are noisy.
lmc_shared_nugget <- function(Y_mat, vgm) {
  P <- ncol(Y_mat)
  nug <- rep(NA_real_, P)
  for (i in seq_len(P)) {
    y <- vapply(vgm$C, function(Cl) Cl[i, i], numeric(1))
    w <- vapply(vgm$N, function(Nl) Nl[i, i], numeric(1))
    vi <- stats::var(Y_mat[, i], na.rm = TRUE)
    if (!is.finite(vi) || sum(w) == 0) next
    fit <- lmc_wls_select(y, w, vgm$h, sill_max = vi)
    if (is.finite(fit$sill)) nug[i] <- vi - fit$sill
  }
  if (all(is.na(nug))) stop("could not estimate a nugget from any process")
  list(nugget = stats::median(nug, na.rm = TRUE), per_process = nug)
}

# Loadings from the eigendecomposition of S - nugget * I, with negative
# eigenvalues set to zero. q is the smallest rank explaining `var_target` of the
# signal variance, restricted to [q_min, q_max].
lmc_loadings <- function(S, nugget, q_max, var_target = 0.80, q_min = 1L) {
  if (length(q_min) != 1L || !is.finite(q_min) || q_min < 1L ||
      q_min != as.integer(q_min) || length(q_max) != 1L ||
      !is.finite(q_max) || q_max < q_min) {
    stop("q_max must allow at least q_min positive-variance factors", call. = FALSE)
  }
  S <- (S + t(S)) / 2
  diag(S) <- diag(S) - nugget
  ee  <- eigen(S, symmetric = TRUE)
  lam <- pmax(ee$values, 0)
  if (sum(lam) <= 0) stop("signal covariance has no positive eigenvalues")
  if (sum(lam > 0) < q_min) {
    stop("signal covariance has fewer positive eigenvalues than q_min",
         call. = FALSE)
  }
  cum <- cumsum(lam) / sum(lam)
  q_rule <- which(cum >= var_target)[1]
  if (is.na(q_rule)) q_rule <- length(lam)
  q <- max(as.integer(q_min), min(as.integer(q_rule), as.integer(q_max)))
  V <- ee$vectors[, seq_len(q), drop = FALSE]
  d <- lam[seq_len(q)]
  list(B = V %*% diag(sqrt(d), nrow = q), V = V, d = d, q = q,
       q_rule = q_rule, frac = cum[q], eigenvalues = lam)
}

# Parameters and likelihoods

# All latent parameters and the nugget as one log-scale vector.
lmc_pack <- function(par_list, nugget) {
  c(unlist(lapply(par_list, function(p) log(as.numeric(p)))), log(nugget))
}

lmc_unpack <- function(theta, families) {
  np  <- vapply(families, lmc_npar, integer(1))
  end <- cumsum(np); start <- end - np + 1L
  list(par = lapply(seq_along(families), function(k) exp(theta[start[k]:end[k]])),
       nugget = exp(theta[length(theta)]))
}

lmc_bounds <- function(families) {
  list(lower = c(unlist(lapply(families, lmc_lower)), log(1e-4)),
       upper = c(unlist(lapply(families, lmc_upper)), log(10)))
}

# Minimise the LMC negative log-likelihood `nll` with L-BFGS-B, using optim()'s
# default controls (at most 100 iterations), and keep the run resumable: the
# best parameters seen are written to `checkpoint_file` every `save_every`
# seconds and when the run ends.
#
# A fresh fit (`resume = FALSE`) starts at `start` and refuses to overwrite an
# existing result or checkpoint. With `resume = TRUE` it warm-restarts from
# whichever of `previous` (the saved result) and the checkpoint has the lower
# nll; optim() does not expose its curvature history, so a restart begins a new
# L-BFGS-B pass. Both must carry the same `identity`, a list naming everything
# the fit depends on (data, structure, bounds). Returns NULL when `previous`
# already converged; otherwise the optim() result plus totals over all passes.
lmc_optimize <- function(nll, start, lower, upper, identity, checkpoint_file,
                         resume = FALSE, previous = NULL, elapsed_before = 0,
                         save_every = 300) {
  initial <- list(theta = start, elapsed_min = elapsed_before, evaluations = 0,
                  passes = 0L, source = "initial values")
  candidates <- list()
  if (!resume) {
    if (!is.null(previous) || file.exists(checkpoint_file))
      stop("an LMC result or checkpoint already exists; resume to continue it",
           call. = FALSE)
    candidates$initial <- initial
  } else {
    if (!is.null(previous)) {
      if (!identical(previous$identity, identity))
        stop("the saved LMC result was fitted to different data or settings",
             call. = FALSE)
      if (identical(previous$convergence, 0L)) return(NULL)
      candidates$result <- c(previous[c("theta", "elapsed_min", "evaluations",
                                        "passes")], source = "saved result")
    }
    if (file.exists(checkpoint_file)) {
      ck <- readRDS(checkpoint_file)
      if (!identical(ck$identity, identity))
        stop("the LMC checkpoint was written for different data or settings",
             call. = FALSE)
      candidates$checkpoint <- c(ck[c("theta", "elapsed_min", "evaluations",
                                      "passes")], source = "checkpoint")
    }
    if (!length(candidates))
      stop("no saved LMC result or checkpoint to resume", call. = FALSE)
  }
  values <- vapply(candidates, function(cand) nll(cand$theta), numeric(1))
  passes <- vapply(candidates, function(cand) as.numeric(cand$passes), numeric(1))
  # Lowest nll; on a tie, the later pass, so the time and evaluations of a pass
  # that crashed without improving are kept.
  best <- order(values, -passes)[1L]
  from <- candidates[[best]]
  if (!is.finite(values[best]) || values[best] >= 1e12)
    stop("invalid likelihood at the ", from$source, call. = FALSE)

  # Track the best evaluation, and save it periodically.
  state <- new.env(parent = emptyenv())
  state$theta <- from$theta
  state$value <- values[best]
  state$evaluations <- 0
  state$last_save <- t_start <- Sys.time()
  save_state <- function() {
    now <- Sys.time()
    ck <- list(identity = identity, theta = state$theta, value = state$value,
               passes = from$passes + 1L,
               evaluations = from$evaluations + state$evaluations,
               elapsed_min = from$elapsed_min +
                 as.numeric(difftime(now, t_start, units = "mins")))
    dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(checkpoint_file, ".tmp-", Sys.getpid())
    saveRDS(ck, tmp)
    if (!file.rename(tmp, checkpoint_file))
      stop("could not save LMC checkpoint: ", checkpoint_file)
    state$last_save <- now
    ck
  }
  objective <- function(theta) {
    value <- nll(theta)
    state$evaluations <- state$evaluations + 1
    if (is.finite(value) && value < state$value) {
      state$theta <- theta
      state$value <- value
    }
    if (as.numeric(difftime(Sys.time(), state$last_save, units = "secs")) >= save_every)
      save_state()
    value
  }

  cat(sprintf("LMC optimization pass %d from the %s, nll %.4f\n",
              from$passes + 1L, from$source, values[best]))
  fit <- optim(from$theta, objective, method = "L-BFGS-B",
               lower = lower, upper = upper, control = list(trace = 1, REPORT = 10))
  if (state$value < fit$value) {     # the best point may not be the last one
    fit$par <- state$theta
    fit$value <- state$value
  }
  total <- save_state()
  if (!is.finite(fit$value) || fit$value >= 1e12)
    stop("LMC optimization ended at an invalid likelihood", call. = FALSE)
  if (fit$convergence == 1L)
    warning("L-BFGS-B stopped at optim()'s 100-iteration default; ",
            "resume to continue from the checkpoint")
  list(fit = fit, identity = identity, theta = fit$par, value = fit$value,
       convergence = fit$convergence, start_value = unname(values[best]),
       start_source = from$source, passes = total$passes,
       evaluations = total$evaluations, elapsed_min = total$elapsed_min)
}

# Negative log-likelihood for colocated data. The columns of B are orthogonal, so
# rotating by the eigenvectors decouples the model: score z_k = V_k' Y has
# covariance d_k R_k + s2 I, and the remaining P - q directions are noise only:
#   -2 logL = sum_k [log|d_k R_k + s2 I| + z_k' (d_k R_k + s2 I)^{-1} z_k]
#             + (P - q) n log s2 + (||Y||^2 - sum_k ||z_k||^2) / s2 + nP log 2pi.
# The term of each latent is cached between evaluations.
lmc_nll_coloc_factory <- function(Y_mat, D, V, d, families) {
  n <- nrow(Y_mat); P <- ncol(Y_mat); q <- length(d)
  stopifnot(!anyNA(Y_mat), length(families) == q, ncol(V) == q)
  Z     <- crossprod(V, t(Y_mat))          # q x n scores
  ss_perp <- sum(Y_mat^2) - sum(Z^2)       # ||V_perp' Y||_F^2
  cache <- new.env(parent = emptyenv())
  cache$par  <- vector("list", q)
  cache$s2   <- rep(NA_real_, q)           # the s2 each cached term was built at
  cache$term <- rep(NA_real_, q)

  function(theta) {
    up <- lmc_unpack(theta, families)
    s2 <- up$nugget
    if (!is.finite(s2) || s2 <= 0) return(1e12)
    tot <- 0
    for (k in seq_len(q)) {
      fresh <- !is.null(cache$par[[k]]) &&
               isTRUE(all.equal(up$par[[k]], cache$par[[k]])) &&
               isTRUE(all.equal(s2, cache$s2[k])) &&
               is.finite(cache$term[k])
      if (!fresh) {
        R <- lmc_cor_safe(families[k], up$par[[k]], D, symmetric = TRUE)
        if (is.null(R)) return(1e12)
        S <- d[k] * R
        diag(S) <- diag(S) + s2
        cs <- tryCatch(chol(S), error = function(e) NULL)
        if (is.null(cs) || any(!is.finite(diag(cs))) || min(diag(cs)) <= 0)
          return(1e12)
        u <- backsolve(cs, Z[k, ], transpose = TRUE)
        cache$term[k]  <- sum(log(diag(cs))) + 0.5 * sum(u^2)
        cache$par[[k]] <- up$par[[k]]
        cache$s2[k]    <- s2
      }
      tot <- tot + cache$term[k]
    }
    v <- tot + 0.5 * ((P - q) * n * log(s2) + ss_perp / s2) +
         0.5 * n * P * log(2 * pi)
    if (is.finite(v)) v else 1e12
  }
}

# Negative log-likelihood for misaligned data, using the Woodbury identity on
# M = I + Z' D^{-1} Z. Block (k, l) of Z' Z is L_k' diag(A_kl) L_l, where
# A_kl(s) sums b_ik b_il over the processes i observed at s, so Z is never
# formed. The blocks do not depend on the nugget and are cached.
lmc_nll_miss_factory <- function(Y_mat, D, B, families) {
  n <- nrow(Y_mat); P <- ncol(Y_mat); q <- ncol(B)
  stopifnot(length(families) == q)
  O  <- !is.na(Y_mat)
  Y0 <- Y_mat; Y0[!O] <- 0
  n_obs <- sum(O)
  ss_y  <- sum(Y0^2)
  YB    <- Y0 %*% B                        # n x q, r_k(s) up to the 1/s2
  Om    <- O * 1
  pair  <- which(upper.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
  A <- lapply(seq_len(nrow(pair)), function(p)
    as.vector(Om %*% (B[, pair[p, 1]] * B[, pair[p, 2]])))

  idx <- function(k) ((k - 1) * n + 1):(k * n)
  cache <- new.env(parent = emptyenv())
  cache$par    <- vector("list", q)
  cache$L      <- vector("list", q)
  cache$blocks <- matrix(0, n * q, n * q)   # Z'Z, i.e. before the 1/s2
  cache$ok     <- FALSE

  function(theta) {
    up <- lmc_unpack(theta, families)
    s2 <- up$nugget
    if (!is.finite(s2) || s2 <= 0) return(1e12)

    changed <- integer(0)
    for (k in seq_len(q)) {
      if (is.null(cache$par[[k]]) ||
          !isTRUE(all.equal(up$par[[k]], cache$par[[k]])) || !cache$ok) {
        R <- lmc_cor_safe(families[k], up$par[[k]], D, symmetric = TRUE)
        if (is.null(R)) { cache$ok <- FALSE; return(1e12) }
        Lk <- tryCatch(lmc_cor_root(R), error = function(e) NULL)
        if (is.null(Lk)) { cache$ok <- FALSE; return(1e12) }
        cache$L[[k]]   <- Lk
        cache$par[[k]] <- up$par[[k]]
        changed <- c(changed, k)
      }
    }
    if (length(changed)) {
      for (k in changed) for (l in seq_len(q)) {
        p <- which(pair[, 1] == min(k, l) & pair[, 2] == max(k, l))
        blk <- crossprod(cache$L[[k]] * A[[p]], cache$L[[l]])
        cache$blocks[idx(k), idx(l)] <- blk
        if (l != k) cache$blocks[idx(l), idx(k)] <- t(blk)
      }
      cache$ok <- TRUE
    }

    M <- cache$blocks / s2
    diag(M) <- diag(M) + 1
    cs <- tryCatch(chol(M), error = function(e) NULL)
    if (is.null(cs) || any(!is.finite(diag(cs))) || min(diag(cs)) <= 0)
      return(1e12)

    z <- unlist(lapply(seq_len(q), function(k) crossprod(cache$L[[k]], YB[, k] / s2)))
    v <- backsolve(cs, z, upper.tri = TRUE, transpose = TRUE)
    quad   <- ss_y / s2 - sum(v^2)
    logdet <- n_obs * log(s2) + 2 * sum(log(diag(cs)))
    val <- 0.5 * (logdet + quad + n_obs * log(2 * pi))
    if (is.finite(val)) val else 1e12
  }
}

# Prediction

# Colocated prediction: krige each score and map back through V. The predictive
# covariance at a site is V diag(var_k) V' + s2 I.
lmc_predict_coloc <- function(Y_mat, D, Dx, V, d, families, par, s2) {
  n <- nrow(Y_mat); P <- ncol(Y_mat); q <- length(d); nt <- nrow(Dx)
  Z  <- Y_mat %*% V                       # n x q scores
  mu <- vr <- matrix(NA_real_, nt, q)
  for (k in seq_len(q)) {
    R <- lmc_cor(families[k], par[[k]], D, symmetric = TRUE)
    S <- d[k] * R; diag(S) <- diag(S) + s2
    cs <- chol(S)
    Rx <- lmc_cor(families[k], par[[k]], Dx)          # nt x n
    a  <- backsolve(cs, backsolve(cs, Z[, k], transpose = TRUE))
    mu[, k] <- d[k] * as.vector(Rx %*% a)
    W <- backsolve(cs, t(Rx), transpose = TRUE)       # n x nt
    vr[, k] <- pmax(d[k] - d[k]^2 * colSums(W^2), 0)
  }
  Ip <- diag(s2, P)
  list(mean = mu %*% t(V),
       cov  = lapply(seq_len(nt), function(i) V %*% (t(V) * vr[i, ]) + Ip))
}

# Misaligned prediction with the same blockwise Woodbury form as the likelihood.
# Returns predictive means and the P x P predictive covariance at each test site.
lmc_predict_miss <- function(Y_mat, D, Dx, B, families, par, s2) {
  n <- nrow(Y_mat); P <- ncol(Y_mat); q <- ncol(B); nt <- nrow(Dx)
  O  <- !is.na(Y_mat); Y0 <- Y_mat; Y0[!O] <- 0
  pair <- which(upper.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
  pidx <- function(k, l) which(pair[, 1] == min(k, l) & pair[, 2] == max(k, l))
  A <- lapply(seq_len(nrow(pair)), function(p)
    as.vector((O * 1) %*% (B[, pair[p, 1]] * B[, pair[p, 2]])))
  idx <- function(k) ((k - 1) * n + 1):(k * n)

  L <- Rx <- vector("list", q)
  for (k in seq_len(q)) {
    L[[k]]  <- lmc_cor_root(lmc_cor(families[k], par[[k]], D, symmetric = TRUE))
    Rx[[k]] <- lmc_cor(families[k], par[[k]], Dx)     # nt x n
  }
  M <- matrix(0, n * q, n * q)
  for (k in seq_len(q)) for (l in k:q) {
    blk <- crossprod(L[[k]] * A[[pidx(k, l)]], L[[l]]) / s2
    M[idx(k), idx(l)] <- blk
    if (l != k) M[idx(l), idx(k)] <- t(blk)
  }
  diag(M) <- diag(M) + 1
  cs <- chol(M)

  # K^{-1} y in n x P form
  YB <- Y0 %*% B
  z  <- unlist(lapply(seq_len(q), function(k) crossprod(L[[k]], YB[, k] / s2)))
  vm <- backsolve(cs, backsolve(cs, z, upper.tri = TRUE, transpose = TRUE))
  G  <- vapply(seq_len(q), function(k) as.vector(L[[k]] %*% vm[idx(k)]), numeric(n))
  Ky <- Y0 - G %*% t(B)
  Ky[!O] <- 0
  Ky <- Ky / s2
  U  <- Ky %*% B                                      # n x q

  Tm <- vapply(seq_len(q), function(k) as.vector(Rx[[k]] %*% U[, k]), numeric(nt))
  pred_mean <- Tm %*% t(B)

  # Psi_k blocks for all test sites at once
  Psi <- vector("list", q)
  for (k in seq_len(q)) {
    Psi[[k]] <- lapply(seq_len(q), function(m)
      crossprod(L[[k]] * A[[pidx(k, m)]], t(Rx[[m]])))  # n x nt
  }
  Gq <- array(NA_real_, c(nt, q, q))
  for (k in seq_len(q)) for (l in k:q) {
    v <- as.vector((Rx[[k]] * Rx[[l]]) %*% A[[pidx(k, l)]])
    Gq[, k, l] <- v; Gq[, l, k] <- v
  }
  Kprior <- tcrossprod(B) + diag(s2, P)
  cov <- vector("list", nt)
  for (t in seq_len(nt)) {
    W <- matrix(0, n * q, P)
    for (k in seq_len(q)) {
      Pk <- vapply(seq_len(q), function(m) Psi[[k]][[m]][, t], numeric(n))
      W[idx(k), ] <- Pk %*% t(B) / s2
    }
    v <- backsolve(cs, W, upper.tri = TRUE, transpose = TRUE)
    part1 <- B %*% Gq[t, , ] %*% t(B) / s2
    cov[[t]] <- Kprior - (part1 - crossprod(v))
  }
  list(mean = pred_mean, cov = cov)
}
