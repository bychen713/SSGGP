## CM3: the spatial parameters, by coordinate-wise ECM -----------------------
##
## Writing h(X_l) = diag(s_l1, ..., s_lP) with s_lj = sqrt(f_j(X_l; iota_j)),
## and u_l = E[W_l W_l^T | Y, theta^(m-1)], the part of the complete-data
## Q-function that depends on iota is
##
##   Q_iota(iota) = - sum_l log|h(X_l)|
##                  - (1/4) sum_l g(X_l) tr[ Q h(X_l)^-1 u_l h(X_l)^-1 ].
##
## Holding iota_k fixed for k != j this reduces to a function of iota_j alone:
##
##   Q_j(iota_j) = - sum_l log s_lj
##                 - (1/4) sum_l g(X_l) [ Q_jj u_l,jj / s_lj^2 + 2 b_lj / s_lj ]
##
## with b_lj = sum_{k != j} Q_jk u_l,jk / s_lk, which does not involve iota_j
## and so is formed once per coordinate; each objective evaluation is then O(L).
##
## The sweep is sequential: b_lj reads s_lk for the neighbors k of j, so
## coordinate j uses the values coordinates 1..j-1 have just produced.

#' Evaluate and validate a spectral density
#'
#' Returns \eqn{f_j(X_\ell;\iota_j)}, rejecting zeros: CM3 divides by the
#' density and takes its logarithm.
#' @keywords internal
#' @noRd
f_checked <- function(fj, X, par) {
  f <- fj(X, par)
  if (anyNA(f) || any(!is.finite(f)) || any(f <= 0)) {
    stop("spectral density must return finite, strictly positive values",
         call. = FALSE)
  }
  f
}

## The importance ratio, and its ceiling ------------------------------------
##
## The model needs the spectral density only through R_lj = g(X_l) / f_j(X_l):
## M = sqrt(R/2), and the CM3 objective is a function of R. Forming R from
## log g - log f lets a light-tailed target survive a heavy-tailed proposal,
## where f underflows long before R overflows.
##
## Under `numerical_guard = "cap"`, R is clamped at R_CAP. That floors f at
## g / R_CAP, so the spectral mass discarded on the clamped set A is
## int_A f = int_A (f/g) g <= P_g(A) / R_CAP.
R_CAP <- 1e300

#' sqrt of the importance ratio, with the ceiling applied
#'
#' @param fam List of `P` resolved registry entries.
#' @param f_par List of `P` named parameter vectors.
#' @param Xmat `L` by `d` frequencies; `r2` their squared norms.
#' @param log_gvals `log g(X_l)`, a constant of the fit.
#' @return A list with `sqrtR` (`L` by `P`) and `capped` (length `P` counts).
#' @keywords internal
#' @noRd
build_sqrt_R <- function(fam, f_par, Xmat, r2, log_gvals) {
  P <- length(fam)
  L <- nrow(Xmat)
  sqrtR <- matrix(NA_real_, L, P)
  capped <- integer(P)
  lcap <- log(R_CAP)
  for (j in seq_len(P)) {
    lfun <- log_density_of(fam[[j]])
    if (is.null(lfun)) {
      stop("numerical_guard = \"cap\" needs an analytic log-density, which ",
           "family '", fam[[j]]$name, "' (process ", j, ") does not have. ",
           "The bound behind the ceiling assumes a valid density and cannot ",
           "be checked for a user-registered one; use the default ",
           "numerical_guard = \"error\".", call. = FALSE)
    }
    lf <- lfun(Xmat, f_par[[j]][fam[[j]]$par_names], r2)
    if (anyNA(lf)) {
      stop("log-density returned NA for process ", j, call. = FALSE)
    }
    lr <- log_gvals - lf
    over <- lr > lcap
    capped[j] <- sum(over)
    lr[over] <- lcap
    sqrtR[, j] <- exp(0.5 * lr)
  }
  if (any(!is.finite(sqrtR))) {
    stop("non-finite importance ratio after the ceiling; the proposal is far ",
         "from the target", call. = FALSE)
  }
  list(sqrtR = sqrtR, capped = capped)
}

#' Row and column indices of a packed lower triangle
#'
#' `feature_moments_cpp` stores each `u_l` as its lower triangle, column-major.
#' These are the maps between that packing and `(i, j)` subscripts, built once
#' per fit.
#'
#' @param P Number of processes.
#' @return A list with `i` and `j` (length `P(P+1)/2` each, `i >= j`) and `idx`,
#'   a `P` by `P` integer matrix of column positions, symmetric.
#' @keywords internal
#' @noRd
tri_pack <- function(P) {
  ii <- integer(0); jj <- integer(0)
  for (j in seq_len(P)) {
    ii <- c(ii, j:P)
    jj <- c(jj, rep.int(j, P - j + 1L))
  }
  idx <- matrix(NA_integer_, P, P)
  idx[cbind(ii, jj)] <- seq_along(ii)
  idx[cbind(jj, ii)] <- seq_along(ii)
  list(i = ii, j = jj, idx = idx)
}

#' The graphical LASSO statistic, recovered from the per-feature moments
#'
#' \eqn{\hat S_{jk} = L^{-1} \sum_\ell M_{\ell j} M_{\ell k} u_{\ell,jk}}.
#'
#' This is the graphical LASSO input of CM1,
#' \eqn{L^{-1}\sum_\ell (g_\ell/2) h_\ell^{-1} u_\ell h_\ell^{-1}}, written in
#' terms of `M`.
#'
#' @param u `L` by `P(P+1)/2` matrix of packed per-feature moments.
#' @param M `L` by `P` matrix from `build_M()`.
#' @param tri Packing maps from `tri_pack()`.
#' @return A `P` by `P` symmetric matrix.
#' @keywords internal
#' @noRd
suffstat_to_s <- function(u, M, tri) {
  suffstat_to_s_cpp(u, M, as.integer(tri$i), as.integer(tri$j))
}

#' The iota-dependent part of the complete-data Q-function
#'
#' The term \eqn{-\sum_\ell \log|h(X_\ell)|} of the complete-data
#' log-likelihood. The trace term is already in
#' \eqn{(L/2)[\log|Q| - \operatorname{tr}(\hat S Q)]}.
#'
#' @param S_mat `L` by `P` matrix of \eqn{s_{\ell j} = \sqrt{f_j(X_\ell)}}, or
#'   of \eqn{\sqrt{R_{\ell j}}} when `log_gvals` is supplied.
#' @param log_gvals `log g(X_\ell)` when `S_mat` holds the ratio, which puts the
#'   value on the same scale as without the ceiling.
#' @return A single number.
#' @keywords internal
#' @noRd
q_iota_logdet <- function(S_mat, log_gvals = NULL) {
  if (is.null(log_gvals)) return(-sum(log(S_mat)))
  # S_mat is sqrt(R) = sqrt(g/f), so -sum(log sqrt(f)) = sum(log S) - P/2 sum(log g)
  sum(log(S_mat)) - 0.5 * ncol(S_mat) * sum(log_gvals)
}

#' One sequential coordinate-ECM sweep over the spatial parameters
#'
#' Cycles `j = 1, ..., P`, maximizing `Q_j` over `iota_j` with the other
#' coordinates held at their current values -- including the values produced
#' earlier in this same sweep.
#'
#' @param f_par List of `P` named parameter vectors, the current \eqn{\iota}.
#' @param fam List of `P` resolved registry entries.
#' @param f_vec List of `P` vectorized densities, from `f_vec_of()`.
#' @param Xmat `L` by `d` matrix of frequencies.
#' @param gvals Length-`L` proposal densities.
#' @param Q `P` by `P` precision, from CM1.
#' @param u `L` by `P(P+1)/2` packed per-feature moments.
#' @param tri Packing maps from `tri_pack()`.
#' @param bounds List of `P` `list(lower =, upper =)` on the log scale.
#' @param S_mat `L` by `P` current \eqn{\sqrt{f_j(X_\ell)}} (or
#'   \eqn{\sqrt{R_{\ell j}}} under the ceiling), refreshed as the sweep goes.
#' @param maxit Optimizer iterations per coordinate.
#' @param verbose Print the optimizer trace.
#' @param r2 Squared frequency norms `rowSums(Xmat^2)`; `NULL` recomputes them.
#' @param gradient Supply `optim()` with the analytic gradient of `Q_j` where
#'   the family has one.
#' @param nthread Threads for the CH density map; `NULL` runs serially.
#' @param guard,log_gvals The numerical guard and `log g(X_l)`; see [fit()].
#' @return A list with `f_par`, the refreshed `S_mat`, and `capped`, the
#'   per-process count of ratios at the ceiling.
#' @keywords internal
#' @noRd
cm3_sweep <- function(f_par, fam, f_vec, Xmat, gvals, Q, u, tri, bounds,
                      S_mat, maxit = 300, verbose = FALSE,
                      r2 = NULL, gradient = FALSE, nthread = 1L,
                      guard = "error", log_gvals = NULL) {
  P <- length(f_par)
  L <- nrow(Xmat)
  d <- ncol(Xmat)
  if (is.null(r2)) r2 <- rowSums(Xmat^2)
  nthread <- if (is.null(nthread)) 1L else max(1L, as.integer(nthread))
  # Under the ceiling S_mat holds sqrt(R) and the objective is written in R.
  # No analytic gradient: where the ceiling binds, R is constant in iota.
  capping <- identical(guard, "cap")
  if (capping) {
    gradient <- FALSE
    log_cap <- log(R_CAP)
  }
  capped <- integer(P)
  for (j in seq_len(P)) {
    par_names <- fam[[j]]$par_names
    fj        <- f_vec[[j]]
    Qjj       <- Q[j, j]
    u_jj      <- u[, tri$idx[j, j]]
    kern      <- cm3_kernel_of(fam[[j]])
    use_grad  <- gradient && !is.null(kern) && kern$grad

    # Neighbor sum b_lj = sum_k Q_jk u_l,jk / s_lk, formed once per
    # coordinate. Under the ceiling it is sum_k Q_jk u_l,jk sqrt(R_lk).
    nb <- which(Q[j, ] != 0)
    nb <- nb[nb != j]
    b  <- numeric(L)
    if (capping) {
      for (k in nb) b <- b + Q[j, k] * u[, tri$idx[j, k]] * S_mat[, k]
    } else {
      for (k in nb) b <- b + Q[j, k] * u[, tri$idx[j, k]] / S_mat[, k]
    }

    # finite but far above any attainable value, so L-BFGS-B walks away from it
    UNEVALUABLE <- 1e10

    # a zero-width box component (lower == upper) holds that parameter fixed
    lower_j <- bounds[[j]]$lower
    upper_j <- bounds[[j]]$upper
    free  <- lower_j < upper_j
    fixed <- !free
    full_start <- log(f_par[[j]][par_names])
    full_start[fixed] <- lower_j[fixed]

    neg_Qj_full <- if (capping) function(lg) {
      par <- exp(lg)
      names(par) <- par_names
      v <- cm3_objective_cpp(kern$id, par, r2, log_gvals, u_jj, b, Qjj, d,
                             FALSE, nthread, log_cap)$value
      if (!is.finite(v)) UNEVALUABLE else v
    } else if (!is.null(kern)) function(lg) {
      v <- cm3_objective_cpp(kern$id, exp(lg), r2, gvals, u_jj, b, Qjj, d,
                             FALSE, nthread)$value
      if (!is.finite(v)) UNEVALUABLE else v
    } else function(lg) {
      par <- exp(lg)
      names(par) <- par_names
      f <- tryCatch(fj(Xmat, par), error = function(e) NULL)
      if (is.null(f) || length(f) != L || anyNA(f) || any(!is.finite(f)) ||
          any(f <= 0)) return(UNEVALUABLE)
      val <- -0.5 * sum(log(f)) -
        0.25 * sum(gvals * (Qjj * u_jj / f + 2 * b / sqrt(f)))
      if (!is.finite(val)) return(UNEVALUABLE)
      -val
    }

    # gradient in the free log-parameters; zero at a rejected point, which the
    # line search rejects on its value
    grad_Qj_full <- if (use_grad) function(lg) {
      g <- cm3_objective_cpp(kern$id, exp(lg), r2, gvals, u_jj, b, Qjj, d,
                             TRUE, nthread)$grad
      if (is.null(g) || anyNA(g) || any(!is.finite(g))) rep(0, length(lg)) else g
    } else NULL

    refresh_state <- function(new) {
      if (capping) {
        one <- build_sqrt_R(fam[j], list(new), Xmat, r2, log_gvals)
        capped[j] <<- one$capped[1]
        one$sqrtR[, 1]
      } else {
        sqrt(f_checked(fj, Xmat, new))
      }
    }

    if (!any(free)) {                 # nothing to optimize
      new <- exp(full_start)
      names(new) <- par_names
      f_par[[j]][par_names] <- new
      S_mat[, j] <- refresh_state(new)
      next
    }

    neg_Qj <- function(lg_free) {
      lg <- full_start
      lg[free] <- lg_free
      neg_Qj_full(lg)
    }
    grad_Qj <- if (use_grad) function(lg_free) {
      lg <- full_start
      lg[free] <- lg_free
      grad_Qj_full(lg)[free]
    } else NULL

    opt <- stats::optim(
      par = full_start[free], fn = neg_Qj, gr = grad_Qj,
      lower = lower_j[free], upper = upper_j[free], method = "L-BFGS-B",
      control = list(maxit = maxit, parscale = rep(1, sum(free)),
                     trace = if (verbose) 100 else 0))
    if (!is.finite(opt$value) || opt$value >= UNEVALUABLE) {
      stop("CM3 optimizer did not return an evaluable solution for coordinate ",
           j, call. = FALSE)
    }
    lg_new <- full_start
    lg_new[free] <- opt$par
    new <- exp(lg_new)
    names(new) <- par_names
    f_par[[j]][par_names] <- new
    S_mat[, j] <- refresh_state(new)
  }

  list(f_par = f_par, S_mat = S_mat, capped = capped)
}

#' The matrix of sqrt(f_j(X_l)) for the current parameters
#'
#' Under the ceiling, the matrix of \eqn{\sqrt{R_{\ell j}}} instead.
#'
#' @inheritParams cm3_sweep
#' @return An `L` by `P` matrix, strictly positive.
#' @keywords internal
#' @noRd
build_S_mat <- function(f_vec, f_par, fam, Xmat, gvals,
                        guard = "error", r2 = NULL, log_gvals = NULL) {
  if (identical(guard, "cap")) {
    if (is.null(r2)) r2 <- rowSums(Xmat^2)
    return(build_sqrt_R(fam, f_par, Xmat, r2, log_gvals)$sqrtR)
  }
  P <- length(f_vec)
  S <- matrix(NA_real_, nrow(Xmat), P)
  for (j in seq_len(P)) {
    S[, j] <- sqrt(f_checked(f_vec[[j]], Xmat,
                             f_par[[j]][fam[[j]]$par_names]))
  }
  S
}
