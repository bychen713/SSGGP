#' Marginal log-likelihood of a fitted model
#'
#' Evaluates the log-likelihood of the observed data after integrating out the
#' latent process, for a fresh draw of `L` random Fourier features. This is the
#' quantity the penalty is selected on: fit across a grid of `rho`, compare by
#' [AIC()], and take the best.
#'
#' @section It is a Monte Carlo estimate:
#' The value depends on the feature draw, so it varies with `seed` and settles
#' as `L` grows. Comparisons across `rho` should use the **same** `L` and
#' `seed`, so that every model sees the same features and the comparison is
#' paired. [ssggp_aic()] does this for you.
#'
#' @param object An `ssggp_fit` from [fit()].
#' @param L Number of random Fourier features. The likelihood is evaluated
#'   exactly for these features, by Woodbury on an \eqn{(Lp) \times (Lp)}
#'   sparse Cholesky.
#' @param seed Optional integer seed for the feature draw.
#' @param g Proposal from [matern_proposal()]; `NULL` inherits the fit's.
#' @param numerical_guard As in [fit()]; `NULL` inherits the fit's.
#' @param threads Thread cap for the duration of the call; see [fit()].
#' @param ... Unused.
#' @return An object of class `logLik`. `df` is the number of **selected
#'   edges**, \eqn{k=\sum_{i<j}\mathbb{1}(\hat\eta_{ij}\neq 0)}: strictly
#'   off-diagonal, so the diagonal and the process-specific parameters -- whose
#'   dimension is constant across the penalty grid -- are excluded. `nobs` is
#'   the number of observations.
#' @export
logLik.ssggp_fit <- function(object, L = 1000, seed = NULL,
                             g = NULL,
                             threads = getOption("SSGGP.threads", 5L),
                             numerical_guard = NULL, ...) {
  with_threads(resolve_threads(threads),
               loglik_impl(object, L, seed, g, numerical_guard))
}

#' @keywords internal
#' @noRd
loglik_impl <- function(object, L, seed, g, numerical_guard = NULL) {
  sp <- object$spec
  P  <- sp$P
  if (is.null(g)) g <- sp$g
  check_proposal(g)

  Y <- object$data$Y
  locs <- object$data$locs
  n <- nrow(Y)

  Y_reordered <- as.vector(t(Y))
  obs_idx <- which(!is.na(Y_reordered))
  Y_obs   <- Y_reordered[obs_idx]
  proc_id <- rep(seq_len(P), n)[obs_idx]

  nug <- object$nugget
  if (length(nug) == 1L) nug <- rep(nug, P)
  d_obs <- nug[proc_id]

  if (!is.null(seed)) set.seed(seed)
  feat <- draw_features(L, sp$d, g)
  ff   <- rff(locs, feat$X, feat$phi, P)
  A_obs <- ff$A[obs_idx, , drop = FALSE]

  gvals <- vapply(feat$X, g$d, numeric(1))
  M <- build_M_for_fit(object, feat$X, gvals, g, numerical_guard)

  val <- ll_cpp(M, as_dgc(object$sigma_inv), as_dgc(A_obs), Y_obs, d_obs)

  # selected edges, each counted once
  si <- object$sigma_inv
  k <- sum(si[upper.tri(si)] != 0)

  structure(val, class = "logLik",
            df = k,
            nobs = length(Y_obs),
            L = L)
}

#' Select the graphical LASSO penalty by AIC
#'
#' Calls [fit()] once for each value of `rho`, evaluates each fit's integrated
#' log-likelihood with [logLik.ssggp_fit()], and selects the penalty with the
#' smallest \eqn{\mathrm{AIC} = -2\,\ell + 2k}, where \eqn{k} is the number of
#' selected edges. Every fit is scored with the same `L_ll` features (drawn with
#' `seed_ll`), so the comparison is not confounded by Monte Carlo variation.
#'
#' @param Y,locs As in [fit()].
#' @param rho Numeric vector of penalties to try.
#' @param L_ll Number of random Fourier features for the likelihood evaluation.
#' @param seed_ll Seed for those features, shared by all fits.
#' @param ... Every other argument of [fit()] -- for example `f_fun`, `f_par`,
#'   `g`, `nugget`, `L`, `MC`, `iter`, `lag` and `seed` -- passed unchanged to
#'   each fit. Do not pass `rho` here.
#' @return A list with `best` (the selected fit), `rho` (the chosen penalty),
#'   `table` (a data frame of `rho`, `logLik`, `df` and `AIC`), and `fits` (all
#'   fits, in the order of `rho`).
#' @examples
#' \donttest{
#' locs <- as.matrix(expand.grid(seq(0, 5, length = 8), seq(0, 5, length = 8)))
#' g    <- matern_proposal(a = 3, nu = 1.35)
#' dat  <- sim(locs, sigma_inv = trimat(3, 2, 1), f_par = c(a = 3, nu = 1.5),
#'             nugget = 0.1, method = "rff", L = 300, g = g, seed = 1)
#' sel <- ssggp_aic(dat$Y, locs, rho = c(0.01, 0.05), L_ll = 200, seed_ll = 1,
#'                  f_par = c(a = 2.5, nu = 1.5), g = g,
#'                  L = 200, MC = 10, iter = 5, seed = 2)
#' sel$table
#' sel$rho
#' }
#' @export
ssggp_aic <- function(Y, locs, rho = c(0, 0.005, 0.01, 0.025, 0.05, 0.1),
                      L_ll = 1000, seed_ll = 1L, ...) {
  fits <- lapply(rho, function(r) fit(Y, locs, rho = r, ...))
  lls  <- lapply(fits, function(f) stats::logLik(f, L = L_ll, seed = seed_ll))
  tab  <- data.frame(
    rho    = rho,
    logLik = vapply(lls, as.numeric, numeric(1)),
    df     = vapply(lls, function(z) attr(z, "df"), numeric(1)))
  tab$AIC <- -2 * tab$logLik + 2 * tab$df
  best <- which.min(tab$AIC)
  list(best = fits[[best]], rho = rho[best], table = tab, fits = fits)
}

#' Select the graphical LASSO penalty by BIC
#'
#' Fits the model across a grid of penalties and compares them by BIC on the
#' marginal likelihood. Every fit is evaluated with the same likelihood
#' features. The caller supplies the sample size `n`; it is used only in the
#' BIC penalty and is not inferred from the dimensions of `Y`.
#'
#' @param Y,locs As in [fit()].
#' @param n Sample size to use in the BIC penalty. Must be one integer greater
#'   than one. For misaligned data, choose the intended definition of sample
#'   size explicitly.
#' @param rho Numeric vector of penalties to try.
#' @param L_ll Features used for likelihood evaluation.
#' @param seed_ll Seed for the likelihood features, shared across all fits.
#' @param ... Passed to [fit()].
#' @return A list with `best` (the selected fit), `rho` (the chosen penalty),
#'   `n` (the supplied sample size), `table` (a data frame of `rho`, `logLik`,
#'   `df` and `BIC`), and `fits` (all fitted models). As in [ssggp_aic()],
#'   `df` counts selected off-diagonal edges; parameters constant across the
#'   penalty grid are omitted from the reported criterion.
#' @export
ssggp_bic <- function(Y, locs, n,
                      rho = c(0, 0.005, 0.01, 0.025, 0.05, 0.1),
                      L_ll = 1000, seed_ll = 1L, ...) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      !is.finite(n) || n < 2 || n != floor(n)) {
    stop("`n` must be one finite integer greater than one", call. = FALSE)
  }

  fits <- lapply(rho, function(r) fit(Y, locs, rho = r, ...))
  lls  <- lapply(fits, function(f) stats::logLik(f, L = L_ll, seed = seed_ll))
  tab  <- data.frame(
    rho    = rho,
    logLik = vapply(lls, as.numeric, numeric(1)),
    df     = vapply(lls, function(z) attr(z, "df"), numeric(1)))
  tab$BIC <- -2 * tab$logLik + log(n) * tab$df
  best <- which.min(tab$BIC)
  list(best = fits[[best]], rho = rho[best], n = n, table = tab, fits = fits)
}
