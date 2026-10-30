#' Matern spectral density
#'
#' The isotropic Matern spectral density on \eqn{\mathbb{R}^d}, normalized to
#' integrate to one. Unit integral is required: the model puts every marginal
#' scale in \eqn{\sigma}, so a spectral density carrying its own scale would
#' make \eqn{\sigma} unidentifiable.
#'
#' The range is parameterized by `a`, the inverse range, so `a = 1/phi` for
#' the usual Matern range `phi`.
#'
#' @param x Numeric vector of length `d`: one frequency.
#' @param par Named numeric with elements `a` (inverse range) and `nu`
#'   (smoothness).
#' @return A single non-negative number, \eqn{f_j(x)}.
#' @export
#' @examples
#' matern_spectral(c(0, 0), c(a = 3, nu = 1.5))
matern_spectral <- function(x, par) {
  a <- par[["a"]]; nu <- par[["nu"]]; d <- length(x)
  a^(2 * nu) * gamma(nu + d / 2) / (gamma(nu) * pi^(d / 2)) *
    (a^2 + sum(x^2))^(-nu - d / 2)
}

#' Vectorized spectral densities
#'
#' `f_vec(X, par)` over the whole `L` by `d` frequency matrix, returning `L`
#' densities: the scalar function with `sum(x^2)` replaced by `rowSums(X^2)`.
#'
#' @param X `L` by `d` matrix of frequencies, one per row.
#' @param par Named parameter vector, as for the scalar form.
#' @return A length-`nrow(X)` numeric vector.
#' @keywords internal
#' @noRd
matern_spectral_vec <- function(X, par) {
  a <- par[["a"]]; nu <- par[["nu"]]; d <- ncol(X)
  a^(2 * nu) * gamma(nu + d / 2) / (gamma(nu) * pi^(d / 2)) *
    (a^2 + rowSums(X^2))^(-nu - d / 2)
}

#' @keywords internal
#' @noRd
se_spectral_vec <- function(X, par) {
  a <- par[["a"]]; d <- ncol(X)
  a^(-d) * (2 * pi)^(-d / 2) * exp(-rowSums(X^2) / (2 * a^2))
}

#' @keywords internal
#' @noRd
ch_spectral_vec <- function(X, par) {
  ch_spectral_cpp(sqrt(rowSums(X^2)), par[["beta"]], par[["nu"]],
                  par[["alpha"]], ncol(X))
}

#' A Matern importance-sampling proposal
#'
#' Bundles the proposal density \eqn{g} with a sampler from it, so the
#' frequencies \eqn{X_\ell \sim g} and the weights \eqn{2f(X_\ell)/g(X_\ell)}
#' use the same \eqn{g}. The proposal should have heavier tails than every
#' marginal spectral density.
#'
#' The sampler is the Gaussian scale mixture \eqn{X = Z/\sqrt{G}},
#' \eqn{G \sim \mathrm{Gamma}(\nu, a^2/2)}, \eqn{Z \sim N(0, I_d)}.
#'
#' @param a Inverse range of the proposal.
#' @param nu Smoothness of the proposal.
#' @return A list with `d(x)`, the density, and `r(L, dim)`, a sampler
#'   returning an `L` by `dim` matrix of frequencies. Also carries `par` and
#'   `name` for reporting.
#' @export
#' @examples
#' g <- matern_proposal(a = 7, nu = 1.25)
#' g$d(c(0, 0))
#' dim(g$r(5, 2))
matern_proposal <- function(a = 7, nu = 1.25) {
  force(a); force(nu)
  stopifnot(is.numeric(a), is.numeric(nu), length(a) == 1L, length(nu) == 1L,
            is.finite(a), is.finite(nu), a > 0, nu > 0)
  list(
    name = "matern",
    par  = c(a = a, nu = nu),
    d    = function(x) matern_spectral(x, c(a = a, nu = nu)),
    r    = function(L, dim) {
      1 / sqrt(stats::rgamma(n = L, shape = nu, rate = (a^2) / 2)) *
        matrix(nrow = L, ncol = dim, stats::rnorm(n = dim * L, sd = 1))
    })
}

#' Squared exponential spectral density
#'
#' The isotropic squared-exponential spectral density on \eqn{\mathbb{R}^d},
#' normalized to integrate to one. It corresponds to the covariance
#' \eqn{\sigma\exp(-\|h\|^2/(2\phi^2))}.
#'
#' As for [matern_spectral()], `a = 1/phi` is the inverse range.
#'
#' @section Choosing a proposal:
#' This density has Gaussian tails, so under a Matern proposal it underflows at
#' frequencies the proposal reaches routinely and \eqn{M = \sqrt{g/(2f)}}
#' overflows. When every process is squared exponential, use
#' [se_proposal()] with `a` larger than every process's `a`. A mixed fit
#' needs the heavy-tailed proposal; use `numerical_guard = "cap"` in [fit()]
#' if the density underflows.
#'
#' @param x Numeric vector of length `d`: one frequency.
#' @param par Named numeric with element `a` (inverse range).
#' @return A single non-negative number, \eqn{f_j(x)}.
#' @export
#' @examples
#' se_spectral(c(0, 0), c(a = 3))
#' check_unit_integral(se_spectral, c(a = 3))
se_spectral <- function(x, par) {
  a <- par[["a"]]; d <- length(x)
  a^(-d) * (2 * pi)^(-d / 2) * exp(-sum(x^2) / (2 * a^2))
}

#' A squared exponential importance-sampling proposal
#'
#' The \eqn{N(0, a^2 I_d)} density and sampler. Use it only when every process
#' is `"se"`, with `a` larger than every process's `a`; Matern and CH
#' densities have polynomial tails that it cannot dominate.
#'
#' @param a Inverse range of the proposal.
#' @return A list with `d(x)`, the density, and `r(L, dim)`, a sampler
#'   returning an `L` by `dim` matrix of frequencies. Also carries `par` and
#'   `name` for reporting.
#' @export
#' @examples
#' g <- se_proposal(a = 5)
#' g$d(c(0, 0))
#' dim(g$r(5, 2))
se_proposal <- function(a = 5) {
  force(a)
  stopifnot(is.numeric(a), length(a) == 1L, is.finite(a), a > 0)
  list(
    name = "se",
    par  = c(a = a),
    d    = function(x) se_spectral(x, c(a = a)),
    r    = function(L, dim) {
      matrix(nrow = L, ncol = dim, stats::rnorm(n = dim * L, sd = a))
    })
}

#' Confluent hypergeometric spectral density
#'
#' The isotropic Confluent Hypergeometric (CH) spectral density on
#' \eqn{\mathbb{R}^d}, normalized to integrate to one:
#'
#' \deqn{f(x) = \frac{\Gamma(\nu + d/2)\,\beta^d}{(2\pi)^{d/2} B(\alpha,\nu)}
#'              \, U\!\left(\nu + \tfrac d2,\, 1-\alpha+\tfrac d2,\,
#'              \tfrac{\beta^2\|x\|^2}{2}\right).}
#'
#' The corresponding covariance is
#' \eqn{\Gamma(\nu+\alpha)/\Gamma(\nu)\; U(\alpha, 1-\nu, \|h\|^2/(2\beta^2))}.
#'
#' @section What the parameters do:
#' \eqn{\nu} controls smoothness as in Matern, \eqn{\alpha} the polynomial
#' covariance tail \eqn{\|h\|^{-2\alpha}}, and \eqn{\beta} is a range (not an
#' inverse range). \eqn{\alpha > d/2} is required for the spectral density to
#' be finite at the origin; the registered bounds assume \eqn{d = 2}, so pass
#' `lower` to [fit()] for other dimensions.
#'
#' @param x Numeric vector of length `d`: one frequency.
#' @param par Named numeric with `beta` (range), `nu` (smoothness) and `alpha`
#'   (tail decay).
#' @return A single non-negative number, \eqn{f_j(x)}.
#' @export
#' @examples
#' ch_spectral(c(0, 0), c(beta = 0.5, nu = 1.5, alpha = 2))
#' check_unit_integral(ch_spectral, c(beta = 0.5, nu = 1.5, alpha = 2))
ch_spectral <- function(x, par) {
  ch_spectral_cpp(sqrt(sum(x^2)), par[["beta"]], par[["nu"]], par[["alpha"]],
                  length(x))[1]
}

#' Validate a proposal object
#' @keywords internal
#' @noRd
check_proposal <- function(g) {
  if (!is.list(g) || !is.function(g$d) || !is.function(g$r)) {
    stop("`g` must be a proposal object with a density `d(x)` and a sampler ",
         "`r(L, dim)`, such as matern_proposal(a = 7, nu = 1.25). Passing a ",
         "bare density is not supported: the sampler must match it, or the ",
         "importance weights do not correspond to what was drawn.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Build the per-feature scaling matrix M
#'
#' `M[l, j] = sqrt( g(X_l) / (2 f_j(X_l)) )`, the diagonal of
#' \eqn{M_\ell} in the simulation algorithm, so that
#' `diag(M[l, ]) %*% sigma_inv %*% diag(M[l, ])` is \eqn{\Sigma_{W_\ell}^{-1}}.
#'
#' @param X List of `L` frequencies, each a length-`d` numeric.
#' @param f_fun List of `P` spectral density functions, or one recycled to `P`.
#' @param f_par List of `P` named parameter vectors, or one recycled to `P`.
#' @param g Proposal object from [matern_proposal()].
#' @param f_vec Optional list of `P` vectorized densities, `function(X, par)`;
#'   `NULL` loops over `f_fun`.
#' @param gvals Optional length-`L` proposal densities; `NULL` evaluates `g$d`.
#' @param s_mat Optional `L` by `P` matrix of \eqn{\sqrt{f_j(X_\ell)}} at the
#'   same `f_par` (from `build_S_mat()`), which skips evaluating the densities.
#'   Under `guard = "cap"` it holds \eqn{\sqrt{R} = \sqrt{g/f}} and is required.
#' @param guard `"error"` or `"cap"`; see [fit()].
#' @return An `L` by `P` numeric matrix.
#' @keywords internal
#' @noRd
build_M <- function(X, f_fun, f_par, g, f_vec = NULL, gvals = NULL,
                    s_mat = NULL, guard = "error") {
  L <- length(X)
  P <- length(f_fun)
  if (is.null(gvals)) gvals <- vapply(X, g$d, numeric(1))

  if (!is.null(s_mat)) {
    # under the ceiling s_mat is sqrt(R) = sqrt(g/f), so M = sqrt(R/2)
    M <- if (identical(guard, "cap")) s_mat / sqrt(2)
         else sqrt(gvals) / (s_mat * sqrt(2))
    dim(M) <- c(L, P)
    return(check_M(M))
  }
  if (identical(guard, "cap")) {
    stop("numerical_guard = \"cap\" requires s_mat = sqrt(R)", call. = FALSE)
  }

  Xmat <- if (is.matrix(X)) X else do.call(rbind, X)
  FD <- vapply(seq_len(P), function(j) {
    if (!is.null(f_vec)) f_vec[[j]](Xmat, f_par[[j]])
    else vapply(X, function(x) f_fun[[j]](x, f_par[[j]]), numeric(1))
  }, numeric(L))
  dim(FD) <- c(L, P)
  if (anyNA(FD) || any(FD < 0)) {
    j <- which(colSums(is.na(FD) | FD < 0) > 0)[1]
    stop("spectral density of process ", j, " returned NA or a negative value",
         call. = FALSE)
  }

  # sqrt(g) / (sqrt(f) * sqrt(2)), the same form as the s_mat branch
  M <- sqrt(gvals) / (sqrt(FD) * sqrt(2))
  dim(M) <- c(L, P)
  check_M(M)
}

#' Reject an importance scaling that is not finite or cannot be squared
#' @keywords internal
#' @noRd
check_M <- function(M) {
  if (any(!is.finite(M))) {
    stop("Non-finite importance scaling in M; the spectral density may have ",
         "underflowed to zero for the chosen proposal", call. = FALSE)
  }
  if (any(M > sqrt(.Machine$double.xmax))) {
    stop("Importance scaling in M is too large to square without overflow; ",
         "choose a proposal closer to the target spectral density",
         call. = FALSE)
  }
  M
}

#' Build M for a fitted object, honouring the guard the fit used
#'
#' `predict()` and `logLik()` rebuild `M` at fresh frequencies. By default they
#' use the guard the model was fitted with.
#'
#' @param object An `ssggp_fit`.
#' @param X List of frequencies; `gvals` the proposal densities.
#' @param g Proposal.
#' @param guard `"error"` or `"cap"`; `NULL` inherits the fit's.
#' @return An `L` by `P` matrix.
#' @keywords internal
#' @noRd
build_M_for_fit <- function(object, X, gvals, g, guard = NULL) {
  P   <- object$spec$P
  fam <- resolve_spectral(object$spec$f_fun, P)
  if (is.null(guard)) {
    guard <- object$spec$numerical_guard
    if (is.null(guard)) guard <- "error"        # fits made before the guard
  }
  guard <- match.arg(guard, c("error", "cap"))
  f_den <- lapply(fam, `[[`, "f")
  if (!identical(guard, "cap")) {
    return(build_M(X, f_den, object$f_par, g, gvals = gvals))
  }
  Xmat <- if (is.matrix(X)) X else do.call(rbind, X)
  sr <- build_sqrt_R(fam, object$f_par, Xmat, rowSums(Xmat^2), log(gvals))
  if (sum(sr$capped)) {
    warning(sum(sr$capped), " of ", length(gvals) * P, " importance ratios hit ",
            "the ceiling (processes ",
            paste(which(sr$capped > 0L), collapse = ", "), ")", call. = FALSE)
  }
  build_M(X, f_den, object$f_par, g, gvals = gvals, s_mat = sr$sqrtR,
          guard = "cap")
}

#' Check that a spectral density integrates to one
#'
#' Radial quadrature of an isotropic spectral density. The model requires each
#' \eqn{f_j} to integrate to one, so that every marginal scale is in
#' \eqn{\sigma}.
#'
#' @param f Spectral density function of `(x, par)`.
#' @param par Named parameter vector.
#' @param d Spatial dimension.
#' @param upper Upper radial limit of integration.
#' @return The approximate integral, which should be near 1.
#' @export
#' @examples
#' check_unit_integral(matern_spectral, c(a = 3, nu = 1.5))
check_unit_integral <- function(f, par, d = 2, upper = Inf) {
  stopifnot(d >= 1)
  # surface area of the unit (d-1)-sphere
  area <- 2 * pi^(d / 2) / gamma(d / 2)
  integrand <- function(r) {
    vapply(r, function(rr) f(c(rr, rep(0, d - 1)), par) * rr^(d - 1), numeric(1))
  }
  area * stats::integrate(integrand, lower = 0, upper = upper,
                          rel.tol = 1e-10)$value
}
