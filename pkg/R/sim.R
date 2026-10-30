#' Simulate a highly-multivariate Gaussian process
#'
#' Draws from a multivariate Gaussian process whose spectral density is
#' \eqn{f(x) = h(x)\sigma h(x)^\top} with
#' \eqn{h(x) = \mathrm{diag}(f_j(x))^{1/2}}, so that zeros in
#' \eqn{\sigma^{-1}} encode process-level conditional independence.
#'
#' Two generators are available and they are **not** interchangeable:
#'
#' * `"exact"` draws from the true covariance function by Cholesky
#'   factorization, so the data are not produced by the same random Fourier
#'   feature approximation being fitted. Cost is \eqn{O(n^3p^3)}, and only some
#'   families have a closed form (see below).
#' * `"rff"` is the simulation algorithm of the paper: an `L`-term random Fourier
#'   feature construction, \eqn{O(Lp^2 + Lpn)}, for any spectral family. The
#'   draw is an approximation that improves with `L`.
#'
#' In both cases the nugget \eqn{\epsilon_i \sim N(0, d)} is added to
#' the simulated process.
#'
#' @section What `"exact"` draws from:
#' It needs a closed form for
#' \eqn{\int e^{i h^\top x}\sigma_{jk}\sqrt{f_j(x)f_k(x)}\,dx}.
#'
#' **Common range.** The geometric mean merges exactly:
#' \eqn{\sqrt{f_jf_k}(x) = \sqrt{c_jc_k}\,(a^2+|x|^2)^{-((\nu_j+\nu_k)/2+d/2)}},
#' a Matern density at the averaged smoothness and the same range. So for any
#' \eqn{\nu_j} this **is** the model -- the parsimonious multivariate Matern of
#' Gneiting, Kleiber and Schlather (2010) -- valid for every positive-definite
#' \eqn{\sigma}, because the spectral matrix factorizes as
#' \eqn{D(x)\sigma D(x)} with \eqn{D=\mathrm{diag}(\sqrt{f_j})}.
#'
#' **Differing ranges.** As a function of \eqn{z=|x|^2} the product has
#' singularities at \eqn{-a_j^2} and \eqn{-a_k^2} while a Matern density has
#' one, so no \eqn{(\tilde a,\tilde\nu)} reproduces it and the familiar
#' \eqn{(a_j+a_k)/2} expression is an ansatz. It is still a valid
#' cross-covariance, so it is permitted -- with a warning, because the
#' generating model is then **not** the one [fit()] estimates. Use
#' `method = "rff"` for the exact model at any \eqn{a_j}.
#'
#' @param locs `n` by `d` matrix of locations.
#' @param sigma `p` by `p` positive-definite cross-process covariance. Supply
#'   this or `sigma_inv`, not both.
#' @param sigma_inv `p` by `p` precision. Zeros here are the conditional
#'   independence structure.
#' @param f_par Named numeric of spectral parameters recycled to all processes,
#'   or a list of `p` such vectors. These are the TRUE values; they are not
#'   passed to [fit()].
#' @param f_fun Spectral family: a registered name, a raw density function, or a
#'   list of length `p`. Unlike [fit()], a raw function is fine here -- no
#'   covariance implementation is needed.
#' @param nugget Scalar or length-`p` measurement error variance.
#' @param method `"exact"` or `"rff"`.
#' @param L Number of random Fourier features; required for `method = "rff"`.
#' @param g Importance-sampling proposal from [matern_proposal()], used by
#'   `method = "rff"` only. May differ from the one used in [fit()].
#' @param numerical_guard As in [fit()]; `"cap"` is needed for a squared exponential (`se`)
#'   process under a Matern proposal.
#' @param seed Optional integer seed.
#' @return A list with `Y` (`n` by `p`) and `locs`.
#' @export
#' @examples
#' locs <- as.matrix(expand.grid(seq(0, 5, length = 8), seq(0, 5, length = 8)))
#' out  <- sim(locs, sigma_inv = trimat(3, 2, 1), f_par = c(a = 3, nu = 1.5),
#'             nugget = 0.1, method = "rff", L = 200, seed = 1)
#' dim(out$Y)
sim <- function(locs, sigma = NULL, sigma_inv = NULL,
                f_par = c(a = 3, nu = 1.5), f_fun = "matern",
                nugget = 0,
                method = c("exact", "rff"), L = NULL,
                g = matern_proposal(a = 7, nu = 1.25),
                numerical_guard = c("error", "cap"),
                seed = NULL) {
  method <- match.arg(method)
  numerical_guard <- match.arg(numerical_guard)
  locs <- as.matrix(locs)
  n <- nrow(locs); d <- ncol(locs)

  if (is.null(sigma) == is.null(sigma_inv)) {
    stop("supply exactly one of `sigma` or `sigma_inv`", call. = FALSE)
  }
  if (is.null(sigma_inv)) sigma_inv <- solve(sigma) else sigma_inv <- as.matrix(sigma_inv)
  if (is.null(sigma))     sigma     <- solve(sigma_inv) else sigma <- as.matrix(sigma)
  P <- ncol(sigma_inv)
  if (nrow(sigma_inv) != P) stop("sigma / sigma_inv must be square", call. = FALSE)
  if (any(eigen(sigma_inv, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("sigma_inv must be positive definite", call. = FALSE)
  }

  fam <- resolve_spectral(f_fun, P)          # raw closures are allowed here
  f_den <- lapply(fam, `[[`, "f")
  if (!is.list(f_par)) f_par <- rep(list(f_par), P)
  if (length(f_par) != P) stop("f_par must be length P = ", P, call. = FALSE)

  nugget <- if (length(nugget) == 1L) rep(nugget, P) else nugget
  if (length(nugget) != P) stop("nugget must be scalar or length P = ", P, call. = FALSE)
  if (any(nugget < 0)) stop("nugget must be non-negative", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)

  Y <- switch(method,
    exact = sim_exact(locs, sigma, f_par, fam, n, P, d),
    rff   = sim_rff(locs, sigma_inv, f_den, f_par, g, L, n, P, d,
                    fam, numerical_guard))

  if (any(nugget > 0)) {
    Y <- Y + matrix(stats::rnorm(n * P), n, P) *
      rep(sqrt(nugget), each = n)
  }
  list(Y = Y, locs = locs)
}

#' Random Fourier feature simulation
#'
#' With \eqn{\sigma^{-1} = \Lambda\Lambda^\top}, solves
#' \eqn{M_\ell W_\ell = Z_\ell} for \eqn{M_\ell = \Lambda^\top
#' \mathrm{diag}(M[\ell, ])}, so \eqn{W_\ell \sim N(0, 2 f(X_\ell)/g(X_\ell))};
#' one Cholesky serves every feature.
#' @keywords internal
#' @noRd
sim_rff <- function(locs, sigma_inv, f_den, f_par, g, L, n, P, d,
                    fam = NULL, numerical_guard = "error") {
  if (is.null(L)) stop("`L` is required for method = \"rff\"", call. = FALSE)
  check_proposal(g)

  feat <- draw_features(L, d, g)
  Xmat <- do.call(rbind, feat$X)
  gvals <- vapply(feat$X, g$d, numeric(1))
  # under the ceiling, features where f has underflowed give that process an
  # essentially zero coefficient
  M <- if (identical(numerical_guard, "cap")) {
    if (is.null(fam)) stop("numerical_guard = \"cap\" needs the resolved ",
                           "families", call. = FALSE)
    sr <- build_sqrt_R(fam, f_par, Xmat, rowSums(Xmat^2), log(gvals))
    if (sum(sr$capped)) {
      warning(sum(sr$capped), " of ", L * P, " importance ratios hit the ",
              "ceiling while simulating (processes ",
              paste(which(sr$capped > 0L), collapse = ", "), ")",
              call. = FALSE)
    }
    build_M(feat$X, f_den, f_par, g, gvals = gvals, s_mat = sr$sqrtR,
            guard = "cap")
  } else build_M(feat$X, f_den, f_par, g)   # L x P

  # sigma_inv = Lambda Lambda^T, so chol() gives R with Lambda^T = R
  R <- chol(sigma_inv)
  Z <- matrix(stats::rnorm(P * L), nrow = P, ncol = L)

  W <- matrix(NA_real_, nrow = P, ncol = L)
  for (l in seq_len(L)) {
    # M_l = Lambda^T diag(M[l, ]); solve M_l W_l = Z_l
    W[, l] <- backsolve(R * rep(M[l, ], each = P), Z[, l])
  }

  # c_{i,l} = L^{-1/2} cos(<X_l, s_i> + psi_l); Y = C W^T
  C <- matrix(NA_real_, nrow = n, ncol = L)
  for (l in seq_len(L)) {
    C[, l] <- 1 / sqrt(L) * cos((locs %*% feat$X[[l]]) + feat$phi[l])
  }
  C %*% t(W)
}

#' Exact draw from the model's cross-covariance
#'
#' Dispatches on the family, which must be the same registered one for every
#' process: the closed forms below are derived per family from
#' \eqn{\sqrt{f_j f_k}}, and there is no such form for a mixture.
#' @keywords internal
#' @noRd
sim_exact <- function(locs, sigma, f_par, fam, n, P, d) {
  fams <- vapply(fam, function(e) if (is.null(e$name)) NA_character_ else e$name,
                 character(1))
  if (anyNA(fams)) {
    stop("method = \"exact\" needs a registered family for every process, but ",
         "a raw density function has no closed-form covariance; use ",
         "method = \"rff\"", call. = FALSE)
  }
  if (length(unique(fams)) != 1L) {
    stop("method = \"exact\" needs one family for every process, got ",
         paste(unique(fams), collapse = " and "),
         "; the cross-covariance of a mixture has no closed form. Use ",
         "method = \"rff\"", call. = FALSE)
  }
  switch(fams[[1]],
    matern   = sim_exact_matern(locs, sigma, f_par, n, P, d),
    se = sim_exact_se(locs, sigma, f_par, n, P, d),
    stop("method = \"exact\" is not implemented for family '", fams[[1]],
         "'; use method = \"rff\"", call. = FALSE))
}

#' Exact draw from the squared-exponential cross-covariance
#'
#' Unlike Matern, this family has an exact closed form for **any** ranges: the
#' geometric mean of two squared-exponential spectral densities is another one,
#' at the averaged squared range phi_jk^2 = (phi_j^2 + phi_k^2)/2, so
#'
#'   C_jk(h) = sigma_jk (phi_j phi_k / phi_jk^2)^(d/2) exp(-|h|^2/(2 phi_jk^2)).
#'
#' There is therefore no unequal-range caveat here.
#' @keywords internal
#' @noRd
sim_exact_se <- function(locs, sigma, f_par, n, P, d) {
  a   <- vapply(f_par, `[[`, numeric(1), "a")
  phi <- 1 / a
  dist_mat <- fields::rdist(locs)
  common_a <- diff(range(a)) <= .Machine$double.eps^0.5

  # Kronecker shortcut, as in the Matern branch: with a common range the
  # prefactor collapses to one and the covariance is exactly sigma %x% C, which
  # avoids the O((nP)^3) dense route. The same rnorm(n * P) is consumed in the
  # same order as the dense branch.
  if (common_a) {
    C  <- exp(-(dist_mat / phi[1])^2 / 2)
    Rc <- chol(C)
    Rs <- chol(sigma)
    Z <- matrix(stats::rnorm(n * P), nrow = n, ncol = P)
    return(crossprod(Rc, Z) %*% Rs)
  }

  cov_mat <- matrix(NA_real_, nrow = n * P, ncol = n * P)
  for (j in seq_len(P)) {
    for (k in j:P) {
      ij <- (n * (j - 1) + 1):(n * j)
      ik <- (n * (k - 1) + 1):(n * k)
      phi2 <- (phi[j]^2 + phi[k]^2) / 2
      block <- sigma[j, k] * (phi[j] * phi[k] / phi2)^(d / 2) *
        exp(-dist_mat^2 / (2 * phi2))
      cov_mat[ik, ij] <- block
      if (j != k) cov_mat[ij, ik] <- t(block)
    }
  }
  ch <- chol(cov_mat)
  matrix(as.vector(t(ch) %*% stats::rnorm(n * P)), nrow = n, ncol = P)
}

#' Exact draw from the parsimonious Matern cross-covariance
#' @keywords internal
#' @noRd
sim_exact_matern <- function(locs, sigma, f_par, n, P, d) {
  a  <- vapply(f_par, `[[`, numeric(1), "a")
  nu <- vapply(f_par, `[[`, numeric(1), "nu")
  common_a  <- diff(range(a))  <= .Machine$double.eps^0.5
  common_nu <- diff(range(nu)) <= .Machine$double.eps^0.5

  # With a COMMON range the two Matern densities merge exactly and the gamma
  # ratio below is the resulting constant, so the covariance built here IS the
  # model for any nu_j. With DIFFERING a_j no (a~, nu~) reproduces the product,
  # and the (a_j + a_k)/2 substitution is an ansatz: a valid cross-covariance,
  # but a different model from the one fit() estimates. See the "What
  # \"exact\" draws from" section of sim().
  if (!common_a) {
    warning("method = \"exact\" with unequal ranges uses the (a_j + a_k)/2 ",
            "cross-covariance, which is NOT the Fourier transform of ",
            "sqrt(f_j f_k); the generating model then differs from the one ",
            "fit() estimates. Use method = \"rff\" for the exact model.",
            call. = FALSE)
  }

  dist_mat <- fields::rdist(locs)

  # Kronecker shortcut, valid only when a AND nu are both common: the gamma
  # prefactor then collapses to 1 and the covariance is exactly sigma %x% C,
  # whose Cholesky factor is chol(sigma) %x% chol(C). Since (A %x% B) vec(X) =
  # vec(B X t(A)), the draw needs only the two small factors -- n x n and P x P
  # instead of nP x nP -- and consumes the same rnorm(n * P) in the same order
  # as the dense branch. Common a with DIFFERING nu no longer factorizes and
  # must fall through.
  if (common_a && common_nu) {
    C  <- fields::Matern(d = dist_mat, smoothness = nu[1], alpha = a[1])
    Rc <- chol(C)
    Rs <- chol(sigma)
    Z <- matrix(stats::rnorm(n * P), nrow = n, ncol = P)
    return(crossprod(Rc, Z) %*% Rs)
  }

  cov_mat  <- matrix(NA_real_, nrow = n * P, ncol = n * P)
  for (j in seq_len(P)) {
    for (k in j:P) {
      ij <- (n * (j - 1) + 1):(n * j)
      ik <- (n * (k - 1) + 1):(n * k)
      block <- sqrt(gamma(nu[j] + d / 2) * gamma(nu[k] + d / 2) /
                    gamma(nu[j]) / gamma(nu[k])) *
        gamma((nu[j] + nu[k]) / 2) / gamma((nu[j] + nu[k] + d) / 2) *
        sigma[j, k] * fields::Matern(d = dist_mat,
                                     smoothness = (nu[j] + nu[k]) / 2,
                                     alpha = (a[j] + a[k]) / 2)
      cov_mat[ik, ij] <- block
      if (j != k) cov_mat[ij, ik] <- t(block)
    }
  }
  ch <- chol(cov_mat)
  # process-major ordering, matching matrix(., n, P, byrow = FALSE)
  matrix(as.vector(t(ch) %*% stats::rnorm(n * P)), nrow = n, ncol = P)
}

#' Tridiagonal precision matrix
#'
#' @param P Dimension.
#' @param diag,offdiag Diagonal and off-diagonal values.
#' @return A `P` by `P` symmetric tridiagonal matrix.
#' @export
#' @examples
#' trimat(4, 2, 1)
trimat <- function(P, diag = 2, offdiag = 1) {
  m <- matrix(0, P, P)
  base::diag(m) <- diag
  if (P > 1) {
    i <- seq_len(P - 1)
    m[cbind(i, i + 1)] <- offdiag
    m[cbind(i + 1, i)] <- offdiag
  }
  m
}
