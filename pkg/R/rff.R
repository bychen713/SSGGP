#' Build random Fourier features
#'
#' Constructs the coefficient matrix \eqn{c_{i\ell} = L^{-1/2}\cos(\langle s_i,
#' X_\ell\rangle + \psi_\ell)} and the design matrix
#' \eqn{A = C \otimes I_P}.
#'
#' `A` has one row per (location, process) pair with the **process index
#' varying fastest**, matching `as.vector(t(Y))`. Every downstream object
#' indexed by observation assumes that ordering.
#'
#' @param locs `n` by `d` matrix of locations.
#' @param X List of `L` frequencies, each a length-`d` numeric.
#' @param phi Length-`L` numeric of phase offsets.
#' @param P Number of processes.
#' @param need_A Build `A`, which costs `O(n L P)` memory. [fit()] and
#'   [predict.ssggp_fit()] work from `c_vals` alone unless
#'   `options(SSGGP.force_general = TRUE)`; the Woodbury calculation in
#'   `logLik()` needs `A`.
#' @return A list with `c_vals` (`n` by `L`, sparse) and `A` (`nP` by `LP`,
#'   sparse), the latter `NULL` when `need_A` is `FALSE`.
#' @keywords internal
#' @noRd
rff <- function(locs, X, phi, P, need_A = TRUE) {
  L <- length(X)
  n <- nrow(locs)
  stopifnot(length(phi) == L)

  c_vals <- matrix(nrow = n, ncol = L)
  for (l in seq_len(L)) {
    c_vals[, l] <- 1 / sqrt(L) * cos((locs %*% X[[l]]) + phi[l])
  }
  c_vals <- as_dgc(c_vals)
  list(c_vals = c_vals,
       A = if (need_A)
             as_dgc(Matrix::kronecker(c_vals,
                                      Matrix::Diagonal(n = P, x = rep(1, P))))
           else NULL)
}

#' Draw the random Fourier feature frequencies and phases
#'
#' Draws `phi` first, then the frequencies from `g$r`; the order fixes the RNG
#' stream.
#'
#' @param L Number of features.
#' @param d Spatial dimension.
#' @param g Proposal object, from [matern_proposal()].
#' @return A list with `phi` (length `L`) and `X` (list of `L` frequencies).
#' @keywords internal
#' @noRd
draw_features <- function(L, d, g) {
  phi <- stats::runif(n = L, 0, 2 * pi)
  prelim <- g$r(L, d)
  list(phi = phi,
       X = lapply(seq_len(nrow(prelim)), function(i) prelim[i, ]))
}
