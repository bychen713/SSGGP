#' Predict at new locations
#'
#' Draws \eqn{W \mid Y, \hat\theta} by Gibbs sampling at fresh random Fourier
#' features and summarizes the latent predictions
#' \eqn{\hat Y(s_*) = L^{-1/2}\sum_\ell W_\ell \cos(\langle s_*, X_\ell\rangle +
#' \psi_\ell)}. The mean and raw draws are for the latent process; the
#' covariance and variance summaries are the sample covariance of the draws
#' plus the fitted nugget. `L`, `MC` and the proposal are independent of those
#' used in [fit()].
#'
#' `summarise = "blockdiag"` returns the \eqn{p \times p} covariance at each
#' location, which is what marginal variances and the variance of a sum across
#' processes need; `"full"` returns the dense
#' \eqn{(n_{\mathrm{test}}p)^2} covariance.
#'
#' @param object An `ssggp_fit` from [fit()].
#' @param locs_new `n_test` by `d` matrix of prediction locations.
#' @param L Number of random Fourier features for prediction.
#' @param MC Number of Gibbs sweeps; the first `burnin` are discarded.
#' @param burnin Sweeps to discard. Defaults to `MC / 2`.
#' @param g Proposal from [matern_proposal()]. `NULL` inherits the fit's.
#' @param keep_chain Return the raw predictive draws (`n_test * p` by
#'   `MC - burnin`) alongside the summaries.
#' @param summarise One of `"blockdiag"` (mean, per-location `p` by `p`
#'   covariance, marginal variances), `"none"`, or `"full"`.
#' @param threads Thread cap for the duration of the call; see [fit()].
#'   `NULL` leaves every setting alone.
#' @param numerical_guard As in [fit()]; `NULL` inherits the fit's.
#' @param seed Optional integer seed.
#' @param ... Unused.
#' @return A list with `locs`, the requested summaries, and -- when
#'   `keep_chain` is `TRUE` -- `draws` (`n_test * p` by `MC - burnin`, process
#'   index varying fastest).
#' @export
predict.ssggp_fit <- function(object, locs_new,
                              L = 1000, MC = 1000, burnin = NULL,
                              g = NULL,
                              keep_chain = FALSE,
                              summarise = c("blockdiag", "none", "full"),
                              threads = getOption("SSGGP.threads", 5L),
                              numerical_guard = NULL,
                              seed = NULL, ...) {
  summarise <- match.arg(summarise)
  with_threads(resolve_threads(threads),
               predict_impl(object, locs_new, L, MC, burnin, g,
                            keep_chain, summarise, seed, numerical_guard))
}

#' @keywords internal
#' @noRd
predict_impl <- function(object, locs_new, L, MC, burnin, g,
                         keep_chain, summarise, seed,
                         numerical_guard = NULL) {
  sp <- object$spec
  P  <- sp$P
  d  <- sp$d
  if (is.null(g)) g <- sp$g
  check_proposal(g)
  if (is.null(burnin)) burnin <- MC %/% 2
  stopifnot(burnin >= 0, burnin < MC, L >= 1)
  if (summarise == "none" && !keep_chain) {
    stop("summarise = \"none\" with keep_chain = FALSE would return neither ",
         "draws nor summaries. Set keep_chain = TRUE to get the raw draws, or ",
         "choose a `summarise` other than \"none\".", call. = FALSE)
  }
  if (summarise != "none" && MC - burnin < 2L) {
    stop("only ", MC - burnin, " draw", if (MC - burnin == 1) "" else "s",
         " would be retained after burn-in, which is too few for a covariance. ",
         "Increase MC, lower burnin, or use summarise = \"none\".",
         call. = FALSE)
  }

  locs_new <- as.matrix(locs_new)
  if (ncol(locs_new) != d) {
    stop("locs_new must have ", d, " columns, matching the fitted locations",
         call. = FALSE)
  }
  Y    <- object$data$Y
  locs <- object$data$locs
  n      <- nrow(locs)
  n_test <- nrow(locs_new)

  # observed training rows, process index varying fastest
  Y_reordered <- as.vector(t(Y))
  obs_idx <- which(!is.na(Y_reordered))
  Y_obs   <- Y_reordered[obs_idx]
  proc_id <- rep(seq_len(P), n)[obs_idx]

  nug <- object$nugget
  if (length(nug) == 1L) nug <- rep(nug, P)
  d_obs  <- nug[proc_id]
  d_test <- nug[rep(seq_len(P), n_test)]

  if (!is.null(seed)) set.seed(seed)

  # fresh features over the combined grid: training rows first, then test.
  # SSGGP.force_general selects the sparse reference branch; see fit().
  feat <- draw_features(L, d, g)
  force_general <- isTRUE(getOption("SSGGP.force_general", FALSE))
  observed <- !is.na(Y)
  use_fast <- all(observed)
  ff   <- rff(rbind(locs, locs_new), feat$X, feat$phi, P, need_A = force_general)
  A_obs  <- if (force_general) ff$A[obs_idx, , drop = FALSE]
            else as_dgc(matrix(0, 0, 0))
  A_test <- if (force_general)
              ff$A[(n * P + 1):((n + n_test) * P), , drop = FALSE]
            else as_dgc(matrix(0, 0, 0))
  A_l_list <- if (force_general) split_col_blocks(A_obs, L, P)
              else vector("list", L)

  # dense branch, as in fit(); only the training design can be misaligned
  c_obs_d  <- if (force_general) matrix(numeric(0), 0, 0)
              else as.matrix(ff$c_vals[seq_len(n), , drop = FALSE])
  c_test_d <- if (force_general) matrix(numeric(0), 0, 0)
              else as.matrix(ff$c_vals[n + seq_len(n_test), , drop = FALSE])
  col_ss   <- if (!force_general && use_fast) Matrix::colSums(
                ff$c_vals[seq_len(n), , drop = FALSE]^2) else numeric(0)
  col_ss_mat <- if (!force_general && !use_fast)
                  masked_col_ss(c_obs_d, observed) else matrix(numeric(0), 0, 0)
  obs_mask <- if (force_general || use_fast) matrix(numeric(0), 0, 0)
              else t(observed) * 1
  Y_gibbs  <- if (force_general || use_fast) Y_obs else {
    v <- Y_reordered; v[is.na(v)] <- 0; v
  }

  gvals <- vapply(feat$X, g$d, numeric(1))
  M <- build_M_for_fit(object, feat$X, gvals, g, numerical_guard)
  updated <- update_matrices_cpp(M, as_dgc(object$sigma_inv), d_obs, A_l_list,
                                 col_ss, col_ss_mat, nug)

  # the sampler keeps only the current and previous state, started at W = 0
  W <- matrix(0, nrow = P * L, ncol = 2L)
  out <- gibbs_W_Y_cpp(W, Y_gibbs, updated, d_obs,
                       as.integer(MC), as.integer(L),
                       A_obs, A_test, A_l_list, c_obs_d, c_test_d,
                       obs_mask, col_ss_mat, nug)

  draws <- out$Y_pred[, (burnin + 1):MC, drop = FALSE]

  res <- list(locs = locs_new, burnin = burnin, MC = MC, L = L)
  if (keep_chain) res$draws <- draws

  if (summarise != "none") {
    res$mean <- matrix(rowMeans(draws), nrow = n_test, ncol = P, byrow = TRUE)
    if (summarise == "blockdiag") {
      idx <- matrix(seq_len(n_test * P), nrow = P)
      cb  <- array(NA_real_, dim = c(n_test, P, P))
      vr  <- matrix(NA_real_, n_test, P)
      for (i in seq_len(n_test)) {
        block <- stats::cov(t(draws[idx[, i], , drop = FALSE]))
        diag(block) <- diag(block) + nug
        cb[i, , ] <- block
        vr[i, ] <- diag(block)
      }
      res$cov_block <- cb
      res$var <- vr
    } else {
      res$cov <- stats::cov(t(draws))
      diag(res$cov) <- diag(res$cov) + d_test
      res$var <- matrix(diag(res$cov), nrow = n_test, ncol = P, byrow = TRUE)
    }
  }
  structure(res, class = "ssggp_prediction")
}

#' Variance of a sum across processes
#'
#' The variance of \eqn{\sum_j Y_j(s)} at each prediction location, obtained
#' from the per-location covariance blocks. This is the quantity that makes
#' cross-process dependence matter for uncertainty quantification.
#'
#' @param x An `ssggp_prediction` made with `summarise = "blockdiag"`.
#' @return A numeric vector of length `n_test`.
#' @export
var_sum <- function(x) {
  if (is.null(x$cov_block)) {
    stop("needs a prediction made with summarise = \"blockdiag\"", call. = FALSE)
  }
  apply(x$cov_block, 1, sum)
}

#' @export
print.ssggp_prediction <- function(x, ...) {
  cat("<ssggp_prediction>\n")
  cat("  locations :", nrow(x$locs), "\n")
  cat("  draws kept:", if (is.null(x$draws)) "none (not retained)"
                       else ncol(x$draws), "of", x$MC,
      paste0("(burn-in ", x$burnin, ")"), "\n")
  cat("  features  :", x$L, "\n")
  if (!is.null(x$mean)) cat("  summaries :",
      paste(intersect(c("mean", "cov_block", "var", "cov"), names(x)),
            collapse = ", "), "\n")
  invisible(x)
}
