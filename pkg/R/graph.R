#' Extract the estimated conditional independence graph
#'
#' Zeros in \eqn{\sigma^{-1}} encode process-level conditional independence,
#' so its sparsity pattern *is* the estimated graph.
#'
#' The pattern is scale invariant: rescaling \eqn{\sigma} by any positive
#' diagonal cannot change which entries are zero.
#'
#' @param object An `ssggp_fit`, or a `p` by `p` precision matrix.
#' @param tol Entries with absolute value at or below this count as zero.
#' @return A logical `p` by `p` adjacency matrix, diagonal `FALSE`.
#' @export
#' @examples
#' graph(trimat(4, 2, 1))
graph <- function(object, tol = 0) {
  si <- if (inherits(object, "ssggp_fit")) object$sigma_inv else as.matrix(object)
  g <- abs(si) > tol
  diag(g) <- FALSE
  g
}

#' Graph recovery metrics
#'
#' Compares an estimated conditional independence graph against the truth,
#' treating edge detection as binary classification over the **off-diagonal**
#' entries. Diagonal entries are excluded: they are non-zero by construction and
#' including them would inflate every rate.
#'
#' Because \eqn{\sigma^{-1}} is symmetric, each edge appears twice; the
#' counts here are over the upper triangle, so an edge is counted once.
#'
#' @param estimate An `ssggp_fit`, or a `p` by `p` precision matrix.
#' @param truth A `p` by `p` precision matrix, or a logical adjacency matrix.
#' @param tol Entries with absolute value at or below this count as zero.
#' @return A named numeric vector: `TP`, `FP`, `FN`, `TN`, `TPR`, `FPR`, `FDR`,
#'   `precision`, `recall`, `F1`, `F2`.
#' @export
#' @examples
#' truth <- trimat(5, 2, 1)
#' graph_metrics(truth, truth)          # perfect recovery
graph_metrics <- function(estimate, truth, tol = 0) {
  est <- graph(estimate, tol)
  tru <- if (is.logical(truth)) truth else graph(truth, tol)
  if (!identical(dim(est), dim(tru))) {
    stop("estimate and truth must have the same dimensions", call. = FALSE)
  }
  ut <- upper.tri(est)                    # each edge once, diagonal excluded
  e <- est[ut]; t <- tru[ut]

  TP <- sum(e & t); FP <- sum(e & !t)
  FN <- sum(!e & t); TN <- sum(!e & !t)

  # a rate is undefined, not zero, when its denominator is empty
  safe <- function(num, den) if (den > 0) num / den else NA_real_
  precision <- safe(TP, TP + FP)
  recall    <- safe(TP, TP + FN)
  f_beta <- function(b) {
    if (is.na(precision) || is.na(recall) || (precision + recall) == 0) return(NA_real_)
    (1 + b^2) * precision * recall / (b^2 * precision + recall)
  }
  c(TP = TP, FP = FP, FN = FN, TN = TN,
    TPR = safe(TP, TP + FN),
    FPR = safe(FP, FP + TN),
    FDR = safe(FP, FP + TP),
    precision = precision, recall = recall,
    F1 = f_beta(1), F2 = f_beta(2))
}

#' Prediction accuracy and calibration
#'
#' Summarises predictive performance against held-out truth, including for the
#' sum across processes -- the quantity that only comes out right if the
#' cross-process dependence is captured.
#'
#' @param prediction An `ssggp_prediction` made with `summarise = "blockdiag"`.
#' @param truth `n_test` by `p` matrix of held-out values. `NA` entries are
#'   excluded; rows with any `NA` are dropped from the summed quantities.
#' @param level Nominal coverage.
#' @return A list with `mse` and `coverage` (per process), `rmse`, `norm` (mean
#'   Euclidean residual), and `sum` (a list of `mse`, `coverage`, `var` for the
#'   sum across processes).
#' @export
prediction_metrics <- function(prediction, truth, level = 0.95) {
  if (is.null(prediction$mean)) {
    stop("needs a prediction with summaries; use summarise = \"blockdiag\"",
         call. = FALSE)
  }
  truth <- as.matrix(truth)
  mu <- prediction$mean
  if (!identical(dim(mu), dim(truth))) {
    stop("truth must be ", nrow(mu), " by ", ncol(mu), call. = FALSE)
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  v <- prediction$var
  resid <- truth - mu

  per_process <- function(f) vapply(seq_len(ncol(truth)), f, numeric(1))
  mse <- per_process(function(j) mean(resid[, j]^2, na.rm = TRUE))
  cov_j <- per_process(function(j) {
    ok <- !is.na(resid[, j])
    mean(abs(resid[ok, j]) <= z * sqrt(v[ok, j]))
  })

  complete <- rowSums(is.na(truth)) == 0
  out <- list(
    mse = mse, coverage = cov_j,
    rmse = sqrt(mean(resid^2, na.rm = TRUE)),
    norm = mean(sqrt(rowSums(resid[complete, , drop = FALSE]^2))))

  if (!is.null(prediction$cov_block) && any(complete)) {
    vs <- var_sum(prediction)[complete]
    rs <- rowSums(truth[complete, , drop = FALSE]) -
          rowSums(mu[complete, , drop = FALSE])
    out$sum <- list(mse = mean(rs^2),
                    coverage = mean(abs(rs) <= z * sqrt(vs)),
                    var = mean(vs))
  }
  out
}
