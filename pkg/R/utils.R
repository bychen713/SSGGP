#' Coerce to a general column-compressed sparse matrix
#'
#' The compiled kernels map sparse arguments with Eigen's `MappedSparseMatrix`,
#' which accepts only `dgCMatrix`; a symmetric matrix would otherwise coerce to
#' `dsCMatrix`.
#'
#' @param m A matrix, or anything coercible by the Matrix package.
#' @return A `dgCMatrix`.
#' @keywords internal
#' @noRd
as_dgc <- function(m) {
  if (!methods::is(m, "Matrix")) m <- Matrix::Matrix(m, sparse = TRUE)
  methods::as(methods::as(methods::as(m, "dMatrix"), "generalMatrix"),
              "CsparseMatrix")
}

#' Split a feature matrix into its per-feature column blocks
#'
#' `lapply(seq_len(L), function(l) A[, ((l-1)*P+1):(l*P)])`, slicing the CSC
#' slots directly, which is much faster than `[` at large `L`.
#'
#' @param A Sparse matrix with `L * P` columns.
#' @param L Number of random Fourier features.
#' @param P Number of processes.
#' @return A list of `L` `dgCMatrix` objects, each with `P` columns.
#' @keywords internal
#' @noRd
split_col_blocks <- function(A, L, P) {
  A  <- as_dgc(A)
  P  <- as.integer(P)
  Ap <- A@p; Ai <- A@i; Ax <- A@x
  nr <- A@Dim[1]
  lapply(seq_len(L), function(l) {
    c1 <- (l - 1L) * P + 1L
    c2 <- l * P
    lo <- Ap[c1]
    hi <- Ap[c2 + 1L]
    sel <- if (hi > lo) (lo + 1L):hi else integer(0)
    methods::new("dgCMatrix",
                 p   = Ap[c1:(c2 + 1L)] - lo,
                 i   = Ai[sel],
                 x   = Ax[sel],
                 Dim = c(nr, P))
  })
}

#' Masked column sums of squares of the feature matrix
#'
#' `out[p, l] = sum_{i observed for process p} c_vals[i, l]^2`, the misaligned
#' generalisation of `colSums(c_vals^2)`: \eqn{A_\ell^\top D^{-1} A_\ell} is
#' diagonal with entries `out[p, l] / d_p`. Computed in column blocks to avoid
#' an `n` by `L` temporary.
#'
#' @param c_dense `n` by `L` dense feature matrix.
#' @param observed `n` by `P` logical, `TRUE` where the process is observed.
#' @param block Columns per block.
#' @return A `P` by `L` numeric matrix.
#' @keywords internal
#' @noRd
masked_col_ss <- function(c_dense, observed, block = 4096L) {
  P <- ncol(observed); L <- ncol(c_dense)
  Z   <- observed * 1
  out <- matrix(0, P, L)
  for (s in seq.int(1L, L, by = block)) {
    e  <- min(s + block - 1L, L)
    cc <- c_dense[, s:e, drop = FALSE]
    out[, s:e] <- crossprod(Z, cc * cc)
  }
  out
}

#' Run an expression under a thread cap
#'
#' Caps the BLAS and OpenMP thread counts, restoring them on exit. The
#' environment variables reach libraries initialised later; `RhpcBLASctl`, when
#' installed, changes the running BLAS and OpenMP runtimes.
#'
#' @param n Thread cap, or `NULL` to leave every setting alone.
#' @param expr Expression to evaluate.
#' @keywords internal
#' @noRd
with_threads <- function(n, expr) {
  if (is.null(n)) return(force(expr))
  n <- as.integer(n)

  vars <- c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS", "GOTO_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
  old_env <- Sys.getenv(vars, names = TRUE, unset = NA_character_)
  do.call(Sys.setenv, stats::setNames(as.list(rep(as.character(n), length(vars))),
                                      vars))
  on.exit({
    keep <- old_env[!is.na(old_env)]
    if (length(keep)) do.call(Sys.setenv, as.list(keep))
    drop <- names(old_env)[is.na(old_env)]
    if (length(drop)) Sys.unsetenv(drop)
  }, add = TRUE)

  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NULL)
    old_omp  <- tryCatch(RhpcBLASctl::omp_get_max_threads(), error = function(e) NULL)
    try(RhpcBLASctl::blas_set_num_threads(n), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(n), silent = TRUE)
    on.exit({
      if (!is.null(old_blas)) try(RhpcBLASctl::blas_set_num_threads(old_blas), silent = TRUE)
      if (!is.null(old_omp))  try(RhpcBLASctl::omp_set_num_threads(old_omp), silent = TRUE)
    }, add = TRUE)
  }
  force(expr)
}

#' Resolve the thread cap
#'
#' `NULL` means leave the BLAS setting untouched. Anything else must be a
#' positive whole number.
#' @keywords internal
#' @noRd
resolve_threads <- function(threads) {
  if (is.null(threads)) return(NULL)
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
      threads < 1 || threads != as.integer(threads)) {
    stop("`threads` must be a single positive whole number, or NULL to leave ",
         "the BLAS setting alone", call. = FALSE)
  }
  as.integer(threads)
}

#' Set the SSGGP thread budget for the session
#'
#' Convenience wrapper for `options(SSGGP.threads = n)`, the default `threads`
#' of every entry point.
#'
#' @param n Threads, or `NULL` to leave thread settings untouched.
#' @return The previous value, invisibly.
#' @export
#' @examples
#' old <- ssggp_threads(2)
#' ssggp_threads(old)
ssggp_threads <- function(n) {
  old <- getOption("SSGGP.threads", 5L)
  options(SSGGP.threads = if (is.null(n)) NULL else resolve_threads(n))
  invisible(old)
}

#' Validate the data and settings before any work is done
#'
#' @param Y,locs,L,MC,burnin,iter,rho,nugget_init,lag,cm_cycles,tol As in
#'   [fit()].
#' @return Invisibly `TRUE`.
#' @keywords internal
#' @noRd
validate_fit_inputs <- function(Y, locs, L, MC, burnin, iter, rho,
                                nugget_init, lag, cm_cycles, tol) {
  n <- nrow(Y); P <- ncol(Y)

  if (n < 1L || P < 1L) stop("Y must have at least one row and one column",
                             call. = FALSE)
  if (nrow(locs) != n) {
    stop("locs has ", nrow(locs), " rows but Y has ", n, call. = FALSE)
  }
  if (!is.numeric(Y)) stop("Y must be numeric", call. = FALSE)
  if (!is.numeric(locs)) stop("locs must be numeric", call. = FALSE)
  if (any(!is.finite(locs))) {
    stop("locs contains NA, NaN or Inf; locations must all be finite",
         call. = FALSE)
  }
  # NA in Y is meaningful (unobserved); Inf and NaN are not
  bad <- is.infinite(Y) | is.nan(Y)
  if (any(bad)) {
    stop("Y contains ", sum(bad), " non-finite value(s). NA marks an ",
         "unobserved process, but Inf and NaN are not supported", call. = FALSE)
  }

  thin <- which(colSums(!is.na(Y)) < 2L)
  if (length(thin)) {
    stop("process", if (length(thin) > 1) "es" else "", " ",
         paste(thin, collapse = ", "), " ha",
         if (length(thin) > 1) "ve" else "s", " fewer than 2 observations",
         call. = FALSE)
  }

  pos_int <- function(x, nm, min = 1L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < min || x != as.integer(x)) {
      stop("`", nm, "` must be a single whole number >= ", min, call. = FALSE)
    }
  }
  pos_int(L, "L"); pos_int(MC, "MC"); pos_int(iter, "iter")
  pos_int(cm_cycles, "cm_cycles")
  if (!is.null(burnin)) {
    pos_int(burnin, "burnin", min = 0L)
    if (burnin >= MC) {
      stop("`burnin` (", burnin, ") must be smaller than `MC` (", MC, ")",
           call. = FALSE)
    }
  }
  if (MC - (if (is.null(burnin)) MC %/% 2 else burnin) < 2L) {
    stop("fewer than 2 draws would be retained after burn-in", call. = FALSE)
  }

  pos_num <- function(x, nm, strict = TRUE) {
    ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
      (if (strict) x > 0 else x >= 0)
    if (!ok) {
      stop("`", nm, "` must be a single finite number ",
           if (strict) "> 0" else ">= 0", call. = FALSE)
    }
  }
  pos_num(nugget_init, "nugget_init")
  pos_num(lag, "lag")
  pos_num(tol, "tol")
  pos_num(rho, "rho", strict = FALSE)
  invisible(TRUE)
}
