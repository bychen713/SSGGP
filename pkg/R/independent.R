#' Fit each process independently
#'
#' The baseline that ignores cross-process dependence: each marginal process is
#' fitted on its own as a univariate Gaussian process with the same spectral
#' family, so \eqn{\sigma} is diagonal by construction. This is the "Ind"
#' comparison in the paper.
#'
#' Unlike [fit()], the likelihood is exact: the full \eqn{n\times n} covariance
#' is formed and factorized per process, at \eqn{O(n^3)} each. Affordable
#' because the processes are independent, but it does not scale the way [fit()]
#' does.
#'
#' @section Two interfaces, one for each parameterization:
#' The legacy interface is Matern-only, with
#' `par_init = (1/variance, a, nu, nugget)` -- a PRECISION first and an INVERSE
#' range second.
#'
#' Supplying `f_par` selects the general interface instead, which mirrors
#' [fit()]: `f_fun` names a family per process and the parameter vector is the
#' registry's `theta` = **variance, covparms..., nugget**, i.e. a variance and a
#' RANGE. Bounds then default to each family's own, `fixed` names parts rather
#' than indices, and processes may carry different families.
#'
#' @param Y `n` by `P` matrix of observations; `NA` where unobserved.
#' @param locs `n` by `d` matrix of locations.
#' @param f_fun Spectral family: a registered name, or a list of length `P`.
#'   Only consulted when `f_par` is supplied.
#' @param f_par Named starting values in the family's own parameters, recycled
#'   to all processes, or a list of `P` such vectors. Supplying this selects the
#'   general interface; leaving it `NULL` keeps the legacy Matern one.
#' @param variance_init,nugget_init Starting variance and nugget for the general
#'   interface. `NULL` variance starts each process at its own `var(Y[, k])`.
#' @param par_init Legacy interface only: length-4 start value
#'   `(1/variance, a, nu, nugget)`. `NA` in the first slot starts each process at
#'   its own `1/var(Y[, k])`.
#' @param fixed Parameters held fixed. Integer indices of `par_init` under the
#'   legacy interface; a character subset of `"variance"`, `"covparms"`,
#'   `"nugget"` under the general one, since an index means different things
#'   for families with different parameter counts.
#' @param lower,upper Log-scale box constraints. Under the legacy interface,
#'   length 4 in the **same** `(1/variance, a, nu, nugget)` order as `par_init`
#'   -- bounds written for `(variance, range, nu, nugget)` must have their first
#'   two inverted and swapped. Under the general interface they are in the
#'   registry's `theta` order and `NULL` takes each family's own, so no
#'   inversion arises.
#' @param threads Thread cap for the duration of the call; see [fit()].
#' @param verbose Print optimizer progress.
#' @return An object of class `ssggp_independent`, shaped like an [fit()]
#'   result so the two can be handled by the same code: `f_par` (a list of named
#'   parameter vectors, one per process), `variance`, `nugget`, `sigma_inv`
#'   (diagonal), `logLik`, `spec` and `call`. `f_par` is a per-process list
#'   rather than a `P` by 4 matrix, since processes may carry different families
#'   with different parameter counts.
#' @export
fit_independent <- function(Y, locs,
                            f_fun = "matern",
                            f_par = NULL,
                            variance_init = NULL,
                            nugget_init = 0.07,
                            par_init = c(1, 2, 1, 0.07),
                            fixed = NULL,
                            lower = NULL,
                            upper = NULL,
                            threads = getOption("SSGGP.threads", 5L),
                            verbose = FALSE) {
  with_threads(resolve_threads(threads),
               # Legacy defaults in (1/variance, a, nu, nugget) order: the Matern
               # bounds variance 0.001-100, range 0.01-5 (a 0.2-100), nu <= 7.
               if (is.null(f_par))
                 fit_independent_matern(Y, locs, par_init, fixed,
                                        if (is.null(lower))
                                          log(c(0.01, 0.2, 0.01, 0.001)) else lower,
                                        if (is.null(upper))
                                          log(c(1000, 100, 7, 100)) else upper,
                                        verbose)
               else
                 fit_independent_family(Y, locs, f_fun, f_par, variance_init,
                                        nugget_init, fixed, lower, upper,
                                        verbose))
}

#' The legacy Matern-only path
#' @keywords internal
#' @noRd
fit_independent_matern <- function(Y, locs, par_init, fixed, lower, upper,
                                   verbose) {
  Y <- as.matrix(Y); locs <- as.matrix(locs)
  P <- ncol(Y); d <- ncol(locs)
  stopifnot(length(par_init) == 4, nrow(locs) == nrow(Y))

  f_par    <- vector("list", P)
  variance <- numeric(P)
  nug      <- numeric(P)
  lls      <- numeric(P)

  for (k in seq_len(P)) {
    keep <- !is.na(Y[, k])
    z    <- Y[keep, k]
    # NA in the precision slot starts at the process's empirical precision
    pk <- par_init
    if (is.na(pk[1])) pk[1] <- 1 / stats::var(z)
    r <- univariate_mle(z, fields::rdist(locs[keep, , drop = FALSE]),
                        pk, fixed, lower, upper, d, verbose)
    f_par[[k]]  <- c(a = unname(r$par[2]), nu = unname(r$par[3]))
    variance[k] <- 1 / r$par[1]          # par_init[1] is a precision
    nug[k]      <- r$par[4]
    lls[k]      <- r$loglik
  }

  structure(list(
    f_par     = f_par,
    variance  = variance,
    nugget    = nug,
    sigma_inv = diag(1 / variance, P),
    logLik    = sum(lls),
    loglik_by_process = lls,
    call      = match.call(),
    spec      = list(P = P, d = d, family = rep("matern", P), fixed = fixed,
                     lower = lower, upper = upper),
    data      = list(Y = Y, locs = locs)),
    class = "ssggp_independent")
}

#' The general path: one registered family per process
#' @keywords internal
#' @noRd
fit_independent_family <- function(Y, locs, f_fun, f_par, variance_init,
                                   nugget_init, fixed, lower, upper, verbose) {
  Y <- as.matrix(Y); locs <- as.matrix(locs)
  P <- ncol(Y); d <- ncol(locs)
  stopifnot(nrow(locs) == nrow(Y))

  fam <- resolve_spectral(f_fun, P)
  require_covariance(fam)

  if (!is.list(f_par)) f_par <- rep(list(f_par), P)
  stopifnot(length(f_par) == P)
  if (!is.null(fixed) && !is.character(fixed)) {
    stop("under the general interface `fixed` must be a character subset of ",
         "'variance', 'covparms', 'nugget'; an index means different things ",
         "to different families", call. = FALSE)
  }

  out_par  <- vector("list", P)
  variance <- numeric(P)
  nug      <- numeric(P)
  lls      <- numeric(P)

  for (k in seq_len(P)) {
    e    <- fam[[k]]
    keep <- !is.na(Y[, k])
    z    <- Y[keep, k]

    need <- e$par_names
    if (!is.null(need) && !all(need %in% names(f_par[[k]]))) {
      stop("f_par[[", k, "]] must be named (", paste(need, collapse = ", "),
           ") for family '", e$name, "'", call. = FALSE)
    }

    v0 <- if (is.null(variance_init)) stats::var(z) else
      variance_init[[min(k, length(variance_init))]]
    theta_ini <- e$to_theta(f_par[[k]], v0, nugget_init)

    kk <- length(e$theta_names)
    lo <- if (is.null(lower)) e$theta_lower else lower
    up <- if (is.null(upper)) e$theta_upper else upper
    if (length(lo) != kk || length(up) != kk) {
      stop("process ", k, " uses family '", e$name, "', whose parameters are (",
           paste(e$theta_names, collapse = ", "), "), so `lower` and ",
           "`upper` must have length ", kk, "; got ", length(lo), " and ",
           length(up), ". Leave them NULL to use each family's own.",
           call. = FALSE)
    }

    r <- univariate_mle_theta(z, fields::rdist(locs[keep, , drop = FALSE]),
                              theta_ini, fixed, lo, up, e, verbose)
    out_par[[k]] <- e$from_theta(r$par)
    variance[k]  <- unname(r$par[1])
    nug[k]       <- unname(r$par[length(r$par)])
    lls[k]       <- r$loglik
  }

  structure(list(
    f_par     = out_par,
    variance  = variance,
    nugget    = nug,
    sigma_inv = diag(1 / variance, P),
    logLik    = sum(lls),
    loglik_by_process = lls,
    call      = match.call(),
    spec      = list(P = P, d = d,
                     family = vapply(fam, function(e) e$name, character(1)),
                     fixed = fixed, lower = lower, upper = upper),
    data      = list(Y = Y, locs = locs)),
    class = "ssggp_independent")
}

#' The marginal covariance of one registered family
#'
#' `theta` is the registry's layout -- variance, covparms..., nugget -- and the
#' family's `cov_fun` supplies the matrix.
#'
#' @param entry A registry entry.
#' @param dist_mat Distance matrix.
#' @param theta Parameter vector in registry order.
#' @param nugget Add the nugget to the diagonal.
#' @keywords internal
#' @noRd
family_cov <- function(entry, dist_mat, theta, nugget = TRUE) {
  cm <- if (!is.null(entry$cov_fun)) {
    entry$cov_fun(dist_mat, theta)
  } else {
    stop("family '", entry$name, "' has no cov_fun and no covariance this ",
         "baseline can build; register it with a `cov_fun`", call. = FALSE)
  }
  if (nugget) diag(cm) <- diag(cm) + theta[length(theta)]
  cm
}

#' Exact univariate MLE for one process, in the registry parameterization
#' @keywords internal
#' @noRd
univariate_mle_theta <- function(z, dist_mat, theta_ini, fixed, lower, upper,
                                 entry, verbose = FALSE) {
  n <- length(z)
  k <- length(theta_ini)
  idx <- fixed_index(fixed, k)

  if (length(idx) == k) {
    names(theta_ini) <- entry$theta_names
    return(list(par = theta_ini, loglik = NA_real_))
  }

  fix_par <- NULL
  if (length(idx)) {
    fix_par   <- theta_ini[idx]
    theta_ini <- theta_ini[-idx]
    lower     <- lower[-idx]
    upper     <- upper[-idx]
  }
  expand <- function(lg) {
    if (length(idx)) {
      th <- rep(1, k); th[idx] <- fix_par; th[-idx] <- exp(lg); th
    } else exp(lg)
  }

  nll <- function(lg) {
    th <- expand(lg)
    val <- tryCatch({
      cm <- family_cov(entry, dist_mat, th)
      ch <- chol(cm)
      -n / 2 * log(2 * pi) - sum(log(diag(ch))) -
        0.5 * sum(backsolve(ch, z, transpose = TRUE)^2)
    }, error = function(e) NaN)
    if (!is.finite(val)) return(5e5)
    -val
  }

  opt <- stats::optim(par = log(theta_ini), fn = nll, method = "L-BFGS-B",
                      lower = lower, upper = upper,
                      control = list(trace = if (verbose) 100 else 0,
                                     maxit = 300,
                                     parscale = rep(1, length(theta_ini))))
  th <- expand(opt$par)
  names(th) <- entry$theta_names
  list(par = th, loglik = -opt$value)
}

#' Exact univariate Matern MLE for one process
#' @keywords internal
#' @noRd
univariate_mle <- function(z, dist_mat, par_init, fixed, lower, upper, d,
                           verbose = FALSE) {
  n <- length(z)
  fix_par <- NULL
  if (!is.null(fixed)) {
    fix_par  <- par_init[fixed]
    par_init <- par_init[-fixed]
    lower    <- lower[-fixed]
    upper    <- upper[-fixed]
  }
  expand <- function(lg) {
    if (!is.null(fixed)) {
      th <- rep(1, 4); th[fixed] <- fix_par; th[-fixed] <- exp(lg); th
    } else exp(lg)
  }

  nll <- function(lg) {
    th <- expand(lg)
    # th[1] is a PRECISION, matching the original parameterization
    cm <- matern_cov(dist_mat, a = th[2], nu = th[3], sigma = 1 / th[1],
                     nugget = th[4], d = d)
    val <- tryCatch({
      ch <- chol(cm)
      -n / 2 * log(2 * pi) - sum(log(diag(ch))) -
        0.5 * sum(backsolve(ch, z, transpose = TRUE)^2)
    }, error = function(e) NaN)
    if (!is.finite(val)) return(5e5)
    -val
  }

  opt <- stats::optim(par = log(par_init), fn = nll, method = "L-BFGS-B",
                      lower = lower, upper = upper,
                      control = list(trace = if (verbose) 100 else 0,
                                     maxit = 300,
                                     parscale = rep(1, length(par_init))))
  list(par = expand(opt$par), loglik = -opt$value)
}

#' Univariate Matern covariance matrix
#'
#' The marginal case of the model's cross-covariance: with a unit-integral
#' spectral density the marginal variance is exactly `sigma`.
#'
#' @param dist_mat Distance matrix.
#' @param a Inverse range.
#' @param nu Smoothness.
#' @param sigma Marginal variance.
#' @param nugget Added to the diagonal.
#' @param d Spatial dimension.
#' @return An `n` by `n` covariance matrix.
#' @keywords internal
#' @noRd
matern_cov <- function(dist_mat, a, nu, sigma, nugget = 0, d = 2) {
  cm <- sigma * fields::Matern(d = dist_mat, smoothness = nu, alpha = a)
  diag(cm) <- diag(cm) + nugget
  cm
}

#' Translate the named parts of `theta` into indices
#'
#' `theta` is laid out as variance, covariance parameters..., nugget.
#'
#' @param fixed Character subset of `"variance"`, `"covparms"`, `"nugget"`.
#' @param k Length of `theta`.
#' @return Sorted integer indices, possibly empty.
#' @keywords internal
#' @noRd
fixed_index <- function(fixed, k) {
  if (is.null(fixed) || !length(fixed)) return(integer(0))
  parts <- c("variance", "covparms", "nugget")
  if (!is.character(fixed)) {
    stop("`fixed` must be a character subset of ",
         paste(sQuote(parts), collapse = ", "), ", not an index vector",
         call. = FALSE)
  }
  bad <- setdiff(fixed, parts)
  if (length(bad)) {
    stop("unknown `fixed` part", if (length(bad) > 1) "s" else "", ": ",
         paste(sQuote(bad), collapse = ", "), ". Valid: ",
         paste(sQuote(parts), collapse = ", "), call. = FALSE)
  }
  idx <- integer(0)
  if ("variance" %in% fixed) idx <- c(idx, 1L)
  if ("covparms" %in% fixed && k > 2L) idx <- c(idx, seq.int(2L, k - 1L))
  if ("nugget"   %in% fixed) idx <- c(idx, k)
  sort(unique(idx))
}

#' @export
print.ssggp_independent <- function(x, ...) {
  cat("<ssggp_independent>\n")
  cat("  processes :", x$spec$P, " locations:", nrow(x$data$Y), "\n")
  cat("  logLik    :", format(x$logLik, digits = 6), "\n")
  # one line per process: with mixed families there is no parameter every
  # process carries
  fam <- x$spec$family
  for (k in seq_len(x$spec$P)) {
    cat(sprintf("  [%d] %-9s %s\n", k,
                if (length(fam) >= k) fam[[k]] else "?",
                paste(names(x$f_par[[k]]),
                      format(unname(x$f_par[[k]]), digits = 4),
                      sep = " = ", collapse = "  ")))
  }
  cat("  variance  :", format(x$variance, digits = 4), "\n")
  cat("  nugget    :", format(x$nugget, digits = 4), "\n")
  invisible(x)
}

#' @export
logLik.ssggp_independent <- function(object, ...) {
  structure(object$logLik, class = "logLik",
            # free parameters per process: its spectral parameters, plus the
            # variance and the nugget
            df = sum(vapply(object$f_par, length, integer(1))) +
                 2L * object$spec$P,
            nobs = sum(!is.na(object$data$Y)))
}
