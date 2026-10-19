## Spectral family registry ------------------------------------------------
##
## Maps a NAME to a spectral density, its vectorised form and the box CM3
## optimises in, plus -- for the independent baseline in independent.R -- a
## covariance function and the map between the spectral parameters and the
## baseline's parameter vector.
##
## Two parameter vectors live here:
##
## * `par_names` / `par_lower` / `par_upper` describe iota_j, the family's own
##   parameters -- (a, nu) for Matern. This is what fit() estimates. There is
##   no variance slot: every marginal scale belongs to sigma, and the density
##   integrates to one.
## * `theta_names` / `theta_lower` / `theta_upper` describe theta, the vector
##   fit_independent() estimates: variance, covariance parameters..., nugget.

.spectral_registry <- new.env(parent = emptyenv())

#' Register a spectral density family
#'
#' Makes a spectral density usable throughout the package.
#'
#' [fit()], [sim()] and [predict.ssggp_fit()] need only `f` (better, `f_vec`),
#' `par_names`, and `par_lower`/`par_upper`. The remaining arguments are needed
#' only to fit the family with [fit_independent()]; leave them `NULL` otherwise.
#'
#' @param name Name to register under.
#' @param f Spectral density, a function of `(x, par)` returning a single
#'   non-negative number. **Must integrate to one** over \eqn{\mathbb{R}^d};
#'   see [check_unit_integral()].
#' @param par_names Character vector of spectral parameter names.
#' @param f_vec Optional vectorised form, `function(X, par)` taking the `L` by
#'   `d` matrix of frequencies and returning `L` densities. Strongly
#'   recommended, since CM3 evaluates the density at every frequency for every
#'   trial point; without it `f` is looped over.
#' @param par_lower,par_upper Box constraints on the **log** scale for the
#'   spectral parameters, in `par_names` order, used by [fit()]. `NULL`
#'   registers a density that only simulation and prediction can use.
#' @param cov_fun `function(dist_mat, theta)` returning the covariance matrix
#'   without the nugget, for [fit_independent()].
#' @param theta_names Names of the [fit_independent()] parameter vector
#'   `theta`, laid out as **variance, covariance parameters..., nugget**.
#' @param to_theta Function `(par, variance, nugget)` returning `theta`.
#' @param from_theta Function `(theta)` returning the named spectral
#'   parameters.
#' @param theta_lower,theta_upper Box constraints on the **log** scale for
#'   `theta`, the same length as `theta_names`.
#' @param overwrite Allow replacing an existing entry.
#' @return Invisibly, the registered entry.
#' @export
#' @examples
#' list_spectral()
register_spectral <- function(name, f, par_names, f_vec = NULL,
                              par_lower = NULL, par_upper = NULL,
                              cov_fun = NULL, theta_names = NULL,
                              to_theta = NULL, from_theta = NULL,
                              theta_lower = NULL, theta_upper = NULL,
                              overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, is.function(f))
  if (!overwrite && exists(name, envir = .spectral_registry, inherits = FALSE)) {
    stop("spectral family '", name, "' is already registered; ",
         "pass overwrite = TRUE to replace it", call. = FALSE)
  }
  if (!is.null(f_vec) && !is.function(f_vec)) {
    stop("`f_vec` must be a function(X, par) or NULL", call. = FALSE)
  }
  if (!is.null(par_lower) || !is.null(par_upper)) {
    np <- length(par_names)
    if (length(par_lower) != np || length(par_upper) != np) {
      stop("`par_lower` and `par_upper` must have length ", np, ", matching ",
           "`par_names` (", paste(par_names, collapse = ", "), ")",
           call. = FALSE)
    }
  }
  if (!is.null(cov_fun)) {
    if (!is.function(cov_fun) || !is.function(to_theta) ||
        !is.function(from_theta)) {
      stop("`cov_fun`, `to_theta` and `from_theta` must all be functions",
           call. = FALSE)
    }
    k <- length(theta_names)
    if (k < 2L) {
      stop("`theta_names` needs at least a variance and a nugget",
           call. = FALSE)
    }
    if (length(theta_lower) != k || length(theta_upper) != k) {
      stop("`theta_lower` and `theta_upper` must have length ", k,
           ", matching `theta_names` (", paste(theta_names, collapse = ", "),
           ")", call. = FALSE)
    }
  }
  entry <- list(name = name, f = f, f_vec = f_vec, par_names = par_names,
                par_lower = par_lower, par_upper = par_upper,
                cov_fun = cov_fun, theta_names = theta_names,
                to_theta = to_theta, from_theta = from_theta,
                theta_lower = theta_lower, theta_upper = theta_upper)
  assign(name, entry, envir = .spectral_registry)
  invisible(entry)
}

#' List registered spectral families
#'
#' @return Character vector of registered names.
#' @export
list_spectral <- function() sort(ls(.spectral_registry))

#' Look up a registered spectral family
#'
#' @param name Registered name.
#' @return The registry entry.
#' @keywords internal
#' @noRd
get_spectral <- function(name) {
  if (!exists(name, envir = .spectral_registry, inherits = FALSE)) {
    stop("unknown spectral family '", name, "'. Registered: ",
         paste(list_spectral(), collapse = ", "),
         ". Use register_spectral() to add one.", call. = FALSE)
  }
  get(name, envir = .spectral_registry)
}

#' Resolve a user `f_fun` specification into P registry-shaped entries
#'
#' Accepts a single name or function recycled to all processes, or a list of
#' length `P` mixing the two. A raw function becomes an entry carrying only the
#' density, which simulation and prediction accept and estimation rejects.
#'
#' @param f_fun Name, function, or list of length `P`.
#' @param P Number of processes.
#' @return A list of `P` entries.
#' @keywords internal
#' @noRd
resolve_spectral <- function(f_fun, P) {
  if (is.character(f_fun) && length(f_fun) > 1L) f_fun <- as.list(f_fun)
  if (!is.list(f_fun)) f_fun <- rep(list(f_fun), P)
  if (length(f_fun) != P) {
    stop("f_fun must be a single family or a list of length P = ", P,
         call. = FALSE)
  }
  lapply(seq_along(f_fun), function(j) {
    el <- f_fun[[j]]
    if (is.character(el)) {
      get_spectral(el)
    } else if (is.function(el)) {
      list(name = NA_character_, f = el)
    } else {
      stop("f_fun[[", j, "]] must be a registered name or a function",
           call. = FALSE)
    }
  })
}

#' Require a covariance implementation for every family
#'
#' @param entries List of resolved entries.
#' @return Invisibly `TRUE`, or an error naming the offending processes.
#' @keywords internal
#' @noRd
require_covariance <- function(entries) {
  bad <- which(vapply(entries, function(e) is.null(e$cov_fun), logical(1)))
  if (length(bad)) {
    stop("f_fun for process", if (length(bad) > 1) "es" else "", " ",
         paste(bad, collapse = ", "),
         " has no covariance function. Register the family with a `cov_fun` ",
         "in register_spectral(), or pass a registered name (", paste(list_spectral(), collapse = ", "), ").", call. = FALSE)
  }
  invisible(TRUE)
}

#' Require parameter names for every family, as CM3 needs
#'
#' @param entries List of resolved entries.
#' @return Invisibly `TRUE`, or an error naming the offending processes.
#' @keywords internal
#' @noRd
require_spectral_estimable <- function(entries) {
  bad <- which(vapply(entries, function(e) is.null(e$par_names), logical(1)))
  if (length(bad)) {
    stop("f_fun for process", if (length(bad) > 1) "es" else "", " ",
         paste(bad, collapse = ", "),
         " supplied as a raw function, which carries no parameter names, so ",
         "CM3 cannot tell which parameters to estimate. Register the family ",
         "with register_spectral(), or pass a registered name (",
         paste(list_spectral(), collapse = ", "), ").", call. = FALSE)
  }
  invisible(TRUE)
}

#' The vectorised density for one entry, looping over `f` if there is none
#'
#' @param e A resolved registry entry.
#' @return A function of `(X, par)` returning a length-`nrow(X)` vector.
#' @keywords internal
#' @noRd
f_vec_of <- function(e) {
  if (!is.null(e$f_vec)) return(e$f_vec)
  fj <- e$f
  force(fj)
  function(X, par) vapply(seq_len(nrow(X)), function(l) fj(X[l, ], par),
                          numeric(1))
}

## Compiled CM3 objectives ---------------------------------------------------
##
## The built-in families have a fused objective in cm3.cpp. Dispatch is by
## name, but `register_spectral(overwrite = TRUE)` can rebind a built-in name,
## so the entry must also still carry the built-in vectorised density and
## parameter order.
.cm3_compiled <- list(
  matern   = list(id = 0L, par_names = c("a", "nu"),           grad = TRUE),
  se       = list(id = 1L, par_names = "a",                    grad = TRUE),
  ch       = list(id = 2L, par_names = c("beta", "nu", "alpha"), grad = FALSE))

#' The compiled CM3 objective for a registry entry, if there is one
#'
#' @param e A resolved registry entry.
#' @return `list(id, par_names, grad)`, or `NULL` for the R evaluator.
#' @keywords internal
#' @noRd
cm3_kernel_of <- function(e) {
  k <- .cm3_compiled[[e$name]]
  if (is.null(k)) return(NULL)
  if (!identical(as.character(e$par_names), k$par_names)) return(NULL)
  builtin <- switch(e$name,
                    matern   = matern_spectral_vec,
                    se       = se_spectral_vec,
                    ch       = ch_spectral_vec)
  if (!identical(e$f_vec, builtin)) return(NULL)
  k
}

#' Analytic log-density for a built-in family
#'
#' Lets the importance ratio be formed as `log g - log f` where `f` itself
#' would underflow. CH goes through its density, whose polynomial tails do not
#' underflow. User-registered families have none, so `numerical_guard = "cap"`
#' is unavailable for them.
#'
#' @param e A resolved registry entry.
#' @return `function(Xmat, par, r2)` returning `log f`, or `NULL`.
#' @keywords internal
#' @noRd
log_density_of <- function(e) {
  if (is.null(cm3_kernel_of(e))) return(NULL)
  switch(e$name,
    matern = function(Xmat, par, r2) {
      a <- par[["a"]]; nu <- par[["nu"]]; d <- ncol(Xmat)
      2 * nu * log(a) + lgamma(nu + d / 2) - lgamma(nu) - (d / 2) * log(pi) -
        (nu + d / 2) * log(a^2 + r2)
    },
    se = function(Xmat, par, r2) {
      a <- par[["a"]]; d <- ncol(Xmat)
      -d * log(a) - (d / 2) * log(2 * pi) - r2 / (2 * a^2)
    },
    ch = function(Xmat, par, r2) log(ch_spectral_vec(Xmat, par)),
    NULL)
}

## built-in families ---------------------------------------------------------

.onLoad <- function(libname, pkgname) {
  register_spectral(
    name        = "matern",
    f           = matern_spectral,
    f_vec       = matern_spectral_vec,
    par_names   = c("a", "nu"),
    # `a` is an inverse range: a in [0.2, 100] is a range in [0.01, 5]
    par_lower   = log(c(a = 0.2, nu = 0.01)),
    par_upper   = log(c(a = 100, nu = 7)),
    cov_fun     = function(dist_mat, theta) {
      theta[1] * fields::Matern(d = dist_mat, smoothness = theta[3],
                                alpha = 1 / theta[2])
    },
    theta_names = c("variance", "range", "smoothness", "nugget"),
    to_theta    = function(par, variance, nugget) {
      c(variance, 1 / par[["a"]], par[["nu"]], nugget)
    },
    from_theta  = function(theta) c(a = 1 / theta[[2]], nu = theta[[3]]),
    theta_lower = log(c(0.001, 0.01, 0.01, 0.001)),
    theta_upper = log(c(100, 5, 7, 100)),
    overwrite   = TRUE)

  register_spectral(
    name        = "se",
    f           = se_spectral,
    f_vec       = se_spectral_vec,
    par_names   = "a",
    par_lower   = log(c(a = 0.2)),
    par_upper   = log(c(a = 100)),
    cov_fun     = function(dist_mat, theta) {
      theta[1] * exp(-(dist_mat / theta[2])^2 / 2)
    },
    theta_names = c("variance", "range", "nugget"),
    to_theta    = function(par, variance, nugget) {
      c(variance, 1 / par[["a"]], nugget)
    },
    from_theta  = function(theta) c(a = 1 / theta[[2]]),
    theta_lower = log(c(0.001, 0.01, 0.001)),
    # a range ceiling of 1 keeps the squared exponential covariance matrix
    # numerically positive definite in the exact baseline likelihood
    theta_upper = log(c(100, 1, 100)),
    overwrite   = TRUE)

  # CH is parameterised by the range beta itself, not an inverse range. The
  # alpha floor is the model's alpha > d/2, assuming d = 2.
  register_spectral(
    name        = "ch",
    f           = ch_spectral,
    f_vec       = ch_spectral_vec,
    par_names   = c("beta", "nu", "alpha"),
    par_lower   = log(c(beta = 0.01, nu = 0.01, alpha = 1.01)),
    par_upper   = log(c(beta = 5,    nu = 5,    alpha = 20)),
    cov_fun     = function(dist_mat, theta) {
      ch_cov_cpp(dist_mat, beta = theta[2], nu = theta[3], alpha = theta[4],
                 variance = theta[1])
    },
    theta_names = c("variance", "range", "smoothness", "tail", "nugget"),
    to_theta    = function(par, variance, nugget) {
      c(variance, par[["beta"]], par[["nu"]], par[["alpha"]], nugget)
    },
    from_theta  = function(theta) {
      c(beta = theta[[2]], nu = theta[[3]], alpha = theta[[4]])
    },
    theta_lower = log(c(0.001, 0.01, 0.01, 1.01, 0.001)),
    theta_upper = log(c(100, 5, 5, 20, 100)),
    overwrite   = TRUE)
}
