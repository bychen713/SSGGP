# Package autostart adapter for the simulation's known-zero mean.
# GpGp::fit_model cannot use our all-zero X column (singular regression).
# Use an intercept ONLY in the preliminary fits, then move that intercept
# into initial W. This preserves initial X*Beta + W without adding an
# intercept to the IOX sampling model. Y, including its NA mask, is untouched.
iox_latent_autostart <- function(Y, X, coords, m = 15L, seed = 2025L,
                                verbose = TRUE) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  coords <- as.matrix(coords)
  stopifnot(nrow(Y) == nrow(X), nrow(Y) == nrow(coords),
            ncol(X) == 1L, all(X == 0),
            all(colSums(is.finite(Y)) >= 3L))
  Y_before <- Y
  set.seed(seed)
  raw <- spiox::autostart(Y = Y, X = matrix(1, nrow(Y), 1L),
                         coords = coords, method = "latent", m = m,
                         nu = 0.5, verbose = verbose)
  stopifnot(identical(Y, Y_before),
            identical(dim(raw$Theta), c(3L, ncol(Y))),
            identical(dim(raw$W), dim(Y)),
            all(is.finite(unlist(raw))),
            all(diag(raw$Sigma) > 0), all(raw$Ddiag > 0),
            all(raw$Theta[2L, ] == .5), all(raw$Theta[3L, ] == 0))
  starting <- raw
  starting$W <- raw$W + matrix(1, nrow(Y), 1L) %*% raw$Beta
  starting$Beta <- matrix(0, 1L, ncol(Y))
  obs <- is.finite(Y)
  stopifnot(isTRUE(all.equal(unname(starting$W[obs]), unname(Y[obs]),
                            tolerance = 1e-10)),
            all(is.finite(starting$W)), all(X %*% starting$Beta == 0))
  # The package's MCMC range parameter is restricted to (.3, 200).
  # Independent GpGp optimization does not enforce these bounds. If needed,
  # project ONLY the start just inside the existing bounds; never change the
  # prior or silently let a bounded-logit proposal start outside its domain.
  phi_raw <- raw$Theta[1L, ]
  phi_start <- pmax(.3 + 1e-4, pmin(200 - 1e-4, phi_raw))
  adjusted <- which(phi_start != phi_raw)
  starting$Theta[1L, ] <- phi_start
  if (length(adjusted)) {
    cat("Autostart range values moved inside package MCMC bounds for outcomes:",
        paste(adjusted, collapse = ","), "\n")
  }
  list(starting = starting,
       details = list(
         method = "spiox::autostart, nu=.5, with zero-mean adapter",
         seed = seed, raw_autostart = raw,
         initialization_intercept = as.numeric(raw$Beta),
         mean_adapter = "Preliminary intercept moved to W; IOX X and Beta start are zero",
         range_bounds = c(.3, 200), range_adjusted_outcomes = adjusted,
         range_raw = phi_raw, range_start = phi_start,
         missing_latent_starts = sum(is.na(Y)),
         observed_Y_unchanged = identical(Y, Y_before),
         note = "Autostart temporarily fills NA locations only to initialize W. Original Y with NA is passed to latent MCMC; missing responses are not imputed by that sampler."
       ))
}
