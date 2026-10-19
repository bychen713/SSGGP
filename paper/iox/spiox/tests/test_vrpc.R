# tests/test_vrpc.R
#
# Smoke / sanity tests for the Vecchia-response preconditioner (VRPC) used in
# the latent_model = 1 (full block Beta-W) sampler.
#
# Run from the package root with:
#   Rscript tests/test_vrpc.R
#
# Exits with non-zero status on first failure.  Uses base R only.

suppressPackageStartupMessages({
  library(spiox)
})

# ---------- helpers ----------

ok <- function(label, expr) {
  cat(sprintf("  [ ] %s ... ", label))
  ans <- tryCatch(isTRUE(expr), error = function(e) {
    cat("ERROR\n"); message(conditionMessage(e)); quit(status = 1L)
  })
  if (!ans) { cat("FAIL\n"); quit(status = 1L) }
  cat("ok\n")
}

simulate_iox <- function(n, q, rho, sd_noise, seed = 1L) {
  set.seed(seed)
  coords     <- cbind(runif(n), runif(n))
  Sigma_true <- matrix(rho, q, q) + diag(1 - rho, q)
  S_true     <- chol(Sigma_true)
  Theta_true <- rbind(phi   = seq(6, 12, length.out = q),
                      nu    = rep(1.5, q),
                      alpha = rep(0,   q))
  W_white <- matrix(rnorm(n*q), n, q)
  W_mix   <- W_white %*% S_true
  W <- matrix(0, n, q)
  for (j in seq_len(q)) {
    K <- spiox::Correlationc(
      coords, coords,
      theta = c(Theta_true["phi", j], 1, Theta_true["nu", j], 0),
      matern = 1, same = TRUE
    )
    L <- chol(K + 1e-6 * diag(n))
    W[, j] <- crossprod(L, W_mix[, j])
  }
  Y <- W + matrix(rnorm(n*q, sd = sd_noise), n, q)
  X <- matrix(1, n, 1)
  list(Y = Y, X = X, coords = coords, W = W,
       Sigma_true = Sigma_true, Theta_true = Theta_true,
       sd_noise = sd_noise)
}

fit_latent_block <- function(d, iter = 40, m = 15) {
  spiox::spiox(
    Y = d$Y, X = d$X, coords = d$coords, m = m,
    method = "latent", fit = "mcmc",
    iter = iter, print_every = 0,
    starting = list(Theta = d$Theta_true,
                    Sigma = d$Sigma_true,
                    W     = d$W,
                    Ddiag = rep(d$sd_noise^2, ncol(d$Y)),
                    Beta  = matrix(0, ncol(d$X), ncol(d$Y))),
    opts = list(update_Theta = c(0, 0, 0),
      cg_preconditioner = "vrpc",   # or "vrpc", same thing
      m_pc = 25, 
      num_threads = 4
    ),
    debug = list(sampling = 1L,
                 sample_Beta = TRUE, sample_Sigma = TRUE, sample_Ddiag = TRUE)
  )
}

# ---------- tests ----------

cat("\n--- spiox VRPC tests (latent_model = 1) ---\n")

cat("\n[1] small problem (n=300, q=2, low cross-corr)\n")
{
  d  <- simulate_iox(n = 500, q = 3, rho = 0.5, sd_noise = 0.3, seed = 11L)
  system.time({
    fit <- fit_latent_block(d, iter = 100, m = 25) })

  ok("fit returns Sigma, W, Ddiag", {
    !is.null(fit$Sigma) && !is.null(fit$W) && !is.null(fit$Ddiag)
  })
  ok("no NA/Inf in Sigma posterior", {
    all(is.finite(fit$Sigma))
  })
  ok("no NA/Inf in W posterior", {
    all(is.finite(fit$W))
  })
  ok("no NA/Inf in Ddiag posterior", {
    all(is.finite(fit$Ddiag))
  })

  # Use last half of chain to give the probe time to lock in
  keep      <- floor(dim(fit$Sigma)[3] / 2) : dim(fit$Sigma)[3]
  Sigma_pm  <- apply(fit$Sigma[, , keep, drop = FALSE], c(1, 2), mean)
  Ddiag_pm  <- rowMeans(fit$Ddiag[, keep, drop = FALSE])

  ok("Sigma diag posterior mean within ±0.6 of truth", {
    max(abs(diag(Sigma_pm) - diag(d$Sigma_true))) < 0.6
  })
  ok("Ddiag posterior mean within ±0.1 of truth", {
    max(abs(Ddiag_pm - d$sd_noise^2)) < 0.1
  })
}

cat("\n[2] moderate problem (n=600, q=3, strong cross-corr)\n")
{
  d  <- simulate_iox(n = 600, q = 3, rho = 0.7, sd_noise = 0.3, seed = 42L)
  fit <- fit_latent_block(d, iter = 40, m = 15)

  ok("fit returns finite Sigma", {
    all(is.finite(fit$Sigma))
  })
  ok("fit returns finite W", {
    all(is.finite(fit$W))
  })

  keep     <- floor(dim(fit$Sigma)[3] / 2) : dim(fit$Sigma)[3]
  Sigma_pm <- apply(fit$Sigma[, , keep, drop = FALSE], c(1, 2), mean)
  Ddiag_pm <- rowMeans(fit$Ddiag[, keep, drop = FALSE])

  ok("Sigma diag within ±0.6 of truth", {
    max(abs(diag(Sigma_pm) - diag(d$Sigma_true))) < 0.6
  })
  ok("Sigma off-diag preserves cross-corr sign", {
    rho_hat <- Sigma_pm[1, 2] / sqrt(Sigma_pm[1, 1] * Sigma_pm[2, 2])
    rho_hat > 0.3
  })
  ok("Ddiag within ±0.1 of truth", {
    max(abs(Ddiag_pm - d$sd_noise^2)) < 0.1
  })
}

cat("\n[3] determinism: same seed -> same chain\n")
{
  d   <- simulate_iox(n = 250, q = 2, rho = 0.3, sd_noise = 0.3, seed = 7L)
  set.seed(99L); fit1 <- fit_latent_block(d, iter = 20, m = 10)
  set.seed(99L); fit2 <- fit_latent_block(d, iter = 20, m = 10)

  ok("Sigma chains identical with same R-side seed", {
    isTRUE(all.equal(fit1$Sigma, fit2$Sigma, tolerance = 1e-12))
  })
  ok("Ddiag chains identical with same R-side seed", {
    isTRUE(all.equal(fit1$Ddiag, fit2$Ddiag, tolerance = 1e-12))
  })
}


# ---------- latent_model = 3 (sequential per-outcome sampler) ----------

fit_latent_seq <- function(d, iter = 60, m = 15, pc = "vrpc") {
  spiox::spiox(
    Y = d$Y, X = d$X, coords = d$coords, m = m,
    method = "latent", fit = "mcmc",
    iter = iter, print_every = 0,
    starting = list(Theta = d$Theta_true,
                    Sigma = d$Sigma_true,
                    W     = d$W,
                    Ddiag = rep(d$sd_noise^2, ncol(d$Y)),
                    Beta  = matrix(0, ncol(d$X), ncol(d$Y))),
    opts = list(update_Theta = c(0, 0, 0),
                cg_preconditioner = pc,
                m_pc = m,
                num_threads = 4),
    debug = list(sampling = 3L,
                 sample_Beta = TRUE, sample_Sigma = TRUE, sample_Ddiag = TRUE)
  )
}

cat("\n--- spiox VRPC tests (latent_model = 3, per-outcome sequential) ---\n")

cat("\n[4] small problem (n=400, q=3, sampling=3, VRPC)\n")
{
  d <- simulate_iox(n = 400, q = 3, rho = 0.5, sd_noise = 0.3, seed = 19L)
  fit <- fit_latent_seq(d, iter = 60, m = 15, pc = "vrpc")

  ok("finite Sigma / W / Ddiag", {
    all(is.finite(fit$Sigma)) && all(is.finite(fit$W)) && all(is.finite(fit$Ddiag))
  })
  ok("vrpc preconditioner code recorded", {
    all(fit$cg_preconditioner == 3L)
  })
  ok("CG iter count is small (< 200 per sweep)", {
    max(fit$cg_iters) < 200
  })
  ok("vrpc_n_builds is 1 (PC built once, frozen thereafter)", {
    fit$vrpc_n_builds == 1L
  })

  keep      <- 31:60
  Sigma_pm  <- apply(fit$Sigma[, , keep, drop = FALSE], c(1, 2), mean)
  Ddiag_pm  <- rowMeans(fit$Ddiag[, keep, drop = FALSE])
  ok("Sigma diag posterior mean within ±0.6 of truth", {
    max(abs(diag(Sigma_pm) - diag(d$Sigma_true))) < 0.6
  })
  ok("Ddiag posterior mean within ±0.1 of truth", {
    max(abs(Ddiag_pm - d$sd_noise^2)) < 0.1
  })
}

cat("\n[5] VRPC vs Jacobi posteriors agree (latent_model=3)\n")
{
  d <- simulate_iox(n = 300, q = 2, rho = 0.4, sd_noise = 0.3, seed = 23L)
  set.seed(7L); fit_jac  <- fit_latent_seq(d, iter = 80, m = 12, pc = "jacobi")
  set.seed(7L); fit_vrpc <- fit_latent_seq(d, iter = 80, m = 12, pc = "vrpc")

  keep      <- 41:80
  Sigma_jac  <- apply(fit_jac$Sigma[, , keep, drop = FALSE], c(1, 2), mean)
  Sigma_vrpc <- apply(fit_vrpc$Sigma[, , keep, drop = FALSE], c(1, 2), mean)
  Dpm_jac    <- rowMeans(fit_jac$Ddiag[, keep, drop = FALSE])
  Dpm_vrpc   <- rowMeans(fit_vrpc$Ddiag[, keep, drop = FALSE])

  ok("Sigma posterior means agree (max abs diff < 0.25)", {
    max(abs(Sigma_jac - Sigma_vrpc)) < 0.25
  })
  ok("Ddiag posterior means agree (max abs diff < 0.05)", {
    max(abs(Dpm_jac - Dpm_vrpc)) < 0.05
  })
}

cat("\nAll VRPC tests passed.\n")
