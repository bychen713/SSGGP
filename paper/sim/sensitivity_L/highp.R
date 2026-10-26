#!/usr/bin/env Rscript

## High-memory p=25/p=50 likelihood-sensitivity calculation. The integrated-W
## and simulated-W stages are deliberately separate: mixed-family integrated-W
## needs triangular process-block accumulators, whereas simulated-W needs one
## observed-space covariance accumulator. Keeping them out of memory at the
## same time makes p=50,n=1000 practical.

options(warn = 1)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "sim", "sensitivity_L", "highp.R")))
  stop("run from paper/", call. = FALSE)

suppressPackageStartupMessages({
  library(SSGGP)
  library(Rcpp)
})

stage <- Sys.getenv("SSGGP_STAGE", "integrated")
family_design <- Sys.getenv("SSGGP_FAMILY", "matern")
precision_design <- Sys.getenv("SSGGP_PRECISION", "negative")
P <- as.integer(Sys.getenv("SSGGP_P", "25"))
replicate_id <- as.integer(Sys.getenv(
  "SSGGP_REPLICATE", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
if (!stage %in% c("integrated", "simulated", "selftest"))
  stop("SSGGP_STAGE must be integrated, simulated, or selftest")
if (!family_design %in% c("matern", "mixed"))
  stop("SSGGP_FAMILY must be matern or mixed")
if (!precision_design %in% c("negative", "positive"))
  stop("SSGGP_PRECISION must be negative or positive")
if (!P %in% c(25L, 50L)) stop("SSGGP_P must be 25 or 50")
if (!is.finite(replicate_id) || replicate_id < 1L || replicate_id > 15L)
  stop("SSGGP_REPLICATE must be from 1 through 15")

n_values <- c(500L, 1000L)
L_values <- c(500L, 1000L, 2500L, 5000L, 10000L, 20000L,
              30000L, 40000L, 50000L)
L_reference_mixed <- 100000L
n_replicates <- 15L
nugget <- 0.1
domain_max <- 5
proposal <- list(a = 3, nu = 1.25)
precision_offdiagonal <- if (precision_design == "positive") 1 else -1
precision_label <- paste0("tridiagonal_", precision_design)
precision_suffix <- if (precision_design == "negative") "" else
  paste0("_", precision_design)

out_dir <- file.path(project_root, "sim", "output", "sensitivity_L",
                     sprintf("%s_p%d%s", family_design, P, precision_suffix))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
integrated_file <- file.path(
  out_dir, sprintf("integrated_replicate_%02d.rds", replicate_id))
simulated_file <- file.path(
  out_dir, sprintf("simulated_replicate_%02d.rds", replicate_id))

Rcpp::sourceCpp(file.path(project_root, "sim", "sensitivity_L", "highp.cpp"),
                rebuild = FALSE, verbose = FALSE)

tridiagonal_precision <- function(P, diagonal = 2,
                                  offdiagonal = precision_offdiagonal) {
  Q <- diag(diagonal, P)
  Q[cbind(seq_len(P - 1L), 2:P)] <- offdiagonal
  Q[cbind(2:P, seq_len(P - 1L))] <- offdiagonal
  Q
}

pairwise_distance <- function(x) {
  ss <- rowSums(x^2)
  sqrt(pmax(outer(ss, ss, "+") - 2 * tcrossprod(x), 0))
}

matern32_correlation <- function(coords, a = 3) {
  ar <- a * pairwise_distance(coords)
  R <- (1 + ar) * exp(-ar)
  diag(R) <- 1
  R
}

log_matern_spectral <- function(X, a, nu) {
  d <- ncol(X)
  2 * nu * log(a) + lgamma(nu + d / 2) - lgamma(nu) -
    (d / 2) * log(pi) -
    (nu + d / 2) * log(a^2 + rowSums(X^2))
}

draw_features <- function(L, d, seed) {
  set.seed(seed)
  phase <- stats::runif(L, 0, 2 * pi)
  scale_mix <- stats::rgamma(L, shape = proposal$nu,
                            rate = proposal$a^2 / 2)
  X <- matrix(stats::rnorm(d * L), nrow = L, ncol = d) /
    sqrt(scale_mix)
  list(X = X, phase = phase)
}

mixed_truth <- function(P) {
  design <- local({
    set.seed(20260811L)
    family <- sample(c(rep("matern", 20), rep("ch", 20),
                       rep("se", 10)))
    parameters <- lapply(family, function(f) switch(f,
      matern = c(nu = sample(c(1.5, 2.5), 1), a = runif(1, 2.5, 3.5)),
      ch = c(beta = runif(1, 0.6, 0.9),
             alpha = runif(1, 1.25, 1.75),
             nu = sample(c(1.5, 2.5), 1)),
      se = c(a = runif(1, 1.75, 2.25))))
    list(family = family, parameters = parameters)
  })
  list(family = design$family[seq_len(P)],
       parameters = design$parameters[seq_len(P)])
}

importance_sqrt <- function(features, family_design) {
  X <- features$X
  log_g <- log_matern_spectral(X, proposal$a, proposal$nu)
  if (identical(family_design, "matern")) {
    log_f <- log_matern_spectral(X, 3, 1.5)
    r <- exp(0.5 * (log(2) + log_f - log_g))
    return(list(r = matrix(r, nrow(X), P), capped = integer(P),
                family = rep("matern", P),
                parameters = rep(list(c(a = 3, nu = 1.5)), P)))
  }
  truth <- mixed_truth(P)
  fam <- SSGGP:::resolve_spectral(as.list(truth$family), P)
  sr <- SSGGP:::build_sqrt_R(
    fam, truth$parameters, X, rowSums(X^2), log_g)
  list(r = sqrt(2) / sr$sqrtR, capped = sr$capped,
       family = truth$family, parameters = truth$parameters)
}

safe_eigen <- function(K) {
  ee <- eigen((K + t(K)) / 2, symmetric = TRUE)
  scale <- max(abs(ee$values), 1)
  if (min(ee$values) < -1e-8 * scale)
    stop("covariance has a materially negative eigenvalue: ", min(ee$values))
  ee$values <- pmax(ee$values, 0)
  ee
}

separable_loglik <- function(Y, spatial_cov, process_eigen) {
  se <- safe_eigen(spatial_cov)
  transformed <- crossprod(se$vectors, Y) %*% process_eigen$vectors
  variance <- outer(se$values, process_eigen$values, "*") + nugget
  -0.5 * (length(Y) * log(2 * pi) + sum(log(variance)) +
            sum(transformed^2 / variance))
}

simulate_exact_matern <- function(coords, process_eigen, seed) {
  se <- safe_eigen(matern32_correlation(coords))
  set.seed(seed)
  Z <- matrix(stats::rnorm(nrow(coords) * P), nrow(coords), P)
  spatial_part <- se$vectors %*% sweep(Z, 1L, sqrt(se$values), "*")
  process_sqrt <- sweep(t(process_eigen$vectors), 1L,
                        sqrt(process_eigen$values), "*")
  signal <- spatial_part %*% process_sqrt
  signal + matrix(stats::rnorm(length(signal), sd = sqrt(nugget)),
                  nrow(coords), P)
}

simulate_mixed <- function(coords, features, r, Sigma, seed,
                           block_size = 500L) {
  set.seed(seed)
  W <- (matrix(stats::rnorm(nrow(r) * P), nrow(r), P) %*% chol(Sigma)) * r
  signal <- matrix(0, nrow(coords), P)
  for (start in seq.int(1L, nrow(r), by = block_size)) {
    end <- min(start + block_size - 1L, nrow(r))
    idx <- start:end
    C <- cos(sweep(coords %*% t(features$X[idx, , drop = FALSE]),
                   2L, features$phase[idx], "+"))
    signal <- signal + C %*% W[idx, , drop = FALSE]
  }
  set.seed(seed + 1000000L)
  signal / sqrt(nrow(r)) +
    matrix(stats::rnorm(length(signal), sd = sqrt(nugget)),
           nrow(coords), P)
}

generate_design <- function(replicate_id, family_design) {
  max_n <- max(n_values)
  max_L <- if (family_design == "mixed") L_reference_mixed else max(L_values)
  Sigma <- solve(tridiagonal_precision(P))
  process_eigen <- eigen(Sigma, symmetric = TRUE)
  set.seed(110000L + replicate_id)
  coords <- matrix(stats::runif(max_n * 2L, 0, domain_max), ncol = 2L)
  features <- draw_features(max_L, 2L, 210000L + replicate_id)
  importance <- importance_sqrt(features, family_design)
  if (family_design == "matern") {
    Y <- simulate_exact_matern(coords, process_eigen, 310000L + replicate_id)
    simulation_capped <- integer(P)
  } else {
    truth_features <- draw_features(
      L_reference_mixed, 2L, 510000L + replicate_id)
    truth_importance <- importance_sqrt(truth_features, family_design)
    Y <- simulate_mixed(coords, truth_features, truth_importance$r, Sigma,
                        310000L + replicate_id)
    simulation_capped <- truth_importance$capped
  }
  list(coords = coords, features = features, r = importance$r, Y = Y,
       Sigma = Sigma, process_eigen = process_eigen,
       importance = importance, simulation_capped = simulation_capped)
}

triangular_pairs <- function(P) {
  out <- do.call(rbind, lapply(seq_len(P), function(j)
    cbind(j = j, k = j:P)))
  storage.mode(out) <- "integer"
  out
}

new_integrated_blocks <- function(pairs, n) {
  lapply(seq_len(nrow(pairs)), function(i) matrix(0, n, n))
}

update_integrated_blocks <- function(blocks, C, r, pairs) {
  weighted <- lapply(seq_len(ncol(r)), function(j)
    sweep(C, 2L, r[, j], "*"))
  for (q in seq_len(nrow(pairs))) {
    j <- pairs[q, "j"]
    k <- pairs[q, "k"]
    blocks[[q]] <- blocks[[q]] + tcrossprod(weighted[[j]], weighted[[k]])
  }
  blocks
}

nested_positions <- function(P, n_small, n_full) {
  n_tail <- n_full - n_small
  lapply(seq_len(P), function(j)
    c((j - 1L) * n_small + seq_len(n_small),
      P * n_small + (j - 1L) * n_tail + seq_len(n_tail)))
}

assemble_integrated_nested <- function(blocks, Sigma, pairs, L,
                                       n_small, n_full) {
  K <- matrix(0, n_full * P, n_full * P)
  pos <- nested_positions(P, n_small, n_full)
  for (q in seq_len(nrow(pairs))) {
    j <- pairs[q, "j"]
    k <- pairs[q, "k"]
    block <- Sigma[j, k] * blocks[[q]] / L
    K[pos[[j]], pos[[k]]] <- block
    if (j != k) K[pos[[k]], pos[[j]]] <- t(block)
  }
  K
}

nested_y <- function(Y, n_small) {
  c(as.vector(Y[seq_len(n_small), , drop = FALSE]),
    as.vector(Y[(n_small + 1L):nrow(Y), , drop = FALSE]))
}

gaussian_loglik <- function(N, logdet, quad) {
  -0.5 * (N * log(2 * pi) + logdet + quad)
}

nested_loglik <- function(K, Y, n_small) {
  y <- nested_y(Y, n_small)
  m <- n_small * ncol(Y)
  st <- chol_nested_stats_inplace(K, y, m, nugget)
  c(`500` = gaussian_loglik(m, st$logdet_small, st$quad_small),
    `1000` = gaussian_loglik(length(y), st$logdet_full, st$quad_full))
}

single_lower_loglik <- function(K, y) {
  st <- chol_stats_lower_inplace(K, y, nugget)
  gaussian_loglik(length(y), st$logdet, st$quad)
}

feature_loglik <- function(C, W, Y) {
  L <- ncol(C)
  C <- C / sqrt(L)
  dvec <- rep(nugget, ncol(Y))
  B <- build_feature_system_highp(C, W, dvec, block_size = 256L)
  CY <- crossprod(C, sweep(Y, 2L, dvec, "/"))
  q <- rowSums(CY * W)
  yquad <- sum(sweep(Y^2, 2L, dvec, "/"))
  st <- chol_feature_stats_inplace(B, q)
  gaussian_loglik(length(Y), nrow(Y) * sum(log(dvec)) + st$logdet,
                  yquad - st$quad)
}

finish_result <- function(rows, reference) {
  ans <- do.call(rbind, rows)
  ans$ll_ref <- unname(reference[as.character(ans$n)])
  ans$abs_error <- abs(ans$ll - ans$ll_ref)
  ans$error_per_n <- ans$abs_error / ans$n
  ans
}

atomic_save <- function(object, path) {
  temp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temp)
  if (!file.rename(temp, path)) stop("could not atomically install ", path)
}

run_integrated <- function() {
  started <- proc.time()[["elapsed"]]
  d <- generate_design(replicate_id, family_design)
  max_n <- max(n_values)
  reference <- setNames(rep(NA_real_, length(n_values)), n_values)
  if (family_design == "matern") {
    for (n in n_values) reference[as.character(n)] <- separable_loglik(
      d$Y[seq_len(n), , drop = FALSE],
      matern32_correlation(d$coords[seq_len(n), , drop = FALSE]),
      d$process_eigen)
  }

  pairs <- triangular_pairs(P)
  blocks <- if (family_design == "mixed")
    new_integrated_blocks(pairs, max_n) else NULL
  spatial_sum <- if (family_design == "matern")
    matrix(0, max_n, max_n) else NULL
  rows <- list()
  row_id <- 0L
  max_L <- if (family_design == "mixed") L_reference_mixed else max(L_values)

  for (start in seq.int(1L, max_L, by = 500L)) {
    end <- min(start + 499L, max_L)
    idx <- start:end
    C <- cos(sweep(d$coords %*% t(d$features$X[idx, , drop = FALSE]),
                   2L, d$features$phase[idx], "+"))
    if (family_design == "matern") {
      weighted <- sweep(C, 2L, d$r[idx, 1L], "*")
      spatial_sum <- spatial_sum + tcrossprod(weighted)
    } else {
      blocks <- update_integrated_blocks(
        blocks, C, d$r[idx, , drop = FALSE], pairs)
    }

    if (end %in% L_values) {
      ll <- if (family_design == "matern") {
        setNames(vapply(n_values, function(n) separable_loglik(
          d$Y[seq_len(n), , drop = FALSE],
          spatial_sum[seq_len(n), seq_len(n), drop = FALSE] / end,
          d$process_eigen), numeric(1)), n_values)
      } else {
        K <- assemble_integrated_nested(
          blocks, d$Sigma, pairs, end, min(n_values), max_n)
        value <- nested_loglik(K, d$Y, min(n_values))
        rm(K); invisible(gc())
        value
      }
      for (n in n_values) {
        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          family = family_design, p = P, n = n, L = end,
          replicate = replicate_id, method = "integrated_W",
          ll = unname(ll[as.character(n)]), stringsAsFactors = FALSE)
      }
      message(sprintf("%s p=%d rep=%d integrated L=%d (%.1f min)",
                      family_design, P, replicate_id, end,
                      (proc.time()[["elapsed"]] - started) / 60))
    }

    if (family_design == "mixed" && end == L_reference_mixed) {
      K <- assemble_integrated_nested(
        blocks, d$Sigma, pairs, end, min(n_values), max_n)
      reference[] <- nested_loglik(K, d$Y, min(n_values))
      rm(K); invisible(gc())
    }
  }

  result <- finish_result(rows, reference)
  metadata <- list(
    stage = "integrated", family = family_design, p = P,
    n_values = n_values, L_values = L_values,
    reference_L = if (family_design == "mixed") L_reference_mixed else Inf,
    replicate = replicate_id, proposal = proposal, nugget = nugget,
    precision = precision_label,
    process_family = d$importance$family,
    process_parameters = d$importance$parameters,
    capped = d$importance$capped,
    simulation_capped = d$simulation_capped,
    reference = reference,
    elapsed_seconds = proc.time()[["elapsed"]] - started)
  atomic_save(list(result = result, metadata = metadata), integrated_file)
  message("wrote ", integrated_file)
}

run_simulated <- function() {
  if (!file.exists(integrated_file))
    stop("integrated stage must finish first: ", integrated_file)
  integrated <- readRDS(integrated_file)
  reference <- integrated$metadata$reference
  started <- proc.time()[["elapsed"]]
  d <- generate_design(replicate_id, family_design)
  max_n <- max(n_values)
  max_L <- max(L_values)
  U_sigma <- chol(d$Sigma)
  set.seed(410000L + replicate_id)
  W <- (matrix(stats::rnorm(max_L * P), max_L, P) %*% U_sigma) *
    d$r[seq_len(max_L), , drop = FALSE]

  N_full <- max_n * P
  N_small <- min(n_values) * P
  Ksum <- matrix(0, N_full, N_full)
  C_all <- matrix(0, max_n, max_L)
  y_nested <- nested_y(d$Y, min(n_values))
  y_small <- y_nested[seq_len(N_small)]
  rows <- list()
  row_id <- 0L

  for (start in seq.int(1L, max_L, by = 500L)) {
    end <- min(start + 499L, max_L)
    idx <- start:end
    C <- cos(sweep(d$coords %*% t(d$features$X[idx, , drop = FALSE]),
                   2L, d$features$phase[idx], "+"))
    C_all[, idx] <- C
    update_nested_covariance_inplace(
      Ksum, C, W[idx, , drop = FALSE], min(n_values))

    if (end %in% L_values) {
      ll <- setNames(rep(NA_real_, length(n_values)), n_values)
      if (end < N_full) {
        ll["1000"] <- feature_loglik(
          C_all[, seq_len(end), drop = FALSE],
          W[seq_len(end), , drop = FALSE], d$Y)
        if (end < N_small) {
          ll["500"] <- feature_loglik(
            C_all[seq_len(min(n_values)), seq_len(end), drop = FALSE],
            W[seq_len(end), , drop = FALSE],
            d$Y[seq_len(min(n_values)), , drop = FALSE])
        } else {
          Ksmall <- Ksum[seq_len(N_small), seq_len(N_small), drop = FALSE] / end
          ll["500"] <- single_lower_loglik(Ksmall, y_small)
          rm(Ksmall); invisible(gc())
        }
      } else {
        Kwork <- Ksum / end
        ll[] <- nested_loglik(Kwork, d$Y, min(n_values))
        rm(Kwork); invisible(gc())
      }

      for (n in n_values) {
        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          family = family_design, p = P, n = n, L = end,
          replicate = replicate_id, method = "simulated_W",
          ll = unname(ll[as.character(n)]), stringsAsFactors = FALSE)
      }
      message(sprintf("%s p=%d rep=%d simulated L=%d (%.1f min)",
                      family_design, P, replicate_id, end,
                      (proc.time()[["elapsed"]] - started) / 60))
    }
  }

  result <- finish_result(rows, reference)
  metadata <- list(
    stage = "simulated", family = family_design, p = P,
    n_values = n_values, L_values = L_values,
    reference_L = integrated$metadata$reference_L,
    replicate = replicate_id, proposal = proposal, nugget = nugget,
    precision = precision_label,
    elapsed_seconds = proc.time()[["elapsed"]] - started)
  atomic_save(list(result = result, metadata = metadata), simulated_file)
  message("wrote ", simulated_file)
}

selftest <- function() {
  set.seed(19)
  n <- 6L; ns <- 3L; p <- 3L; b <- 5L
  C <- matrix(rnorm(n * b), n, b)
  W <- matrix(rnorm(b * p), b, p)
  K <- matrix(0, n * p, n * p)
  update_nested_covariance_inplace(K, C, W, ns)
  Phi_pm <- do.call(rbind, lapply(seq_len(p), function(j)
    sweep(C, 2L, W[, j], "*")))
  ord <- c(unlist(lapply(0:(p - 1L), function(j) j * n + seq_len(ns))),
           unlist(lapply(0:(p - 1L), function(j) j * n + (ns + 1L):n)))
  want <- tcrossprod(Phi_pm[ord, , drop = FALSE])
  stopifnot(isTRUE(all.equal(K[lower.tri(K, diag = TRUE)],
                            want[lower.tri(want, diag = TRUE)],
                            tolerance = 1e-10)))
  y <- rnorm(n * p)
  Kwork <- K / b
  st <- chol_nested_stats_inplace(Kwork, y, ns * p, nugget)
  ll_full <- gaussian_loglik(length(y), st$logdet_full, st$quad_full)
  ll_small <- gaussian_loglik(ns * p, st$logdet_small, st$quad_small)
  direct <- function(A, z) {
    U <- chol(A + diag(nugget, nrow(A)))
    gaussian_loglik(length(z), 2 * sum(log(diag(U))),
                    sum(backsolve(U, z, transpose = TRUE)^2))
  }
  stopifnot(isTRUE(all.equal(ll_full, direct(want / b, y), tolerance = 1e-9)),
            isTRUE(all.equal(ll_small,
                             direct(want[seq_len(ns * p), seq_len(ns * p)] / b,
                                    y[seq_len(ns * p)]), tolerance = 1e-9)))
  Y <- matrix(rnorm(n * p), n, p)
  ll_feature <- feature_loglik(C, W, Y)
  Phi_lm <- do.call(rbind, lapply(seq_len(n), function(i)
    t(W) * rep(C[i, ] / sqrt(b), each = p)))
  ll_feature_direct <- direct(tcrossprod(Phi_lm), as.vector(t(Y)))
  stopifnot(isTRUE(all.equal(ll_feature, ll_feature_direct, tolerance = 1e-9)))
  message("high-p covariance, nested Cholesky, and feature-space self-tests passed")
}

switch(stage,
       integrated = run_integrated(),
       simulated = run_simulated(),
       selftest = selftest())
