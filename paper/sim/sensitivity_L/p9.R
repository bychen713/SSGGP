#!/usr/bin/env Rscript

## Compare two finite-L likelihood calculations for p=9:
##   integrated_W  -- W is integrated analytically;
##   simulated_W   -- the supplement's Phi Phi' covariance from simulated W.
##
## The Matern experiment uses the exact covariance likelihood as its reference.
## The mixed-family experiment uses the integrated-W likelihood at L=100,000.

options(warn = 1)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "sim", "sensitivity_L", "p9.R")))
  stop("run from paper/", call. = FALSE)

suppressPackageStartupMessages(library(SSGGP))

action <- Sys.getenv("SSGGP_ACTION", "run")
family_design <- Sys.getenv("SSGGP_FAMILY", "matern")
precision_design <- Sys.getenv("SSGGP_PRECISION", "negative")
replicate_id <- as.integer(Sys.getenv(
  "SSGGP_REPLICATE", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
if (!family_design %in% c("matern", "mixed"))
  stop("SSGGP_FAMILY must be 'matern' or 'mixed'")
if (!precision_design %in% c("negative", "positive"))
  stop("SSGGP_PRECISION must be 'negative' or 'positive'")
if (!is.finite(replicate_id) || replicate_id < 1L)
  stop("SSGGP_REPLICATE must be a positive integer")

P <- 9L
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
                     paste0(family_design, "_p", P, precision_suffix))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tridiagonal_precision <- function(P, diagonal = 2,
                                  offdiagonal = precision_offdiagonal) {
  Q <- diag(diagonal, P)
  if (P > 1L) {
    Q[cbind(seq_len(P - 1L), 2:P)] <- offdiagonal
    Q[cbind(2:P, seq_len(P - 1L))] <- offdiagonal
  }
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

direct_loglik <- function(K, Y) {
  diag(K) <- diag(K) + nugget
  ## Both covariance builders fill a symmetric matrix exactly.  Avoid making
  ## several additional np-by-np copies merely to symmetrise it again.
  U <- chol(K)
  y <- as.vector(Y)
  -0.5 * length(y) * log(2 * pi) - sum(log(diag(U))) -
    0.5 * sum(backsolve(U, y, transpose = TRUE)^2)
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

assemble_integrated <- function(blocks, Sigma, pairs, L, n, max_n) {
  N <- n * P
  K <- matrix(0, N, N)
  for (q in seq_len(nrow(pairs))) {
    j <- pairs[q, "j"]
    k <- pairs[q, "k"]
    rows <- (j - 1L) * n + seq_len(n)
    cols <- (k - 1L) * n + seq_len(n)
    block <- Sigma[j, k] * blocks[[q]][seq_len(n), seq_len(n), drop = FALSE] / L
    K[rows, cols] <- block
    if (j != k) K[cols, rows] <- t(block)
  }
  K
}

update_simulated_covariance <- function(accumulator, C, W) {
  Phi <- do.call(rbind, lapply(seq_len(P), function(j)
    sweep(C, 2L, W[, j], "*")))
  accumulator + tcrossprod(Phi)
}

subset_process_major <- function(n, max_n) {
  unlist(lapply(0:(P - 1L), function(j) j * max_n + seq_len(n)),
         use.names = FALSE)
}

simulate_mixed_reference <- function(coords, features, r, Sigma, seed,
                                     block_size = 500L) {
  U_sigma <- chol(Sigma)
  set.seed(seed)
  Z <- matrix(stats::rnorm(nrow(r) * P), nrow(r), P)
  W <- (Z %*% U_sigma) * r
  signal <- matrix(0, nrow(coords), P)
  for (start in seq.int(1L, nrow(r), by = block_size)) {
    end <- min(start + block_size - 1L, nrow(r))
    idx <- start:end
    C <- cos(sweep(coords %*% t(features$X[idx, , drop = FALSE]),
                   2L, features$phase[idx], "+"))
    signal <- signal + C %*% W[idx, , drop = FALSE]
  }
  signal <- signal / sqrt(nrow(r))
  set.seed(seed + 1000000L)
  signal + matrix(stats::rnorm(length(signal), sd = sqrt(nugget)),
                  nrow(coords), P)
}

run_replicate <- function(replicate_id, family_design) {
  if (replicate_id > n_replicates)
    stop("replicate exceeds configured total of ", n_replicates)
  started <- proc.time()[["elapsed"]]
  max_n <- max(n_values)
  max_L <- if (family_design == "mixed") L_reference_mixed else max(L_values)

  precision <- tridiagonal_precision(P)
  Sigma <- solve(precision)
  process_eigen <- eigen(Sigma, symmetric = TRUE)
  U_sigma <- chol(Sigma)

  set.seed(110000L + replicate_id)
  coords <- matrix(stats::runif(max_n * 2L, 0, domain_max), ncol = 2L)
  features <- draw_features(max_L, 2L, 210000L + replicate_id)
  importance <- importance_sqrt(features, family_design)
  r <- importance$r

  Y <- if (family_design == "matern") {
    simulate_exact_matern(coords, process_eigen, 310000L + replicate_id)
  } else {
    ## Match the simulation design with a 100k-feature mixed process, but use
    ## a feature draw independent of the one used to evaluate every likelihood.
    ## Otherwise the lower-L likelihoods would be privileged prefixes of the
    ## very same random covariance that generated Y.
    truth_features <- draw_features(
      L_reference_mixed, 2L, 510000L + replicate_id)
    truth_importance <- importance_sqrt(truth_features, family_design)
    simulate_mixed_reference(coords, truth_features, truth_importance$r, Sigma,
                             310000L + replicate_id)
  }

  reference <- setNames(rep(NA_real_, length(n_values)), n_values)
  if (family_design == "matern") {
    for (n in n_values) {
      reference[as.character(n)] <- separable_loglik(
        Y[seq_len(n), , drop = FALSE],
        matern32_correlation(coords[seq_len(n), , drop = FALSE]),
        process_eigen)
    }
  }

  pairs <- triangular_pairs(P)
  integrated_blocks <- if (family_design == "mixed")
    new_integrated_blocks(pairs, max_n) else NULL
  integrated_spatial <- if (family_design == "matern")
    matrix(0, max_n, max_n) else NULL
  simulated_sum <- matrix(0, max_n * P, max_n * P)

  ## Independent W draw for the supplement approximation, paired across L.
  set.seed(410000L + replicate_id)
  Z_approx <- matrix(stats::rnorm(max(L_values) * P), max(L_values), P)
  W_approx <- (Z_approx %*% U_sigma) * r[seq_len(max(L_values)), , drop = FALSE]

  rows <- list()
  row_id <- 0L
  block_size <- 500L
  for (start in seq.int(1L, max_L, by = block_size)) {
    end <- min(start + block_size - 1L, max_L)
    idx <- start:end
    C <- cos(sweep(coords %*% t(features$X[idx, , drop = FALSE]),
                   2L, features$phase[idx], "+"))

    if (family_design == "matern") {
      weighted <- sweep(C, 2L, r[idx, 1L], "*")
      integrated_spatial <- integrated_spatial + tcrossprod(weighted)
    } else {
      integrated_blocks <- update_integrated_blocks(
        integrated_blocks, C, r[idx, , drop = FALSE], pairs)
    }

    if (end <= max(L_values)) {
      simulated_sum <- update_simulated_covariance(
        simulated_sum, C, W_approx[idx, , drop = FALSE])
    }

    if (end %in% L_values) {
      for (n in n_values) {
        Yi <- Y[seq_len(n), , drop = FALSE]
        ll_integrated <- if (family_design == "matern") {
          separable_loglik(
            Yi,
            integrated_spatial[seq_len(n), seq_len(n), drop = FALSE] / end,
            process_eigen)
        } else {
          direct_loglik(
            assemble_integrated(integrated_blocks, Sigma, pairs, end, n, max_n),
            Yi)
        }

        obs <- subset_process_major(n, max_n)
        ll_simulated <- direct_loglik(
          simulated_sum[obs, obs, drop = FALSE] / end, Yi)

        for (method in c("integrated_W", "simulated_W")) {
          row_id <- row_id + 1L
          rows[[row_id]] <- data.frame(
            family = family_design, p = P, n = n, L = end,
            replicate = replicate_id, method = method,
            ll = if (method == "integrated_W") ll_integrated else ll_simulated,
            stringsAsFactors = FALSE)
        }
      }
      message(sprintf("%s replicate %d: completed L=%d (%.1f min)",
                      family_design, replicate_id, end,
                      (proc.time()[["elapsed"]] - started) / 60))
    }

    if (family_design == "mixed" && end == L_reference_mixed) {
      for (n in n_values) {
        reference[as.character(n)] <- direct_loglik(
          assemble_integrated(integrated_blocks, Sigma, pairs, end, n, max_n),
          Y[seq_len(n), , drop = FALSE])
      }
    }
  }

  result <- do.call(rbind, rows)
  result$ll_ref <- unname(reference[as.character(result$n)])
  result$abs_error <- abs(result$ll - result$ll_ref)
  result$error_per_n <- result$abs_error / result$n

  metadata <- list(
    family = family_design, p = P, n_values = n_values, L_values = L_values,
    reference_L = if (family_design == "mixed") L_reference_mixed else Inf,
    replicate = replicate_id, proposal = proposal, nugget = nugget,
    precision = precision_label, process_family = importance$family,
    process_parameters = importance$parameters,
    capped = importance$capped,
    simulation_capped = if (family_design == "mixed")
      truth_importance$capped else integer(P),
    elapsed_seconds = proc.time()[["elapsed"]] - started)
  result_file <- file.path(out_dir, sprintf("replicate_%02d.rds", replicate_id))
  temp_file <- paste0(result_file, ".tmp-", Sys.getpid())
  saveRDS(list(result = result, metadata = metadata), temp_file)
  if (!file.rename(temp_file, result_file))
    stop("could not atomically install ", result_file)
  message("wrote ", result_file)
  invisible(result)
}

selftest <- function() {
  p_old <- P
  n <- 7L
  L <- 5L
  set.seed(9)
  C <- matrix(stats::rnorm(n * L), n, L)
  r <- matrix(stats::runif(L * p_old, 0.2, 1.5), L, p_old)
  Sigma <- solve(tridiagonal_precision(p_old))
  pairs <- triangular_pairs(p_old)
  blocks <- update_integrated_blocks(new_integrated_blocks(pairs, n), C, r, pairs)
  got <- assemble_integrated(blocks, Sigma, pairs, L, n, n)
  design <- do.call(cbind, lapply(seq_len(L), function(l)
    kronecker(diag(r[l, ], p_old), C[, l, drop = FALSE])))
  want <- design %*% kronecker(diag(L), Sigma) %*% t(design) / L
  stopifnot(isTRUE(all.equal(got, want, tolerance = 1e-10)))

  U <- chol(Sigma)
  W <- matrix(stats::rnorm(L * p_old), L, p_old) %*% U
  app <- update_simulated_covariance(matrix(0, n * p_old, n * p_old), C, W)
  Phi <- do.call(rbind, lapply(seq_len(p_old), function(j)
    sweep(C, 2L, W[, j], "*")))
  stopifnot(isTRUE(all.equal(app, tcrossprod(Phi), tolerance = 1e-12)))

  ## Check the process-major stacking convention used by direct_loglik()
  ## against the separable eigensystem calculation.
  spatial <- tcrossprod(matrix(stats::rnorm(n * 4L), n, 4L)) / 4
  pe <- eigen(Sigma, symmetric = TRUE)
  Y <- matrix(stats::rnorm(n * p_old), n, p_old)
  explicit <- kronecker(Sigma, spatial)
  stopifnot(isTRUE(all.equal(
    direct_loglik(explicit, Y), separable_loglik(Y, spatial, pe),
    tolerance = 1e-9)))
  message("self-test passed: integrated and simulated-W covariance builders agree")
}

switch(action,
       run = run_replicate(replicate_id, family_design),
       selftest = selftest(),
       stop("unknown SSGGP_ACTION: ", action))
