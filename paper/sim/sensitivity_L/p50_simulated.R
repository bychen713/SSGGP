#!/usr/bin/env Rscript

## Optimized p=50 simulated-W stage. For n=1000, L never exceeds np=50,000,
## so the supplement's L-dimensional feature system is used (including equality
## at L=50,000). Only the n=500 observed covariance is accumulated, reducing
## covariance-formation work and memory by a factor of four.

if (Sys.getenv("SSGGP_P", "") != "50")
  stop("this entry point requires SSGGP_P=50")
Sys.setenv(SSGGP_STAGE = "selftest")
source("sim/sensitivity_L/highp.R", local = FALSE)
Rcpp::sourceCpp(file.path(project_root, "sim", "sensitivity_L", "p50_sim.cpp"),
                rebuild = FALSE, verbose = FALSE)

run_p50_simulated <- function() {
  if (!file.exists(integrated_file))
    stop("integrated stage must finish first: ", integrated_file)
  integrated <- readRDS(integrated_file)
  reference <- integrated$metadata$reference
  started <- proc.time()[["elapsed"]]
  d <- generate_design(replicate_id, family_design)
  max_n <- max(n_values)
  n_small <- min(n_values)
  max_L <- max(L_values)
  set.seed(410000L + replicate_id)
  W <- (matrix(stats::rnorm(max_L * P), max_L, P) %*% chol(d$Sigma)) *
    d$r[seq_len(max_L), , drop = FALSE]

  N_small <- n_small * P
  Ksmall_sum <- matrix(0, N_small, N_small)
  C_all <- matrix(0, max_n, max_L)
  y_small <- as.vector(d$Y[seq_len(n_small), , drop = FALSE])
  rows <- list()
  row_id <- 0L

  for (start in seq.int(1L, max_L, by = 500L)) {
    end <- min(start + 499L, max_L)
    idx <- start:end
    C <- cos(sweep(d$coords %*% t(d$features$X[idx, , drop = FALSE]),
                   2L, d$features$phase[idx], "+"))
    C_all[, idx] <- C
    update_process_covariance_inplace(
      Ksmall_sum, C[seq_len(n_small), , drop = FALSE],
      W[idx, , drop = FALSE])

    if (end %in% L_values) {
      ll <- setNames(rep(NA_real_, length(n_values)), n_values)
      ## Here L <= np for n=1000 at every requested point.
      ll["1000"] <- feature_loglik(
        C_all[, seq_len(end), drop = FALSE],
        W[seq_len(end), , drop = FALSE], d$Y)
      if (end < N_small) {
        ll["500"] <- feature_loglik(
          C_all[seq_len(n_small), seq_len(end), drop = FALSE],
          W[seq_len(end), , drop = FALSE],
          d$Y[seq_len(n_small), , drop = FALSE])
      } else {
        Kwork <- Ksmall_sum / end
        ll["500"] <- single_lower_loglik(Kwork, y_small)
        rm(Kwork); invisible(gc())
      }

      for (n in n_values) {
        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          family = family_design, p = P, n = n, L = end,
          replicate = replicate_id, method = "simulated_W",
          ll = unname(ll[as.character(n)]), stringsAsFactors = FALSE)
      }
      message(sprintf("%s p=50 rep=%d simulated L=%d (%.1f min)",
                      family_design, replicate_id, end,
                      (proc.time()[["elapsed"]] - started) / 60))
    }
  }

  result <- finish_result(rows, reference)
  metadata <- list(
    stage = "simulated", orientation = "feature_n1000_observed_n500",
    family = family_design, p = P, n_values = n_values, L_values = L_values,
    reference_L = integrated$metadata$reference_L,
    replicate = replicate_id, proposal = proposal, nugget = nugget,
    precision = precision_label,
    elapsed_seconds = proc.time()[["elapsed"]] - started)
  atomic_save(list(result = result, metadata = metadata), simulated_file)
  message("wrote ", simulated_file)
}

run_p50_simulated()
