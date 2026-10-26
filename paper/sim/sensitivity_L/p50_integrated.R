#!/usr/bin/env Rscript

## Split p=50 mixed-family integrated-W work into scheduler-safe stages:
##   curves    -- nine plotted L values through 50,000;
##   reference -- stream once to L=100,000 and factor only there;
##   merge     -- attach the reference and install the standard integrated RDS.

if (Sys.getenv("SSGGP_P", "") != "50" ||
    Sys.getenv("SSGGP_FAMILY", "") != "mixed")
  stop("this entry point requires SSGGP_P=50 and SSGGP_FAMILY=mixed")
p50_stage <- Sys.getenv("SSGGP_P50_STAGE", "curves")
if (!p50_stage %in% c("curves", "reference", "merge"))
  stop("SSGGP_P50_STAGE must be curves, reference, or merge")
Sys.setenv(SSGGP_STAGE = "selftest")
source("sim/sensitivity_L/highp.R", local = FALSE)

curves_file <- file.path(
  out_dir, sprintf("integrated_curves_replicate_%02d.rds", replicate_id))
reference_file <- file.path(
  out_dir, sprintf("reference_replicate_%02d.rds", replicate_id))

run_curves <- function() {
  started <- proc.time()[["elapsed"]]
  d <- generate_design(replicate_id, family_design)
  pairs <- triangular_pairs(P)
  blocks <- new_integrated_blocks(pairs, max(n_values))
  rows <- list(); row_id <- 0L
  for (start in seq.int(1L, max(L_values), by = 500L)) {
    end <- min(start + 499L, max(L_values))
    idx <- start:end
    C <- cos(sweep(d$coords %*% t(d$features$X[idx, , drop = FALSE]),
                   2L, d$features$phase[idx], "+"))
    blocks <- update_integrated_blocks(
      blocks, C, d$r[idx, , drop = FALSE], pairs)
    if (end %in% L_values) {
      K <- assemble_integrated_nested(
        blocks, d$Sigma, pairs, end, min(n_values), max(n_values))
      ll <- nested_loglik(K, d$Y, min(n_values))
      rm(K); invisible(gc())
      for (n in n_values) {
        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          family = family_design, p = P, n = n, L = end,
          replicate = replicate_id, method = "integrated_W",
          ll = unname(ll[as.character(n)]), stringsAsFactors = FALSE)
      }
      message(sprintf("mixed p=50 rep=%d curves L=%d (%.1f min)",
                      replicate_id, end,
                      (proc.time()[["elapsed"]] - started) / 60))
    }
  }
  metadata <- list(
    stage = "integrated_curves", family = family_design, p = P,
    n_values = n_values, L_values = L_values, replicate = replicate_id,
    proposal = proposal, nugget = nugget,
    precision = precision_label,
    process_family = d$importance$family,
    process_parameters = d$importance$parameters,
    capped = d$importance$capped,
    simulation_capped = d$simulation_capped,
    elapsed_seconds = proc.time()[["elapsed"]] - started)
  atomic_save(list(result = do.call(rbind, rows), metadata = metadata),
              curves_file)
  message("wrote ", curves_file)
}

run_reference <- function() {
  started <- proc.time()[["elapsed"]]
  d <- generate_design(replicate_id, family_design)
  pairs <- triangular_pairs(P)
  blocks <- new_integrated_blocks(pairs, max(n_values))
  for (start in seq.int(1L, L_reference_mixed, by = 500L)) {
    end <- min(start + 499L, L_reference_mixed)
    idx <- start:end
    C <- cos(sweep(d$coords %*% t(d$features$X[idx, , drop = FALSE]),
                   2L, d$features$phase[idx], "+"))
    blocks <- update_integrated_blocks(
      blocks, C, d$r[idx, , drop = FALSE], pairs)
    if (end %% 10000L == 0L)
      message(sprintf("mixed p=50 rep=%d reference accumulation L=%d (%.1f min)",
                      replicate_id, end,
                      (proc.time()[["elapsed"]] - started) / 60))
  }
  K <- assemble_integrated_nested(
    blocks, d$Sigma, pairs, L_reference_mixed,
    min(n_values), max(n_values))
  reference <- nested_loglik(K, d$Y, min(n_values))
  metadata <- list(
    stage = "reference", family = family_design, p = P,
    reference_L = L_reference_mixed, replicate = replicate_id,
    elapsed_seconds = proc.time()[["elapsed"]] - started)
  atomic_save(list(reference = reference, metadata = metadata), reference_file)
  message("wrote ", reference_file)
}

run_merge <- function() {
  if (!file.exists(curves_file) || !file.exists(reference_file))
    stop("curves and reference stages must both finish before merge")
  curves <- readRDS(curves_file)
  ref <- readRDS(reference_file)
  result <- curves$result
  result$ll_ref <- unname(ref$reference[as.character(result$n)])
  result$abs_error <- abs(result$ll - result$ll_ref)
  result$error_per_n <- result$abs_error / result$n
  metadata <- curves$metadata
  metadata$stage <- "integrated"
  metadata$reference_L <- L_reference_mixed
  metadata$reference <- ref$reference
  metadata$reference_elapsed_seconds <- ref$metadata$elapsed_seconds
  atomic_save(list(result = result, metadata = metadata), integrated_file)
  message("wrote ", integrated_file)
}

switch(p50_stage,
       curves = run_curves(),
       reference = run_reference(),
       merge = run_merge())
