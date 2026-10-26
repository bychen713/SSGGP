#!/usr/bin/env Rscript
# Estimation/prediction tables, fitted precision, empirical correlation, and maps.
# Run from paper/ once all fits have finished:  Rscript MERRA/result/summarize.R
#
#   MERRA/output/prediction_by_fold.csv   one row per design x method x fold
#   MERRA/output/prediction_summary.csv   Table 4: mean (SD) over folds
#
# Methods: IND (independent processes), SSGGP dense and BIC-selected, LMC, and
# IOX. Columns a method does not report are NA.
source("MERRA/settings.R")

read_metrics <- function(path, design) {
  if (!file.exists(path)) stop("missing result: ", path, call. = FALSE)
  m <- readRDS(path)$metrics
  cbind(design = design, m)
}

# IOX metrics from MERRA/iox/run.R's prediction file, in the columns used above.
# Intervals are the 2.5% and 97.5% posterior predictive quantiles (type 7), for
# each process and for the sum across processes within each draw.
read_iox <- function(out, fold, design) {
  path <- file.path(out, "iox", sprintf("pred_fold%d.rds", fold))
  if (!file.exists(path)) stop("missing result: ", path, call. = FALSE)
  p <- readRDS(path); m <- p$metrics
  Y_out <- as.matrix(readRDS(file.path(merra_data_root, design,
                                       sprintf("fold%d.rds", fold)))$Y_mat_test)
  q <- function(x, prob) quantile(x, prob, type = 7, names = FALSE)
  lo <- apply(p$draws, 1:2, q, 0.025); hi <- apply(p$draws, 1:2, q, 0.975)
  sums <- apply(p$draws, c(1L, 3L), sum)
  lo_sum <- apply(sums, 1L, q, 0.025); hi_sum <- apply(sums, 1L, q, 0.975)
  truth_sum <- rowSums(Y_out)
  data.frame(design = design, fold = fold, method = "IOX",
             mse = mean((Y_out - p$mean)^2), norm = m$norm,
             coverage = mean(Y_out >= lo & Y_out <= hi),
             pi_length = mean(hi - lo), sum_mse = m$rmse_sum^2,
             sum_coverage = mean(truth_sum >= lo_sum & truth_sum <= hi_sum),
             sum_pi_length = mean(hi_sum - lo_sum),
             estimation_min = m$total_estimate_minutes,
             prediction_min = m$predict_minutes)
}

rows <- list()
for (design in c("misaligned", "colocated")) {
  out <- file.path(merra_output_base, design)
  for (fold in seq_len(merra_n_folds)) {
    rows[[length(rows) + 1L]] <- read_metrics(
      file.path(out, "ssggp", sprintf("ind_%d.rds", fold)), design)
    for (method in c("dense", "BIC"))
      rows[[length(rows) + 1L]] <- read_metrics(
        file.path(out, "ssggp",
                  sprintf("pred_%s_fold%d%s.rds", method, fold, merra_pred_suffix)),
        design)
    rows[[length(rows) + 1L]] <- read_metrics(
      file.path(out, "lmc", sprintf("fit_fold%d.rds", fold)), design)
    rows[[length(rows) + 1L]] <- read_iox(out, fold, design)
  }
}

all_names <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(x) {
  for (nm in setdiff(all_names, names(x))) x[[nm]] <- NA
  x[all_names]
})
detail <- do.call(rbind, rows)
detail$sum_rmse <- sqrt(detail$sum_mse)   # per fold, before averaging

# Table 4: mean (SD) over folds.
metrics <- c(norm = "Norm", coverage = "Coverage", sum_rmse = "Sum RMSE",
             sum_coverage = "Sum coverage", pi_length = "PI", sum_pi_length = "PI sum")
summary <- do.call(rbind, lapply(split(detail, detail[c("design", "method")], drop = TRUE),
  function(x) data.frame(design = x$design[1], method = x$method[1],
    metric = unname(metrics),
    value = sprintf("%.3f (%.3f)", vapply(names(metrics), function(m) mean(x[[m]]), 1),
                    vapply(names(metrics), function(m) stats::sd(x[[m]]), 1)))))
summary <- reshape(summary, idvar = c("design", "metric"), timevar = "method",
                   direction = "wide")
names(summary) <- sub("^value\\.", "", names(summary))
summary <- summary[order(match(summary$design, c("colocated", "misaligned")),
                         match(summary$metric, metrics)), ]

write.csv(detail, file.path(merra_output_base, "prediction_by_fold.csv"), row.names = FALSE)
write.csv(summary, file.path(merra_output_base, "prediction_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)

source("MERRA/result/report.R")
write_merra_report()
