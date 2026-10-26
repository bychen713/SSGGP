#!/usr/bin/env Rscript
# Precision-recall curve (Figure 2) and a simulated-field example (Figure S1).
# Run from paper/ after data generation and the penalty-path fits.
if (!file.exists("sim/settings.R")) stop("run from paper/", call. = FALSE)
source("sim/settings.R")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
truth <- true_precision()
edges <- upper.tri(truth)
positive <- truth[edges] != 0
if (!any(positive)) stop("true graph has no edges", call. = FALSE)

curve <- do.call(rbind, lapply(seq_len(nrow(grid_params)), function(i) {
  path <- file.path(out_dir, sprintf("est%d.RData", i))
  if (!file.exists(path)) stop("missing fit: ", path, call. = FALSE)
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  selected <- abs(e$result_est$sigma_inv_est[edges]) > 1e-8
  data.frame(sim = grid_params$sim[i], rho = grid_params$rho[i],
             precision = if (any(selected)) sum(selected & positive) / sum(selected) else 1,
             recall = sum(selected & positive) / sum(positive))
}))
write.csv(curve, file.path(fig_dir, "precision_recall_by_replicate.csv"),
          row.names = FALSE)
mean_curve <- aggregate(cbind(precision, recall) ~ rho, curve, mean)
plot <- ggplot2::ggplot(mean_curve, ggplot2::aes(recall, precision)) +
  ggplot2::geom_path(colour = "#145A86") +
  ggplot2::geom_point(ggplot2::aes(colour = rho == 0), size = 2) +
  ggplot2::scale_colour_manual(values = c("FALSE" = "#145A86", "TRUE" = "#B4413C"),
                               guide = "none") +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(x = "Recall", y = "Precision") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(fig_dir, "precision_recall.png"), plot,
                width = 5, height = 5, dpi = 300)

data_path <- file.path(data_dir, "1.RData")
if (!file.exists(data_path)) stop("missing simulated data: ", data_path,
                                  call. = FALSE)
e <- new.env(parent = emptyenv())
load(data_path, envir = e)
show <- seq_len(min(4L, ncol(e$Y_mat)))
field <- do.call(rbind, lapply(show, function(j) {
  data.frame(x = c(e$grid[, 1], e$grid_test[, 1]),
             y = c(e$grid[, 2], e$grid_test[, 2]),
             value = c(e$Y_mat[, j], e$Y_mat_test[, j]),
             process = paste("Process", j))
}))
plot <- ggplot2::ggplot(field, ggplot2::aes(x, y, colour = value)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::facet_wrap(~ process, ncol = 2) +
  ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "white",
                                  high = "#B2182B", midpoint = 0, na.value = "grey85") +
  ggplot2::coord_equal() + ggplot2::theme_bw() +
  ggplot2::labs(x = NULL, y = NULL, colour = "Simulated value")
ggplot2::ggsave(file.path(fig_dir, "simulated_fields.png"), plot,
                width = 8, height = 7, dpi = 300)
