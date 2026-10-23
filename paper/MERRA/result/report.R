# Fitted-parameter tables and descriptive figures from saved MERRA results.

# Parameter estimates in each model's native parameterization; no true MERRA
# parameters are available for an estimation-error calculation.
write_estimation_table <- function() {
  rows <- list()
  add <- function(design, fold, method, component, parameter, value) {
    rows[[length(rows) + 1L]] <<- data.frame(
      design, fold, method, component, parameter, estimate = as.numeric(value))
  }
  for (design in c("misaligned", "colocated")) {
    out <- file.path(merra_output_base, design)
    for (fold in seq_len(merra_n_folds)) {
      dat <- readRDS(file.path(merra_data_root, design, sprintf("fold%d.rds", fold)))
      ind <- readRDS(file.path(out, "ssggp", sprintf("ind_%d.rds", fold)))$fit
      fits <- list(IND = ind)
      for (method in c("dense", "BIC")) {
        pred <- readRDS(file.path(out, "ssggp",
          sprintf("pred_%s_fold%d%s.rds", method, fold, merra_pred_suffix)))
        ri <- if (method == "dense") 1L else
          which.min(replace(pred$bic, 1L, Inf))
        fits[[method]] <- readRDS(file.path(out, "ssggp",
          sprintf("est_fold%d_rho%d.rds", fold, ri)))$fit
      }
      for (method in names(fits)) {
        fit <- fits[[method]]
        for (j in seq_len(merra_p)) {
          par <- c(fit$f_par[[j]], nugget = fit$nugget[j],
                   precision_diagonal = fit$sigma_inv[j, j])
          add(design, fold, method, dat$process_names[j], names(par), par)
        }
      }
      lmc <- readRDS(file.path(out, "lmc_process_nugget", sprintf("fit_fold%d.rds", fold)))
      for (k in seq_along(lmc$est$par)) {
        par <- setNames(lmc$est$par[[k]], c("a", "nu"))
        add(design, fold, "LMC_Matern", paste0("factor", k), names(par), par)
        add(design, fold, "LMC_Matern", dat$process_names,
            paste0("loading_factor", k), lmc$pilot$ld$B[, k])
      }
      add(design, fold, "LMC_Matern", dat$process_names, "nugget", lmc$est$nugget)
      # IOX: posterior means over the iterations used for prediction
      keep <- readRDS(file.path(out, "iox", sprintf("pred_fold%d.rds", fold)))$keep
      fit <- readRDS(file.path(out, "iox", sprintf("fit_fold%d.rds", fold)))$fit
      par <- rbind(phi = rowMeans(fit$Theta[1L, , keep, drop = FALSE][1L, , ]),
                   nu = rowMeans(fit$Theta[2L, , keep, drop = FALSE][1L, , ]),
                   nugget = rowMeans(fit$Ddiag[, keep, drop = FALSE]))
      rm(fit)
      for (j in seq_len(merra_p))
        add(design, fold, "IOX", dat$process_names[j], rownames(par), par[, j])
    }
  }
  write.csv(do.call(rbind, rows),
            file.path(merra_output_base, "estimation_by_fold.csv"), row.names = FALSE)
}

# Full-grid empirical correlation and first-fold fitted precision matrices.
write_matrix_figures <- function(full, figure_dir) {
  heatmap <- function(x, title, limits) {
    p <- nrow(x)
    d <- expand.grid(row = seq_len(p), col = seq_len(p))
    d$value <- as.vector(x)
    ggplot2::ggplot(d, ggplot2::aes(col, row, fill = value)) +
      ggplot2::geom_raster() + ggplot2::coord_equal() +
      ggplot2::scale_y_reverse(breaks = seq_len(p), labels = full$process_names) +
      ggplot2::scale_x_continuous(breaks = seq_len(p), labels = full$process_names) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                    limits = limits, midpoint = 0) +
      ggplot2::labs(title = title, x = NULL, y = NULL, fill = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       size = if (p > 10) 5 else 9),
                     axis.text.y = ggplot2::element_text(size = if (p > 10) 5 else 9),
                     panel.grid = ggplot2::element_blank())
  }
  ggplot2::ggsave(file.path(figure_dir, "correlation.png"),
    heatmap(cor(full$normalized), "Empirical full-grid residual correlation", c(-1, 1)),
    width = 9, height = 8, dpi = 300)
  matrices <- list()
  for (design in c("misaligned", "colocated")) {
    out <- file.path(merra_output_base, design, "ssggp")
    for (method in c("dense", "BIC")) {
      pred <- readRDS(file.path(out,
        sprintf("pred_%s_fold1%s.rds", method, merra_pred_suffix)))
      ri <- if (method == "dense") 1L else which.min(replace(pred$bic, 1L, Inf))
      est <- readRDS(file.path(out, sprintf("est_fold1_rho%d.rds", ri)))
      title <- sprintf("%s: %s (fold 1, rho=%g)", design, method, est$rho)
      matrices[[title]] <- est$fit$sigma_inv
    }
  }
  limit <- max(abs(unlist(matrices)))
  panels <- Map(function(x, title) heatmap(x, title, c(-limit, limit)),
                matrices, names(matrices))
  figure <- patchwork::wrap_plots(panels, ncol = 2, guides = "collect")
  ggplot2::ggsave(file.path(figure_dir, "precision.png"), figure,
                  width = 15, height = 13, dpi = 300)
}

# Compare raw concentrations on the full grid with the processed values of
# colocated fold 1 at its training and test locations.
write_data_maps <- function(full, figure_dir) {
  fold <- readRDS(file.path(merra_data_root, "colocated", "fold1.rds"))
  processes <- unique(c(1L, ceiling(merra_p / 2), merra_p))
  fold_index <- c(fold$train_index, fold$test_index)
  processed <- rbind(fold$Y_mat, fold$Y_mat_test)
  processed_limit <- max(abs(processed[, processes, drop = FALSE]))
  world <- ggplot2::map_data("world")
  make_panel <- function(j, raw) {
    rows <- if (raw) seq_len(nrow(full$coords)) else fold_index
    d <- data.frame(lon = full$coords[rows, 1L], lat = full$coords[rows, 2L],
                    value = if (raw) full$raw[, j] else processed[, j])
    plot <- ggplot2::ggplot(d, ggplot2::aes(lon, lat, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::geom_path(data = world, ggplot2::aes(long, lat, group = group),
                         inherit.aes = FALSE, colour = "grey45", linewidth = 0.25) +
      ggplot2::coord_quickmap(xlim = range(full$coords[, 1L]),
                              ylim = range(full$coords[, 2L]), expand = FALSE) +
      ggplot2::labs(title = full$process_names[j], x = "Longitude", y = "Latitude",
                     subtitle = if (raw) "Raw concentration" else "Processed (colocated fold 1)",
                     fill = if (raw) "kg/kg" else "Standardized\nlog residual") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     # grey marks locations outside the fold in processed maps
                     panel.background = ggplot2::element_rect(
                       fill = if (raw) "white" else "grey80", colour = NA))
    if (raw) {
      trans <- scales::pseudo_log_trans(sigma = full$offset[j], base = 10)
      breaks <- trans$inverse(seq(trans$transform(min(d$value)),
                                   trans$transform(max(d$value)), length.out = 3L))
      plot + ggplot2::scale_fill_viridis_c(option = "C",
        trans = trans, breaks = breaks, labels = scales::label_scientific(digits = 1))
    } else {
      plot + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
        high = "#B2182B", limits = c(-processed_limit, processed_limit),
        breaks = c(-processed_limit, 0, processed_limit),
        labels = scales::label_number(accuracy = 0.1), midpoint = 0)
    }
  }
  panels <- c(lapply(processes, make_panel, raw = TRUE),
              lapply(processes, make_panel, raw = FALSE))
  ggplot2::ggsave(file.path(figure_dir, "data_maps.png"),
    patchwork::wrap_plots(panels, ncol = length(processes)),
    width = 6 * length(processes), height = 8, dpi = 300)
}

write_merra_report <- function() {
  write_estimation_table()
  full <- readRDS(merra_source)
  stopifnot(ncol(full$normalized) == merra_p)
  figure_dir <- file.path(merra_output_base, "figure")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  write_matrix_figures(full, figure_dir)
  write_data_maps(full, figure_dir)
}
