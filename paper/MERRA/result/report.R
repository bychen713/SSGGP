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
        q <- as.matrix(fit$sigma_inv)   # share of zero off-diagonal entries
        add(design, fold, method, "all", "sparsity", 100 * mean(q[upper.tri(q)] == 0))
      }
      lmc <- readRDS(file.path(out, "lmc", sprintf("fit_fold%d.rds", fold)))
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
  detail <- do.call(rbind, rows)
  write.csv(detail, file.path(merra_output_base, "estimation_by_fold.csv"),
            row.names = FALSE)
  median_rows <- aggregate(estimate ~ design + method + component + parameter,
                           detail, stats::median)
  write.csv(median_rows, file.path(merra_output_base, "estimation_median.csv"),
            row.names = FALSE)

  # Table S8: median across folds of the within-fold mean over processes
  # (LMC: over latent factors for a and nu).
  s8 <- detail[detail$parameter %in% c("a", "phi", "nu", "nugget", "sparsity"), ]
  s8$parameter[s8$parameter == "phi"] <- "a"
  s8 <- aggregate(estimate ~ design + method + fold + parameter, s8, mean)
  s8 <- aggregate(estimate ~ design + method + parameter, s8, stats::median)
  s8 <- reshape(s8, idvar = c("design", "parameter"), timevar = "method",
                direction = "wide")
  names(s8) <- sub("^estimate\\.", "", names(s8))
  s8 <- s8[order(match(s8$design, c("colocated", "misaligned")),
                 match(s8$parameter, c("a", "nu", "nugget", "sparsity"))), ]
  write.csv(s8, file.path(merra_output_base, "estimation_summary.csv"), row.names = FALSE)
}

# Full-grid empirical correlation and fold-median fitted precision matrices.
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
  write_precision_figures(full, figure_dir)
}

# Figure 4: colocated sample precision and the BIC-selected and dense SAECM
# estimates, each an entrywise median over folds, one image per matrix. Each
# colour scale is linear and clipped at the 99th percentile of the matrix's
# absolute off-diagonal entries.
write_precision_figures <- function(full, figure_dir) {
  out <- file.path(merra_output_base, "colocated", "ssggp")
  fold_median <- function(x) apply(simplify2array(x), c(1L, 2L), stats::median)
  sample_q <- dense_q <- bic_q <- vector("list", merra_n_folds)
  for (fold in seq_len(merra_n_folds)) {
    y <- readRDS(file.path(merra_data_root, "colocated", sprintf("fold%d.rds", fold)))$Y_mat
    sample_q[[fold]] <- solve(stats::cov(y))
    bic <- vapply(seq_along(merra_rho), function(ri) {
      ll <- readRDS(file.path(out, sprintf("ll_fold%d_rho%d.rds", fold, ri)))
      -2 * ll$logLik + log(merra_n_train) * ll$edges   # colocated: n = 500
    }, numeric(1))
    est <- function(ri) as.matrix(readRDS(file.path(out,
      sprintf("est_fold%d_rho%d.rds", fold, ri)))$fit$sigma_inv)
    dense_q[[fold]] <- est(1L)
    bic_q[[fold]] <- est(which.min(replace(bic, 1L, Inf)))
  }
  panels <- list(sample = fold_median(sample_q), bic = fold_median(bic_q),
                 dense = fold_median(dense_q))
  labels <- sub("\\(~.*$", "", full$process_names)
  palette <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(201)
  p <- length(labels)
  for (name in names(panels)) {
    q <- panels[[name]]
    limit <- as.numeric(stats::quantile(abs(q[upper.tri(q)]), 0.99))
    path <- file.path(figure_dir, sprintf("precision_colocated_%s.png", name))
    grDevices::png(path, width = 3200, height = 2750, res = 300)
    graphics::layout(matrix(1:2, nrow = 1L), widths = c(1, 0.18))
    graphics::par(mar = c(8.2, 8.8, 1.5, 0.5), pty = "s")
    graphics::image(seq_len(p), seq_len(p), t(pmin(pmax(q, -limit), limit)[p:1L, ]),
                    col = palette, zlim = c(-limit, limit), axes = FALSE,
                    xlab = "", ylab = "", asp = 1)
    graphics::axis(1, at = seq_len(p), labels = labels, las = 2, cex.axis = 1.08,
                   tick = FALSE, line = -0.35, gap.axis = 0)
    graphics::axis(2, at = p - seq_len(p) + 1L, labels = labels, las = 2,
                   cex.axis = 1.08, tick = FALSE, line = -0.35, gap.axis = 0)
    graphics::box()
    graphics::par(mar = c(8.2, 0.5, 1.5, 2.2), pty = "m")   # colour bar
    values <- seq(-limit, limit, length.out = length(palette))
    graphics::image(1, values, matrix(values, nrow = 1L), col = palette,
                    zlim = c(-limit, limit), axes = FALSE, xlab = "", ylab = "")
    graphics::axis(4, at = pretty(c(-limit, limit), n = 6), las = 1, cex.axis = 1.22)
    graphics::mtext(expression(eta), side = 3, line = 0.25, cex = 1.35)
    graphics::box()
    grDevices::dev.off()
  }
}

# Compare raw concentrations on the full grid with the processed values of
# colocated fold 1 at its training and test locations.
write_data_maps <- function(full, figure_dir) {
  fold <- readRDS(file.path(merra_data_root, "colocated", "fold1.rds"))
  processes <- unique(c(1L, ceiling(merra_p / 2), merra_p))
  fold_index <- c(fold$train_index, fold$test_index)
  processed <- rbind(fold$Y_mat, fold$Y_mat_test)
  if (anyNA(processed[, processes, drop = FALSE]))
    stop("Fold contains masked training values; use a colocated fold for this map")
  world <- ggplot2::map_data("world")
  make_panel <- function(j, raw) {
    rows <- if (raw) seq_len(nrow(full$coords)) else fold_index
    d <- data.frame(lon = full$coords[rows, 1L], lat = full$coords[rows, 2L],
                    value = if (raw) full$raw[, j] else processed[, j])
    plot <- ggplot2::ggplot(d, ggplot2::aes(lon, lat, fill = value)) +
      ggplot2::geom_raster(interpolate = FALSE) +
      ggplot2::geom_path(data = world, ggplot2::aes(long, lat, group = group),
                         inherit.aes = FALSE, colour = "grey45", linewidth = 0.25) +
      ggplot2::scale_x_continuous(breaks = c(-95, -85, -75, -65)) +
      ggplot2::scale_y_continuous(breaks = c(25, 35, 45, 50)) +
      ggplot2::coord_quickmap(xlim = range(d$lon), ylim = range(d$lat),
                              expand = FALSE) +
      ggplot2::labs(x = if (raw) NULL else "Longitude", y = "Latitude",
                     fill = if (raw) "kg/kg" else "Detrended\nlog ratio") +
      ggplot2::theme_bw(base_size = 17) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     # grey marks locations outside the fold in processed maps
                     panel.background = ggplot2::element_rect(
                       fill = if (raw) "white" else "#B8B8B8", colour = NA),
                     axis.title = ggplot2::element_text(size = 16),
                     axis.text = ggplot2::element_text(size = 14),
                     legend.position = "right",
                     legend.title = ggplot2::element_text(size = 13),
                     legend.text = ggplot2::element_text(size = 11))
    legend_guide <- ggplot2::guide_colorbar(direction = "vertical",
      title.position = "top", barwidth = grid::unit(0.14, "in"),
      barheight = grid::unit(1.8, "in"))
    if (raw) {
      data_range <- range(d$value, na.rm = TRUE)
      breaks <- 10^seq(log10(data_range[1L]), log10(data_range[2L]), length.out = 3L)
      plot + ggplot2::scale_fill_viridis_c(option = "C",
        trans = "log10", breaks = breaks, labels = scales::label_scientific(digits = 1),
        guide = legend_guide)
    } else {
      limit <- max(abs(d$value), na.rm = TRUE)
      plot + ggplot2::scale_fill_gradientn(
        colors = c("#67001F", "#B2182B", "#E66B5B", "#FFFFFF",
                   "#6798D0", "#2166AC", "#053061"),
        values = c(0, 0.25, 0.375, 0.5, 0.625, 0.75, 1),
        limits = c(-limit, limit), breaks = c(-limit, 0, limit),
        labels = scales::label_number(accuracy = 0.1), guide = legend_guide)
    }
  }
  panels <- c(lapply(processes, make_panel, raw = TRUE),
              lapply(processes, make_panel, raw = FALSE))
  ggplot2::ggsave(file.path(figure_dir, "data_maps.png"),
    patchwork::wrap_plots(panels, ncol = length(processes), guides = "keep"),
    width = 4.8 * length(processes), height = 6.6, dpi = 300)
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
