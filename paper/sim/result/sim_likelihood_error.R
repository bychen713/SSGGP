#!/usr/bin/env Rscript
# Figure S2: random-feature likelihood error for the first BIC-selected
# simulation fit. A larger-feature Monte Carlo median is the reference, not an
# analytic likelihood. Run from paper/ after the selected fit exists.
if (!file.exists("sim/settings.R")) stop("run from paper/", call. = FALSE)
source("sim/settings.R")
setup_compute()

sim <- 1L
path <- file.path(out_dir, sprintf("BIC_%s_%d.RData", setting_name, sim))
if (!file.exists(path)) stop("missing selected fit: ", path, call. = FALSE)
e <- new.env(parent = emptyenv())
load(path, envir = e)
rows <- which(grid_params$sim == sim & grid_params$rho == e$result_final$rho)
if (length(rows) != 1L) stop("selected rho is not in the fit grid", call. = FALSE)
path <- file.path(out_dir, sprintf("est%d.RData", rows))
if (!file.exists(path)) stop("missing estimate: ", path, call. = FALSE)
load(path, envir = e)
fit <- e$result_est$fit
if (!inherits(fit, "ssggp_fit")) stop("invalid saved fit", call. = FALSE)

feature_counts <- c(100L, 250L, 500L, 1000L, 2000L)
reference_features <- 5000L
seeds <- seq_len(5L)
evaluate <- function(L, seed) as.numeric(stats::logLik(
  fit, L = L, seed = seed, g = proposal, threads = threads))
reference <- stats::median(vapply(seeds, function(seed)
  evaluate(reference_features, seed), numeric(1)))
n_observed <- sum(!is.na(fit$data$Y))
rows <- expand.grid(L = feature_counts, seed = seeds)
rows$logLik <- mapply(evaluate, rows$L, rows$seed)
rows$absolute_error_per_observation <- abs(rows$logLik - reference) / n_observed
rows$reference_features <- reference_features
rows$reference_logLik <- reference
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(rows, file.path(fig_dir, "likelihood_error.csv"), row.names = FALSE)

plot <- ggplot2::ggplot(rows,
  ggplot2::aes(factor(L), absolute_error_per_observation)) +
  ggplot2::geom_boxplot(fill = "#B9D5E9", outlier.shape = NA) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1),
                      alpha = 0.5, size = 1) +
  ggplot2::labs(x = "Random Fourier features",
                y = "Absolute log-likelihood error per observation") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(fig_dir, "likelihood_error.png"), plot,
                width = 6, height = 4, dpi = 300)
