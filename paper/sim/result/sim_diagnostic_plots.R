# Diagnostic traces for the BIC-selected simulation fits. This file is sourced
# by sim_result.R after the ordinary result tables and figures are produced.

.read_rdata_object <- function(path, name) {
  if (!file.exists(path)) stop("missing file: ", path, call. = FALSE)
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists(name, envir = e, inherits = FALSE))
    stop("missing `", name, "` in ", path, call. = FALSE)
  e[[name]]
}

.unpack_precision <- function(packed, P) {
  Q <- matrix(0, P, P)
  upper <- upper.tri(Q, diag = TRUE)
  if (length(packed) != sum(upper))
    stop("invalid packed precision length", call. = FALSE)
  Q[upper] <- packed
  Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
  Q
}

.microergodic_trace <- function(fit, family) {
  chain <- fit$chain
  P <- fit$spec$P
  process <- which(fit$spec$family == family)
  if (!length(process)) return(NULL)
  n_iter <- length(chain$spectral)
  out <- matrix(NA_real_, nrow = length(process), ncol = n_iter,
                dimnames = list(as.character(process), NULL))
  for (i in seq_len(n_iter)) {
    Q <- .unpack_precision(chain$sigma_inv[, i], P)
    marginal_variance <- diag(chol2inv(chol(Q)))[process]
    pars <- chain$spectral[[i]][process]
    if (family == "matern") {
      a <- vapply(pars, `[[`, numeric(1), "a")
      nu <- vapply(pars, `[[`, numeric(1), "nu")
      out[, i] <- marginal_variance * a^(2 * nu)
    } else if (family == "ch") {
      beta <- vapply(pars, `[[`, numeric(1), "beta")
      alpha <- vapply(pars, `[[`, numeric(1), "alpha")
      nu <- vapply(pars, `[[`, numeric(1), "nu")
      # High-frequency tail under the package's beta^2 |omega|^2 / 2
      # parameterisation.
      out[, i] <- exp(log(marginal_variance) + nu * log(2) +
                        lgamma(alpha + nu) - lgamma(alpha) -
                        2 * nu * log(beta))
    } else {
      stop("no polynomial-tail microergodic scale for family `", family, "`",
           call. = FALSE)
    }
  }
  if (any(!is.finite(out)) || any(out <= 0))
    stop("non-finite microergodic trace", call. = FALSE)
  out
}

.precision_frobenius_trace <- function(fit) {
  P <- fit$spec$P
  pairs <- which(upper.tri(matrix(0, P, P), diag = TRUE), arr.ind = TRUE)
  weights <- ifelse(pairs[, 1] == pairs[, 2], 1, 2)
  values <- sqrt(colSums(fit$chain$sigma_inv^2 * weights))
  terminal <- sqrt(sum(as.matrix(fit$sigma_inv)^2))
  if (any(!is.finite(values)) ||
      !isTRUE(all.equal(tail(values, 1), terminal, tolerance = 1e-10)))
    stop("invalid precision Frobenius-norm trace", call. = FALSE)
  values
}

.draw_trace_plot <- function(path, values, ylab,
                             iterations = seq_len(nrow(values)), title = NULL) {
  values <- as.matrix(values)
  if (nrow(values) != length(iterations) || any(!is.finite(values)))
    stop("invalid values supplied for ", basename(path), call. = FALSE)
  center <- rowMeans(values)
  spread <- if (ncol(values) > 1L) apply(values, 1L, stats::sd) else
    rep(0, nrow(values))
  limits <- range(values, center - spread, center + spread, finite = TRUE)
  pad <- max(diff(limits) * 0.04, 0.002)
  limits <- limits + c(-pad, pad)
  if (!is.null(path)) {
    grDevices::png(path, width = 1950, height = 1500, res = 300,
                   type = "cairo")
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  # Place the y label just beyond the widest tick label.
  label_line <- 1.6 + 0.85 * max(nchar(format(pretty(limits))))
  graphics::par(mar = c(5.0, label_line + 2.2, 0.7, 2.5),
                mgp = c(3.5, 0.75, 0), las = 1)
  graphics::plot(NA_real_, xlim = range(iterations), ylim = limits,
                 xaxs = "i", xlab = "Iteration", ylab = "",
                 cex.axis = if (is.null(path)) 0.8 else 1.5,
                 cex.lab = if (is.null(path)) 0.8 else 1.4,
                 main = title)
  graphics::mtext(ylab, side = 2, line = label_line, las = 0,
                  cex = if (is.null(path)) 0.8 else
                    min(1.4, 1.4 * 36 / nchar(ylab)))
  graphics::abline(h = pretty(limits, 4), col = "#ECECEC", lwd = 0.6)
  graphics::polygon(c(iterations, rev(iterations)),
                    c(center - spread, rev(center + spread)),
                    col = "#B9D5E9", border = NA)
  for (j in seq_len(ncol(values)))
    graphics::lines(iterations, values[, j],
                    col = grDevices::adjustcolor("#5C7E98", alpha.f = 0.22),
                    lwd = 0.6)
  graphics::lines(iterations, center, col = "#145A86", lwd = 2.1)
  if (length(iterations) == 1L)   # a single point draws no line
    graphics::points(iterations, center, pch = 19, col = "#145A86")
  graphics::box(col = "#777777")
  invisible(path)
}

.nugget_trace <- function(fit) {
  nugget <- as.matrix(fit$chain$nugget)
  colMeans(nugget)
}

.draw_diagnostic_grid <- function(fig_dir, setting_name, traces) {
  trace_path <- file.path(fig_dir, "diagnostic_traces.rds")
  saveRDS(traces, trace_path)
  other <- if (grepl("_miss", setting_name)) sub("_miss", "_coloc", setting_name) else
    sub("_coloc", "_miss", setting_name)
  # fig_dir is sim/output/figure/<setting>/diagnostics.
  other_path <- file.path(dirname(dirname(fig_dir)), other, "diagnostics",
                          "diagnostic_traces.rds")
  if (!file.exists(other_path)) return(invisible(NULL))
  pair <- list(miss = if (grepl("_miss", setting_name)) traces else readRDS(other_path),
               coloc = if (grepl("_coloc", setting_name)) traces else readRDS(other_path))
  path <- file.path(dirname(dirname(fig_dir)),
                    sub("_(miss|coloc)", "", setting_name), "diagnostics_2x4.png")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = 2400, height = 1200, res = 200, type = "cairo")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 4), oma = c(0, 0, 0, 0))
  for (design in c("coloc", "miss")) {
    for (name in c("matern", "ch", "precision", "nugget")) {
      values <- pair[[design]][[name]]
      if (is.null(values)) {
        graphics::plot.new()
        graphics::title(main = paste(design, name))
      } else {
        .draw_trace_plot(NULL, values, ylab = name,
                         title = paste(design, name))
      }
    }
  }
  invisible(path)
}

.fit_at_iteration <- function(fit, iteration) {
  chain <- fit$chain
  P <- fit$spec$P
  fit$sigma_inv <- .unpack_precision(chain$sigma_inv[, iteration], P)
  fit$nugget <- if (identical(fit$spec$nugget, "shared"))
    chain$nugget[1L, iteration] else chain$nugget[, iteration]
  fit$f_par <- chain$spectral[[iteration]]
  fit$chain <- NULL
  fit
}

.penalized_observed_ll <- function(fit, iteration, L_ll, seed_ll, threads) {
  state <- .fit_at_iteration(fit, iteration)
  observed_ll <- as.numeric(stats::logLik(state, L = L_ll, seed = seed_ll,
                                          threads = threads))
  off_diagonal <- row(state$sigma_inv) != col(state$sigma_inv)
  penalty <- state$spec$L * state$spec$rho / 2 *
    sum(abs(state$sigma_inv[off_diagonal]))
  c(observed_ll = observed_ll, penalty = penalty,
    penalized_ll = observed_ll - penalty)
}

.observed_ll_identity <- function(estimate_md5, selection_md5, rho, iteration,
                                  L_ll, seed_ll) {
  list(version = 1L, estimate_md5 = estimate_md5,
       selection_md5 = selection_md5, rho = rho,
       iteration = as.integer(iteration), L_ll = as.integer(L_ll),
       seed_ll = as.integer(seed_ll))
}

.collect_penalized_observed_ll <- function(selected, diagnostic_dir,
                                           iterations, L_ll, seed_ll, threads,
                                           compute) {
  cache_dir <- file.path(diagnostic_dir,
                         sprintf("penalized_observed_ll_L%d_seed%d",
                                 L_ll, seed_ll))
  paths <- unlist(lapply(seq_along(selected), function(r)
    file.path(cache_dir, sprintf("rep%02d_iter%04d.rds", r, iterations))))
  if (compute) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  k <- 0L
  for (r in seq_along(selected)) for (iteration in iterations) {
    k <- k + 1L
    z <- selected[[r]]
    identity <- .observed_ll_identity(z$estimate_md5, z$selection_md5, z$rho,
                                      iteration, L_ll, seed_ll)
    valid <- FALSE
    if (file.exists(paths[[k]])) {
      saved <- readRDS(paths[[k]])
      valid <- is.list(saved) && identical(saved$identity, identity) &&
        is.numeric(saved$value) && length(saved$value) == 3L &&
        all(is.finite(saved$value))
      if (!valid && !compute)
        stop("stale diagnostic cache: ", paths[[k]], call. = FALSE)
    }
    if (!valid && compute) {
      value <- .penalized_observed_ll(z$fit, iteration, L_ll, seed_ll,
                                      threads)
      tmp <- paste0(paths[[k]], ".tmp-", Sys.getpid())
      saveRDS(list(identity = identity, value = value), tmp)
      if (!file.rename(tmp, paths[[k]])) {
        unlink(tmp)
        stop("could not save diagnostic cache: ", paths[[k]], call. = FALSE)
      }
      message("computed penalized observed likelihood: replicate ", r,
              ", iteration ", iteration)
    }
  }
  if (any(!file.exists(paths))) return(NULL)
  saved <- lapply(paths, readRDS)
  matrix(vapply(saved, function(x) unname(x$value[["penalized_ll"]]),
                numeric(1)), nrow = length(iterations))
}

#' Plot simulation diagnostics from the selected SAECM chains
#'
#' Microergodic scales and precision norms use the saved chains and are cheap.
#' Penalized observed log-likelihoods require a fresh marginal-likelihood
#' factorization at every plotted iteration. They are cached one point at a time
#' and are computed only when `compute_observed_ll` is true.
plot_simulation_diagnostics <- function(out_dir, fig_dir, setting_name,
                                        grid_params, n_sim, P, family_names,
                                        n_iter, L_ll,
                                        seed_ll, threads,
                                        compute_observed_ll = FALSE,
                                        observed_ll_every = 50L) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  selected <- lapply(seq_len(n_sim), function(r) {
    selected_path <- file.path(
      out_dir, sprintf("BIC_%s_%d.RData", setting_name, r))
    result_final <- .read_rdata_object(selected_path, "result_final")
    rows <- which(grid_params$sim == r)
    selected_row <- rows[match(result_final$rho, grid_params$rho[rows])]
    if (length(selected_row) != 1L || is.na(selected_row))
      stop("selected rho is not in grid_params for replicate ", r,
           call. = FALSE)
    est_path <- file.path(out_dir, sprintf("est%d.RData", selected_row))
    result_est <- .read_rdata_object(est_path, "result_est")
    fit <- result_est$fit
    if (!inherits(fit, "ssggp_fit") || fit$spec$P != P ||
        !identical(fit$spec$family, family_names) ||
        !isTRUE(all.equal(result_est$rho, result_final$rho, tolerance = 0)) ||
        !isTRUE(all.equal(fit$sigma_inv, result_final$sigma_inv_est,
                          tolerance = 1e-10, check.attributes = FALSE)))
      stop("selected result does not match estimate for replicate ", r,
           call. = FALSE)
    list(fit = fit, rho = result_final$rho,
         estimate_md5 = unname(tools::md5sum(est_path)),
         selection_md5 = unname(tools::md5sum(selected_path)))
  })
  fits <- lapply(selected, `[[`, "fit")
  if (any(vapply(fits, function(x) length(x$chain$spectral), integer(1)) != n_iter))
    stop("expected ", n_iter, " iterations in every selected fit", call. = FALSE)

  written <- character()
  trace_grid <- list()
  for (family in intersect(c("matern", "ch"), unique(family_names))) {
    traces <- lapply(fits, .microergodic_trace, family = family)
    # Average replicate-specific log10 scales for each process, retaining the
    # between-process trajectories in the plot.
    values <- t(Reduce(`+`, lapply(traces, log10)) / length(traces))
    path <- file.path(fig_dir, paste0(family, "_microergodic.png"))
    .draw_trace_plot(path, values, ylab = "log10 microergodic scale")
    trace_grid[[family]] <- values
    written <- c(written, path)
  }

  precision <- vapply(fits, .precision_frobenius_trace, numeric(n_iter))
  path <- file.path(fig_dir, "precision_frobenius.png")
  .draw_trace_plot(path, precision, ylab = "Precision Frobenius norm")
  trace_grid$precision <- precision
  written <- c(written, path)

  nugget <- vapply(fits, .nugget_trace, numeric(n_iter))
  path <- file.path(fig_dir, "nugget.png")
  .draw_trace_plot(path, nugget, ylab = "Nugget variance")
  trace_grid$nugget <- nugget
  written <- c(written, path)
  .draw_diagnostic_grid(fig_dir, setting_name, trace_grid)

  # Runs shorter than observed_ll_every are evaluated at their last iteration.
  observed_iterations <- if (n_iter >= observed_ll_every)
    seq.int(observed_ll_every, n_iter, by = observed_ll_every) else n_iter
  values <- .collect_penalized_observed_ll(
    selected, file.path(out_dir, "diagnostics"),
    observed_iterations, L_ll, seed_ll, threads, compute_observed_ll)
  if (is.null(values)) {
    message("penalized observed-likelihood plot skipped: cached evaluations ",
            "are incomplete. Set SSGGP_COMPUTE_DIAGNOSTIC_LOGLIK=true and ",
            "rerun sim_result.R to compute them.")
  } else {
    path <- file.path(fig_dir, "penalized_observed_loglik.png")
    .draw_trace_plot(path, values / 1000,
                     "Penalized observed log-likelihood (thousands)",
                     observed_iterations)
    written <- c(written, path)
  }
  invisible(written)
}
