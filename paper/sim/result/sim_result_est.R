# Parameter-estimation tables by covariance family, for mixed-family designs.
#   Rscript sim/result/sim_result_est.R       (also sourced by sim_result.R)
# Reads fit$f_par from the joint and independent fits and the BIC files for the
# selected penalty. Tables go to sim/output/tables/<setting>/.
# Family errors average process-specific biases and SDs across replicates,
# including the shared or process-specific nugget.

if (sys.nframe() == 0L) rm(list = ls())   # only when run as a script
if (!file.exists("settings.R")) {  # work from sim/, where settings.R lives
  if (file.exists(file.path("sim", "settings.R"))) setwd("sim") else
  if (file.exists(file.path("..", "settings.R"))) setwd("..") else
    stop("run this script from paper/")
}
source("settings.R")

# Settings to summarize (the active one by default).
scan_settings <- data.frame(prec_structure = prec_structure,
                            offdiag_sign   = offdiag_sign,
                            missing_TF     = missing_TF,
                            stringsAsFactors = FALSE)

# TRUE restricts the summaries to converged fits.
use_converged_only <- FALSE

# Collect

# True parameters, one row per (process, parameter).
true_long <- do.call(rbind, lapply(seq_len(P), function(j) {
  fam <- if (is.character(f_fun_true[[j]])) f_fun_true[[j]] else NA_character_
  par <- c(f_par_true[[j]], nugget = rep_len(nugget_true, P)[j])
  # par_order keeps each family's own parameter order in every table.
  data.frame(process = j, family_true = fam, parameter = names(par),
             par_order = seq_along(par),
             truth = as.numeric(par), stringsAsFactors = FALSE)
}))

cat("true parameters:\n")
print(true_long[, setdiff(names(true_long), "par_order")], row.names = FALSE)
cat("\n")

collect_one <- function(prec_structure, offdiag_sign, missing_TF) {
  Sys.setenv(SSGGP_PREC = prec_structure,
             SSGGP_SIGN = offdiag_sign,
             SSGGP_MISS = missing_TF)
  source("settings.R", local = FALSE)

  if (!dir.exists(out_dir)) {
    warning("no output directory for ", setting_name, ": ", out_dir,
            call. = FALSE)
    return(NULL)
  }

  # One row per (process, parameter) of a fit.
  unpack <- function(f_par, fam_fit, nugget) {
    if (!length(nugget) %in% c(1L, P))
      stop("Expected one shared nugget or one nugget per process")
    nugget <- rep_len(nugget, P)
    do.call(rbind, lapply(seq_len(P), function(j) {
      est <- f_par[[j]]
      if (is.null(est) || is.null(names(est))) return(NULL)
      est <- c(est, nugget = as.numeric(nugget[j]))
      data.frame(process    = j,
                 family_fit = fam_fit[[j]],
                 parameter  = names(est),
                 estimate   = as.numeric(est),
                 stringsAsFactors = FALSE)
    }))
  }

  rows <- seq_len(nrow(grid_params))
  out <- list()
  n_missing <- 0L

  # Joint fits, one per (replicate, penalty)
  for (i in rows) {
    f <- file.path(out_dir, paste0("est", i, ".RData"))
    if (!file.exists(f)) { n_missing <- n_missing + 1L; next }

    result_est <- NULL
    load(f)                                   # result_est, time_approx
    fit_obj <- result_est$fit
    if (is.null(fit_obj) || is.null(fit_obj$f_par)) {
      warning("est", i, ".RData in ", setting_name,
              " has no fit$f_par (pre-package file?); skipped", call. = FALSE)
      next
    }

    fam_fit <- fit_obj$spec$family
    if (is.null(fam_fit)) fam_fit <- rep(NA_character_, P)

    d <- unpack(fit_obj$f_par, fam_fit, fit_obj$nugget)
    if (is.null(d)) next
    d$setting   <- setting_name
    d$method    <- "saem"
    d$sim       <- grid_params$sim[i]
    d$rho       <- grid_params$rho[i]
    d$converged <- isTRUE(result_est$converged)
    out[[length(out) + 1L]] <- d
  }

  # Independent-process fits, one per replicate.
  n_ind <- 0L
  for (j in seq_len(n_sim)) {
    f <- file.path(out_dir, paste0("ind_", j, ".RData"))
    if (!file.exists(f)) next
    result_ind <- NULL
    load(f)                                   # result_ind, time_true
    if (is.null(result_ind$f_par)) {
      warning("ind_", j, ".RData in ", setting_name, " has no f_par; skipped",
              call. = FALSE)
      next
    }
    fam_fit <- result_ind$fit$spec$family
    if (is.null(fam_fit)) fam_fit <- rep("matern", P)

    nugget <- result_ind$nugget
    if (is.null(nugget)) nugget <- result_ind$fit$nugget
    d <- unpack(result_ind$f_par, fam_fit, nugget)
    if (is.null(d)) next
    d$setting   <- setting_name
    d$method    <- "ind"
    d$sim       <- j
    d$rho       <- NA_real_
    d$converged <- NA   # fit_independent reports no convergence flag
    out[[length(out) + 1L]] <- d
    n_ind <- n_ind + 1L
  }

  if (n_missing) {
    message("  ", setting_name, ": ", n_missing, " of ", length(rows),
            " est files absent")
  }
  message("  ", setting_name, ": ", n_ind, " of ", n_sim,
          " ind_*.RData found")
  do.call(rbind, out)
}

cat("collecting:\n")
long <- do.call(rbind, lapply(seq_len(nrow(scan_settings)), function(k) {
  s <- scan_settings[k, ]
  cat("  ", paste(s$prec_structure, s$offdiag_sign, s$missing_TF, sep = "_"), "\n")
  collect_one(s$prec_structure, s$offdiag_sign, s$missing_TF)
}))

if (is.null(long) || !nrow(long)) stop("no est<row>.RData files found to score")

# Attach the true values; estimates without a matching truth keep truth = NA.
long <- merge(long, true_long, by = c("process", "parameter"), all.x = TRUE)

# Only estimates from the generating family are compared with the truth.
long$comparable <- !is.na(long$family_fit) & !is.na(long$family_true) &
  long$family_fit == long$family_true

mismatch <- unique(long[!long$comparable & !is.na(long$family_fit),
                        c("method", "process", "family_true", "family_fit")])
if (nrow(mismatch)) {
  message("\nfitted family differs from the generating family, so `truth` is ",
          "not the same quantity as the estimate for these rows:")
  message(paste(sprintf("  %-5s process %d: generated '%s', fitted '%s'",
                        mismatch$method, mismatch$process,
                        mismatch$family_true, mismatch$family_fit),
                collapse = "\n"))
}

if (use_converged_only) {
  n_drop <- sum(!long$converged)
  long <- long[long$converged, , drop = FALSE]
  message("dropped ", n_drop, " rows from fits that did not converge")
}

# Summarize

# Mean and SD of the estimates across replicates, for one grouping.
summarise_est <- function(d, by) {
  key <- interaction(d[by], drop = TRUE, sep = "\r")
  parts <- split(seq_len(nrow(d)), key)
  out <- do.call(rbind, lapply(parts, function(ix) {
    v  <- d$estimate[ix]
    ok <- is.finite(v)
    cbind(d[ix[1], by, drop = FALSE],
          data.frame(par_order = d$par_order[ix[1]],
                     truth    = d$truth[ix[1]],
                     n        = sum(ok),
                     n_conv   = sum(d$converged[ix] & ok),
                     mean_est = if (any(ok)) mean(v[ok]) else NA_real_,
                     sd       = if (sum(ok) > 1L) stats::sd(v[ok]) else NA_real_,
                     comparable = d$comparable[ix[1]],
                     stringsAsFactors = FALSE))
  }))
  rownames(out) <- NULL
  keys <- intersect(c("setting", "process", "par_order", "method", "rho"),
                    names(out))
  out <- out[do.call(order, unname(out[keys])), ]
  out$par_order <- NULL
  out
}

# Estimates by penalty (joint fits only).
saem <- long[long$method == "saem", , drop = FALSE]
by_rho <- summarise_est(saem, c("setting", "process", "family_fit",
                                "parameter", "rho"))

# Penalty of the dense fit (the first) and of the BIC-selected fit (from
# BIC_<setting>_<sim>.RData) for each replicate.
selected_rows <- function(setting_row) {
  Sys.setenv(SSGGP_PREC = setting_row$prec_structure,
             SSGGP_SIGN = setting_row$offdiag_sign,
             SSGGP_MISS = setting_row$missing_TF)
  source("settings.R", local = FALSE)

  dense <- data.frame(setting = setting_name, sim = seq_len(n_sim),
                      method = "dense", rho = rho_vals[1],
                      stringsAsFactors = FALSE)

  bic <- do.call(rbind, lapply(seq_len(n_sim), function(j) {
    f <- file.path(out_dir, paste0("BIC_", setting_name, "_", j, ".RData"))
    if (!file.exists(f)) return(NULL)
    result_final <- NULL
    load(f)
    if (is.null(result_final$rho)) return(NULL)
    data.frame(setting = setting_name, sim = j, method = "BIC",
               rho = as.numeric(result_final$rho), stringsAsFactors = FALSE)
  }))
  if (is.null(bic)) {
    message("  ", setting_name,
            ": no BIC_*.RData -- BIC rows omitted (run sim_BIC.R first)")
  }
  rbind(dense, bic)
}

sel <- do.call(rbind, lapply(seq_len(nrow(scan_settings)),
                             function(k) selected_rows(scan_settings[k, ])))

# Match the selected penalty to the nearest grid value rather than with ==.
sel$rho_key <- NA_real_
saem$rho_key <- saem$rho
for (i in seq_len(nrow(sel))) {
  cand <- unique(saem$rho[saem$setting == sel$setting[i]])
  if (!length(cand)) next
  sel$rho_key[i] <- cand[which.min(abs(cand - sel$rho[i]))]
}

# Dense, BIC, and independent estimates side by side.
long_sel <- merge(saem[, setdiff(names(saem), "method")],
                  sel[, c("setting", "sim", "method", "rho_key")],
                  by = c("setting", "sim", "rho_key"))
long_sel$rho_key <- NULL
ind_rows <- long[long$method == "ind", , drop = FALSE]
if (nrow(ind_rows)) long_sel <- rbind(long_sel, ind_rows[names(long_sel)])

by_method <- if (nrow(long_sel))
  summarise_est(long_sel, c("setting", "process", "family_fit", "parameter",
                            "method")) else NULL

# Estimates by process and method as "mean (SD)".
method_cols <- c("BIC", "dense", "ind")

widen <- function(d, cols = method_cols) {
  if (is.null(d) || !nrow(d)) return(NULL)
  key  <- paste(d$setting, d$process, d$parameter, sep = "\r")
  ord  <- unique(key)                   # d is already in family-declared order
  base <- d[match(ord, key), c("setting", "process", "family_fit", "parameter",
                               "truth")]
  names(base)[names(base) == "family_fit"] <- "family"
  names(base)[names(base) == "parameter"]  <- "par"

  cell <- function(m, s) {
    out <- ifelse(is.na(m), NA_character_,
                  ifelse(is.na(s), sprintf("%.3f", m),
                         sprintf("%.3f (%.3f)", m, s)))
    out
  }
  for (mth in cols) {
    ix <- match(ord, key[d$method == mth])
    sub <- d[d$method == mth, , drop = FALSE]
    base[[mth]] <- if (!nrow(sub)) NA_character_ else
      cell(sub$mean_est[ix], sub$sd[ix])
  }
  rownames(base) <- NULL
  base
}

by_method_wide <- widen(by_method)

# Calculate bias and SD across replicates separately for each process, then
# average those biases and SDs within each family and parameter. Repeating a
# shared nugget over processes leaves its across-replicate bias and SD intact.
err_by_family <- function(d, cols = method_cols) {
  if (is.null(d) || !nrow(d)) return(NULL)
  keep <- d$comparable & is.finite(d$estimate) & is.finite(d$truth)
  if (sum(!keep))
    message(sum(!keep), " of ", length(keep), " rows left out of the error ",
            "table: fitted family not the generating one, or no finite estimate")
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) return(NULL)

  parts <- split(seq_len(nrow(d)),
                 interaction(d$setting, d$family_fit, d$parameter,
                             drop = TRUE, sep = "\r"))
  out <- do.call(rbind, lapply(parts, function(ix) {
    x   <- d[ix, , drop = FALSE]
    row <- data.frame(setting   = x$setting[1],
                      family    = x$family_fit[1],
                      par       = x$parameter[1],
                      par_order = x$par_order[1],   # dropped after ordering
                      P         = length(unique(x$process)),
                      truth_rng = sprintf("%.2f-%.2f", min(x$truth),
                                          max(x$truth)),
                      stringsAsFactors = FALSE)
    for (mth in cols) {
      xm <- x[x$method == mth, , drop = FALSE]
      process_errors <- split(xm$estimate - xm$truth, xm$process)
      row[[paste0("bias_", mth)]] <- if (length(process_errors))
        mean(vapply(process_errors, mean, numeric(1))) else NA_real_
      row[[paste0("sd_", mth)]] <- if (length(process_errors))
        mean(vapply(process_errors, stats::sd, numeric(1))) else NA_real_
      row[[paste0("n_rep_", mth)]] <- length(unique(xm$sim))
    }
    row
  }))
  out <- out[order(out$setting, out$family, out$par_order), ]
  out$par_order <- NULL
  rownames(out) <- NULL
  out
}

err_wide <- err_by_family(long_sel)

# Report

# Print one block per setting.
report <- function(d, title, stat_cols = c("truth", "n", "n_conv", "mean_est",
                                           "sd", "comparable")) {
  cat("\n================ ", title, " ================\n", sep = "")
  for (s in unique(d$setting)) {
    x <- d[d$setting == s, , drop = FALSE]
    x$setting <- NULL
    names(x)[names(x) == "family_fit"] <- "family"
    x <- x[, c(setdiff(names(x), stat_cols), intersect(stat_cols, names(x)))]
    num <- vapply(x, is.numeric, logical(1))
    x[num] <- lapply(x[num], function(z) round(z, 4))
    cat("\n[", s, "]\n", sep = "")
    print(x, row.names = FALSE)
  }
}

old_width <- options(width = 200)

print_by_rho <- FALSE
if (print_by_rho)
  report(by_rho, "estimates by penalty (SAEM): mean and sd over replicates")

if (!is.null(by_method_wide)) {
  cat("\n================ estimates by method: mean (sd) over replicates ",
      "================\n", sep = "")
  for (s in unique(by_method_wide$setting)) {
    x <- by_method_wide[by_method_wide$setting == s, , drop = FALSE]
    x$setting <- NULL
    cat("\n[", s, "]\n", sep = "")
    print(x, row.names = FALSE, na.print = "-")
  }
}

if (!is.null(err_wide)) {
  cat("\n================ error by family: bias and SD over replicates per process, ",
      "then averaged over processes ================\n", sep = "")
  for (s in unique(err_wide$setting)) {
    x <- err_wide[err_wide$setting == s, , drop = FALSE]
    x$setting <- NULL
    num <- vapply(x, is.numeric, logical(1))
    x[num] <- lapply(x[num], function(z) round(z, 4))
    cat("\n[", s, "]\n", sep = "")
    print(x, row.names = FALSE, na.print = "-")
  }
}

options(old_width)

# Save the tables
for (setting in unique(by_method_wide$setting)) {
  table_dir <- here("sim", "output", "tables", setting)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(by_method_wide[by_method_wide$setting == setting, , drop = FALSE],
            file.path(table_dir, "estimation.csv"), row.names = FALSE)
  if (!is.null(err_wide))
    write.csv(err_wide[err_wide$setting == setting, , drop = FALSE],
              file.path(table_dir, "covariance_estimation_error.csv"), row.names = FALSE)
}
