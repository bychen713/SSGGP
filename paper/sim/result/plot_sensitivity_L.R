# Figure S2: absolute log-likelihood error |l_L - l_ref| / (n p) against the
# number of features L, for p = 9, 25, 50 (mixed families), plotted from saved
# results of the L-sensitivity experiment (no fitting).
#   Rscript sim/result/plot_sensitivity_L.R [RESULTS_DIR]
# RESULTS_DIR (default sim/output/sensitivity_L) holds mixed_p9/, mixed_p25/ and
# mixed_p50/: replicate_01.rds ... for p = 9, and integrated_replicate_01.rds ...
# plus simulated_replicate_01.rds ... for p = 25 and 50, each with a `result`
# data frame (columns n, p, L, method, abs_error).
# Writes sensitivity_L_mixed_p<P>.png to sim/output/figure/.

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args)) args[1L] else file.path("sim", "output", "sensitivity_L")
fig_dir <- file.path("sim", "output", "figure")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

n_values <- c(500L, 1000L)
L_values <- c(500L, 1000L, 2500L, 5000L, 10000L, 20000L, 30000L, 40000L, 50000L)
colours <- c(`500` = "#0072B2", `1000` = "#D55E00")
line_types <- c(integrated_W = 1, simulated_W = 2)
point_shapes <- c(integrated_W = 21, simulated_W = 24)
method_labels <- c(integrated_W = "W integrated", simulated_W = "W simulated")

read_results <- function(P) {
  dir <- file.path(results_dir, sprintf("mixed_p%d", P))
  files <- if (P == 9L) sprintf("replicate_%02d.rds", 1:15) else
    c(sprintf("integrated_replicate_%02d.rds", 1:15),
      sprintf("simulated_replicate_%02d.rds", 1:15))
  files <- file.path(dir, files)
  if (any(!file.exists(files)))
    stop("missing results: ", paste(files[!file.exists(files)], collapse = ", "))
  do.call(rbind, lapply(files, function(f) readRDS(f)$result))
}

# Mean and SD over replicates of the error per n p, by n, L and method.
summarise <- function(P) {
  dat <- read_results(P)
  dat$error <- dat$abs_error / (dat$n * dat$p)
  out <- aggregate(error ~ n + L + method, dat, function(x) c(mean = mean(x), sd = sd(x)))
  data.frame(out[c("n", "L", "method")], mean = out$error[, "mean"], sd = out$error[, "sd"])
}

for (P in c(9L, 25L, 50L)) {
  s <- summarise(P)
  path <- file.path(fig_dir, sprintf("sensitivity_L_mixed_p%d.png", P))
  png(path, width = 1800, height = 1200, res = 180)
  par(mar = c(5.5, 7.2, 4.0, 1.5))
  plot(range(L_values) * c(1 / 1.10, 1.15), c(0, max(s$mean + s$sd) * 1.06),
       type = "n", log = "x", xaxt = "n", xaxs = "i", yaxs = "i", las = 1,
       xlab = "", ylab = expression(abs(ell[L] - ell[ref]) / (n * p)),
       main = sprintf("p = %d", P), cex.axis = 1.452, cex.lab = 1.672, cex.main = 1.958)
  labels <- format(L_values, big.mark = ",", scientific = FALSE, trim = TRUE)
  labels[L_values %in% c(30000L, 40000L)] <- ""
  axis(1, at = L_values, labels = labels, cex.axis = 1.364, gap.axis = 0.15)
  grid(nx = NA, ny = NULL, col = "grey90", lty = 1)
  for (n in n_values) for (method in names(line_types)) {
    z <- s[s$n == n & s$method == method, ]
    z <- z[order(z$L), ]
    col <- colours[as.character(n)]; lty <- line_types[method]
    lo <- pmax(z$mean - z$sd, 0); hi <- z$mean + z$sd
    segments(z$L, lo, z$L, hi, col = col, lty = lty, lwd = 1.4)        # +/- 1 SD
    segments(z$L / 1.018, lo, z$L * 1.018, lo, col = col, lty = lty)
    segments(z$L / 1.018, hi, z$L * 1.018, hi, col = col, lty = lty)
    lines(z$L, z$mean, col = col, lty = lty, lwd = 2)
    points(z$L, z$mean, col = col, bg = "white", pch = point_shapes[method],
           lwd = 1.5, cex = 1.05)
  }
  legend("topright", bty = "n", cex = 1.474, pt.bg = "white", lwd = 2,
         y.intersp = 1.12, text.font = c(2, 1, 1, 1, 2, 1, 1),
         legend = c("Sample size", paste0("n = ", n_values), "", "Likelihood",
                    unname(method_labels)),
         col = c(NA, colours, NA, NA, "black", "black"),
         lty = c(NA, 1, 1, NA, NA, unname(line_types)),
         pch = c(NA, NA, NA, NA, NA, unname(point_shapes)))
  mtext("L", side = 1, line = 3.5, cex = 1.628)
  invisible(dev.off())
  cat("Saved", path, "\n")
}
