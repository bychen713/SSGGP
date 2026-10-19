# Latent IOX baseline for one MERRA-2 fold: spiox autostart() starts and the
# single-site sampler (sampling = 2). Run IOX_STAGE=estimate, then
# IOX_STAGE=predict (slurm/MERRA/iox/run.sh does both). Writes init_, fit_ and
# pred_fold<k>.rds under MERRA/output/<layout>/iox/, read by MERRA/result/.
# Run from paper/ (the launchers do).
root <- normalizePath(Sys.getenv("SLURM_SUBMIT_DIR", getwd()), mustWork = TRUE)
if (!file.exists(file.path(root, "iox", "latent_autostart.R")))
  stop("run from paper/", call. = FALSE)
layout <- Sys.getenv("IOX_LAYOUT", "colocated")
stopifnot(layout %in% c("colocated", "misaligned"))
setting <- paste0("merra_", layout)
source(file.path(root, "iox/latent_autostart.R"))
stage <- Sys.getenv("IOX_STAGE", "estimate")
fold <- as.integer(Sys.getenv("IOX_REP", "1"))
mcmc <- as.integer(Sys.getenv("IOX_MCMC", "10000"))
m_nn <- as.integer(Sys.getenv("IOX_NEIGHBORS", "15"))
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
n_keep <- as.integer(Sys.getenv("IOX_PRED_DRAWS", as.character(mcmc %/% 4L)))
stopifnot(stage %in% c("estimate", "predict", "validate"),
          !anyNA(c(fold, mcmc, m_nn, threads, n_keep)), fold %in% 1:3,
          mcmc >= 8L, m_nn >= 1L, threads >= 1L, n_keep >= 2L, n_keep <= mcmc)
if (!requireNamespace("spiox", quietly = TRUE) ||
    utils::packageVersion("spiox") < "0.8.1") stop("Local spiox 0.8.1 is required")
package_path <- normalizePath(find.package("spiox"))
package_version <- as.character(utils::packageVersion("spiox"))
package_dll <- system.file("libs", paste0("spiox", .Platform$dynlib.ext), package = "spiox")
stopifnot(nzchar(package_dll), file.exists(package_dll))
package_md5 <- unname(tools::md5sum(package_dll))
data_file <- normalizePath(file.path(root, "MERRA/data", layout,
                                     sprintf("fold%d.rds", fold)), mustWork = TRUE)
data_md5 <- unname(tools::md5sum(data_file))
data <- readRDS(data_file)
Y_in <- as.matrix(data$Y_mat)
Y_out <- as.matrix(data$Y_mat_test)
cx_in <- as.matrix(data$grid)
cx_out <- as.matrix(data$grid_test)
q <- ncol(Y_in)
stopifnot(is.numeric(Y_in), all(is.finite(Y_in[!is.na(Y_in)])),
          all(rowSums(!is.na(Y_in)) > 0L), all(colSums(!is.na(Y_in)) >= 3L),
          identical(anyNA(Y_in), layout == "misaligned"), all(is.finite(Y_out)),
          identical(ncol(Y_out), q),
          identical(dim(cx_in), c(nrow(Y_in), 2L)),
          identical(dim(cx_out), c(nrow(Y_out), 2L)),
          all(is.finite(cx_in)), all(is.finite(cx_out)), m_nn < nrow(Y_in),
          identical(dim(data$masked_train), dim(Y_in)),
          all(is.na(Y_in) == data$masked_train),
          length(data$train_index) == nrow(Y_in), length(data$test_index) == nrow(Y_out),
          !any(data$train_index %in% data$test_index))
X_in <- matrix(0, nrow(Y_in), 1L)
X_out <- matrix(0, nrow(Y_out), 1L)
out_dir <- file.path(root, "MERRA/output", layout, "iox")
init_file <- file.path(out_dir, paste0("init_fold", fold, ".rds"))
fit_file <- file.path(out_dir, paste0("fit_fold", fold, ".rds"))
pred_file <- file.path(out_dir, paste0("pred_fold", fold, ".rds"))
metadata <- list(data_file = data_file, data_md5 = data_md5, setting = setting,
                 dataset = "MERRA-2", preprocessing = data$preprocessing,
                 process_names = data$process_names,
                 train_index = data$train_index, test_index = data$test_index,
                 observed_per_process = colSums(!is.na(Y_in)),
                 fold = fold, model = "latent", implementation = "github", layout = layout,
                 initialization = "autostart_zero_mean", initial_nu = .5,
                 cg_preconditioner = "not_used",
                 update_Theta = c(1L, 1L, 0L),
                 runner_md5 = unname(tools::md5sum(file.path(root, "MERRA/iox/run.R"))),
                 initializer_md5 = unname(tools::md5sum(file.path(root, "iox/latent_autostart.R"))),
                 mcmc = mcmc, m_nn = m_nn, package_version = package_version,
                 package_path = package_path, package_md5 = package_md5,
                 threads = threads, sampling = 2L, fit_seed = 2024L + fold,
                 job_id = Sys.getenv("SLURM_JOB_ID"), node = Sys.getenv("SLURMD_NODENAME"))

if (stage == "validate") {
  stopifnot(!file.exists(fit_file), !file.exists(pred_file), !file.exists(init_file))
  cat("Validated", setting, "replicate", fold, "data MD5", data_md5, "\n")
} else if (stage == "estimate") {
  if (file.exists(fit_file) || file.exists(pred_file) || file.exists(init_file))
    stop("Refusing to overwrite an existing fit or prediction in ", out_dir)
  cat("Initializing with package autostart (nu=.5); IOX mean remains zero.\n")
  RhpcBLASctl::blas_set_num_threads(1L)
  initialization_elapsed <- system.time({
    init <- iox_latent_autostart(Y_in, X_in, cx_in, m = m_nn,
                                seed = metadata$fit_seed)
  })
  starting <- init$starting
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(c(list(starting = starting, details = init$details,
                 initialization_elapsed = initialization_elapsed), metadata),
          init_file, compress = FALSE)
  cat("Saved initialization:", init_file, "\n")
  set.seed(metadata$fit_seed)
  elapsed <- system.time({
    fit <- spiox::spiox(
      Y = Y_in, X = X_in, coords = cx_in, m = m_nn,
      method = "latent", fit = "mcmc", iter = mcmc, print_every = 100L,
      starting = starting,
      opts = list(matern = 1L, num_threads = threads, update_Theta = c(1L, 1L, 0L)),
      debug = list(sample_Beta = FALSE, sampling = 2L))
  })
  stopifnot(identical(dim(fit$W), c(nrow(Y_in), q, mcmc)),
            identical(dim(fit$Ddiag), c(q, mcmc)),
            identical(dim(fit$Theta), c(3L, q, mcmc)),
            all(is.finite(fit$W)), all(is.finite(fit$Sigma)),
            all(is.finite(fit$Ddiag)), all(fit$Ddiag > 0),
            all(is.finite(fit$Theta)), all(fit$Theta[3L, , ] == 0))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(c(list(fit = fit, elapsed = elapsed, initialization_elapsed = initialization_elapsed,
                 initialization_file = init_file), metadata), fit_file, compress = FALSE)
  cat("Saved fit:", fit_file, "\n")
} else {
  if (file.exists(pred_file)) stop("Refusing to overwrite prediction: ", pred_file)
  if (!file.exists(fit_file)) stop("Missing fit: ", fit_file)
  saved <- readRDS(fit_file)
  for (key in c("data_file", "data_md5", "setting", "fold", "model", "implementation",
                "mcmc", "m_nn", "package_version", "package_path", "package_md5", "sampling",
                "layout", "initialization", "initial_nu", "update_Theta", "cg_preconditioner",
                "runner_md5", "initializer_md5")) {
    if (!identical(saved[[key]], metadata[[key]])) stop("Fit metadata mismatch: ", key)
  }
  fit <- saved$fit
  keep <- tail(seq_len(mcmc), n_keep)
  stopifnot(identical(dim(fit$W), c(nrow(Y_in), q, mcmc)),
            identical(dim(fit$Ddiag), c(q, mcmc)),
            identical(dim(fit$Theta), c(3L, q, mcmc)),
            identical(dim(fit$Beta), c(1L, q, mcmc)),
            identical(dim(fit$Sigma), c(q, q, mcmc)))
  # Match W, Sigma, noise, and spatial parameters at the same retained iteration.
  # The predictor expects four theta rows, whereas spiox() returns three.
  theta <- array(0, c(4L, q, n_keep))
  theta[1L, , ] <- fit$Theta[1L, , keep]
  theta[2L, , ] <- 1
  theta[3L, , ] <- fit$Theta[2L, , keep]
  theta[4L, , ] <- fit$Theta[3L, , keep]
  dag <- spiox::dag_vecchia_predict(cx_in, cx_out, m_nn)
  set.seed(52024L + fold)
  elapsed <- system.time({
    draws <- spiox::spiox_latent_predict(
      X_new = X_out, coords_new = cx_out, coords = cx_in, dag = dag,
      W = fit$W[, , keep, drop = FALSE], B = fit$Beta[, , keep, drop = FALSE],
      Sigma = fit$Sigma[, , keep, drop = FALSE], Dvec = fit$Ddiag[, keep, drop = FALSE],
      theta = theta, matern = 1L, num_threads = threads)$Y
  })
  stopifnot(identical(dim(draws), c(nrow(Y_out), q, n_keep)), all(is.finite(draws)))
  pred_mean <- apply(draws, c(1L, 2L), mean)
  pred_sd <- apply(draws, c(1L, 2L), stats::sd)
  sum_draws <- apply(draws, c(1L, 3L), sum)
  sum_mean <- rowMeans(sum_draws)
  sum_sd <- apply(sum_draws, 1L, stats::sd)
  metrics <- list(
    norm = mean(sqrt(rowSums((Y_out - pred_mean)^2))),
    coverage = mean(abs(Y_out - pred_mean) <= 1.96 * pred_sd),
    coverage_per_process = colMeans(abs(Y_out - pred_mean) <= 1.96 * pred_sd),
    rmse_sum = sqrt(mean((rowSums(Y_out) - sum_mean)^2)),
    coverage_sum = mean(abs(rowSums(Y_out) - sum_mean) <= 1.96 * sum_sd),
    # Arithmetic mean pointwise width, matching draft.tex. RMS widths are
    # retained under explicit names for comparison with older stored metrics.
    pi_length = 3.92 * mean(pred_sd), pi_length_per_process = 3.92 * colMeans(pred_sd),
    pi_length_sum = 3.92 * mean(sum_sd), pi_length_rms = 3.92 * sqrt(mean(pred_sd^2)),
    pi_length_sum_rms = 3.92 * sqrt(mean(sum_sd^2)),
    initialization_minutes = unname(saved$initialization_elapsed["elapsed"]) / 60,
    estimate_minutes = unname(saved$elapsed["elapsed"]) / 60,
    total_estimate_minutes = (unname(saved$initialization_elapsed["elapsed"]) + unname(saved$elapsed["elapsed"])) / 60,
    predict_minutes = unname(elapsed["elapsed"]) / 60)
  saveRDS(c(list(draws = draws, mean = pred_mean, sd = pred_sd, metrics = metrics,
                 elapsed = elapsed, fit_file = fit_file, keep = keep,
                 pi_length_definition = "3.92 times arithmetic mean of pointwise posterior SD"),
            metadata), pred_file, compress = FALSE)
  cat("Saved predictions:", pred_file, "\n")
  print(metrics)
}
