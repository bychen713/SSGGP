# Joint SSGGP estimation, one (replicate, penalty) pair per array task.
#   sbatch slurm/sim/main/sim_main.sh        (one task per row of grid_params)
# Reads the dataset and, with use_independent_init, the independent fit of the
# same replicate. Writes est<row>.RData with `result_est` and `time_approx`.

rm(list = ls())
if (!file.exists("settings.R")) {  # work from sim/, where settings.R lives
  if (file.exists(file.path("sim", "settings.R"))) setwd("sim") else
  if (file.exists(file.path("..", "settings.R"))) setwd("..") else
    stop("run this script from paper/")
}
source("settings.R")
setup_compute()

array_id <- job_index(1L)
if (array_id < 1L || array_id > nrow(grid_params)) {
  stop("THISJOBVALUE selects array row ", array_id, ", but grid_params has ",
       nrow(grid_params), " rows", call. = FALSE)
}
pars <- grid_params[array_id, , drop = FALSE]
seed <- pars$sim[[1]]

cat(sprintf("setting %s | array %d | sim %d | rho %g | L %d | MC %d | lag %d\n",
            setting_name, array_id, seed, pars$rho[[1]], pars$L[[1]],
            pars$MC[[1]], pars$lag[[1]]))

data_file <- file.path(data_dir, paste0(seed, ".RData"))
if (!file.exists(data_file)) {
  stop("dataset for simulation ", seed, " does not exist: ", data_file,
       ". Run sim/data_generation.R first.", call. = FALSE)
}
load(data_file)  # Y_mat, grid, Y_mat_test_obs, Y_mat_test, grid_test

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

independent_init <- NULL
fit_f_par <- f_par_init
fit_nugget_init <- nugget_init
if (isTRUE(use_independent_init)) {
  ind_file <- file.path(out_dir, paste0("ind_", seed, ".RData"))
  if (!file.exists(ind_file)) {
    stop("independent initialization is enabled but the required file is ",
         "missing: ", ind_file, ". Run sim/main/sim_ind.R (or slurm/sim/main/sim_ind.sh) before ",
         "the joint estimation array.", call. = FALSE)
  }
  load(ind_file)  # result_ind, time_true
  if (!isTRUE(all.equal(result_ind$Y, Y_mat, tolerance = 0,
                        check.attributes = FALSE)) ||
      !isTRUE(all.equal(result_ind$grid, grid, tolerance = 0,
                        check.attributes = FALSE))) {
    stop("the independent-initialization file does not match the current ",
         "dataset: ", ind_file, ". Rerun sim_ind.sh before sim_main.sh.",
         call. = FALSE)
  }
  fit_f_par <- result_ind$f_par
  fit_nugget_init <- mean(result_ind$nugget)
  independent_init <- list(
    f_par = result_ind$f_par,
    variance = result_ind$variance,
    nugget = result_ind$nugget,
    sigma_inv = result_ind$sigma_inv,
    logLik = result_ind$logLik,
    elapsed_minutes = as.numeric(time_true))
  # Families have different parameters, so summarize the starts by name.
  init_nms <- unique(unlist(lapply(fit_f_par, names)))
  init_summ <- vapply(init_nms, function(nm) {
    v <- unlist(lapply(fit_f_par, function(z)
      if (nm %in% names(z)) unname(z[[nm]]) else NULL))
    sprintf("mean %s %.4f", nm, mean(v))
  }, character(1))
  cat(sprintf("independent initialization: %s | mean nugget %.4f\n",
              paste(init_summ, collapse = " | "), fit_nugget_init))
}

a_time <- Sys.time()
checkpoint_file <- if (identical(Sys.getenv("SSGGP_CHECKPOINT", "0"), "1"))
  file.path(out_dir, "checkpoints", paste0("est", array_id, ".rds")) else NULL
f <- fit(
  Y     = Y_mat,
  locs  = grid,
  f_fun = f_fun_init,
  f_par = fit_f_par,
  g     = proposal,
  rho   = pars$rho[[1]],
  nugget      = nugget_type,
  nugget_init = fit_nugget_init,
  L = pars$L[[1]], MC = pars$MC[[1]], iter = iter, lag = pars$lag[[1]],
  cm_cycles = cm_cycles,
  tol = tol, diag_pen = diag_pen,
  lower = NULL, upper = NULL,
  threads = threads, seed = seed, verbose = TRUE,
  # light-tailed (e.g. Gaussian) targets need the importance-ratio cap
  numerical_guard = if (family_mix) "cap" else "error",
  checkpoint_file = checkpoint_file, checkpoint_every = 10L)
b_time <- Sys.time()
time_approx <- as.numeric(difftime(b_time, a_time, units = "secs"))

# Flattened a and nu estimates for the result tables (NA where a family has no
# such parameter); the full estimates are in result_est$fit$f_par.
extract_par <- function(x, nm) {
  if (nm %in% names(x)) unname(x[[nm]]) else NA_real_
}
a_est  <- vapply(f$f_par, extract_par, numeric(1), nm = "a")
nu_est <- vapply(f$f_par, extract_par, numeric(1), nm = "nu")

result_est <- list(
  Y          = as.vector(t(Y_mat)),
  Y_mat      = Y_mat,
  grid       = grid,
  Y_test     = as.vector(t(Y_mat_test)),
  Y_mat_test = Y_mat_test,
  grid_test  = grid_test,
  a_est      = a_est,
  nu_est     = nu_est,
  sigma_inv_est      = f$sigma_inv,
  sigma_inv_est_vals = f$chain$sigma_inv,
  D_est_vals = f$nugget,
  rho        = pars$rho[[1]],
  ll_prec    = f$chain$loglik[, "prec"],
  ll_nugget  = f$chain$loglik[, "nugget"],
  iterations = f$iterations,
  converged  = f$converged,
  fit        = f
)

save(result_est, time_approx, independent_init,
     file = file.path(out_dir, paste0("est", array_id, ".RData")))
cat(sprintf("array %d finished in %.1f seconds after %d iterations (converged: %s)\n",
            array_id, time_approx, f$iterations, f$converged))
