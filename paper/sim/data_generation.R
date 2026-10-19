# Generate the simulated datasets of the active setting.
#   Rscript sim/data_generation.R            (or: sbatch slurm/sim/data_generation.sh)
# Writes sim/data/<setting>/<sim>.RData with Y_mat and grid (training
# data, NA where unobserved), Y_mat_test and grid_test (held-out data), and
# Y_mat_test_obs and new_obs_ind.

rm(list = ls())
if (!file.exists("settings.R")) {
  if (file.exists(file.path("sim", "settings.R"))) setwd("sim") else
    stop("run this script from paper/ or paper/sim/")
}
source("settings.R")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

Sigma_inv <- true_precision()
cat("setting:", setting_name, "\n")
cat("true precision:\n"); print(round(Sigma_inv, 3))
cat("edges:", sum(Sigma_inv[upper.tri(Sigma_inv)] != 0), "of", choose(P, 2), "\n\n")

for (array_id in seq_len(n_sim)) {
  set.seed(array_id)

  grid_all <- matrix(stats::runif(2 * dim_grid^2, min = 0, max = domain_max),
                     ncol = 2)
  n_all    <- nrow(grid_all)

  out <- sim(grid_all, sigma_inv = Sigma_inv,
             f_fun = f_fun_true, f_par = f_par_true,
             nugget = nugget_true,
             method = sim_method, L = sim_L, g = proposal_sim,
             # mixed families need the same importance-ratio cap as the fit
             numerical_guard = if (family_mix) "cap" else "error",
             seed = NULL)
  Y_mat_all <- out$Y

  # Locations are iid uniform, so the first train_prop of them form a random split.
  idx_train <- seq_len(floor(n_all * train_prop))
  idx_test  <- setdiff(seq_len(n_all), idx_train)

  grid       <- grid_all[idx_train, , drop = FALSE]
  Y_mat      <- Y_mat_all[idx_train, , drop = FALSE]
  grid_test  <- grid_all[idx_test,  , drop = FALSE]
  Y_mat_test <- Y_mat_all[idx_test, , drop = FALSE]

  # Misaligned design: each process loses its own random subset of training locations.
  miss_obs_prop     <- missing_prop$miss_obs_prop
  miss_new_obs_prop <- missing_prop$miss_new_obs_prop
  stopifnot(length(miss_obs_prop) == P, length(miss_new_obs_prop) == P)

  n <- nrow(Y_mat)
  for (j in seq_len(P)) {
    if (miss_obs_prop[j] != 0)
      Y_mat[sample(seq_len(n), miss_obs_prop[j] * n), j] <- NA
  }

  # Y_mat_test_obs holds test values revealed to the predictor and Y_mat_test the
  # values held out for scoring. With miss_new_obs_prop = 1 nothing is revealed.
  n_test         <- nrow(Y_mat_test)
  Y_mat_test_obs <- Y_mat_test
  new_obs_ind    <- sample(seq_len(n_test),
                           (1 - utils::tail(miss_new_obs_prop, 1)) * n_test)
  for (j in seq_len(P)) {
    if (miss_new_obs_prop[j] == 0) {
      Y_mat_test[, j] <- NA
    } else if (miss_new_obs_prop[j] == 1) {
      Y_mat_test_obs[, j] <- NA
    } else {
      Y_mat_test_obs[-new_obs_ind, j] <- NA
      Y_mat_test[new_obs_ind, j]      <- NA
    }
  }

  save(Y_mat, grid, Y_mat_test_obs, new_obs_ind, Y_mat_test, grid_test,
       file = file.path(data_dir, paste0(array_id, ".RData")))
  cat(sprintf("  sim %2d: %d train x %d processes (%.0f%% observed), %d test\n",
              array_id, nrow(Y_mat), P, 100 * mean(!is.na(Y_mat)), nrow(grid_test)))
}

cat("\nwrote", n_sim, "datasets to", data_dir, "\n")
