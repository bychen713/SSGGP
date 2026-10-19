# Independent-process baseline ("Ind"): each process is fitted separately, so the
# precision matrix is diagonal. Its estimates also initialize the joint fit.
#   sbatch slurm/sim/main/sim_ind.sh         (one task per replicate)
# Writes ind_<sim>.RData with `result_ind` and `time_true`.

rm(list = ls())
if (!file.exists("settings.R")) {  # work from sim/, where settings.R lives
  if (file.exists(file.path("sim", "settings.R"))) setwd("sim") else
  if (file.exists(file.path("..", "settings.R"))) setwd("..") else
    stop("run this script from paper/")
}
source("settings.R")
setup_compute()

j <- job_index(1L)
stopifnot(j >= 1, j <= n_sim)

load(file.path(data_dir, paste0(j, ".RData")))   # Y_mat, grid, Y_mat_test, grid_test
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("setting %s | replicate %d | independent baseline\n", setting_name, j))

time_a <- Sys.time()
ind <- fit_independent(
  Y = Y_mat, locs = grid,
  f_fun = f_fun_init,
  f_par = f_par_init,
  variance_init = NULL,               # start at each process's sample variance
  nugget_init = nugget_init,
  fixed = NULL,
  lower = NULL, upper = NULL,
  threads = threads)
time_b <- Sys.time()
time_true <- difftime(time_b, time_a, units = "mins")

print(ind)

# Univariate kriging of each process with its fitted covariance.
n_test <- nrow(grid_test)
pred_mean <- matrix(NA_real_, n_test, P)
pred_var  <- matrix(NA_real_, n_test, P)
entries <- SSGGP:::resolve_spectral(f_fun_init, P)
for (k in seq_len(P)) {
  keep <- !is.na(Y_mat[, k])
  o    <- grid[keep, , drop = FALSE]
  z    <- Y_mat[keep, k]
  e    <- entries[[k]]
  theta_k <- e$to_theta(ind$f_par[[k]], ind$variance[k], ind$nugget[k])

  nt   <- nrow(grid_test); no <- nrow(o)
  Call <- SSGGP:::family_cov(e, fields::rdist(rbind(grid_test, o)), theta_k,
                             nugget = FALSE)
  Cno <- Call[seq_len(nt), nt + seq_len(no), drop = FALSE]
  Coo <- Call[nt + seq_len(no), nt + seq_len(no), drop = FALSE]
  diag(Coo) <- diag(Coo) + ind$nugget[k]
  ch  <- chol(Coo)
  pred_mean[, k] <- as.vector(Cno %*% backsolve(ch, backsolve(ch, z, transpose = TRUE)))
  # predictive variance includes the nugget
  v <- backsolve(ch, t(Cno), transpose = TRUE)
  pred_var[, k]  <- ind$variance[k] + ind$nugget[k] - colSums(v^2)
}

# Per-process predictions and parameters in the layout read by result/helpers.R.
slot <- function(p, nm) if (nm %in% names(p)) unname(p[[nm]]) else NA_real_
opt_true <- lapply(seq_len(P), function(k)
  list(par = c(1 / ind$variance[k], slot(ind$f_par[[k]], "a"),
               slot(ind$f_par[[k]], "nu"), ind$nugget[k])))
names(opt_true) <- paste0("process", seq_len(P))

result_ind <- list(
  Y          = Y_mat,
  grid       = grid,
  Y_test     = Y_mat_test,
  grid_test  = grid_test,
  prediction = split(pred_mean, col(pred_mean)),
  pred_var   = split(pred_var,  col(pred_var)),
  opt_true   = opt_true,
  f_par      = ind$f_par,
  variance   = ind$variance,
  nugget     = ind$nugget,
  sigma_inv  = ind$sigma_inv,
  logLik     = ind$logLik,
  fit        = ind)

save(time_true, result_ind,
     file = file.path(out_dir, paste0("ind_", j, ".RData")))
cat("wrote ind_", j, ".RData\n", sep = "")
