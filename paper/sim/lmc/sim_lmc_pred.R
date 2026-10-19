# LMC stage 2: joint maximum-likelihood fit with the structure from stage 1, then
# prediction at the test locations.
#   sbatch slurm/sim/lmc/sim_lmc_pred.sh        (one task per replicate)
# Writes lmcscov_<sim>.RData, read by sim/result/sim_result.R.

rm(list = ls())

if (file.exists("sim/settings.R")) {
  source("sim/settings.R")
} else if (file.exists("../settings.R")) {
  source("../settings.R")
} else {
  stop("run from paper/ or paper/sim/lmc/", call. = FALSE)
}
source(here("sim", "lmc", "lmc_model.R"))
invisible(list2env(lmc_dataset_config(project_root, setting_name), envir = environment()))
setup_compute()

array_id <- job_index()
cat("setting:", setting_name, " replicate:", array_id, "\n")

# SSGGP_LMC_CONTINUE=true resumes an unfinished fit from its result/checkpoint.
resume <- tolower(Sys.getenv("SSGGP_LMC_CONTINUE", "false")) %in% c("true", "1", "yes")
data_file <- file.path(data_dir, paste0(array_id, ".RData"))
result_file <- file.path(out_dir, paste0("lmcscov_", array_id, ".RData"))
checkpoint_file <- file.path(out_dir, paste0("lmc_checkpoint_", array_id, ".rds"))
load(data_file)
est_file <- file.path(out_dir, paste0("lmc_structure_", array_id, ".RData"))
if (!file.exists(est_file)) {
  stop("stage 1 output not found: ", est_file,
       "\nRun sim/lmc/sim_lmc_structure.R first (submit_sim_lmc_all.sh chains the stages).",
       call. = FALSE)
}
load(est_file)
if (!exists("fam", inherits = FALSE) ||
    !exists("selection_method", inherits = FALSE) ||
    !exists("selection_aic", inherits = FALSE) ||
    length(fam) != q_fit) {
  stop("stage 1 output lacks valid family information; rerun ",
       "sim/lmc/sim_lmc_structure.R before prediction", call. = FALSE)
}
if (family_mix) {
  if (!exists("require_all_families", inherits = FALSE)) require_all_families <- TRUE
  if (!exists("family_selection", inherits = FALSE) ||
      (require_all_families && !all(LMC_FAMILIES %in% fam)) ||
      !identical(fam, family_selection$family) ||
      is.null(selection_aic)) {
    stop("mixed-family stage 1 output lacks valid AIC selection; rerun ",
         "sim/lmc/sim_lmc_structure.R before prediction", call. = FALSE)
  }
} else if (!identical(fam, rep("matern", q_fit)) ||
           !is.null(selection_aic) || !is.null(family_selection) ||
           !grepl("^Fixed Matern family;", selection_method)) {
  stop("all-Matern stage 1 output is not a fixed-Matern fit; rerun ",
       "sim/lmc/sim_lmc_structure.R before prediction", call. = FALSE)
}

P <- ncol(Y_mat); n <- nrow(grid)
q_pad <- q                       # number of factors reported (settings.R)
stopifnot(q_fit <= q_pad)
cat(sprintf("q = %d of %d retained (%.1f%% of signal variance); families: %s\n",
            q_fit, q_pad, 100 * ld$frac, paste(fam, collapse = " ")))

D  <- fields::rdist(grid)
Dx <- fields::rdist(grid_test, grid)

# 1. Joint maximum likelihood
bnd <- lmc_bounds(fam)
# The nugget is on the log scale, so a zero pilot nugget starts at the lower bound.
theta0 <- lmc_pack(par_raw, max(nugget_hat, exp(tail(bnd$lower, 1L))))
nll <- if (identical(missing_TF, "coloc")) {
  lmc_nll_coloc_factory(Y_mat, D, V, d, fam)
} else {
  lmc_nll_miss_factory(Y_mat, D, B, fam)
}

# Everything the fit depends on; a saved result or checkpoint is reused only
# when this matches exactly.
identity <- list(setting = setting_name, replicate = array_id,
                 data_md5 = unname(tools::md5sum(data_file)),
                 structure_md5 = unname(tools::md5sum(est_file)),
                 family = fam, q_pad = q_pad, lower = bnd$lower, upper = bnd$upper)
previous <- NULL
if (file.exists(result_file)) {
  old <- new.env()
  load(result_file, envir = old)
  previous <- old$result_lmc$opt_lmc$optimization
  if (is.null(previous))
    stop("existing result predates resumable fits; delete it to refit: ",
         result_file, call. = FALSE)
}
opt <- lmc_optimize(nll, theta0, bnd$lower, bnd$upper, identity, checkpoint_file,
                    resume = resume, previous = previous,
                    elapsed_before = as.numeric(time_lmc_est),
                    save_every = as.numeric(Sys.getenv("SSGGP_LMC_CHECKPOINT_SECONDS", "300")))
if (is.null(opt)) {
  cat("converged LMC result already exists:", result_file, "\n")
  quit(status = 0L)
}
fit <- opt$fit
cat(sprintf("optim: nll %.4f, convergence %d, %d evaluations over %d pass(es)\n",
            fit$value, fit$convergence, opt$evaluations, opt$passes))

est <- lmc_unpack(fit$par, fam)
est_par <- est$par
nugget_mle <- est$nugget
cat(sprintf("nugget: stage 1 (covariogram) %.4f -> stage 2 (MLE) %.4f\n",
            nugget_hat, nugget_mle))
time_lmc <- opt$elapsed_min        # stage 1 and all optimization passes

# 2. Prediction
t_pred <- Sys.time()
pred <- if (identical(missing_TF, "coloc")) {
  lmc_predict_coloc(Y_mat, D, Dx, V, d, fam, est_par, nugget_mle)
} else {
  lmc_predict_miss(Y_mat, D, Dx, B, fam, est_par, nugget_mle)
}
if (!all(is.finite(pred$mean)) ||
    !all(vapply(pred$cov, function(x) all(is.finite(x)) &&
                min(diag(x)) >= -1e-8, logical(1)))) {
  stop("nonfinite predictions or invalid predictive variances; refusing to save results")
}
Y_mat_pred    <- pred$mean
pred_cov_list <- pred$cov
cat(sprintf("prediction: %.2f min\n",
            as.numeric(difftime(Sys.time(), t_pred, units = "mins"))))

# Parameter summaries padded to q_pad factors.
est_B <- matrix(NA_real_, P, q_pad); est_B[, seq_len(q_fit)] <- B
canon <- vapply(seq_len(q_fit), function(k) lmc_canonical(fam[k], est_par[[k]]),
                numeric(2))
est_a  <- rep(NA_real_, q_pad); est_a[seq_len(q_fit)]  <- canon["a", ]
est_nu <- rep(NA_real_, q_pad); est_nu[seq_len(q_fit)] <- canon["nu", ]
est_nuggets <- rep(nugget_mle, P)

result <- list(est_B = est_B, est_a = est_a, est_nu = est_nu,
               est_nuggets = est_nuggets,
               family = fam, par_raw = est_par, par_init = par_raw,
               q_fit = q_fit, var_frac = ld$frac, d = d,
               nugget_covariogram = nugget_hat, nugget_mle = nugget_mle,
               selection_method = selection_method,
               family_selection = family_selection, selection_aic = selection_aic,
               optimization = opt[names(opt) != "fit"])

# 3. Prediction summary
pred_var_mat <- t(vapply(pred_cov_list, diag, numeric(P)))
miss_new <- lapply(seq_len(P), function(j) is.na(Y_mat_test_obs[, j]))
cov_j <- mse_j <- numeric(P)
for (k in seq_len(P)) {
  ix <- miss_new[[k]] & !is.na(Y_mat_test[, k])
  se <- 1.96 * sqrt(pred_var_mat[ix, k])
  cov_j[k] <- mean(Y_mat_pred[ix, k] - se <= Y_mat_test[ix, k] &
                   Y_mat_test[ix, k] <= Y_mat_pred[ix, k] + se)
  mse_j[k] <- mean((Y_mat_pred[ix, k] - Y_mat_test[ix, k])^2)
}
cat(sprintf("coverage: mean %.3f [%.3f, %.3f];  rmse: mean %.3f\n",
            mean(cov_j), min(cov_j), max(cov_j), sqrt(mean(mse_j))))

result_lmc <- list(Y = Y_mat, grid = grid,
                   Y_test = Y_mat_test, grid_test = grid_test,
                   opt_lmc = result,
                   prediction = Y_mat_pred,
                   pred_var = pred_cov_list,
                   new_obs_ind = new_obs_ind)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save(time_lmc, result_lmc,
     file = result_file)
cat("wrote", result_file,
    sprintf("(estimation %.2f min)\n", time_lmc))
