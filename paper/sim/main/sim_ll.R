# Marginal log-likelihood of each penalized fit, for BIC selection.
#   sbatch slurm/sim/main/sim_ll.sh          (one task per fit)
# All penalties use the same features and seed (seed_ll), so their likelihoods
# are compared on equal terms. n_par is the number of selected edges; parameters
# whose number does not change with the penalty are left out of the BIC.
# Writes ll_<row>.RData with `ll_k` and `n_par`.

rm(list = ls())
if (!file.exists("settings.R")) {  # work from sim/, where settings.R lives
  if (file.exists(file.path("sim", "settings.R"))) setwd("sim") else
  if (file.exists(file.path("..", "settings.R"))) setwd("..") else
    stop("run this script from paper/")
}
source("settings.R")
setup_compute()

array_id <- job_index(1L)
est_file <- file.path(out_dir, paste0("est", array_id, ".RData"))

if (!file.exists(est_file)) stop("missing fit: ", est_file, call. = FALSE)

load(est_file)                        # result_est, time_approx

ll <- logLik(result_est$fit, L = L_ll, seed = seed_ll, g = proposal,
             threads = threads)

ll_k  <- as.numeric(ll)
n_par <- attr(ll, "df")               # = sum_{i<j} 1(eta_ij != 0)

cat(sprintf("array %d | rho %g | logLik %.4f | edges %d\n",
            array_id, result_est$rho, ll_k, n_par))

save(ll_k, n_par, file = file.path(out_dir, paste0("ll_", array_id, ".RData")))
