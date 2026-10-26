# Simulation settings, sourced by every stage. Run scripts from paper/.

library(SSGGP)
library(Matrix)

# paper/, located from a working directory of paper/, paper/sim/, or below.
.find_project_root <- function() {
  candidates <- unique(normalizePath(c(".", "..", "../.."), mustWork = TRUE))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "..", "pkg", "DESCRIPTION")) &&
        file.exists(file.path(candidate, "sim", "settings.R"))) {
      return(candidate)
    }
  }
  stop("could not locate the paper/ directory; run from paper/, paper/sim/, ",
       "or a folder under paper/sim/", call. = FALSE)
}

project_root <- .find_project_root()
here <- function(...) file.path(project_root, ...)

# Setting: four binary factors, taken from the environment so that one launcher
# can run any of the sixteen combinations.
prec_structure <- Sys.getenv("SSGGP_PREC",   "random")      # tridiagonal | random
offdiag_sign   <- Sys.getenv("SSGGP_SIGN",   "neg")         # pos | neg
missing_TF     <- Sys.getenv("SSGGP_MISS",   "miss")       # coloc | miss
family_design  <- Sys.getenv("SSGGP_FAMILY", "mixed")       # mixed | matern
if (!family_design %in% c("mixed", "matern"))
  stop("SSGGP_FAMILY must be 'mixed' or 'matern'", call. = FALSE)

# Data-generating design
P          <- as.integer(Sys.getenv("SSGGP_P", "50"))
dim_grid   <- 25          # dim_grid^2 locations in total
train_prop <- 0.8         # the rest is held out for prediction
n_train    <- floor(dim_grid^2 * train_prop)   # training locations
n_sim      <- 15          # replicate datasets per setting
d          <- 2           # spatial dimension
domain_max <- 5           # locations are uniform on [0, domain_max]^2

# Process families (in random order) and their true parameters. local() keeps
# the caller's RNG state unchanged. The all-Matern design uses nu = 1.5, a = 3.
truth_seed <- 20260811L
truth_design <- if (family_design == "matern") {
  list(family = rep("matern", P),
       parameters = replicate(P, c(nu = 1.5, a = 3), simplify = FALSE))
} else local({
  set.seed(truth_seed)
  family <- sample(c(rep("matern", 20), rep("ch", 20),
                     rep("se", 10)))
  parameters <- lapply(family, function(f) switch(f,
    matern = c(nu=sample(c(1.5, 2.5),1), a=runif(1, 2.5, 3.5)),
    # CH alpha must exceed d/2.
    ch = c(beta=runif(1, 0.6, 0.9), alpha=runif(1, 1.25, 1.75),
           nu=sample(c(1.5, 2.5),1)),
    se = c(a=runif(1, 1.75, 2.25))))
  list(family=family[seq_len(P)], parameters=parameters[seq_len(P)])
})


f_fun_true <- as.list(truth_design$family)
f_par_true <- truth_design$parameters
nugget_true <- 0.1

# Mixed-family cross-covariances have no closed form, so data are simulated by RFF.
sim_method <- "rff"
sim_L      <- 100000       # only used when sim_method = "rff"

# Estimation
L_vals   <- as.integer(Sys.getenv("SSGGP_L", "20000")) # random Fourier features
MC_vals  <- as.integer(Sys.getenv("SSGGP_MC", "30")) # Monte Carlo draws per E step
lag_vals <- c(200)         # gamma_m = lag / (m - 1 + lag)

# Penalty grid. The first value is the dense fit (no selection); BIC selects
# among the others. The three smallest positive values refine the interval
# between the dense fit and rho = 0.003.
rho_vals <- c(0, 0.0003, 0.0006, 0.001, 0.003, 0.006, 0.012, 0.025, 0.05, 0.1)

iter     <- 500           # maximum SAEM iterations
cm_cycles <- 1L           # CM cycles per expensive E step
tol      <- 1e-3          # convergence tolerance on the objective
diag_pen <- FALSE         # whether to penalise the diagonal
nugget_type <- "shared"   # "shared" or "by_process"
use_independent_init <- TRUE


# Starting values for the independent fits, which in turn initialize the joint
# fit (use_independent_init).
f_fun_init <- f_fun_true
f_par_init <- lapply(unlist(f_fun_init), function(f) switch(f,
  matern = c(a=2, nu=1),
  ch = c(beta=0.75, nu=2, alpha=1.5),
  se = c(a=2)))
nugget_init <- 0.05


# Importance-sampling proposals; their tails must be heavier than every f_j.
proposal_sim  <- matern_proposal(a = 3, nu = 1.25)
proposal      <- matern_proposal(a = 3, nu = 1.25)
proposal_pred <- matern_proposal(a = 3, nu = 1.25)   # may differ from the above

# Likelihood and prediction
L_ll      <- 1000         # features for the marginal likelihood
L_ll_diag <- 10000        # features for the diagnostic likelihood traces
L_pred    <- 50000        # features for prediction
MC_pred   <- 10000        # Gibbs sweeps for prediction (the first half is burn-in)
seed_ll   <- 1L           # shared across penalties so the BIC comparison is paired

# Compute: threads within each fit
threads   <- 32L

# LMC baseline: saved results are padded to q factor columns. The fits use 10
# factors (random precision) or 6 (tridiagonal); see sim/lmc/sim_lmc_structure.R.
q <- as.integer(Sys.getenv("SSGGP_LMC_Q", as.character(floor(0.25 * P))))

# Derived quantities

# One row per estimation job: THISJOBVALUE = i runs row i, so the --array range
# of slurm/sim/main/sim_main.sh is 0..(nrow(grid_params) - 1). `sim` varies fastest.
grid_params <- expand.grid(sim = seq_len(n_sim), rho = rho_vals,
                           L = L_vals, lag = lag_vals, MC = MC_vals)

# Mixed-family and all-Matern designs are stored in separate directories.
family_mix   <- !all(unlist(f_fun_true) == "matern")

setting_name <- paste(prec_structure, offdiag_sign, missing_TF, sep = "_")
if (family_mix) setting_name <- paste0(setting_name, "_mix")
data_dir     <- here("sim/data", setting_name)
out_dir      <- here("sim/output/simulation", setting_name)
fig_dir      <- here("sim/output/figure", setting_name)

# True precision matrix of the active setting.
true_precision <- function() {
  if (prec_structure == "tridiagonal") {
    trimat(P, diag = 2, offdiag = if (offdiag_sign == "pos") 1 else -1)
  } else {
    if (offdiag_sign == "pos") {
      n_non_zero <- as.integer(round(0.8 * P))
      prec_seed <- 91772;   offdiag_val <-  0.8
    } else {
      n_non_zero <- 50; prec_seed <- 1118335; offdiag_val <- -0.8
    }
    set.seed(prec_seed)
    Sigma_inv <- matrix(0, P, P)
    # Pairs are drawn with replacement, symmetrized, and the diagonal reset.
    random_ij <- cbind(sample(1:P, n_non_zero, replace = TRUE),
                       sample(1:P, n_non_zero, replace = TRUE))
    random_ij <- rbind(random_ij, cbind(random_ij[, 2], random_ij[, 1]))
    Sigma_inv[random_ij] <- offdiag_val
    diag(Sigma_inv) <- 2
    Sigma_inv
  }
}

# Missingness: miss_obs_prop[j] is the share of training locations removed for
# process j; miss_new_obs_prop[j] = 1 holds out every test value of process j.
missing_prop <- if (missing_TF == "miss") {
  list(miss_obs_prop = rep(0.4, P), miss_new_obs_prop = rep(1, P))
} else {
  list(miss_obs_prop = rep(0,   P), miss_new_obs_prop = rep(1, P))
}

# Thread setup, called at the start of every stage.
setup_compute <- function() {
  options(SSGGP.threads = threads)
  invisible(TRUE)
}

# Array index from the launcher; 1 when run interactively.
job_index <- function(default = 1L) {
  v <- as.integer(Sys.getenv("THISJOBVALUE"))
  if (is.na(v)) default else v
}
