# A short demonstration of the SSGGP package with p = 3
# processes, few features and iterations. It runs in about a minute, and the
# estimates are not meant to be accurate.
#
#   Rscript demo.R
#
# Sections:
#   1. simulate a colocated multivariate process
#   2. fit with a shared nugget
#   3. make the data misaligned and fit again
#   4. predict, and score the predictions
#   5. compare against the independent-process baseline
#   6. select the penalty by BIC
#   7. define a spectral family and fit a mixed-family process

library(SSGGP)
set.seed(1)

## 1. Simulate
# Zeros in sigma^-1 are the conditional independence graph. A tridiagonal
# precision means each process is conditionally independent of all but its
# neighbours.
P <- 3
Sigma_inv <- trimat(P, diag = 2, offdiag = 1)
cat("true precision:\n"); print(Sigma_inv)
cat("true graph (off-diagonal non-zeros):\n"); print(graph(Sigma_inv))

locs_all <- as.matrix(expand.grid(seq(0, 5, length = 12), seq(0, 5, length = 12)))
train <- 1:100
test  <- setdiff(seq_len(nrow(locs_all)), train)

# All three processes are Matern with a = 3 and nu = 1.5. The Matern proposal
# has a = 3 and nu = 1.35; its smaller nu gives it heavier tails than the target.
f_par_true <- c(a = 3, nu = 1.5)
proposal   <- matern_proposal(a = 3, nu = 1.35)

# Simulate with random Fourier features.
sim_out <- sim(locs_all, sigma_inv = Sigma_inv, f_par = f_par_true,
               f_fun = "matern", nugget = 0.1, method = "rff", L = 600,
               g = proposal, seed = 1)
Y_all <- sim_out$Y
cat("\nsimulated", nrow(Y_all), "locations x", P, "processes\n")
cat("empirical marginal variances:", format(apply(Y_all, 2, var), digits = 3), "\n")

Y  <- Y_all[train, ];  locs  <- locs_all[train, ]
Yt <- Y_all[test, ];   locst <- locs_all[test, ]

## 2. Shared nugget
cat("\n== shared nugget ==\n")
ctl <- list(
  f_fun = "matern",
  f_par = c(a = 2.8, nu = 1.6),
  g = proposal,
  L = 1000,
  MC = 50,
  iter = 30,
  lag = 10,
  rho = 0.05,
  seed = 2
)
fit_shared <- do.call(fit, c(list(Y = Y, locs = locs, nugget = "shared"), ctl))
print(fit_shared)
cat("\ngraph recovery:\n")
print(graph_metrics(fit_shared, Sigma_inv)[c("TPR", "FPR", "F1")])

## 3. Misaligned observations
# NA marks a process unobserved at a location.
Y_miss <- Y
for (j in seq_len(P)) Y_miss[sample(nrow(Y), 40), j] <- NA
cat("observed fraction:", format(mean(!is.na(Y_miss)), digits = 3), "\n")

fit_miss <- do.call(fit, c(list(Y = Y_miss, locs = locs, nugget = "by_process"), ctl))
print(fit_miss)

## 4. Predict and score
cat("\n== prediction ==\n")
# summarise = "blockdiag" returns the per-location P x P covariance rather than
# the dense cross-location covariance matrix.
pred <- predict(fit_shared, locst, L = 400, MC = 200, seed = 3)
print(pred)

m <- prediction_metrics(pred, Yt)
cat("\nper-process MSE     :", format(m$mse, digits = 3), "\n")
cat("per-process coverage:", format(m$coverage, digits = 3), "  (nominal 0.95)\n")
cat("sum across processes: MSE", format(m$sum$mse, digits = 3),
    " coverage", format(m$sum$coverage, digits = 3), "\n")
# The variance of the sum includes the cross-process covariances, so it differs
# from the sum of the marginal variances.
cat("var of the sum      :", format(mean(var_sum(pred)), digits = 4),
    " vs sum of marginals", format(sum(colMeans(pred$var)), digits = 4), "\n")

## 5. Independent baseline
# Each process fitted alone as a Matern process, so sigma is diagonal.
cat("\n== independent Matern baseline ==\n")
fit_ind <- fit_independent(Y, locs, par_init = c(1, 2, 1, 0.07))
print(fit_ind)
cat("its graph is empty by construction:", sum(graph(fit_ind$sigma_inv)) == 0, "\n")

## 6. Penalty selection by BIC
# ssggp_bic() fits the model once per value of `rho`, then scores every fit by
# BIC = -2 logLik + log(n) k on the integrated likelihood, where k is the
# number of selected edges. `n`, `rho`, `L_ll` and `seed_ll` belong to
# ssggp_bic(); every other argument is passed unchanged to each fit() call.
# All fits are scored with the same L_ll features (seed_ll), so their
# likelihoods are directly comparable. ssggp_aic() works the same way, without `n`.
cat("\n== BIC over the penalty grid ==\n")
sel <- ssggp_bic(
  Y, locs,
  n       = nrow(Y),                   # sample size: number of training locations
  rho     = c(0.01, 0.05, 0.2),        # penalty grid
  L_ll    = 400,                       # features for the likelihood evaluation
  seed_ll = 1,                         # same features for every penalty
  # passed to fit():
  f_fun = "matern", f_par = c(a = 2.8, nu = 1.6), g = proposal,
  L = 1000, MC = 50, iter = 30, lag = 10, seed = 2
)
print(sel$table, digits = 5)
cat("selected rho:", sel$rho, "\n")
cat("its graph vs truth:\n"); print(graph_metrics(sel$best, Sigma_inv)[c("TPR","FPR","F1")])

## 7. A self-defined family, mixed with Matern and CH
cat("\n== a self-defined family, mixed with Matern and CH ==\n")

# gcauchy: an isotropic (hence stationary) spectral density for d = 2,
#   f(x; b) = exp(-|x| / b) / (2 pi b^2),
# with covariance C(h) = (1 + b^2 |h|^2)^(-3/2). Like every family it must
# integrate to one, since all marginal scale is in sigma.
gcauchy_vec <- function(X, par) {          # X: one frequency per row
  b <- par[["b"]]
  exp(-sqrt(rowSums(X^2)) / b) / (2 * pi * b^2)
}
gcauchy <- function(x, par) gcauchy_vec(matrix(x, nrow = 1), par)
cat("integral of f:", check_unit_integral(gcauchy, c(b = 0.5)), "\n")

# fit() needs only the density, its parameter names and log-scale bounds.
register_spectral(
  name      = "gcauchy",
  f         = gcauchy,
  f_vec     = gcauchy_vec,
  par_names = "b",
  par_lower = log(c(b = 0.01)),
  par_upper = log(c(b = 10)),
  overwrite = TRUE)
cat("registered families:", paste(list_spectral(), collapse = ", "), "\n")

# One family per process, observed at 400 random locations.
f_fun_mix <- list("matern", "ch", "gcauchy")
f_par_mix <- list(c(a = 3, nu = 1.5),
                  c(beta = 0.6, nu = 2, alpha = 1.5),   # CH needs alpha > d/2
                  c(b = 0.5))
locs_mix <- matrix(stats::runif(800, 0, 5), ncol = 2)
mix <- sim(locs_mix, sigma_inv = Sigma_inv, f_fun = f_fun_mix, f_par = f_par_mix,
           nugget = 0.1, method = "rff", L = 5000, g = proposal, seed = 4)

# SAEM moves range parameters slowly, so start near the truth here; in
# practice, take starting values from fit_independent().
fit_mix <- fit(mix$Y, locs_mix, f_fun = f_fun_mix,
               f_par = list(c(a = 2.5, nu = 1.2),
                            c(beta = 0.8, nu = 1.5, alpha = 2),
                            c(b = 0.4)),
               g = proposal, rho = 0.01, L = 1000, MC = 30, iter = 30,
               lag = 10, seed = 2, trace_every = 0)
print(rbind(truth    = unlist(f_par_mix),
            estimate = unlist(fit_mix$f_par)), digits = 3)

cat("\ndemo complete.\n")
