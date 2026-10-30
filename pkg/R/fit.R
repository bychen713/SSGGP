#' Fit a sparse spectral graphical Gaussian process
#'
#' Estimates a sparse inverse cross-spectral density by stochastic
#' approximation EM with a graphical LASSO penalty, using random Fourier
#' features for the latent process and a coordinate-wise ECM update for the
#' process-specific spectral parameters.
#'
#' Handles both colocated and misaligned designs: `NA` entries in `Y` mark
#' locations where a process is unobserved. The nugget may be shared across
#' processes or estimated per process; `"shared"` is the constrained case
#' \eqn{d_1 = \dots = d_P}.
#'
#' @param Y `n` by `P` matrix of observations; `NA` where unobserved.
#' @param locs `n` by `d` matrix of locations.
#' @param f_par Named numeric `c(a =, nu =)` of starting values, recycled to all
#'   processes, or a list of `P` such vectors.
#' @param f_fun Spectral family: a registered name (see [list_spectral()]), or a
#'   list of length `P` of names. Raw density functions must be registered with
#'   [register_spectral()] first. Densities must integrate to one; see
#'   [check_unit_integral()].
#' @param g Importance-sampling proposal, from [matern_proposal()].
#' @param rho Graphical LASSO penalty, passed to `glasso::glasso()` unchanged.
#'   CM1 maximizes
#'   \eqn{(L/2)[\log|Q| - \operatorname{tr}(\hat S Q)] - (L/2)\rho\|Q\|_1}
#'   over the off-diagonal (or all of `Q` when `diag_pen = TRUE`), so the
#'   penalty coefficient on the observed-data log-likelihood is
#'   \eqn{\lambda = L\rho/2}, i.e. \eqn{\rho = 2\lambda/L}. Keep `L * rho`
#'   fixed to hold \eqn{\lambda} fixed across `L`; the estimates can still
#'   differ, because the random-feature likelihood itself changes with `L`.
#' @param nugget One of `"shared"` or `"by_process"`.
#' @param nugget_init Starting nugget variance.
#' @param sigma_inv_init Optional `P` by `P` starting value for
#'   \eqn{\sigma^{-1}}, symmetric positive definite. `NULL` uses the default
#'   described under Details.
#' @param L Number of random Fourier features.
#' @param MC Gibbs sweeps per E step.
#' @param burnin Sweeps discarded at the start of each E step. Defaults to
#'   `MC / 2`; the remaining `MC - burnin` draws form the Monte Carlo average.
#' @param iter Maximum EM iterations.
#' @param lag SAEM step-size lag: at iteration `m` the step size is
#'   \eqn{\gamma_m = \mathrm{lag} / (m - 1 + \mathrm{lag})}, so
#'   \eqn{\gamma_1 = 1}.
#' @param cm_cycles Number of CM1/CM3 cycles per E step. CM2 and the
#'   stochastic-approximation update run once per iteration; further cycles
#'   repeat CM1 and CM3 with the same statistic.
#' @param numerical_guard What to do when the importance ratio \eqn{R = g/f}
#'   is not representable, which happens when a light-tailed target (a
#'   squared exponential (`se`) process under a Matern proposal) has essentially no spectral
#'   mass at a drawn frequency. `"error"` (default) stops. `"cap"` clamps
#'   \eqn{R} at `1e300`, i.e. floors \eqn{f} at \eqn{g/10^{300}}; the target
#'   mass discarded on the clamped set \eqn{A} is at most \eqn{P_g(A)/10^{300}}.
#'   Under `"cap"` CM3 is evaluated in \eqn{R} and `cm3_gradient` is ignored.
#'   Available only for the built-in families. Counts are recorded per
#'   iteration in `chain$capped`.
#' @param cm3_gradient Give `optim()` the analytic gradient of \eqn{Q_j} in
#'   CM3 (Matern and squared exponential only). It reduces the number of objective
#'   evaluations and moves the optimum within finite-difference noise, so it is
#'   off by default.
#' @param checkpoint_file Optional path to an RDS checkpoint. When supplied,
#'   the fit saves its resumable state and resumes from a compatible existing
#'   checkpoint.
#' @param checkpoint_every Save every this many completed iterations. A final
#'   checkpoint is always written when the call stops.
#' @param run_until Optional iteration at which this call stops, for staged
#'   runs; `iter` remains the overall target.
#' @param tol Convergence tolerance on the change in the complete-data
#'   penalized criterion between iterations.
#' @param diag_pen Penalize the diagonal in the graphical LASSO.
#' @param lower,upper Log-scale box constraints on the spectral parameters
#'   \eqn{\iota_j}, in each family's `par_names` order -- `(a, nu)` for Matern.
#'   A component with `lower == upper` is held fixed. `NULL` uses each family's
#'   registered bounds, which is required when families with different
#'   parameter counts are mixed.
#' @param threads BLAS/OpenMP threads, capped for the duration of the call and
#'   restored afterwards. `NULL` leaves the setting untouched. Defaults to
#'   `getOption("SSGGP.threads", 5)`.
#' @param seed Integer seed.
#' @param keep_W What to retain of the latent draws: `"none"` (default),
#'   `"last"` for the final iteration's, or `"all"` for one entry per iteration
#'   (`P * L * (MC - burnin)` doubles each).
#' @param verbose Print progress and per-stage timings.
#' @param trace_every Print the current parameter estimates every this many
#'   iterations (and at iteration 1); `0` disables it.
#' @section Algorithm:
#' Each iteration runs, in order: an E step of `MC` Gibbs sweeps started from
#' the previous iteration's final state; the stochastic-approximation update of
#' the per-feature second moments \eqn{u_\ell} and of the nugget; CM1 (graphical
#' LASSO for \eqn{\sigma^{-1}}); CM2 (nugget); and CM3, a sequential
#' coordinate-wise update of each process's spectral parameters with the others
#' held at their current values. Starting values: unless `sigma_inv_init` is
#' given, \eqn{\sigma^{-1}} is the inverse sample covariance when the data
#' are colocated and the identity otherwise; [fit_independent()] is a good
#' source of `f_par`, `nugget_init` and a diagonal `sigma_inv_init`.
#' @return An object of class `ssggp_fit`.
#' @export
fit <- function(Y, locs,
                f_par = c(a = 2, nu = 1),
                f_fun = "matern",
                g = matern_proposal(a = 7, nu = 1.25),
                rho = 0.05,
                nugget = c("shared", "by_process"),
                nugget_init = 0.07, sigma_inv_init = NULL,
                L = 300, MC = 50, burnin = NULL, iter = 100, lag = 15,
                tol = 1e-1, diag_pen = FALSE,
                lower = NULL, upper = NULL,
                threads = getOption("SSGGP.threads", 5L),
                seed = NULL,
                keep_W = c("none", "last", "all"),
                verbose = FALSE, trace_every = 10L,
                cm_cycles = 1L, cm3_gradient = FALSE,
                numerical_guard = c("error", "cap"),
                checkpoint_file = NULL, checkpoint_every = 25L,
                run_until = NULL) {

  nugget <- match.arg(nugget)
  keep_W <- match.arg(keep_W)
  numerical_guard <- match.arg(numerical_guard)
  threads <- resolve_threads(threads)
  with_threads(threads, fit_impl(
    Y, locs, f_par, f_fun, g, rho, nugget, nugget_init,
    L, MC, burnin, iter, lag, cm_cycles, tol, diag_pen,
    lower, upper, threads, seed, keep_W, verbose, trace_every,
    checkpoint_file, checkpoint_every, run_until, cm3_gradient,
    numerical_guard, sigma_inv_init))
}

# SAEM step size gamma_m, shared by the sufficient-statistic and nugget updates
.saem_gain <- function(i, lag) lag / (i - 1 + lag)

#' @keywords internal
#' @noRd
fit_impl <- function(Y, locs, f_par, f_fun, g, rho, nugget,
                     nugget_init, L, MC, burnin, iter, lag, cm_cycles, tol,
                     diag_pen,
                     lower, upper, threads, seed,
                     keep_W, verbose, trace_every = 10L,
                     checkpoint_file = NULL, checkpoint_every = 25L,
                     run_until = NULL, cm3_gradient = FALSE,
                     numerical_guard = "error", sigma_inv_init = NULL) {
  Y <- as.matrix(Y); locs <- as.matrix(locs)
  n <- nrow(Y); P <- ncol(Y); d <- ncol(locs)
  validate_fit_inputs(Y, locs, L, MC, burnin, iter, rho, nugget_init,
                      lag, cm_cycles, tol)
  if (!is.null(checkpoint_file) &&
      (!is.character(checkpoint_file) || length(checkpoint_file) != 1L ||
       is.na(checkpoint_file) || !nzchar(checkpoint_file)))
    stop("`checkpoint_file` must be NULL or one non-empty path", call. = FALSE)
  if (!is.logical(cm3_gradient) || length(cm3_gradient) != 1L ||
      is.na(cm3_gradient))
    stop("`cm3_gradient` must be TRUE or FALSE", call. = FALSE)
  checkpoint_every <- as.integer(checkpoint_every)
  if (length(checkpoint_every) != 1L || is.na(checkpoint_every) ||
      checkpoint_every < 1L)
    stop("`checkpoint_every` must be a positive integer", call. = FALSE)
  if (is.null(run_until)) run_until <- iter
  run_until <- as.integer(run_until)
  if (length(run_until) != 1L || is.na(run_until) || run_until < 1L ||
      run_until > iter)
    stop("`run_until` must be between 1 and `iter`", call. = FALSE)

  check_proposal(g)
  fam <- resolve_spectral(f_fun, P)
  require_spectral_estimable(fam)
  f_den <- lapply(fam, `[[`, "f")
  f_vec <- lapply(fam, f_vec_of)

  if (!is.list(f_par)) f_par <- rep(list(f_par), P)
  stopifnot(length(f_par) == P)
  for (j in seq_len(P)) {
    need <- fam[[j]]$par_names
    if (!all(need %in% names(f_par[[j]]))) {
      stop("f_par[[", j, "]] must be named (",
           paste(need, collapse = ", "), ") for family '",
           fam[[j]]$name, "'", call. = FALSE)
    }
  }

  # observations in process-fastest order, matching as.vector(t(Y))
  observed    <- !is.na(Y)
  Y_reordered <- as.vector(t(Y))
  obs_idx     <- which(!is.na(Y_reordered))
  Y_obs       <- Y_reordered[obs_idx]
  proc_id     <- rep(seq_len(P), n)[obs_idx]

  # box for each iota_j: the family's own bounds unless overridden
  bounds <- lapply(seq_along(fam), function(j) {
    e  <- fam[[j]]
    lo <- if (is.null(lower)) e$par_lower else lower
    up <- if (is.null(upper)) e$par_upper else upper
    k  <- length(e$par_names)
    if (is.null(lo) || is.null(up)) {
      stop("family '", e$name, "' (process ", j, ") was registered without ",
           "`par_lower`/`par_upper`; re-register it with bounds, or pass ",
           "`lower` and `upper` in the order (",
           paste(e$par_names, collapse = ", "), ").", call. = FALSE)
    }
    if (length(lo) != k || length(up) != k) {
      stop("process ", j, " uses family '", e$name, "' with spectral parameters (",
           paste(e$par_names, collapse = ", "), "), so `lower` and `upper` ",
           "must have length ", k, "; leave them NULL when mixing families ",
           "with different parameter counts.", call. = FALSE)
    }
    if (any(lo > up)) stop("`lower` exceeds `upper` for process ", j,
                           call. = FALSE)
    list(lower = unname(lo), upper = unname(up))
  })

  if (!is.null(seed)) set.seed(seed)

  # --- initialization --------------------------------------------------------
  BurnIn   <- if (is.null(burnin)) MC %/% 2 else burnin
  n_keep   <- MC - BurnIn
  complete <- rowSums(!observed) == 0
  sigma_est <- if (any(!observed)) diag(1, P) else stats::var(Y[complete, , drop = FALSE])
  sigma_inv_est <- solve(sigma_est)
  if (!is.null(sigma_inv_init)) {
    sigma_inv_init <- as.matrix(sigma_inv_init)
    if (!is.numeric(sigma_inv_init) || !identical(dim(sigma_inv_init), c(P, P)) ||
        any(!is.finite(sigma_inv_init)) || !isSymmetric(unname(sigma_inv_init)) ||
        inherits(try(chol(sigma_inv_init), silent = TRUE), "try-error"))
      stop("`sigma_inv_init` must be a symmetric positive definite P x P matrix",
           call. = FALSE)
    sigma_inv_est <- unname(sigma_inv_init)
  }

  d_cur <- rep(nugget_init, P)
  feat  <- draw_features(L, d, g)

  # `use_fast` (all processes observed everywhere) is part of the checkpoint
  # spec. SSGGP.force_general selects the sparse reference branch of the
  # kernels, which forms A = c_vals (x) I_P: a slower reference against which
  # the dense branches can be checked.
  use_fast      <- all(observed)
  force_general <- isTRUE(getOption("SSGGP.force_general", FALSE))
  masked        <- !use_fast && !force_general
  capping       <- identical(numerical_guard, "cap")
  tri           <- tri_pack(P)

  # P x n observation mask, in the residual's layout, and Y on the full grid
  # with unobserved entries at zero (misaligned dense branch only)
  obs_mask <- if (masked) t(observed) * 1 else matrix(numeric(0), 0, 0)
  Y_gibbs  <- if (masked) {
    v <- Y_reordered; v[is.na(v)] <- 0; v
  } else Y_obs

  chain <- list(
    sigma_inv = matrix(NA_real_, nrow = P * (P + 1) / 2, ncol = iter),
    nugget    = matrix(NA_real_, nrow = P, ncol = iter),
    # a list: with mixed families the parameter names differ between processes
    spectral  = vector("list", iter),
    loglik    = matrix(NA_real_, nrow = iter, ncol = 2,
                       dimnames = list(NULL, c("prec", "nugget"))),
    # importance ratios at the ceiling, per iteration and process
    capped    = matrix(0L, nrow = iter, ncol = P))

  u_sa <- NULL              # stochastic-approximation state, L x P(P+1)/2
  prev_W <- NULL
  warned_capped <- FALSE
  timing <- matrix(NA_real_, iter, 6,
                   dimnames = list(NULL, c("estep", "moments", "cm1", "cm2",
                                           "cm3", "refresh")))
  .tic <- Sys.time()
  tock <- function() {
    now <- Sys.time(); d <- as.numeric(difftime(now, .tic, units = "secs"))
    .tic <<- now; d
  }
  W_keep <- if (keep_W == "all") vector("list", iter) else NULL
  converged <- FALSE

  # The checkpoint holds only the stochastic and mutable state; everything
  # derived from the feature draw is rebuilt on resume.
  checkpoint_spec <- list(
    version = 3L, Y = Y, locs = locs, P = P, d = d, L = L, MC = MC,
    BurnIn = BurnIn, iter = iter, lag = lag, cm_cycles = cm_cycles,
    rho = rho, nugget = nugget, nugget_init = nugget_init, tol = tol,
    diag_pen = diag_pen, keep_W = keep_W, seed = seed,
    f_par_init = f_par, lower = lower, upper = upper, use_fast = use_fast,
    family = vapply(fam, function(e) e$name, character(1)), g_par = g$par,
    numerical_guard = numerical_guard, cm3_gradient = cm3_gradient,
    sigma_inv_init = sigma_inv_init)
  completed <- 0L
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    ck <- readRDS(checkpoint_file)
    if (!is.list(ck) || !identical(ck$spec, checkpoint_spec))
      stop("checkpoint is incompatible with this fit: ", checkpoint_file,
           call. = FALSE)
    completed <- ck$completed
    if (!is.numeric(completed) || length(completed) != 1L ||
        completed < 0L || completed > iter)
      stop("checkpoint has an invalid completed iteration", call. = FALSE)
    feat <- ck$state$feat
    f_par <- ck$state$f_par
    sigma_inv_est <- ck$state$sigma_inv_est
    d_cur <- ck$state$d_cur
    u_sa <- ck$state$u_sa
    prev_W <- ck$state$prev_W
    chain <- ck$state$chain
    timing <- ck$state$timing
    W_keep <- ck$state$W_keep
    converged <- isTRUE(ck$state$converged)
    assign(".Random.seed", ck$rng_state, envir = .GlobalEnv)
    message("Resuming ", checkpoint_file, " after iteration ", completed)
  }

  save_checkpoint <- function(i, converged) {
    if (is.null(checkpoint_file)) return(invisible(NULL))
    dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
    payload <- list(
      spec = checkpoint_spec, completed = as.integer(i),
      state = list(feat = feat, f_par = f_par,
                   sigma_inv_est = as.matrix(sigma_inv_est), d_cur = d_cur,
                   u_sa = u_sa, prev_W = prev_W, chain = chain,
                   timing = timing, W_keep = W_keep,
                   converged = converged),
      rng_state = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    tmp <- paste0(checkpoint_file, ".tmp-", Sys.getpid())
    saveRDS(payload, tmp, compress = FALSE)
    if (!file.rename(tmp, checkpoint_file)) {
      unlink(tmp)
      stop("could not atomically replace checkpoint: ", checkpoint_file,
           call. = FALSE)
    }
    invisible(NULL)
  }

  # --- constants of the feature draw -----------------------------------------
  Xmat      <- do.call(rbind, feat$X)
  gvals     <- vapply(feat$X, g$d, numeric(1))
  log_gvals <- log(gvals)
  r2        <- rowSums(Xmat^2)           # squared frequency norms, read by CM3

  ff <- rff(locs, feat$X, feat$phi, P, need_A = force_general)
  A  <- if (force_general) ff$A[obs_idx, , drop = FALSE] else as_dgc(matrix(0, 0, 0))
  A_l_list <- if (force_general) split_col_blocks(A, L, P) else vector("list", L)
  # dense n x L cosine matrix; empty selects the sparse reference branch
  c_dense <- if (force_general) matrix(numeric(0), 0, 0) else as.matrix(ff$c_vals)
  # A_l^T D^-1 A_l is diagonal: ||c_l||^2 / d_p when colocated, and the column
  # sums of squares over each process's observed locations when misaligned
  col_ss <- if (!force_general && use_fast) Matrix::colSums(ff$c_vals^2)
            else numeric(0)
  col_ss_mat <- if (masked) masked_col_ss(c_dense, observed)
                else matrix(numeric(0), 0, 0)
  rm(ff)

  # M = sqrt(g / 2f) from S_mat = sqrt(f) (or sqrt(R) under the ceiling), and
  # the per-feature matrices the Gibbs sampler reads
  scaling <- function(S_mat) {
    build_M(feat$X, f_den, f_par, g, f_vec = f_vec, gvals = gvals,
            s_mat = S_mat, guard = numerical_guard)
  }
  sampler_matrices <- function(M, sigma_inv, d_cur) {
    update_matrices_cpp(M, as_dgc(sigma_inv), d_cur[proc_id], A_l_list,
                        col_ss, col_ss_mat, d_cur)
  }

  S_mat   <- build_S_mat(f_vec, f_par, fam, Xmat, gvals,
                         numerical_guard, r2, log_gvals)
  M       <- scaling(S_mat)
  updated <- sampler_matrices(M, sigma_inv_est, d_cur)

  # A W for every retained draw, as an n_obs x T matrix. W is (P L) x T with
  # the process index fastest, so per draw A W = vec(Wmat c^T); all draws are
  # done in one (n x L) (L x P T) product.
  latent_fit <- function(W) {
    if (force_general) return(as.matrix(A %*% W))
    n_draw <- ncol(W)
    rhs  <- aperm(array(W, c(P, L, n_draw)), c(2, 1, 3))
    cw   <- c_dense %*% matrix(rhs, L, P * n_draw)
    full <- matrix(aperm(array(cw, c(n, P, n_draw)), c(2, 1, 3)), n * P, n_draw)
    if (use_fast) full else full[obs_idx, , drop = FALSE]
  }

  iteration_ids <- if (!converged && completed < run_until)
    seq.int(completed + 1L, run_until) else integer(0)
  for (i in iteration_ids) {
    if (verbose) message("  EM iteration ", i)
    g_k <- .saem_gain(i, lag)
    timing[i, ] <- 0

    # --- E step: MC Gibbs sweeps from the previous iteration's final state ---
    W <- matrix(0, nrow = P * L, ncol = MC)
    if (i != 1) W[, 1] <- prev_W
    W <- gibbs_W_cpp(W, updated, d_cur[proc_id], as.integer(MC),
                     as.integer(L), Y_gibbs, A, A_l_list, c_dense,
                     obs_mask, col_ss_mat, d_cur)$W_vals_sample
    prev_W <- W[, MC]
    W <- W[, (BurnIn + 1):MC, drop = FALSE]
    timing[i, "estep"] <- tock()
    if (keep_W == "last") W_keep <- W else if (keep_W == "all") W_keep[[i]] <- W

    # --- stochastic approximation of u_l = E[W_l W_l^T | Y] ------------------
    # in the raw W scale, so that both CM1 and CM3 read the same statistic
    u_hat <- feature_moments_cpp(W, as.integer(L))
    u_sa  <- if (i == 1) u_hat else (1 - g_k) * u_sa + g_k * u_hat
    timing[i, "moments"] <- tock()

    capped_now <- integer(P)
    for (cycle in seq_len(cm_cycles)) {
      # --- CM1: graphical LASSO on S = L^-1 sum_l (g_l/2) h^-1 u_l h^-1 -------
      s  <- suffstat_to_s(u_sa, M, tri)
      gl <- glasso::glasso(s, rho = rho, nobs = L * n_keep,
                           penalize.diagonal = diag_pen)
      sigma_inv_est <- Matrix::symmpart(gl$wi)
      timing[i, "cm1"] <- timing[i, "cm1"] + tock()

      # --- CM2: nugget, from the residuals of this E step's draws ------------
      # It depends only on W, so it runs once per iteration.
      if (cycle == 1L) {
        AW <- latent_fit(W)
        resid2 <- (Y_obs - AW)^2
        v <- if (nugget == "by_process") {
          vapply(seq_len(P), function(p) mean(resid2[proc_id == p, ]), numeric(1))
        } else {
          rep(mean(resid2), P)
        }
        d_cur <- if (i == 1) v else (1 - g_k) * d_cur + g_k * v
        chain$nugget[, i] <- d_cur
        timing[i, "cm2"] <- tock()
      }

      # --- CM3: spatial parameters, one process at a time --------------------
      sw <- cm3_sweep(f_par, fam, f_vec, Xmat, gvals,
                      as.matrix(sigma_inv_est), u_sa, tri, bounds,
                      S_mat, verbose = verbose, r2 = r2,
                      gradient = cm3_gradient, nthread = threads,
                      guard = numerical_guard, log_gvals = log_gvals)
      f_par <- sw$f_par
      S_mat <- sw$S_mat
      capped_now <- pmax(capped_now, sw$capped)
      timing[i, "cm3"] <- timing[i, "cm3"] + tock()

      M <- scaling(S_mat)
      updated <- sampler_matrices(M, sigma_inv_est, d_cur)
      timing[i, "refresh"] <- timing[i, "refresh"] + tock()
    }

    chain$capped[i, ] <- capped_now
    if (capping && !warned_capped && any(capped_now > 0L)) {
      warned_capped <- TRUE
      warning(sum(capped_now), " of ", L * P, " importance ratios hit the ",
              "ceiling at iteration ", i, " (processes ",
              paste(which(capped_now > 0L), collapse = ", "),
              "); see chain$capped.", call. = FALSE)
    }
    chain$spectral[[i]] <- f_par
    chain$sigma_inv[, i] <-
      as.matrix(sigma_inv_est)[upper.tri(diag(P), diag = TRUE)]

    # Complete-data penalized criterion at this iteration's (sigma, d, iota),
    # up to constants. The W part is glasso's objective
    # (L/2)[log|Q| - tr(SQ) - rho ||Q||_1], with S re-formed at the updated
    # iota, plus -sum_l log|h(X_l)|; the Y part is the Gaussian log-density of
    # the residuals, averaged over draws.
    Q_now  <- as.matrix(sigma_inv_est)
    s_post <- suffstat_to_s(u_sa, M, tri)
    pen <- if (diag_pen) sum(abs(Q_now))
           else sum(abs(Q_now[row(Q_now) != col(Q_now)]))
    chain$loglik[i, "prec"] <-
      (L / 2) * (as.numeric(determinant(Q_now, logarithm = TRUE)$modulus) -
                 sum(s_post * Q_now) - rho * pen) +
      q_iota_logdet(S_mat, if (capping) log_gvals else NULL)
    chain$loglik[i, "nugget"] <-
      sum(stats::dnorm(rep(Y_obs, n_keep),
                       mean = as.vector(AW),
                       sd = sqrt(rep(d_cur[proc_id], n_keep)),
                       log = TRUE)) / n_keep

    timing[i, "refresh"] <- timing[i, "refresh"] + tock()
    if (verbose) {
      message("    ", paste(sprintf("%s %.1fs", colnames(timing), timing[i, ]),
                            collapse = "  "))
    }

    # parameter summaries, by parameter name across the processes carrying it
    if (trace_every > 0L && (i %% trace_every == 0L || i == 1L)) {
      nms <- unique(unlist(lapply(f_par, names)))
      summ <- vapply(nms, function(nm) {
        v <- unlist(lapply(f_par, function(z)
          if (nm %in% names(z)) unname(z[[nm]]) else NULL))
        sprintf("%s %.3f [%.3f, %.3f]", nm, mean(v), min(v), max(v))
      }, character(1))
      message(sprintf(
        "  [iter %d] %s | nugget %.4f | diag(Q) %.3f | crit %.3f",
        i, paste(summ, collapse = "  "), mean(d_cur), mean(diag(Q_now)),
        sum(chain$loglik[i, ])))
    }

    if (i > 1 && abs(sum(chain$loglik[i, ]) - sum(chain$loglik[i - 1, ])) < tol)
      converged <- TRUE
    if (!is.null(checkpoint_file) &&
        (i %% checkpoint_every == 0L || i == run_until || converged))
      save_checkpoint(i, converged)
    if (converged) break
  }

  n_done <- if (length(iteration_ids)) i else completed
  if (!is.null(checkpoint_file) && n_done > 0L &&
      (n_done %% checkpoint_every != 0L || n_done == completed))
    save_checkpoint(n_done, converged)
  if (n_done < 1L)
    stop("fit stopped before completing an iteration", call. = FALSE)
  keep <- seq_len(n_done)
  structure(list(
    sigma_inv = as.matrix(sigma_inv_est),
    nugget    = if (nugget == "shared") d_cur[[1]] else d_cur,
    f_par     = f_par,
    # the graphical LASSO statistic at the final spectral parameters
    s         = suffstat_to_s(u_sa, M, tri),
    chain     = list(
      sigma_inv = chain$sigma_inv[, keep, drop = FALSE],
      nugget    = chain$nugget[, keep, drop = FALSE],
      spectral  = chain$spectral[keep],
      loglik    = chain$loglik[keep, , drop = FALSE],
      capped    = chain$capped[keep, , drop = FALSE],
      timing    = timing[keep, , drop = FALSE]  # seconds per stage
    ),
    W          = if (keep_W == "all") W_keep[keep] else W_keep,
    iterations = n_done,
    converged  = converged,
    call       = match.call(),
    spec = list(f_fun = f_fun, family = vapply(fam, function(e) e$name, character(1)),
                g = g, P = P, d = d,
                L = L, MC = MC, lag = lag, cm_cycles = cm_cycles,
                cm3_gradient = cm3_gradient,
                rho = rho, nugget = nugget,
                diag_pen = diag_pen, cm3 = "coordinate_ecm",
                numerical_guard = numerical_guard,
                lower = lower, upper = upper, seed = seed,
                threads = threads),
    data = list(Y = Y, locs = locs)
  ), class = "ssggp_fit")
}

#' @export
print.ssggp_fit <- function(x, ...) {
  cat("<ssggp_fit>\n")
  cat("  processes      :", x$spec$P, " locations:", nrow(x$data$Y), "\n")
  cat("  features (L)   :", x$spec$L, " MC:", x$spec$MC, "\n")
  cat("  CM cycles      :",
      if (is.null(x$spec$cm_cycles)) 1L else x$spec$cm_cycles,
      "per E-step\n")
  cat("  iterations     :", x$iterations,
      if (x$converged) "(converged)" else "(max reached)", "\n")
  cat("  nugget         :", x$spec$nugget, "=",
      format(x$nugget, digits = 4), "\n")
  cat("  penalty (rho)  :", x$spec$rho, "\n")
  cp <- x$chain$capped
  if (!is.null(cp) && any(cp > 0L)) {
    L <- x$spec$L; P <- x$spec$P
    cat(sprintf("  ratio ceiling   : hit on %d of %d iterations, max %.4f%% of ",
                sum(rowSums(cp) > 0L), nrow(cp), 100 * max(rowSums(cp)) / (L * P)))
    cat(sprintf("entries\n                    processes affected: %s\n",
                paste(which(colSums(cp) > 0L), collapse = ", ")))
  }
  cat("  nonzero off-diagonals in sigma_inv:",
      sum(x$sigma_inv[upper.tri(x$sigma_inv)] != 0), "of",
      choose(x$spec$P, 2), "\n")

  nms <- unique(unlist(lapply(x$f_par, names)))
  if (length(nms)) {
    cat("  spectral parameters:\n")
    for (nm in nms) {
      v <- vapply(x$f_par, function(z)
        if (nm %in% names(z)) z[[nm]] else NA_real_, numeric(1))
      v <- v[!is.na(v)]
      if (length(v) <= 6L) {
        cat(sprintf("    %-6s %s\n", nm,
                    paste(format(v, digits = 4), collapse = " ")))
      } else {
        q <- stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
        cat(sprintf("    %-6s n=%d  min %s  q1 %s  med %s  q3 %s  max %s\n",
                    nm, length(v), format(q[1], digits = 4),
                    format(q[2], digits = 4), format(q[3], digits = 4),
                    format(q[4], digits = 4), format(q[5], digits = 4)))
      }
    }
  }
  invisible(x)
}
