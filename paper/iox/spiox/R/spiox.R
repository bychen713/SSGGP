#' Fit GP-IOX models for multivariate spatial outcomes
#'
#' @description
#' Fits Gaussian Process Independent Output (GP-IOX) models for multivariate spatial data. 
#' The function supports two primary model formulations (Response or Latent) and two 
#' inference frameworks (Markov Chain Monte Carlo or Variational Inference). 
#' 
#' Internally, the function is highly automated: it detects regular spatial grids to 
#' optimize Vecchia Directed Acyclic Graph (DAG) construction, manages OpenMP and BLAS 
#' threading to prevent thread contention, and handles missing data appropriately 
#' based on the chosen inference method.
#'
#' @details
#' \strong{Key Capabilities:}
#' \itemize{
#'   \item \strong{Model Formulations:} Choose `method = "response"` for direct spatial modeling of the outcomes, or `method = "latent"` for spatial modeling of underlying latent factors.
#'   \item \strong{Inference Engines:} Choose `fit = "mcmc"` for full Bayesian inference via Markov Chain Monte Carlo, or `fit = "vi"` for Variational Inference.
#'   \item \strong{Missing Data Handling:} 
#'     \itemize{
#'       \item For `response` + `mcmc`: Automatically imputes missing values in `Y` during sampling and returns the posterior samples for missing indices.
#'       \item For `response` + `vi`: Automatically filters `Y`, `X`, and `coords` to complete cases (dropping rows with NAs) and issues a warning.
#'       \item For `latent` + `mcmc`: Automatically samples the latent effects at all locations. Returns all samples (including at NA locations).
#'       \item For `latent` + `vi`: Automatically finds the posterior mean at all locations. Optionally returns posterior predictive samples for the latent effects by using the `opts$vi_pred_smp` argument.
#'     }
#'   \item \strong{Grid Detection & DAG Optimization:} The function automatically checks if `coords` fall on a regular grid (within a machine-precision tolerance). If true, it optimizes the Vecchia DAG construction.
#'   \item \strong{Thread Management:} Safely sets BLAS threads to 1 when OpenMP multi-threading (`opts$num_threads`) is engaged to prevent explosive CPU contention.
#' }
#'
#' @param Y Numeric matrix of responses with dimensions `n x q` (rows = locations, columns = outcomes). Can contain `NA` values (behavior depends on `method` and `fit`).
#' @param X Numeric design matrix with dimensions `n x p` (must have the same number of rows as `Y`).
#' @param coords Numeric matrix of spatial coordinates with dimensions `n x 2`.
#' @param m Integer. Vecchia neighbor size used to construct the DAG. Defaults to 15.
#' @param method Character. The model formulation: `"response"` or `"latent"`.
#' @param fit Character. The inference method: `"mcmc"` or `"vi"`.
#' @param iter Integer. Number of MCMC iterations if `fit="mcmc"`. If `fit="vi"`, this serves as the maximum number of iterations allowed. Defaults to 1000.
#' @param print_every Integer. Print progress every `print_every` iterations. Set to `0` for silent execution.
#' @param starting Character (`"auto"`) or named list of starting values. If `"auto"`, parameters are initialized internally using `autostart()`. If a list, omitted elements are initialized to standard defaults. 
#'   \itemize{
#'     \item `Beta`: numeric matrix `p x q` (starting regression coefficients, defaults to all `0`).
#'     \item `Sigma`: numeric matrix `q x q` (starting outcome covariance / coregionalization, defaults to `diag(q)`).
#'     \item `Theta`: numeric matrix `3 x q` with rows corresponding to spatial parameters `phi`, `nu`, `alpha`. Must be supplied when `fit="vi"` or when MCMC updates for Theta are disabled.
#'     \item `W`: numeric matrix `n x q` (starting latent factors, required/used only for `method="latent"`).
#'     \item `Ddiag`: numeric vector of length `q` (diagonal noise variance, used only for `method="latent"`; defaults to `rep(1, q)`).
#'   }
#' @param opts Optional named list of algorithmic and model options:
#'   \itemize{
#'     \item `matern`: integer (0 = squared exponential where `nu` acts as the exponent; 1 = Matern where `nu` is smoothness). Defaults to 1.
#'     \item `nu`: numeric. Baseline smoothness parameter used for `"auto"` starting values. Defaults to 0.5.
#'     \item `num_threads`: integer. Number of OpenMP threads to use. Automatically defaults to available cores via `RhpcBLASctl`.
#'     \item `update_Theta`: integer vector of length 3 indicating which of (`phi`, `nu`, `alpha`) to update via MCMC (1 = update, 0 = fixed). Defaults to `c(0,0,0)` for latent, `c(1,0,1)` for response. \emph{Note: `alpha` is forced to 0 for `method="latent"`}.
#'     \item `tol`: numeric. Convergence tolerance for Variational Inference. Defaults to `1e-2`.
#'     \item `vi_pred_smp`: integer. Number of predictive samples to draw when using VI. Defaults to 0.
#'     \item `cg_preconditioner`: character or integer. CG preconditioner for
#'       the joint BW block sampler (`method = "latent"`, `fit = "mcmc"`,
#'       `debug$sampling = 1L`). Choices:
#'       \itemize{
#'         \item `"auto"` (default): probe POSTERIOR for the first 5 Gibbs sweeps,
#'           then JACOBI for the next 5 with a budget cap of `2 · max(POSTERIOR
#'           CG iters)`. Lock in whichever converged in fewer CG iterations on
#'           average (POSTERIOR by default if JACOBI hits the cap without
#'           converging on any probe sweep).
#'         \item `"jacobi"`: diagonal of the joint precision operator. Cheap per
#'           apply; usually weak.
#'         \item `"posterior"`: block-diagonal-on-(B,W) preconditioner. Exact
#'           dense Cholesky on the B block (per-outcome `p×p`), Σ-mixed
#'           Vecchia-precision factors on the W block (per-outcome `H_A,j`
#'           sandwiched with the correlation matrix of Σ). Built once per
#'           Gibbs sweep; tends to use far fewer CG iters than JACOBI when Σ
#'           has substantial off-diagonal cross-outcome structure.
#'       }
#'       Ignored for `debug$sampling != 1L` and for `method = "response"`.
#'   }
#' @param debug Optional named list of MCMC controls, primarily for targeted development and sampling restriction:
#'   \itemize{
#'     \item `sampling`: integer. Specifies the latent-model sampler variant. Used only when `method="latent"` and `fit="mcmc"`. Options are `1L` for full block, `2L` for single-site sequential, `3L` for single-outcome sequential.
#'     \item `sample_Beta`: logical. Whether to sample regression coefficients (`Beta`). Defaults to `TRUE`.
#'     \item `sample_Sigma`: logical. Whether to sample `Sigma`. Defaults to `TRUE`.
#'     \item `sample_Ddiag`: logical. Whether to sample `Ddiag` in latent-model MCMC. Defaults to `TRUE`.
#'   }
#'
#' @return A comprehensively formatted list containing:
#' \itemize{
#'   \item \strong{Model Parameters}: Fitted matrices/arrays for `Theta`, `Beta`, `Sigma`. 
#'   \item \strong{Latent States}: If `method="latent"`, includes `W` (reordered to match the original input `Y` coordinates) and `Ddiag`.
#'   \item \strong{Imputed Data}: If `method="response"`, `fit="mcmc"`, and `Y` contained `NA`s, returns `Y_missing_indices`, `Y_missing_col`, and `Y_missing_samples` mapped back to the original spatial ordering.
#'   \item \strong{Metadata}: `call`, `method`, `fit`, `gridded` (boolean indicating if a regular grid was detected), `dag_opts`, `dag_info`, and original `coords`.
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Fit a GP-IOX response model using Variational Inference (VI)
#' fit_vi <- spiox(
#'   Y = Y, 
#'   X = X, 
#'   coords = coords, 
#'   m = 15, 
#'   method = "response", 
#'   fit = "vi", 
#'   iter = 200, 
#'   print_every = 100,
#'   opts = list(num_threads = 4)
#' )
#' 
#' # Example 2: Fit a GP-IOX latent model using MCMC (with automatic imputation if Y has NAs)
#' fit_mcmc <- spiox(
#'   Y = Y, 
#'   X = X, 
#'   coords = coords, 
#'   method = "latent", 
#'   fit = "mcmc", 
#'   iter = 1000, 
#'   starting = "auto"
#' )
#' }
#' @export
spiox <- function(Y, X, coords, m = 15, 
                  method = c("response", "latent"),
                  fit = c("mcmc", "vi"),
                  iter = 1000, print_every = 100, 
                  starting = "auto", opts = NULL, debug = NULL) {
  
  # ---------------------------------------------------------------------------
  # 1. Argument Matching & Data Validation
  # ---------------------------------------------------------------------------
  method <- match.arg(method)
  fit    <- match.arg(fit)
  
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  coords <- as.matrix(coords)
  
  if( (fit == "vi") & (method == "response") ){
    cc <- complete.cases(Y)
    
    Y <- Y[cc,,drop=FALSE]
    X <- X[cc,,drop=FALSE]
    coords <- coords[cc,,drop=FALSE]
    
    nava <- sum(cc)
    message(paste0("Fitting on ", nava, " points (others are partly missing)"))
  }
  
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  d <- ncol(coords)
  
  stopifnot(
    "X must have same number of rows as Y" = nrow(X) == n, 
    "coords must have same number of rows as Y" = nrow(coords) == n, 
    "coords must be 2-dimensional" = d == 2L
  )
  
  # ---------------------------------------------------------------------------
  # 2. Spatial Grid Evaluation & Vecchia DAG Construction
  # ---------------------------------------------------------------------------
  # Helper to check if coordinates fall on a regular grid
  is_gridded <- function(xy, tol = sqrt(.Machine$double.eps)) {
    x <- sort(unique(xy[, 1]))
    y <- sort(unique(xy[, 2]))
    dx <- diff(x)
    dy <- diff(y)
    
    ok_step <- function(d) {
      length(d) <= 1L || 
        max(d, na.rm = TRUE) <= tol || 
        (max(d) - min(d) <= tol * max(1, max(d)))
    }
    
    nx <- length(x)
    ny <- length(y)
    
    (nrow(unique(xy)) == nx * ny) && ok_step(dx) && ok_step(dy)
  }
  
  gridded <- is_gridded(coords)
  dag <- dag_vecchia_o(coords, m, gridded)
  dag_opts <- if (gridded) -1L else 0L
  
  # ---------------------------------------------------------------------------
  # 3. Model & Debug Options Setup
  # ---------------------------------------------------------------------------
  # Set method-specific defaults for Theta updates
  default_update_theta <- if (method == "latent") c(0L, 0L, 0L) else c(1L, 0L, 1L)
  
  opts_defaults <- list(
    update_Theta      = default_update_theta,
    num_threads       = RhpcBLASctl::get_num_cores(),
    tol               = 1e-2,
    matern            = 1,
    nu                = 0.5,
    vi_pred_smp       = 0,
    # CG preconditioner choice for the joint BW block sampler (sampling = 1).
    # One of: "auto" (probe POSTERIOR vs JACOBI), "jacobi", "posterior".
    # Ignored for other samplers and for method = "response".
    cg_preconditioner = "auto"
  )
  opts <- modifyList(opts_defaults, if (is.null(opts)) list() else opts)
  opts$vi_pred_smp <- as.integer(opts$vi_pred_smp)
  stopifnot("opts$update_Theta must be length 3" = length(opts$update_Theta) == 3L)

  # Translate opts$cg_preconditioner from string to integer code expected by
  # the C++ side.  Codes (sampling = 1, the joint BW block sampler):
  #   0 = auto       (probe POSTERIOR for 5 sweeps, then JACOBI for 5 with a
  #                   2·max(POSTERIOR iters) budget cap, lock in the winner)
  #   1 = jacobi     (diagonal of the joint precision operator)
  #   2 = posterior  (block-diagonal-on-(B,W) PC with exact dense B half and
  #                   Σ-mixed Vecchia-precision W half)
  cg_pc_codes <- c(auto = 0L, jacobi = 1L, posterior = 2L)
  cg_pc_key <- if (is.numeric(opts$cg_preconditioner)) {
    as.integer(opts$cg_preconditioner)
  } else {
    code <- cg_pc_codes[tolower(as.character(opts$cg_preconditioner))]
    if (is.na(code)) {
      stop("Invalid opts$cg_preconditioner: '",
           opts$cg_preconditioner,
           "'. Use one of: ", paste(names(cg_pc_codes), collapse = ", "), ".")
    }
    code
  }
  opts$cg_preconditioner_int <- cg_pc_key


  # Thread management
  if (opts$num_threads > 1) {
    RhpcBLASctl::blas_set_num_threads(1)
    message(paste0("BLAS threads set to 1, OMP threads set to ", opts$num_threads))
  } else {
    message("OMP threads set to 1.")
  }
  
  # Debug controls (primarily for MCMC development)
  debug_defaults <- list(
    sampling     = 1L,  # 1 = full block
    sample_Beta  = TRUE,
    sample_Sigma = TRUE,
    sample_Ddiag = TRUE
  )
  debug <- modifyList(debug_defaults, if (is.null(debug)) list() else debug)
  
  # ---------------------------------------------------------------------------
  # 4. Starting Values Initialization
  # ---------------------------------------------------------------------------
  auto_do <- is.character(starting) && length(starting) == 1 && starting == "auto"
  
  if (!auto_do && !is.list(starting)) {
    stop("Invalid input for 'starting'. Specify 'auto' or provide a list.")
  }
  
  if (auto_do) {
    # Generate automatic starting values
    starting <- autostart(Y, X, coords, method, m = m, opts$nu)
    Beta_start  <- starting$Beta
    Sigma_start <- starting$Sigma
    Theta       <- starting$Theta
    
    if (method == "latent") {
      W_start     <- starting$W
      Ddiag_start <- starting$Ddiag
    }
  } else {
    # Process user-provided starting values
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    
    # Common variables
    Beta_start  <- starting$Beta  %||% matrix(0, nrow = p, ncol = q)
    Sigma_start <- starting$Sigma %||% diag(q)
    
    # Latent-only variables
    if (method == "latent") {
      W_start     <- starting$W     %||% matrix(0, nrow = n, ncol = q)
      Ddiag_start <- starting$Ddiag %||% diag(q)
    } 
    
    # Validate Theta based on user's fit/update options
    if (all(opts$update_Theta == 0L) && is.null(starting$Theta)) {
      stop("starting$Theta must be supplied when update_Theta = c(0,0,0).")
    }
    if (fit == "vi" && is.null(starting$Theta)) {
      stop("starting$Theta (3 x q with rows phi, nu, alpha) must be supplied for fit='vi'.")
    }
    
    Theta <- starting$Theta %||% NULL
    
    if (!is.null(Theta)) {
      Theta <- as.matrix(Theta)
      stopifnot("Theta must be 3 x q" = nrow(Theta) == 3L, ncol(Theta) == q)
    } else {
      # Fallback defaults if Theta is partially updated
      alpha_default <- ifelse(method == "latent", 0, 0.1)
      Theta <- rbind(10, 0.5, alpha_default)
    }
    rownames(Theta) <- c("phi", "nu", "alpha")
    
    # Alert user if they passed invalid alpha for latent method
    if (method == "latent" && !is.null(Theta)) {
      alpha <- as.numeric(Theta[3, ])
      if (any(abs(alpha) > 0)) {
        warning("In method='latent', alpha (last row of starting$Theta) should be 0; nonzero values were supplied.")
      }
    }
  }
  
  # Construct the full Theta start matrix (adding placeholder sigmasq)
  Theta_start <- rbind(
    phi     = Theta[1, , drop = FALSE],
    sigmasq = rep(1, q),                 # Internal tracking, discarded later
    nu      = Theta[2, , drop = FALSE],
    alpha   = Theta[3, , drop = FALSE]
  )
  
  # Format the update bools (phi, sigmasq (fixed 0), nu, alpha)
  update_Theta_full <- as.integer(
    c(opts$update_Theta[1], 0L, opts$update_Theta[2], opts$update_Theta[3])
  )
  if(method == "latent"){
    if(opts$update_Theta[3] == 1){
      update_Theta_full[4] <- 0
      warning("MCMC for alpha disabled because method = 'latent' (opts$update_Theta[3] forced to 0)")
    }
  }
  
  # ---------------------------------------------------------------------------
  # 5. Model Dispatch
  # ---------------------------------------------------------------------------
  out <- switch(
    paste(method, fit, sep = ":"),
    
    "response:mcmc" = spiox_response(
      Y[dag$order,,drop=F], X[dag$order,,drop=F], coords[dag$order,,drop=F], 
      dag$dag,
      Beta_start, Sigma_start, Theta_start,
      mcmc         = iter, 
      print_every  = print_every, 
      matern       = opts$matern, 
      dag_opts     = dag_opts,
      sample_Beta  = debug$sample_Beta, 
      sample_Sigma = debug$sample_Sigma,
      update_Theta = update_Theta_full,
      num_threads  = as.integer(opts$num_threads)
    ),
    
    "latent:mcmc" = spiox_latent(
      Y[dag$order,,drop=F], X[dag$order,,drop=F], coords[dag$order,,drop=F],
      dag$dag,
      Beta_start,
      W_start[dag$order,,drop=F],
      Sigma_start, Theta_start, Ddiag_start,
      mcmc              = iter,
      print_every       = print_every,
      matern            = opts$matern,
      dag_opts          = dag_opts,
      sample_Beta       = debug$sample_Beta,
      sample_Sigma      = debug$sample_Sigma,
      sample_Ddiag      = debug$sample_Ddiag,
      update_Theta      = update_Theta_full,
      num_threads       = as.integer(opts$num_threads),
      sampling          = as.integer(debug$sampling),
      cg_preconditioner = opts$cg_preconditioner_int
    ),
    
    "response:vi" = spiox_response_vi(
      Y[dag$order,,drop=F], X[dag$order,,drop=F], coords[dag$order,,drop=F],  
      dag$dag, dag_opts,
      Theta        = Theta_start, 
      Sigma_start  = Sigma_start, 
      Beta_start   = Beta_start,
      print_every  = print_every, 
      matern       = opts$matern,
      num_threads  = as.integer(opts$num_threads)
    ),
    
    "latent:vi" = spiox_latent_vi(
      Y[dag$order,,drop=F], X[dag$order,,drop=F], coords[dag$order,,drop=F], 
      dag$dag, dag_opts,
      Theta        = Theta_start, 
      Sigma_start  = Sigma_start, 
      Beta_start   = Beta_start,
      W_start      = W_start[dag$order,,drop=F],
      Ddiag_start  = Ddiag_start,
      matern       = opts$matern, 
      num_threads  = as.integer(opts$num_threads),
      print_every  = print_every, 
      tol          = opts$tol, 
      max_iter     = iter,
      vi_pred_smp = opts$vi_pred_smp
    ),
    
    stop("Unknown method/fit combination.")
  )
  
  # ---------------------------------------------------------------------------
  # 6. Format and Return Output
  # ---------------------------------------------------------------------------
  # Remove the internal 'sigmasq' placeholder from Theta for consistency
  if ("Theta" %in% names(out)) {
    out$Theta <- out$Theta[-2, , , drop = FALSE]
  } else {
    out$Theta <- Theta_start[-2, , drop = FALSE]
  }
  
  # reorder W
  if ("W" %in% names(out)) {
    if(fit == "mcmc"){
      # an array with samples in the third dimension
      out$W <- out$W[order(dag$order),,,drop=FALSE]  
    } else {
      # not an array
      out$W <- out$W[order(dag$order),,drop=FALSE]
    }
  } 
  
  # indices of imputed, response:mcmc
  if ( (method == "response") & (fit == "mcmc") & any(is.na(Y)) ) {
    ord <- dag$order
    cur <- which(is.na(Y[ord, ]))
    tgt <- as.numeric(which(is.na(Y)))
    
    # map current indices to original-space indices
    n <- nrow(Y)
    row_reord <- ((cur - 1) %% n) + 1
    col_idx   <- ((cur - 1) %/% n) + 1
    row_orig  <- ord[row_reord]
    orig_idx  <- row_orig + (col_idx - 1) * n
    
    # orig_idx is a permutation of tgt; find the reordering
    reorder <- match(tgt, orig_idx)
    
    out$Y_missing_indices <- tgt
    out$Y_missing_samples <- out$Y_missing_samples[reorder, , drop=FALSE]
    out$Y_missing_col <- ((tgt - 1) %/% n) + 1
  }
  
  # Append metadata to the output object
  out$call     <- match.call()
  out$method   <- method
  out$fit      <- fit
  out$gridded  <- gridded
  out$dag_opts <- dag_opts
  out$dag_info <- dag
  out$coords   <- coords
  
  return(out)
}