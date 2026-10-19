#' Simulate GP-IOX data 
#'
#' @description
#' Simulates GP-IOX data over a set of coordinates.
#'
#' @param coords A numeric matrix of spatial coordinates. 
#'   Must not contain missing values (`NA`).
#' @param Sigma A numeric, square covariance matrix of dimension `q x q`.
#' @param Theta A numeric matrix of spatial parameters of dimension `3 x q`. 
#'   The rows correspond to `phi`, `nu`, and `alpha`, respectively.
#' @param m A single positive integer specifying the number of nearest neighbors 
#'   to condition on for the Vecchia approximation. Default is `15`.
#' @param matern A single numeric value indicating the covariance 
#'   specification. Default is `1` for Matern.
#' @param silent Logical. If `FALSE`, prints progress messages to the console. 
#'   Default is `TRUE`.
#' @param num_threads Integer. Number of OpenMP threads to use for the simulation. 
#'   Defaults to the number of available cores via `RhpcBLASctl::get_num_cores()`.
#'
#' @return A numeric matrix of simulated values with dimensions `n x q`, ordered 
#'   to match the original input `coords`.
#'
#' @export
rgpiox <- function(coords, Sigma, Theta, m=15, matern=1, silent=TRUE, num_threads  = RhpcBLASctl::get_num_cores()){

  if (any(is.na(coords))) {
    stop("Input 'coords' must not contain any missing values (NA).")
  }
  coords <- as.matrix(coords)
  if (!is.numeric(coords)) {
    stop("Input 'coords' must be numeric.")
  }
  
  # Check 'Sigma'
  if (!is.matrix(Sigma) || !is.numeric(Sigma)) {
    stop("Input 'Sigma' must be a numeric matrix.")
  }
  q <- nrow(Sigma)
  if (ncol(Sigma) != q) {
    stop(sprintf("Input 'Sigma' must be a square matrix. Current dimensions: %d x %d.", q, ncol(Sigma)))
  }
  
  # Check 'Theta'
  if (!is.matrix(Theta) || !is.numeric(Theta)) {
    stop("Input 'Theta' must be a numeric matrix.")
  }
  if (nrow(Theta) != 3) {
    stop(sprintf("Input 'Theta' must have exactly 3 rows. Current rows: %d.", nrow(Theta)))
  }
  if (ncol(Theta) != q) {
    stop(sprintf("The number of columns in 'Theta' (%d) must match the dimensions of 'Sigma' (%d).", ncol(Theta), q))
  }
  
  # Check 'm'
  if (!is.numeric(m) || length(m) != 1 || m <= 0 || m %% 1 != 0) {
    stop("Input 'm' must be a single positive integer.")
  }
  
  # Check 'matern'
  if (!is.numeric(matern) || length(matern) != 1) {
    stop("Input 'matern' must be a single numeric value.")
  }
  
  if(!silent) print("Building DAG.")
  dag <- dag_vecchia_o(coords, m)
  
  # Construct the full Theta start matrix (adding placeholder sigmasq)
  Theta_internal <- rbind(
    phi     = Theta[1, , drop = FALSE],
    sigmasq = rep(1, q),                 # Internal tracking, discarded later
    nu      = Theta[2, , drop = FALSE],
    alpha   = Theta[3, , drop = FALSE]
  )
  
  if(!silent) print("Simulating...")
  W <- spiox_simulate(coords[dag$order,,drop=F], dag$dag, Sigma, Theta_internal, matern, num_threads)
  
  if(!silent) print("Done.")
  return(W[order(dag$order),,drop=FALSE])
}