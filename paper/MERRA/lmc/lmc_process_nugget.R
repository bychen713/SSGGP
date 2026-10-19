## Exact LMC likelihood and prediction with a diagonal, process-specific
## observation nugget. Kept separate from LMC/lmc_common.R so existing shared-
## nugget fits and simulation code remain unchanged.

lmc_pn_pack <- function(par, nugget) {
  stopifnot(all(is.finite(nugget)), all(nugget > 0))
  c(unlist(lapply(par, function(x) log(as.numeric(x)))), log(nugget))
}

lmc_pn_unpack <- function(theta, families, P) {
  np <- vapply(families, lmc_npar, integer(1))
  m <- sum(np)
  stopifnot(length(theta) == m + P)
  end <- cumsum(np)
  start <- end - np + 1L
  list(par = lapply(seq_along(families), function(k)
         exp(theta[start[k]:end[k]])),
       nugget = exp(theta[m + seq_len(P)]))
}

lmc_pn_bounds <- function(families, P) {
  list(lower = c(unlist(lapply(families, lmc_lower)), rep(log(1e-4), P)),
       upper = c(unlist(lapply(families, lmc_upper)), rep(log(10), P)))
}

# K = sum_k R_k (x) b_k b_k' + I_n (x) diag(tau_1^2,...,tau_P^2).
# For missing responses, restrict K to observed entries. Woodbury reduces the
# Cholesky from n_obs x n_obs to nq x nq; inverse-nugget weights enter both the
# inner blocks and the quadratic term. This also handles colocated data, where
# the shared-nugget rotation is no longer valid after tau_i^2 varies with i.
lmc_pn_nll_factory <- function(Y_mat, D, B, families) {
  n <- nrow(Y_mat); P <- ncol(Y_mat); q <- ncol(B)
  stopifnot(length(families) == q, identical(dim(D), c(n, n)),
            all(is.finite(B)), all(colSums(!is.na(Y_mat)) > 0))
  O <- !is.na(Y_mat)
  Y0 <- Y_mat; Y0[!O] <- 0
  n_obs <- sum(O)
  count <- colSums(O)
  pair <- which(upper.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
  pair_index <- matrix(NA_integer_, q, q)
  for (p in seq_len(nrow(pair)))
    pair_index[pair[p, 1L], pair[p, 2L]] <-
      pair_index[pair[p, 2L], pair[p, 1L]] <- p
  bprod <- vapply(seq_len(nrow(pair)), function(p)
    B[, pair[p, 1L]] * B[, pair[p, 2L]], numeric(P))
  idx <- lapply(seq_len(q), function(k) ((k - 1L) * n + 1L):(k * n))
  cache <- new.env(parent = emptyenv())
  cache$par <- vector("list", q)
  cache$L <- vector("list", q)
  cache$nugget <- NULL
  cache$blocks <- matrix(0, n * q, n * q)
  cache$ok <- FALSE

  function(theta) {
    up <- lmc_pn_unpack(theta, families, P)
    tau2 <- up$nugget
    if (any(!is.finite(tau2)) || any(tau2 <= 0)) return(1e12)
    nugget_changed <- !identical(tau2, cache$nugget)
    changed <- integer(0)
    for (k in seq_len(q)) {
      if (!cache$ok || !identical(up$par[[k]], cache$par[[k]])) {
        R <- lmc_cor_safe(families[k], up$par[[k]], D, symmetric = TRUE)
        if (is.null(R)) { cache$ok <- FALSE; return(1e12) }
        Lk <- tryCatch(lmc_cor_root(R), error = function(e) NULL)
        if (is.null(Lk)) { cache$ok <- FALSE; return(1e12) }
        cache$L[[k]] <- Lk
        cache$par[[k]] <- up$par[[k]]
        changed <- c(changed, k)
      }
    }
    weights <- sweep(O * 1, 2L, tau2, "/")
    A <- weights %*% bprod
    YB <- (Y0 * weights) %*% B
    if (nugget_changed || !cache$ok) {
      for (k in seq_len(q)) for (l in k:q) {
        blk <- crossprod(cache$L[[k]] * A[, pair_index[k, l]], cache$L[[l]])
        cache$blocks[idx[[k]], idx[[l]]] <- blk
        if (l != k) cache$blocks[idx[[l]], idx[[k]]] <- t(blk)
      }
    } else if (length(changed)) {
      for (k in changed) for (l in seq_len(q)) {
        blk <- crossprod(cache$L[[k]] * A[, pair_index[k, l]], cache$L[[l]])
        cache$blocks[idx[[k]], idx[[l]]] <- blk
        if (l != k) cache$blocks[idx[[l]], idx[[k]]] <- t(blk)
      }
    }
    cache$nugget <- tau2
    cache$ok <- TRUE

    M <- cache$blocks
    diag(M) <- diag(M) + 1
    cs <- tryCatch(chol(M), error = function(e) NULL)
    if (is.null(cs) || any(!is.finite(diag(cs))) || min(diag(cs)) <= 0)
      return(1e12)
    z <- unlist(lapply(seq_len(q), function(k)
      crossprod(cache$L[[k]], YB[, k])))
    v <- backsolve(cs, z, upper.tri = TRUE, transpose = TRUE)
    quad <- sum(Y0^2 * weights) - sum(v^2)
    logdet <- sum(count * log(tau2)) + 2 * sum(log(diag(cs)))
    value <- 0.5 * (logdet + quad + n_obs * log(2 * pi))
    if (is.finite(value)) value else 1e12
  }
}

lmc_pn_predict <- function(Y_mat, D, Dx, B, families, par, nugget) {
  n <- nrow(Y_mat); P <- ncol(Y_mat); q <- ncol(B); nt <- nrow(Dx)
  stopifnot(length(nugget) == P, all(nugget > 0), length(par) == q,
            length(families) == q)
  O <- !is.na(Y_mat)
  Y0 <- Y_mat; Y0[!O] <- 0
  weights <- sweep(O * 1, 2L, nugget, "/")
  pair <- which(upper.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
  pair_index <- matrix(NA_integer_, q, q)
  for (p in seq_len(nrow(pair)))
    pair_index[pair[p, 1L], pair[p, 2L]] <-
      pair_index[pair[p, 2L], pair[p, 1L]] <- p
  bprod <- vapply(seq_len(nrow(pair)), function(p)
    B[, pair[p, 1L]] * B[, pair[p, 2L]], numeric(P))
  A <- weights %*% bprod
  idx <- lapply(seq_len(q), function(k) ((k - 1L) * n + 1L):(k * n))
  L <- Rx <- vector("list", q)
  for (k in seq_len(q)) {
    L[[k]] <- lmc_cor_root(lmc_cor(families[k], par[[k]], D,
                                   symmetric = TRUE))
    Rx[[k]] <- lmc_cor(families[k], par[[k]], Dx)
  }
  M <- matrix(0, n * q, n * q)
  for (k in seq_len(q)) for (l in k:q) {
    blk <- crossprod(L[[k]] * A[, pair_index[k, l]], L[[l]])
    M[idx[[k]], idx[[l]]] <- blk
    if (l != k) M[idx[[l]], idx[[k]]] <- t(blk)
  }
  diag(M) <- diag(M) + 1
  cs <- chol(M)

  YB <- (Y0 * weights) %*% B
  z <- unlist(lapply(seq_len(q), function(k) crossprod(L[[k]], YB[, k])))
  vm <- backsolve(cs, backsolve(cs, z, upper.tri = TRUE, transpose = TRUE))
  G <- vapply(seq_len(q), function(k)
    as.vector(L[[k]] %*% vm[idx[[k]]]), numeric(n))
  Ky <- weights * (Y0 - G %*% t(B))
  U <- Ky %*% B
  Tm <- vapply(seq_len(q), function(k)
    as.vector(Rx[[k]] %*% U[, k]), numeric(nt))
  pred_mean <- Tm %*% t(B)

  Psi <- lapply(seq_len(q), function(k)
    lapply(seq_len(q), function(m)
      crossprod(L[[k]] * A[, pair_index[k, m]], t(Rx[[m]]))))
  Gq <- array(NA_real_, c(nt, q, q))
  for (k in seq_len(q)) for (l in k:q) {
    value <- as.vector((Rx[[k]] * Rx[[l]]) %*% A[, pair_index[k, l]])
    Gq[, k, l] <- value
    Gq[, l, k] <- value
  }
  Kprior <- tcrossprod(B) + diag(nugget, P)
  cov <- vector("list", nt)
  for (t in seq_len(nt)) {
    W <- matrix(0, n * q, P)
    for (k in seq_len(q)) {
      Pk <- vapply(seq_len(q), function(m) Psi[[k]][[m]][, t], numeric(n))
      W[idx[[k]], ] <- Pk %*% t(B)
    }
    v <- backsolve(cs, W, upper.tri = TRUE, transpose = TRUE)
    part1 <- B %*% Gq[t, , ] %*% t(B)
    Ct <- Kprior - part1 + crossprod(v)
    cov[[t]] <- (Ct + t(Ct)) / 2
  }
  list(mean = pred_mean, cov = cov)
}
