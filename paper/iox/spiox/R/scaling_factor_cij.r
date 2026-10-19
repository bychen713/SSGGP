scaling_factor_at_zero <- function(S, thetamat, m=20, nmax=1000, matern=TRUE, n_threads=1){
  nr <- nrow(S)
  
  subs <- round(seq(1, nr, length.out=min(c(nmax, nr))))
  S_sub <- S[subs,]
  daginfo <- spiox:::dag_vecchia_o(S_sub, m)
  
  R <- sfact(daginfo$dag, S_sub[daginfo$order,], thetamat, matern, n_threads)
  C <- R
  C[upper.tri(C)] <- t(R)[upper.tri(R)]
  
  return(
    C
  )
}

