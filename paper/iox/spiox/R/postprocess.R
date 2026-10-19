iox_fij <- function(iox_obj, i, j, 
                    n_threads=1, matern=1, tail_mcmc=NULL, max_dist=1000, max_range=1){
  
  if(!is.null(tail_mcmc)){
    iox_obj$Theta <- tail(iox_obj$Theta, c(NA, NA, tail_mcmc))
  }
  
  n <- nrow(iox_obj$coords)
  q <- dim(iox_obj$Theta)[2]
  mcmc <- dim(iox_obj$Theta)[3]
  
  c_o <- iox_obj$coords[iox_obj$dag_info$order, ]
  
  Theta4 <- array(1, dim=c(4, q, mcmc))
  Theta4[1,,] <- iox_obj$Theta[1,,]
  Theta4[3:4,,] <- iox_obj$Theta[2:3,,]
  
  result <- spiox:::iox_make_fij(i-1, j-1, c_o, 
                                iox_obj$dag_info$dag, 
                                iox_obj$dag_opts, 
                                Theta4, matern, n_threads, max_dist, max_range)

  return(result)
}
  
  
iox_fij0 <- function(iox_obj, 
                    n_threads=1, matern=1, tail_mcmc=NULL){
  
  if(!is.null(tail_mcmc)){
    iox_obj$Theta <- tail(iox_obj$Theta, c(NA, NA, tail_mcmc))
  }
  
  n <- nrow(iox_obj$coords)
  q <- dim(iox_obj$Theta)[2]
  mcmc <- dim(iox_obj$Theta)[3]
  
  c_o <- iox_obj$coords[iox_obj$dag_info$order, ]
  
  Theta4 <- array(1, dim=c(4, q, mcmc))
  Theta4[1,,] <- iox_obj$Theta[1,,]
  Theta4[3:4,,] <- iox_obj$Theta[2:3,,]
  
  result <- spiox:::iox_make_fij0(c_o, 
                                 iox_obj$dag_info$dag, 
                                 iox_obj$dag_opts, 
                                 Theta4, matern, n_threads)
  
  return(result)
}


