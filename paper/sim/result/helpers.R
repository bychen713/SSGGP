# Helpers for the simulation result tables: true and estimated precision
# matrices, graph-recovery rates, and prediction summaries for each method.

trimat <- function(P,diag,offdiag){
  Sigma_inv <- matrix(nrow = P, ncol = P, 0)
  diag(Sigma_inv) <- diag
  if (P > 1) {
    r_entries1 <- 1:(P - 1)
    r_entries2 <- (1:(P - 1)) + 1
    Sigma_inv[cbind(r_entries1,r_entries2)] <- offdiag
    Sigma_inv[cbind(r_entries2,r_entries1)] <- offdiag
  }
  return(Sigma_inv)
}

build_matrix <- function(row, P, type_prefix) {
  mat <- matrix(0, P, P)
  for (i in 1:P) {
    for (j in i:P) {
      key <- paste0(type_prefix, i, j)
      if (key %in% names(row)) {
        val <- row[[key]]
        mat[i, j] <- val
        mat[j, i] <- val
      }
    }
  }
  return(mat)
}

# False discovery, false positive, and true positive rates of the estimated graph,
# and the F2 score.
compute_sparsity <- function(row_df, P, prefix, true_pos_mask, true_neg_mask, n_true_pos, n_true_neg, tol = 1e-8) {
  est_mat     <- build_matrix(row_df, P, prefix)
  est_nonzero <- abs(est_mat) > tol
  
  TP <- sum(est_nonzero & true_pos_mask)
  FP <- sum(est_nonzero & true_neg_mask)
  FN <- sum(!est_nonzero & true_pos_mask)
  
  precision = TP/(TP+FP)
  recall = TP/(TP+FN)
  
  FPR <- FP / n_true_neg
  TPR <- 1 - (FN / n_true_pos)
  FDR <- if ((FP + TP) > 0) FP / (FP + TP) else 0
  F1 <- 2*TP/(2*TP+FP+FN)
  F2 <- 5*precision*recall / (4*precision+recall)
  
  return(c(FDR = FDR, FPR = FPR, TPR = TPR, F2=F2))
}

# Prediction and parameter summaries of one fit. approach is "ind", "lmc", or a
# SSGGP fit ("dense" or "BIC"); the result is list(data_plot = ...).
get_result <- function(result, approach){
  if(approach=="ind"){
    
    Y_obs_ind <- c(t(result$Y))
    Y_mat_test_ind <- result$Y_test
    Y_mat_test_ind_reord <- as.vector(t(Y_mat_test_ind))
    Y_mat_pred_ind <- result$prediction
    grid_ind <- result$grid
    pred_var_ind <- result$pred_var
    
    # prediction error
    mask <- rowSums(is.na(Y_mat_test_ind))==0  # test locations with every process held out
    Y_mat_pred_ind <- coverage_var_ind <- Y_mat_test_ind
    for(k in 1:P){
      mask_vec <- !is.na(Y_mat_pred_ind[,k])
      Y_mat_pred_ind[mask_vec,k] <- result$prediction[[k]]
      coverage_var_ind[mask_vec,k] <- pred_var_ind[[k]]
    }
    Y_mat_test_ind <- Y_mat_test_ind[mask,]
    Y_mat_pred_ind <- Y_mat_pred_ind[mask,]
    coverage_var_ind <- coverage_var_ind[mask,]
    
    n <- nrow(result$grid)
    n_test <- nrow(Y_mat_pred_ind)
    n_all <- n + n_test
    
    norm_ind <- mean(apply((Y_mat_test_ind-Y_mat_pred_ind),1,function(x) sqrt(sum(x^2))))
    mse_ind <- apply((Y_mat_test_ind-Y_mat_pred_ind)**2,2,mean)
    var_ind <- colMeans(sqrt(coverage_var_ind))   # mean predictive SD
    
    # inference
    coverage_ind <- vector("list",P)
    for(k in 1:P){
      coverage_ind[[k]] <- Y_mat_pred_ind[,k]-1.96*sqrt(coverage_var_ind[,k])<=Y_mat_test_ind[,k] & 
        Y_mat_test_ind[,k]<=Y_mat_pred_ind[,k]+1.96*sqrt(coverage_var_ind[,k])
    }
    coverage_ind <- unlist(lapply(coverage_ind,mean))
    
    # sum of the processes
    mse_sum_ind <- mean((rowSums(Y_mat_test_ind)-rowSums(Y_mat_pred_ind))**2)
    coverage_sum_var_ind <- rowSums(coverage_var_ind)
    varsum_ind <- mean(sqrt(coverage_sum_var_ind))
    
    coveragesum_ind <- rowSums(Y_mat_pred_ind)-1.96*sqrt(coverage_sum_var_ind)<=rowSums(Y_mat_test_ind) & 
      rowSums(Y_mat_test_ind)<=rowSums(Y_mat_pred_ind)+1.96*sqrt(coverage_sum_var_ind)
    coveragesum_ind <- sum(coveragesum_ind)/n_test
    
    data_plot <- list(sigma_inv_est = sapply(result$opt_true,"[[","par")[1,],
                      a_est = sapply(result$opt_true,"[[","par")[2,],
                      nu_est = sapply(result$opt_true,"[[","par")[3,],
                      D_est = sapply(result$opt_true,"[[","par")[4,],
                      norm = norm_ind,
                      mse = mse_ind, 
                      coverage = coverage_ind, 
                      var_width = var_ind,
                      mse_sum_ind=mse_sum_ind,
                      coveragesum_ind=coveragesum_ind,
                      varsum_ind=varsum_ind
    )
    output <- list(data_plot=data_plot)
    return(output)
  }else if(approach=="lmc"){
    
    Y_mat_test_lmc <- result$Y_test
    Y_mat_pred_lmc <- result$prediction
    pred_var_lmc <- pred_var_mat <- t(sapply(result$pred_var, diag))
    
    # prediction error
    mask <- rowSums(is.na(Y_mat_test_lmc))==0  # test locations with every process held out
    Y_mat_pred_lmc <- coverage_var_lmc <- Y_mat_test_lmc
    for(k in 1:P){
      mask_vec <- !is.na(Y_mat_pred_lmc[,k])
      Y_mat_pred_lmc[mask_vec,k] <- result$prediction[,k]
      coverage_var_lmc[mask_vec,k] <- pred_var_lmc[,k]
    }
    Y_mat_test_lmc <- Y_mat_test_lmc[mask,]
    Y_mat_pred_lmc <- Y_mat_pred_lmc[mask,]
    coverage_var_lmc <- coverage_var_lmc[mask,]
    
    n <- nrow(result$grid)
    n_test <- nrow(Y_mat_pred_lmc)
    n_all <- n + n_test
    
    norm_lmc <- mean(apply((Y_mat_test_lmc-Y_mat_pred_lmc),1,function(x) sqrt(sum(x^2))))
    mse_lmc <- apply((Y_mat_test_lmc-Y_mat_pred_lmc)**2,2,mean)
    var_lmc <- colMeans(sqrt(coverage_var_lmc))   # mean predictive SD
    
    # inference
    coverage_lmc <- vector("list",P)
    for(k in 1:P){
      coverage_lmc[[k]] <- Y_mat_pred_lmc[,k]-1.96*sqrt(coverage_var_lmc[,k])<=Y_mat_test_lmc[,k] & 
        Y_mat_test_lmc[,k]<=Y_mat_pred_lmc[,k]+1.96*sqrt(coverage_var_lmc[,k])
    }
    coverage_lmc <- unlist(lapply(coverage_lmc,mean))
    
    # sum of the processes
    mse_sum_lmc <- mean((rowSums(Y_mat_test_lmc)-rowSums(Y_mat_pred_lmc))**2)
    # Variance of the sum is 1' Sigma 1 over the full predictive covariance.
    coverage_sum_var_lmc <- vapply(result$pred_var[mask], sum, numeric(1))
    varsum_lmc <- mean(sqrt(coverage_sum_var_lmc))
    
    coveragesum_lmc <- rowSums(Y_mat_pred_lmc)-1.96*sqrt(coverage_sum_var_lmc)<=rowSums(Y_mat_test_lmc) & 
      rowSums(Y_mat_test_lmc)<=rowSums(Y_mat_pred_lmc)+1.96*sqrt(coverage_sum_var_lmc)
    coveragesum_lmc <- sum(coveragesum_lmc)/n_test
    
    data_plot <- list(sigma_est = result$opt_lmc[[1]],   # est_B
                      a_est = result$opt_lmc[[2]],
                      nu_est = result$opt_lmc[[3]],
                      D_est = result$opt_lmc[[4]],       # est_nuggets
                      norm = norm_lmc,
                      mse = mse_lmc, 
                      coverage = coverage_lmc,
                      var_width = var_lmc,
                      mse_sum_lmc=mse_sum_lmc,
                      coveragesum_lmc=coveragesum_lmc,
                      varsum_lmc=varsum_lmc
    )
    output <- list(data_plot=data_plot)
  }else{
    Y_obs <- result$Y
    Y_mat_test <- result$Y_mat_test
    grid <- result$grid
    
    tmp <- t(Y_mat_test)
    mask <- which(!is.na(tmp))
    
    tmp[mask] <- as.numeric(result$prediction_gibbs)
    Y_mat_pred_approx_gibbs <- t(tmp)
    
    pred_var_approx <- t(tmp)
    tmp[mask] <- as.numeric(diag(result$pred_var_gibbs))
    pred_var_approx_gibbs <- t(tmp)
    
    mask <- rowSums(is.na(Y_mat_test))==0
    Y_mat_test <- Y_mat_test[mask,]
    Y_mat_pred_approx_gibbs <- Y_mat_pred_approx_gibbs[mask,]
    pred_var_approx_gibbs <- pred_var_approx_gibbs[mask,]
    n_test <- nrow(Y_mat_pred_approx_gibbs)
    
    # parameter estimation
    D_est_vals_approx <- result$D_est_vals
    sigma_inv_est_vals_approx <- result$sigma_inv_est_vals
    a_approx <- result$a_est
    nu_approx <- result$nu_est
    sigma_inv_est_vals_approx <- sigma_inv_est_vals_approx[,!is.na(sigma_inv_est_vals_approx[1,]),drop=F]
    sigma_inv_est_approx <- sigma_inv_est_vals_approx[,ncol(sigma_inv_est_vals_approx)]
    ind <- which(upper.tri(matrix(nrow=P, ncol=P),diag=T),arr.ind=T)
    params_approx_approx <- sigma_inv_est_approx[order(ind[,2]-ind[,1])]
    
    # prediction error
    norm_approx_gibbs <- mean(apply((Y_mat_test-Y_mat_pred_approx_gibbs),1,function(x) sqrt(sum(x^2))))
    mse_approx_gibbs <- apply((Y_mat_test-Y_mat_pred_approx_gibbs)**2, 2, mean)
    
    # inference
    coverage_approx_gibbs <- vector('list',P)
    for(k in 1:P){
      coverage_approx_gibbs[[k]] <- Y_mat_pred_approx_gibbs[,k]-1.96*sqrt(pred_var_approx_gibbs[,k])<=Y_mat_test[,k] &
        Y_mat_test[,k]<=Y_mat_pred_approx_gibbs[,k]+1.96*sqrt(pred_var_approx_gibbs[,k])
    }
    coverage_approx_gibbs <- unlist(lapply(coverage_approx_gibbs,mean))
    var_approx_gibbs <- colMeans(sqrt(pred_var_approx_gibbs))   # mean predictive SD
    
    # sum of the processes
    mse_sum_approx_gibbs <- mean((rowSums(Y_mat_test)-rowSums(Y_mat_pred_approx_gibbs))**2)
    coverage_sum_var_approx_gibbs <- unlist(lapply(1:n_test,function(x){as.vector(t(rep(1,P))%*%result$pred_var_gibbs[((1:P)+P*(x-1)),((1:P)+P*(x-1))]%*%rep(1,P))}))
    varsum_approx_gibbs <- mean(sqrt(coverage_sum_var_approx_gibbs))
    coveragesum_approx_gibbs <- rowSums(Y_mat_pred_approx_gibbs)-1.96*sqrt(coverage_sum_var_approx_gibbs)<=rowSums(Y_mat_test) &
      rowSums(Y_mat_test)<=rowSums(Y_mat_pred_approx_gibbs)+1.96*sqrt(coverage_sum_var_approx_gibbs)
    coveragesum_approx_gibbs <- sum(coveragesum_approx_gibbs)/n_test
    
    data_plot <- list(params_approx_approx = params_approx_approx,
                      D_est = tail(D_est_vals_approx[!is.na(D_est_vals_approx)],1),
                      a_est = a_approx,
                      nu_est = nu_approx,
                      rho = result$rho,
                      ll = result$ll_MCEM,
                      norm_gibbs = norm_approx_gibbs,
                      mse_gibbs = mse_approx_gibbs,
                      coverage_gibbs = coverage_approx_gibbs,
                      var_width_gibbs = var_approx_gibbs,
                      mse_sum_gibbs = mse_sum_approx_gibbs, 
                      coveragesum_gibbs = coveragesum_approx_gibbs,
                      varsum_gibbs=varsum_approx_gibbs 
    )
    output <- list(data_plot=data_plot)
    output
  }
}

