
dag_vecchia_orig_order <- function(coords, m=15, gridded=FALSE){
  # returns the dag for the original coordinate system
  # this makes the H matrices not lower triangular : not a problem if we don't need them to be
  if(gridded){
    dag_for_gridded_cols(coords, m)
  } else {
    ixmm <- MaxMincpp(coords)
    ixmm_order <- order(ixmm)
    
    coords_mm <- coords[ixmm,]
    nn_dag_mat <- findOrderedNN_kdtree2(coords_mm, m)
    nn_dag <- apply(nn_dag_mat, 1, function(x){ x[!is.na(x)][-1]-1 })
    nn_dag_reorder <- (nn_dag %>% lapply(\(dag) ixmm[dag+1]-1 ))[ixmm_order] 
    return(list(dag=nn_dag_reorder, order=ixmm_order))
  }
}


dag_vecchia_o <- function(coords, m=15, gridded=FALSE){
  # returns the dag on the new order
  # this makes the H matrices lower triangular : good for using sparse triangular solver
  if(gridded){
    d <- ncol(coords)
    ord <- do.call(order, lapply(d:1, function(k) coords[, k]))
    nn_dag <- dag_for_gridded_cols(coords[ord, , drop = FALSE], m)
    return(list(dag = nn_dag, order = ord))
  } else {
    ixmm <- MaxMincpp(coords)
    ixmm_order <- order(ixmm)
    
    coords_mm <- coords[ixmm,]
    nn_dag_mat <- findOrderedNN_kdtree2(coords_mm, m)
    nn_dag <- apply(nn_dag_mat, 1, function(x){ x[!is.na(x)][-1]-1 })#[ixmm]
    return(list(dag=nn_dag, order=ixmm))
  }
}


dag_vecchia_predict <- function(coords, newcoords, m=15){
  nn.found <- FNN::get.knnx(coords, newcoords, k=m)
  
  pred_dag <- as.list(as.data.frame(t(nn.found$nn.index-1)))
  return(pred_dag)
}