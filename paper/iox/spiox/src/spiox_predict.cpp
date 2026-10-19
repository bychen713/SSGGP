#include "spiox.h"

//[[Rcpp::export]]
Rcpp::List spiox_predict(
                   const arma::mat& X_new,
                   const arma::mat& coords_new,
                   
                   const arma::mat& Y, 
                   const arma::mat& X, 
                   const arma::mat& coords,
                   
                   const arma::field<arma::uvec>& dag,
                   
                   const arma::cube& B,
                   const arma::cube& Sigma,
                   const arma::cube& theta,
                   int matern = 1,
                   int num_threads = 1
){
  int q = Y.n_cols;
  int p = X.n_cols;
  int mcmc = B.n_slices;
  
  arma::uvec oneuv = arma::ones<arma::uvec>(1);
  
  int bessel_ws_inc = MAT_NU_MAX;//see bessel_k.c for working space needs
  double * bessel_ws = (double *) R_alloc(num_threads*bessel_ws_inc, sizeof(double));

  int ntrain = coords.n_rows;
  int ntest = coords_new.n_rows;
  
  arma::mat cxall = arma::join_vert(coords, coords_new);
  
  arma::cube Y_out_mcmc = arma::zeros(ntest, q, mcmc);
  arma::cube random_stdnormal = arma::randn(mcmc, q, ntest);
  
  // loop over test locations
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(num_threads)
  #endif
  for(int i=0; i<ntest; i++){
    
    arma::uvec ix = oneuv * i + ntrain;
    arma::uvec px = dag(i);
    
    arma::mat theta_current = arma::zeros(4, q);
    
    arma::field<arma::mat> CC(q);
    arma::field<arma::mat> CPt(q);
    arma::field<arma::mat> PPi(q);
    arma::field<arma::mat> ht(q);
    arma::field<arma::mat> CPh(q);
    arma::vec Dj = arma::zeros(q);
    
    for(int m=0; m<mcmc; m++){
      
      arma::vec xb_new = arma::trans( X_new.row(i) * B.slice(m) );
      arma::mat yxb_old = Y - X * B.slice(m);
      
      arma::mat theta_m = theta.slice(m);
      double theta_diff = abs(arma::accu(theta_m - theta_current));
      
      arma::mat Sigmam = Sigma.slice(m);
      arma::mat Sigmauchol = arma::chol(Sigmam, "upper");
      
      //arma::mat W_out = arma::zeros(ntest, q);
      arma::vec Y_out = arma::zeros(q);
      arma::mat rndnorm_m = random_stdnormal.row(m);
    
      
      arma::vec Yspat = arma::zeros(q);
      
      // loop over outcomes
      for(int j=0; j<q; j++){
        // predict spatially
        arma::uvec outjx = oneuv * j;

        if(theta_diff > 1e-15){
          CC(j) = Correlationf(cxall, ix, ix, theta_m.col(j), bessel_ws, matern, true); // 1 for matern 
          CPt(j) = Correlationf(cxall, px, ix, theta_m.col(j), bessel_ws, matern, false);
          PPi(j) = arma::inv_sympd( Correlationf(cxall, px, px, theta_m.col(j), bessel_ws, matern, true) );  
          ht(j) = PPi(j) * CPt(j);
          CPh(j) = CPt(j).t() * ht(j);
          Dj(j) = sqrt( abs(arma::conv_to<double>::from(
            CC(j) - CPh(j) )) ); // abs for numerical zeros
        }
        
        Yspat(j) = arma::conv_to<double>::from(ht(j).t() * yxb_old(px, outjx)); 
        Y_out(j) = xb_new(j) + Yspat(j); // post mean only
      }
      
      Y_out += arma::diagmat(Dj) * Sigmauchol.t() * rndnorm_m.col(i) ; // pred uncert
      Y_out_mcmc.subcube(i, 0, m, i, q-1, m) = Y_out;
      
      theta_current = theta_m;
    }
    
    //Rcpp::Rcout << "saving m" << endl;
    //Y_out_mcmc.slice(m) = Y_out;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("Y") = Y_out_mcmc,
    Rcpp::Named("dag") = dag,
    Rcpp::Named("seedmvn") = random_stdnormal
  );
}


//[[Rcpp::export]]
Rcpp::List spiox_predict_part(
    const arma::mat& Y_new, 
    const arma::mat& X_new,
    const arma::mat& coords_new,
    
    const arma::mat& Y, 
    const arma::mat& X, 
    const arma::mat& coords,
    
    const arma::field<arma::uvec>& dag,
    
    const arma::cube& B,
    const arma::cube& Sigma,
    const arma::cube& theta,
    int matern = 1,
    int num_threads = 1
){
  // Y_new is a matrix with as many rows as X_new and coords_new
  // and q columns. NA in Y_new are to be predicted
  int q = Y.n_cols;
  int p = X.n_cols;
  int mcmc = B.n_slices;
  
  arma::uvec oneuv = arma::ones<arma::uvec>(1);
  
  int bessel_ws_inc = MAT_NU_MAX;//see bessel_k.c for working space needs
  double * bessel_ws = (double *) R_alloc(num_threads*bessel_ws_inc, sizeof(double));
  
  int ntrain = coords.n_rows;
  int ntest = coords_new.n_rows;
  
  arma::mat cxall = arma::join_vert(coords, coords_new);
  
  arma::cube Y_out_mcmc = arma::zeros(ntest, q, mcmc);
  arma::cube random_stdnormal = arma::randn(mcmc, q, ntest);

  // loop over test locations
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
  for(int i=0; i<ntest; i++){
      //Rcpp::Rcout << "i: " << i << endl;
      
      arma::uvec missing = arma::find_nan(Y_new.row(i));
      arma::uvec observd = arma::find_finite(Y_new.row(i));
      int n_miss = missing.n_elem;
      int n_obsv = observd.n_elem; 
      
      arma::uvec ix = oneuv * i + ntrain;
      arma::uvec px = dag(i);
      
      arma::mat theta_current = arma::zeros(4, q);
      
      arma::field<arma::mat> CC(q);
      arma::field<arma::mat> CPt(q);
      arma::field<arma::mat> PPi(q);
      
      for(int m=0; m<mcmc; m++){
        
        //Rcpp::Rcout << "              m: " << m << endl;
        
        arma::vec xb_new = arma::trans( X_new.row(i) * B.slice(m) );
        arma::mat yxb_old = Y - X * B.slice(m);
        
        arma::mat theta_m = theta.slice(m);
        double theta_diff = abs(arma::accu(theta_m - theta_current));
        
        arma::mat Sigmam = Sigma.slice(m);
        
        //arma::mat W_out = arma::zeros(ntest, q);
        arma::vec Y_out = arma::zeros(q);
        arma::mat rndnorm_m = random_stdnormal.row(m);
        
        arma::vec Djm = arma::zeros(n_miss);
        arma::vec Djo = arma::zeros(n_obsv);
        
        arma::vec Ym_spatl = arma::zeros(n_miss);
        arma::vec Yo_resid = arma::zeros(n_obsv);
        
        //Rcpp::Rcout << "building comps " << endl;
        for(int j=0; j<q; j++){
          if(theta_diff > 1e-15){
            CC(j) = Correlationf(cxall, ix, ix, theta_m.col(j), bessel_ws, matern, true); // 1 for matern 
            CPt(j) = Correlationf(cxall, px, ix, theta_m.col(j), bessel_ws, matern, false);
            PPi(j) = arma::inv_sympd( Correlationf(cxall, px, px, theta_m.col(j), bessel_ws, matern, true) );  
          }
        }
        
        //Rcpp::Rcout << "building obvs" << endl;
        for(int jx=0; jx<n_obsv; jx++){
          int j = observd(jx);
          arma::uvec outjx = oneuv * j;
          arma::vec ht = PPi(j) * CPt(j);
          double sqrtR = sqrt( abs(arma::conv_to<double>::from(
            CC(j) - CPt(j).t() * ht )) ); // abs for numerical zeros
          Yo_resid(jx) = Y_new(i, j) - xb_new(j) - 
            arma::conv_to<double>::from(ht.t() * yxb_old(px, outjx)); // post mean only
          Djo(jx) = sqrtR;
        }
        
        //Rcpp::Rcout << "building misss" << endl;
        for(int jx=0; jx<n_miss; jx++){
          int j = missing(jx);
          arma::uvec outjx = oneuv * j;
          arma::vec ht = PPi(j) * CPt(j);
          double sqrtR = sqrt( abs(arma::conv_to<double>::from(
            CC(j) - CPt(j).t() * ht )) ); // abs for numerical zeros
          Ym_spatl(jx) = xb_new(j) + 
            arma::conv_to<double>::from(ht.t() * yxb_old(px, outjx)); // post mean only
          Djm(jx) = sqrtR;
        }
        
        //Rcpp::Rcout << "finally" << endl;
        
        arma::mat S_mo = Sigmam(missing, observd);
        arma::mat S_o_inv = arma::inv_sympd(Sigmam(observd, observd));
        arma::mat S_m = Sigmam(missing, missing);
        
        arma::mat Hmo = arma::diagmat(Djm) * S_mo * S_o_inv * arma::diagmat(1.0/Djo);
        arma::mat Rmo = arma::diagmat(Djm) * (S_m - S_mo * S_o_inv * S_mo.t()) * arma::diagmat(Djm);
        
        arma::vec missing_post_mean = Hmo * Yo_resid;
        arma::mat missing_post_cholvar = arma::chol(arma::symmatu(Rmo), "lower");
        
        //Rcpp::Rcout << "finally 3" << endl;
        arma::mat rnorm_all = random_stdnormal.subcube(m, 0, i, m, q-1, i) ;
        arma::vec rnorm_miss = arma::trans( rnorm_all.cols(missing) );
        arma::vec Ym_pred = Ym_spatl + missing_post_mean + missing_post_cholvar * rnorm_miss;
        
        for(int jx=0; jx<n_miss; jx++){
          int j = missing(jx);
          Y_out_mcmc(i, j, m) = Ym_pred(jx);
        }
        for(int jx=0; jx<n_obsv; jx++){
          int j = observd(jx);
          Y_out_mcmc(i, j, m) = Y_new(i, j);
        }
        
        theta_current = theta_m;
      }
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("Y") = Y_out_mcmc,
    Rcpp::Named("dag") = dag
  );
}


//[[Rcpp::export]]
Rcpp::List spiox_latent_predict(
    const arma::mat& X_new,
    const arma::mat& coords_new,
    
    const arma::mat& coords,
    
    const arma::field<arma::uvec>& dag,
    
    const arma::cube& W,
    const arma::cube& B,
    const arma::cube& Sigma,
    const arma::mat& Dvec,
    const arma::cube& theta,
    int matern = 1,
    int num_threads = 1
){
  int q = W.n_cols;
  int p = X_new.n_cols;
  int mcmc = B.n_slices;
  
  arma::uvec oneuv = arma::ones<arma::uvec>(1);
  
  int bessel_ws_inc = MAT_NU_MAX;//see bessel_k.c for working space needs
  double * bessel_ws = (double *) R_alloc(num_threads*bessel_ws_inc, sizeof(double));
  
  int ntrain = coords.n_rows;
  int ntest = coords_new.n_rows;
  
  arma::mat cxall = arma::join_vert(coords, coords_new);
  
  arma::cube Y_out_mcmc = arma::zeros(ntest, q, mcmc);
  arma::cube random_stdnormal_w = arma::randn(mcmc, q, ntest);
  arma::cube random_stdnormal_y = arma::randn(mcmc, q, ntest);
  
  // loop over test locations
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
  for(int i=0; i<ntest; i++){
    
    arma::uvec ix = oneuv * i + ntrain;
    arma::uvec px = dag(i);
    
    arma::mat theta_current = arma::zeros(4, q);
    
    arma::field<arma::mat> CC(q);
    arma::field<arma::mat> CPt(q);
    arma::field<arma::mat> PPi(q);
    
    for(int m=0; m<mcmc; m++){
      
      arma::mat xb_new = X_new * B.slice(m);
      arma::mat W_mcmc = W.slice(m);
      
      arma::mat theta_m = theta.slice(m);
      double theta_diff = abs(arma::accu(theta_m - theta_current));
      
      arma::mat Sigmam = Sigma.slice(m);
      arma::mat Sigmauchol = arma::chol(Sigmam, "upper");
      
      arma::vec deltsqrt = sqrt(Dvec.col(m));
      
      arma::mat W_out = arma::zeros(ntest, q);
      arma::mat Y_out = arma::zeros(ntest, q);
      arma::mat rndnorm_m = random_stdnormal_w.row(m);
      arma::mat rndnorm_y = random_stdnormal_y.row(m);
      arma::vec Dj = arma::zeros(q);
      
      // loop over outcomes
      for(int j=0; j<q; j++){
        // predict spatially
        arma::uvec outjx = oneuv * j;
        
        if(theta_diff > 1e-15){
          CC(j) = Correlationf(cxall, ix, ix, theta_m.col(j), bessel_ws, matern, true); // 1 for matern 
          CPt(j) = Correlationf(cxall, px, ix, theta_m.col(j), bessel_ws, matern, false);
          PPi(j) = arma::inv_sympd( Correlationf(cxall, px, px, theta_m.col(j), bessel_ws, matern, true) );  
        }
        
        arma::vec ht = PPi(j) * CPt(j);
        double sqrtR = sqrt( abs(arma::conv_to<double>::from(
          CC(j) - CPt(j).t() * ht )) ); // abs for numerical zeros
        
        W_out(i, j) = 
          arma::conv_to<double>::from(ht.t() * W_mcmc(px, outjx)); // post mean only
        Dj(j) = sqrtR;
      }
      
      W_out.row(i) += arma::trans( arma::diagmat(Dj) * Sigmauchol.t() * rndnorm_m.col(i) ); // pred uncert
      
      Y_out_mcmc.subcube(i, 0, m, i, q-1, m) = 
        xb_new.row(i) + W_out.row(i) + 
        arma::trans(arma::diagmat(deltsqrt) * rndnorm_y.col(i));
  
      theta_current = theta_m;
    }
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("Y") = Y_out_mcmc,
    Rcpp::Named("dag") = dag
  );
}

