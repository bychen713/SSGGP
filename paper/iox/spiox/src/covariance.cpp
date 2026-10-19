#include "covariance.h"

using namespace std;

// matern
void matern_inplace(arma::mat& res, 
                             const arma::mat& coords,
                             const arma::uvec& ix, const arma::uvec& iy, 
                             const double& phi, const double& nu, 
                             const double& sigmasq, 
                             double * bessel_ws,
                             const double& alpha=0,  bool same=false){
  
  int threadid = 0;
#ifdef _OPENMP
  threadid = omp_get_thread_num();
#endif
  
  // 1-alpha is the proportion of variance explained by spatial component
  // alpha is is the proportion explained by nugget effect
  
  double pow2_nu1_gammanu_sigmasq = (1-alpha) * sigmasq * pow(2.0, 1.0-nu) / R::gammafn(nu);
  
  if(same){
    for(unsigned int i=0; i<ix.n_rows; i++){
      arma::rowvec cri = coords.row(ix(i));//x.row(i);
      for(unsigned int j=i; j<iy.n_rows; j++){
        arma::rowvec delta = cri - coords.row(iy(j));//y.row(j);
        double hphi = arma::norm(delta) * phi;
        if(hphi > 0.0){
          res(i, j) = pow(hphi, nu) * pow2_nu1_gammanu_sigmasq *
            R::bessel_k_ex(hphi, nu, 1.0, &bessel_ws[threadid*MAT_NU_MAX]);
        } else {
          res(i, j) = sigmasq;
        }
      }
    }
    res = arma::symmatu(res);
  } else {
    for(unsigned int i=0; i<ix.n_rows; i++){
      arma::rowvec cri = coords.row(ix(i));//x.row(i);
      for(unsigned int j=0; j<iy.n_rows; j++){
        arma::rowvec delta = cri - coords.row(iy(j));//y.row(j);
        double hphi = arma::norm(delta) * phi;
        if(hphi > 0.0){
          res(i, j) = pow(hphi, nu) * pow2_nu1_gammanu_sigmasq *
            R::bessel_k_ex(hphi, nu, 1.0, &bessel_ws[threadid*MAT_NU_MAX]);
        } else {
          res(i, j) = sigmasq;
        }
      }
    }
  }
  //return res;
}


// powered exponential nu<2
void powerexp_inplace(arma::mat& res, 
                      const arma::mat& coords,
                      const arma::uvec& ix, const arma::uvec& iy, 
                      const double& phi, const double& nu, 
                      const double& sigmasq, const double& alpha,
                      bool same){
  
  // 1-alpha is the proportion of variance explained by spatial component
  // alpha is is the proportion explained by nugget effect
  
  
  if(same){
    for(unsigned int i=0; i<ix.n_rows; i++){
      arma::rowvec cri = coords.row(ix(i));
      for(unsigned int j=i; j<iy.n_rows; j++){
        arma::rowvec delta = cri - coords.row(iy(j));
        double hh = arma::norm(delta);
        if(hh==0){
          res(i, j) = sigmasq;
        } else {
          double hnuphi = pow(hh, nu) * phi;
          res(i, j) = exp(-hnuphi) * sigmasq * (1-alpha);
        }
      }
    }
    res = arma::symmatu(res);
  } else {
    for(unsigned int i=0; i<ix.n_rows; i++){
      arma::rowvec cri = coords.row(ix(i));
      for(unsigned int j=0; j<iy.n_rows; j++){
        arma::rowvec delta = cri - coords.row(iy(j));
        double hh = arma::norm(delta);
        if(hh==0){
          res(i, j) = sigmasq;
        } else {
          double hnuphi = pow(hh, nu) * phi;
          res(i, j) = exp(-hnuphi) * sigmasq * (1-alpha);
        }
      }
    }
  }
}

void wave_inplace(arma::mat& res, 
                      const arma::mat& coords,
                      const arma::uvec& ix, const arma::uvec& iy, 
                      const double& phi, const double& nu, 
                      const double& sigmasq, const double& alpha,
                      bool same){
  
  // exp(-phi * d) cos(nu * d),  phi >= nu

  // 1-alpha is the proportion of variance explained by spatial component
  // alpha is is the proportion explained by nugget effect
  
  if(same){
    for(unsigned int i=0; i<ix.n_rows; i++){
      arma::rowvec cri = coords.row(ix(i));
      for(unsigned int j=i; j<iy.n_rows; j++){
        arma::rowvec delta = cri - coords.row(iy(j));
        double hh = arma::norm(delta);
        if(hh==0){
          res(i, j) = sigmasq;
        } else {
          double hphi = hh * phi;
          double hnu = hh * nu;
          res(i, j) = exp(-hphi) * cos(hnu) * sigmasq * (1-alpha);
        }
      }
    }
    res = arma::symmatu(res);
  } else {
    for(unsigned int i=0; i<ix.n_rows; i++){
      arma::rowvec cri = coords.row(ix(i));
      for(unsigned int j=0; j<iy.n_rows; j++){
        arma::rowvec delta = cri - coords.row(iy(j));
        double hh = arma::norm(delta);
        if(hh==0){
          res(i, j) = sigmasq;
        } else {
          double hphi = hh * phi;
          double hnu = hh * nu;
          res(i, j) = exp(-hphi) * cos(hnu) * sigmasq * (1-alpha);
        }
      }
    }
  }
}


arma::mat Correlationf(
    const arma::mat& coords,
    const arma::uvec& ix, const arma::uvec& iy,
    const arma::vec& theta,
    double * bessel_ws,
    int matern,
    bool same){
  arma::mat res = arma::zeros(ix.n_rows, iy.n_rows);
  
  double phi = theta(0);
  double sigmasq = theta(1);
  double nu = theta(2); // exponent in the power exponential case
  
  double alpha = theta(3);
  // 1-alpha is the proportion of variance explained by spatial component
  // alpha is is the proportion explained by nugget effect
  
  if(matern == 1){
    matern_inplace(res, coords, ix, iy, phi, nu, sigmasq, bessel_ws, alpha, same);   
  } else if(matern == 0) {
    powerexp_inplace(res, coords, ix, iy, phi, nu, sigmasq, alpha, same); 
  } else if(matern == 2) { 
    wave_inplace(res, coords, ix, iy, phi, nu, sigmasq, alpha, same);
  }
  return res;
}

//[[Rcpp::export]]
arma::mat Correlationc(
    const arma::mat& coordsx,
    const arma::mat& coordsy,
    const arma::vec& theta,
    int matern,
    bool same){

  int nthreads = 1;
  //int bessel_ws_inc = 3;//see bessel_k.c for working space needs
  double *bessel_ws = (double *) R_alloc(nthreads*MAT_NU_MAX, sizeof(double));
  
  if(same){
    arma::uvec ix = arma::regspace<arma::uvec>(0, coordsx.n_rows-1);
    
    return Correlationf(coordsx, ix, ix, theta, bessel_ws, matern, same);
  } else {
    arma::mat coords = arma::join_vert(coordsx, coordsy);
    arma::uvec ix = arma::regspace<arma::uvec>(0, coordsx.n_rows-1);
    arma::uvec iy = arma::regspace<arma::uvec>(coordsx.n_rows, coords.n_rows-1);
    
    return Correlationf(coords, ix, iy, theta, bessel_ws, matern, same);
  }
  
}

