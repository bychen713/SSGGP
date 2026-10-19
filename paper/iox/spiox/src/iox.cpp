#include <RcppArmadillo.h>
#include "covariance.h"
#include "omp_import.h"
#include "daggp.h"
#include "interrupt.h"

using namespace std;

//[[Rcpp::export]]
arma::mat iox(const arma::mat& x, const arma::mat& y, int i, int j,
                   const arma::mat& S, 
                   const arma::mat& theta, 
                   bool matern = true,
                   bool diag_only=false, bool at_limit=false){

  arma::vec thetai = theta.col(i-1);
  arma::mat Ki = Correlationc(S, S, thetai, matern, true);
  arma::mat Li = arma::chol(Ki, "lower");
  arma::mat Li_inv = arma::inv(arma::trimatl(Li));
  arma::mat Ki_inv = Li_inv.t() * Li_inv;
  arma::mat rhoi_xS = Correlationc(x, S, thetai, matern, false);
  
  arma::vec thetaj = theta.col(j-1);
  arma::mat Kj = Correlationc(S, S, thetaj, matern, true);
  arma::mat Lj = arma::chol(Kj, "lower");
  arma::mat Lj_inv = arma::inv(arma::trimatl(Lj));
  arma::mat Kj_inv = Lj_inv.t() * Lj_inv;
  arma::mat rhoj_yS = Correlationc(y, S, thetaj, matern, false);
  
  arma::mat R = arma::zeros(x.n_rows, y.n_rows);
  arma::vec Rdiag = arma::zeros(x.n_rows);
  
  if(!at_limit){
    arma::mat fixcoord = arma::zeros(1, 2);
    double atZerodisti = Correlationc(fixcoord, fixcoord, thetai, matern, true)(0,0);
    double atZerodistj = Correlationc(fixcoord, fixcoord, thetaj, matern, true)(0,0);
    
    for(unsigned int r=0; r<x.n_rows; r++){
      if(!diag_only){
        bool Ri_computed = false; // avoid computing until we actually need 
        double Ri = 0;
        for(unsigned int c=0; c<y.n_rows; c++){
          double dist = arma::accu(abs( x.row(r) - y.row(c) ));
          if(dist == 0){
            if(!Ri_computed){ // compute Ri 
              Ri = atZerodisti - arma::conv_to<double>::from(
                rhoi_xS.row(r) * Ki_inv * arma::trans(rhoi_xS.row(r)) );
              if(Ri < 0){ Ri = 0; } // numerical zero, coord is in S
            }
            double Rj = atZerodistj - arma::conv_to<double>::from(
              rhoj_yS.row(c) * Kj_inv * arma::trans(rhoj_yS.row(c)) );
            if(Rj < 0){ Rj = 0; } // numerical zero, coord is in S
            R(r,c) = sqrt(Ri * Rj);
          }
        }
      } else {
        bool Ri_computed = false; // avoid computing until we actually need 
        double Ri = 0;
        int c=r; 
        double dist = arma::accu(abs( x.row(r) - y.row(c) ));
        if(dist == 0){
          if(!Ri_computed){ // compute Ri 
            Ri = atZerodisti - arma::conv_to<double>::from(
              rhoi_xS.row(r) * Ki_inv * arma::trans(rhoi_xS.row(r)) );
            if(Ri < 0){ Ri = 0; } // numerical zero, coord is in S
          }
          double Rj = atZerodistj - arma::conv_to<double>::from(
            rhoj_yS.row(c) * Kj_inv * arma::trans(rhoj_yS.row(c)) );
          if(Rj < 0){ Rj = 0; } // numerical zero, coord is in S
          Rdiag(r) = sqrt(Ri * Rj);
        } 
      }
      
    } 
  }
  
  if(diag_only){
    arma::mat maincomp = arma::sum((rhoi_xS * Li_inv.t())%(rhoj_yS * Lj_inv.t()), 1);
    return maincomp + Rdiag;//R.diag();
  } else {
    return rhoi_xS * Li_inv.t() * Lj_inv * rhoj_yS.t() + R;
  }
}


//' @title Scaling factor of IOX.
//' @description This function computes the scaling factor used to compute the cross-covariance at zero distance for IOX.
//'
//' @param dag an object returned from `spiox::dag_vecchia` or similar
//' @param S a matrix of reference coordinates
//' @param theta a matrix of dimension (4,q) where each column are the marginal covariance parameters of the corresponding variable
//' @param matern an integer. Options are 0: power exponential, 1: matern, 2: wave
//' @param n_threads integer number of threads for multi-threaded operations
//'
//' @return A matrix of dimension (q,q) filled in its lower-triangular portion. The (i,j) element is the scaling factor for \eqn{C_{ij}}.
//' The IOX cross-covariance at zero distance will be \eqn{Sigma_{ij}} multiplied by this scaling factor.
//'
//' @details The function is designed for computing the IOX scaling factor used to compute the cross-covariances \eqn{C_{ij}}.
//'
//' @export
// [[Rcpp::export]]
arma::mat sfact(const arma::field<arma::uvec>& dag, const arma::mat& S, 
                const arma::mat& theta, int matern=1, int n_threads=1){
  int n = S.n_rows;
  int q = theta.n_cols;
  arma::mat I1 = arma::eye(n, n);
  arma::mat result = arma::zeros(q,q);
  
  int dag_opts = 0;

  std::vector<DagGP> daggps(q);
  arma::cube Hinvs(n,n,q);
  for(int i=0; i<q; i++){
    daggps[i] = DagGP(S, theta.col(i), dag, dag_opts, matern, n_threads);
    Hinvs.slice(i) = daggps[i].H_solve_A(I1);
    for(int j=0; j<=i; j++){
      result(i,j) = arma::accu(Hinvs.slice(i)%Hinvs.slice(j))/(n+.0);
    }
  }
  
  return(result);
}

// [[Rcpp::export]]
arma::cube Sigma_x_sfact_cpp(const arma::field<arma::uvec>& dag, const arma::mat& S, 
                const arma::cube& Sigma, const arma::cube& theta, int matern=1, int n_threads=1){
  int n = S.n_rows;
  int q = theta.n_cols;
  arma::mat I1 = arma::eye(n, n);
  int mcmc = theta.n_slices;
  arma::cube result_cube = arma::zeros(q,q, mcmc);
  
  std::vector<DagGP> daggps(q);
  arma::cube Hinvs(n,n,q);
  
  int dag_opts=0;

  for(int i=0; i<q; i++){
    arma::mat theta0 = theta.slice(0);
    arma::vec thetai = theta0.col(i);
    daggps[i] = DagGP(S, thetai, dag, dag_opts, matern, n_threads);
    Hinvs.slice(i) = daggps[i].H_solve_A(I1);
  }
  arma::mat theta_curr = theta.slice(0);
  
  for(int m=0; m<mcmc; m++){
    arma::mat result = arma::zeros(q,q);
    const arma::mat thetam = theta.slice(m);
    
    // update theta-dependent state once per m if needed
    bool theta_changed = !arma::approx_equal(theta_curr, thetam, "absdiff", 1e-10);
    if (theta_changed) theta_curr = thetam;
    
    for(int i=0; i<q; i++){
      arma::mat thetam = theta.slice(m);
      if(theta_changed){
        // update
        theta_curr = thetam;
        arma::vec thetai = theta_curr.col(i);
        daggps[i] = DagGP(S, thetai, dag, dag_opts, matern, n_threads);
        Hinvs.slice(i) = daggps[i].H_solve_A(I1);
      }
      for(int j=0; j<=i; j++){
        double v = arma::accu( Sigma(i,j,m) * (Hinvs.slice(i) % Hinvs.slice(j)) ) / double(n);
        result(i,j) = v;
        result(j,i) = v;  
      }
    }
    result_cube.slice(m) = result;
  }
  return(result_cube);
}


arma::mat rvec(const arma::mat& x, int i, 
                      const arma::mat& S, 
                      const arma::mat& theta, 
                      bool matern = true){
  
  arma::vec thetai = theta.col(i-1);
  arma::mat rhoi_xS = Correlationc(x, S, thetai, matern, false);
  
  arma::mat Ki = Correlationc(S, S, thetai, matern, true);
  arma::mat Li = arma::chol(Ki, "lower");
  arma::mat Li_inv = arma::inv(arma::trimatl(Li));
  arma::mat Ki_inv = Li_inv.t() * Li_inv;
  
  arma::vec R = arma::zeros(x.n_rows);
  arma::mat fixcoord = arma::zeros(1, 2);
  double atZerodisti = Correlationc(fixcoord, fixcoord, thetai, matern, true)(0,0);
  
  for(unsigned int r=0; r<x.n_rows; r++){
    double Ri = atZerodisti - arma::conv_to<double>::from(
      rhoi_xS.row(r) * Ki_inv * arma::trans(rhoi_xS.row(r)) );
    if(Ri < 0){ Ri = 0; } // numerical zero, coord is in S
    R(r) = sqrt(Ri);
  } 
  return R;
}


arma::mat iox_precomp(const arma::mat& x, const arma::mat& y, int i, int j,
              const arma::mat& S, 
              const arma::mat& theta, 
              const arma::field<arma::mat>& L_invs,
              bool matern = true,
              bool diag_only=false,
              bool D_only=false){
  
  arma::vec thetai = theta.col(i-1);
  arma::mat rhoi_xS = Correlationc(x, S, thetai, matern, false);
  
  arma::vec thetaj = theta.col(j-1);
  arma::mat rhoj_yS = Correlationc(y, S, thetaj, matern, false);
  
  arma::mat Ki_inv = L_invs(i-1).t() * L_invs(i-1);
  arma::mat Kj_inv = L_invs(j-1).t() * L_invs(j-1);
  
  arma::mat R = arma::zeros(x.n_rows, y.n_rows);
  arma::mat fixcoord = arma::zeros(1, 2);
  double atZerodisti = Correlationc(fixcoord, fixcoord, thetai, matern, true)(0,0);
  double atZerodistj = Correlationc(fixcoord, fixcoord, thetaj, matern, true)(0,0);
  for(unsigned int r=0; r<x.n_rows; r++){
    double Ri = atZerodisti - arma::conv_to<double>::from(
      rhoi_xS.row(r) * Ki_inv * arma::trans(rhoi_xS.row(r)) );
    if(Ri < 0){ Ri = 0; } // numerical zero, coord is in S
    for(unsigned int c=0; c<y.n_rows; c++){
      double dist = arma::accu(abs( x.row(r) - y.row(c) ));
      if(dist == 0){
        double Rj = atZerodistj - arma::conv_to<double>::from(
          rhoj_yS.row(c) * Kj_inv * arma::trans(rhoj_yS.row(c)) );
        if(Rj < 0){ Rj = 0; } // numerical zero, coord is in S
        R(r,c) = sqrt(Ri * Rj);
      }
    }
  } 
  if(D_only){
    return R;
  }
  
  if(diag_only){
    return arma::sum((rhoi_xS * L_invs(i-1).t())%(rhoj_yS * L_invs(j-1).t()), 1) + R.diag();
  } else {
    return rhoi_xS * L_invs(i-1).t() * L_invs(j-1) * rhoj_yS.t() + R;
  }
}

//[[Rcpp::export]]
arma::uvec make_ix(int q, int n){
  return arma::regspace<arma::uvec>(0, n, n*q-1);
}

// IOX cross-covariance matrix function core
//[[Rcpp::export]]
arma::mat iox_mat(const arma::mat& x, const arma::mat& y, 
                  const arma::mat& S, const arma::mat& theta, 
                  bool matern = true, bool D_only=false){
  
  int nx = x.n_rows;
  int ny = y.n_rows;
  int q = theta.n_cols;
  
  arma::field<arma::mat> L_invs(q);
  for(unsigned int i=0; i<q; i++){
    arma::vec thetai = theta.col(i);
    arma::mat Ki = Correlationc(S, S, thetai, matern, true);
    arma::mat Li = arma::chol(Ki, "lower");
    L_invs(i) = arma::inv(arma::trimatl(Li));
  }
  
  arma::mat result = arma::zeros(nx * q, ny * q);
  for(unsigned int ni=0; ni<nx; ni++){
    for(unsigned int nj=0; nj<ny; nj++){
      for(unsigned int i=0; i<q; i++){
        for(unsigned int j=0; j<q; j++){
          int ir = i+1;
          int jr = j+1;
          arma::mat ijres = iox_precomp(x.row(ni), y.row(nj), ir, jr, S, theta, L_invs, matern, false, D_only);
          result(i * nx + ni, j * ny + nj) = ijres(0, 0);
        }
      }
    }
  }
  return(result);
}

