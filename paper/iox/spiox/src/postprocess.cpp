#include <RcppArmadillo.h>
#include "daggp.h"

using namespace std;

//[[Rcpp::export]]
arma::cube S_to_Sigma(const arma::cube& S){
  // S is upper cholesky factor of sigma
  arma::cube Sigma(arma::size(S));
  int mcmc = S.n_slices;
  for(int i=0; i<mcmc; i++){
    Sigma.slice(i) = S.slice(i).t() * S.slice(i);
  }
  return Sigma;
}

//[[Rcpp::export]]
arma::cube S_to_Q(const arma::cube& S){
  // S is upper cholesky factor of sigma
  arma::cube Q(arma::size(S));
  int mcmc = S.n_slices;
  for(int i=0; i<mcmc; i++){
    arma::mat Si = arma::inv(arma::trimatu(S.slice(i)));
    Q.slice(i) = Si * Si.t();
  }
  return Q;
}

//[[Rcpp::export]]
arma::cube Sigma_to_correl(const arma::cube& Sigma){
  int q = Sigma.n_rows;
  int k = Sigma.n_cols;
  int m = Sigma.n_slices;
  
  arma::cube Omega = arma::zeros(arma::size(Sigma));
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(int i=0; i<m; i++){
    arma::mat dllt = arma::diagmat(1.0/sqrt( Sigma.slice(i).diag() ));
    Omega.slice(i) = dllt * Sigma.slice(i) * dllt;
  }
  return Omega;
}

//[[Rcpp::export]]
arma::cube Sigma_identify(const arma::cube& Sigma, const arma::cube& theta){
  
  arma::cube Sigma_id = arma::cube(arma::size(Sigma));
  for(unsigned int m = 0; m<Sigma.n_slices; m++){
    
    arma::mat theta_local = theta.slice(m);
    arma::mat Ider = arma::diagmat(sqrt(theta_local.row(1)));
    
    Sigma_id.slice(m) = Ider * Sigma.slice(m) * Ider;
  }
  return Sigma_id;
}


// iox cross covariance computed at pairs of locs of reference set
//[[Rcpp::export]]
Rcpp::List iox_make_fij(int i, int j,
                        const arma::mat& coords,
                        const arma::field<arma::uvec>& custom_dag,
                        int dag_opts, 
                        const arma::cube& theta, 
                        int cov_model_matern,
                        int num_threads,
                        int n_bins = 1000,
                        double max_range = 0.5,
                        bool silent=true){
  
  int n = coords.n_rows;
  int q = theta.n_cols;
  int mcmc = theta.n_slices;
  
  // Assign each pair to a uniform bin in [0, max_range]
  double bin_width = max_range / n_bins;
  arma::imat bin_id(n, n);
  arma::ivec pair_count(n_bins, arma::fill::zeros);
  
  for (int s1 = 0; s1 < n; s1++) {
    for (int s2 = s1; s2 < n; s2++) {
      double d = arma::norm(coords.row(s1) - coords.row(s2));
      if (d > max_range) {
        bin_id(s1, s2) = -1;
        bin_id(s2, s1) = -1;
      } else {
        int b = (int)(d / bin_width);
        if (b >= n_bins) b = n_bins - 1;  // edge case: d == max_range
        bin_id(s1, s2) = b;
        bin_id(s2, s1) = b;
        pair_count(b) += (s1 == s2) ? 1 : 2;
      }
    }
  }
  
  // Bin midpoints as reported distances; drop empty bins
  std::vector<int> kept_bins;
  for (int b = 0; b < n_bins; b++) {
    if (pair_count(b) > 0) kept_bins.push_back(b);
  }
  int n_kept = kept_bins.size();
  
  // Map original bin index -> compact index; -1 for empty
  arma::ivec bin_remap(n_bins);
  bin_remap.fill(-1);
  arma::vec unique_dists(n_kept);
  arma::ivec kept_pair_count(n_kept);
  for (int k = 0; k < n_kept; k++) {
    bin_remap(kept_bins[k]) = k;
    unique_dists(k) = (kept_bins[k] + 0.5) * bin_width;
    kept_pair_count(k) = pair_count(kept_bins[k]);
  }
  
  // Remap bin_id to compact indices
  for (int s1 = 0; s1 < n; s1++) {
    for (int s2 = s1; s2 < n; s2++) {
      int b = bin_id(s1, s2);
      int new_b = (b >= 0) ? bin_remap(b) : -1;
      bin_id(s1, s2) = new_b;
      bin_id(s2, s1) = new_b;
    }
  }
  
  arma::mat f_ij(n_kept, mcmc, arma::fill::zeros);
  
  int n_batches = 10;
  int batch_size = mcmc / n_batches;
  int remainder = mcmc % n_batches;
  
  for (int b = 0; b < n_batches; b++) {
    int m_start = b * batch_size + std::min(b, remainder);
    int m_end = m_start + batch_size + (b < remainder ? 1 : 0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (int m = m_start; m < m_end; m++) {
      
      arma::mat theta_m = theta.slice(m);
      
      std::vector<arma::mat> L(q);
      for (int ii = 0; ii < q; ii++) {
        arma::vec theta_m_j = theta_m.col(ii);
        DagGP daggp(coords, theta_m_j, custom_dag,
                    dag_opts,
                    cov_model_matern,
                    1);

        L[ii] = daggp.H_solve_A(arma::eye(n, n));
      }
      
      arma::mat M = L[i] * L[j].t();
      
      arma::vec bin_sum(n_kept, arma::fill::zeros);
      for (int s1 = 0; s1 < n; s1++) {
        for (int s2 = s1; s2 < n; s2++) {
          int idx = bin_id(s1, s2);
          if (idx < 0) continue;
          if (s1 == s2) {
            bin_sum(idx) += M(s1, s2);
          } else {
            bin_sum(idx) += M(s1, s2) + M(s2, s1);
          }
        }
      }
      
      for (int h = 0; h < n_kept; h++) {
        f_ij(h, m) = bin_sum(h) / kept_pair_count(h);
      }
    }
    
    if(!silent){
      Rcpp::Rcout << "Batch " << (b + 1) << "/" << n_batches
                  << " done (" << m_end << "/" << mcmc << " iterations)\n";
    }
    
    Rcpp::checkUserInterrupt();
  }
  
  return Rcpp::List::create(
    Rcpp::Named("dists") = unique_dists,
    Rcpp::Named("f_ij") = f_ij
  );
}

//[[Rcpp::export]]
arma::cube iox_make_fij0(const arma::mat& coords,
                         const arma::field<arma::uvec>& custom_dag,
                         int dag_opts,
                         const arma::cube& theta,
                         int cov_model_matern,
                         int num_threads){
  
  int n = coords.n_rows;
  int q = theta.n_cols;
  int mcmc = theta.n_slices;
  
  arma::cube result(q, q, mcmc);
  arma::mat prev_theta;
  arma::mat F(q, q);
  
  for (int m = 0; m < mcmc; m++) {
    arma::mat theta_m = theta.slice(m);
    
    if (m == 0 || !arma::approx_equal(theta_m, prev_theta, "absdiff", 1e-12)) {
      std::vector<arma::mat> L(q);
      for (int j = 0; j < q; j++) {
        arma::vec theta_j = theta_m.col(j);
        DagGP daggp(coords, theta_j, custom_dag,
                    dag_opts, cov_model_matern,
                    num_threads);
        L[j] = daggp.H_solve_A(arma::eye(n, n));
      }
      
      for (int i = 0; i < q; i++) {
        for (int j = i; j < q; j++) {
          double s = 0.0;
          for (int s1 = 0; s1 < n; s1++) {
            s += arma::dot(L[i].row(s1), L[j].row(s1));
          }
          F(i, j) = s / n;
          F(j, i) = F(i, j);
        }
      }
      
      prev_theta = theta_m;
    }
    
    result.slice(m) = F;
  }
  
  return result;
}

