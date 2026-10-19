
#include <RcppArmadillo.h>
#include "omp_import.h"


static inline arma::vec autocorr_fft(const arma::vec& x) {
  const arma::uword n = x.n_elem;
  arma::vec y = x - arma::mean(x);
  
  arma::uword m = 1;
  while (m < 2 * n) m <<= 1;
  
  arma::vec y_pad(m, arma::fill::zeros);
  y_pad.head(n) = y;
  
  arma::cx_vec Y = arma::fft(y_pad);
  arma::vec acf = arma::real(arma::ifft(Y % arma::conj(Y)));
  
  acf = acf.head(n);
  if (acf(0) > 0.0) acf /= acf(0);
  return acf;
}

static inline double ess_single(const arma::vec& x) {
  const arma::uword n = x.n_elem;
  if (n < 4) return static_cast<double>(n);
  
  // Guard against a constant chain (variance 0 -> ESS undefined; return n).
  if (arma::var(x) <= 0.0) return static_cast<double>(n);
  
  arma::vec rho = autocorr_fft(x);
  
  // Geyer's initial positive sequence: pair lags into P_k = rho[2k] + rho[2k+1]
  // and stop at the first non-positive pair.
  std::vector<double> P;
  P.reserve(n / 2);
  const arma::uword max_pairs = n / 2;
  for (arma::uword k = 0; k < max_pairs; ++k) {
    const double pk = rho(2 * k) + rho(2 * k + 1);
    if (pk <= 0.0) break;
    P.push_back(pk);
  }
  if (P.empty()) return static_cast<double>(n);
  
  // Geyer's initial monotone sequence: enforce non-increasing pair sums.
  for (std::size_t k = 1; k < P.size(); ++k) {
    if (P[k] > P[k - 1]) P[k] = P[k - 1];
  }
  
  double sumP = 0.0;
  for (double v : P) sumP += v;
  
  // tau = 1 + 2 * sum_{t>=1} rho_t  =  2 * sum_k P_k - 1
  double tau = 2.0 * sumP - 1.0;
  if (tau < 1.0) tau = 1.0;
  
  double ess = static_cast<double>(n) / tau;
  if (ess > static_cast<double>(n)) ess = static_cast<double>(n);
  return ess;
}

//' Compute ESS for every entry of an MCMC sample of N x K matrices.
//'
//' @param draws Cube of dimension N x K x T (slice t = sample t).
//' @return N x K matrix of effective sample sizes.
// [[Rcpp::export]]
arma::mat effective_sample_size(const arma::cube& draws) {
 const arma::uword N = draws.n_rows;
 const arma::uword K = draws.n_cols;
 const arma::uword T = draws.n_slices;
 
 if (T < 4) {
   Rcpp::stop("Need at least 4 MCMC iterations to estimate ESS.");
 }
 
 arma::mat out(N, K);
 
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
 for (arma::uword i = 0; i < N; ++i) {
   for (arma::uword j = 0; j < K; ++j) {
     arma::vec chain = arma::vectorise(draws.tube(i, j));
     out(i, j) = ess_single(chain);
   }
 }
 return out;
}