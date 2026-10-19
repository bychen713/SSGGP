// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using Eigen::MatrixXd;

// Per-feature second moments of the latent process,
//
//   u_l = (1/T) sum_tau W_l^(tau) W_l^(tau)^T,
//
// in the raw W scale: CM3 varies h(X_l; iota), so the statistic must not carry
// h. The graphical lasso input is recovered from it as
// S_jk = (1/L) sum_l M[l,j] M[l,k] u_l,jk (suffstat_to_s_cpp).
//
// W is (P * L) x T with the process index fastest inside each feature block.
// Each u_l is stored as its packed lower triangle, column-major:
// (0,0), (1,0), ..., (P-1,0), (1,1), ..., (P-1,P-1).

// [[Rcpp::export]]
Rcpp::NumericMatrix feature_moments_cpp(const MatrixXd& W, const int L) {
  const int rows = W.rows();
  const int T    = W.cols();

  const int P    = rows / L;
  const int ntri = P * (P + 1) / 2;

  Rcpp::NumericMatrix out(L, ntri);
  const double inv_T = 1.0 / static_cast<double>(T);

  // copy the feature's rows into a contiguous block so the Gram matrix is a
  // dense rank update
  MatrixXd blk(P, T);
  MatrixXd G(P, P);
  for (int l = 0; l < L; ++l) {
    blk = W.middleRows(l * P, P);
    G.setZero();
    G.selfadjointView<Eigen::Lower>().rankUpdate(blk, inv_T);
    // pack the lower triangle, column major within the feature
    int k = 0;
    for (int j = 0; j < P; ++j) {
      for (int i = j; i < P; ++i, ++k) {
        out(l, k) = G(i, j);
      }
    }
  }
  return out;
}
