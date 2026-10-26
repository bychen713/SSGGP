// p=50 simulated-W accumulator for the n=500 observed-space orientation.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using Rcpp::NumericMatrix;
using Eigen::Map;
using Eigen::MatrixXd;

// Add Phi Phi' to the LOWER triangle of K. Rows are process-major and C has
// exactly the locations needed for this covariance (n=500 here).
// [[Rcpp::export]]
void update_process_covariance_inplace(NumericMatrix K,
                                       const NumericMatrix& C,
                                       const NumericMatrix& W) {
  const int n = C.nrow(), b = C.ncol(), P = W.ncol(), N = n * P;
  if (W.nrow() != b || K.nrow() != N || K.ncol() != N)
    Rcpp::stop("incompatible dimensions in update_process_covariance_inplace");
  Map<MatrixXd> Km(REAL(K), N, N);
  Map<const MatrixXd> Cm(REAL(C), n, b);
  Map<const MatrixXd> Wm(REAL(W), b, P);
  MatrixXd Phi(N, b);
  for (int j = 0; j < P; ++j)
    for (int l = 0; l < b; ++l)
      Phi.block(j * n, l, n, 1) = Cm.col(l) * Wm(l, j);
  Km.selfadjointView<Eigen::Lower>().rankUpdate(Phi, 1.0);
}
