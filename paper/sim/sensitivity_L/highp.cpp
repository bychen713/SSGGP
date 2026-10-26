// Memory-conscious kernels for the p=25 and p=50 likelihood-sensitivity runs.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <R_ext/Lapack.h>

using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Eigen::Map;
using Eigen::MatrixXd;

// Add Phi Phi' to the LOWER triangle of K without allocating another N x N
// matrix. Rows are ordered as: every process at the first n_small locations,
// then every process at the remaining locations. This makes the n_small
// covariance a leading principal submatrix of the n_full covariance.
// [[Rcpp::export]]
void update_nested_covariance_inplace(NumericMatrix K,
                                      const NumericMatrix& C,
                                      const NumericMatrix& W,
                                      int n_small) {
  const int n = C.nrow(), b = C.ncol(), P = W.ncol(), N = n * P;
  if (W.nrow() != b || K.nrow() != N || K.ncol() != N ||
      n_small < 1 || n_small >= n)
    Rcpp::stop("incompatible dimensions in update_nested_covariance_inplace");

  Map<MatrixXd> Km(REAL(K), N, N);
  Map<const MatrixXd> Cm(REAL(C), n, b);
  Map<const MatrixXd> Wm(REAL(W), b, P);
  MatrixXd Phi(N, b);
  const int n_tail = n - n_small;

  for (int j = 0; j < P; ++j) {
    for (int l = 0; l < b; ++l) {
      const double w = Wm(l, j);
      for (int i = 0; i < n_small; ++i)
        Phi(j * n_small + i, l) = Cm(i, l) * w;
      for (int i = 0; i < n_tail; ++i)
        Phi(P * n_small + j * n_tail + i, l) = Cm(n_small + i, l) * w;
    }
  }
  Km.selfadjointView<Eigen::Lower>().rankUpdate(Phi, 1.0);
}

// Factor K + nugget I in place using its lower triangle. Return Gaussian
// log-determinant and quadratic-form components for both the leading n_small
// variables and the complete vector. K is intentionally destroyed.
// [[Rcpp::export]]
List chol_nested_stats_inplace(NumericMatrix K, NumericVector y,
                               int n_small, double nugget) {
  const int N = K.nrow();
  if (K.ncol() != N || y.size() != N || n_small < 1 || n_small >= N ||
      !R_finite(nugget) || nugget <= 0)
    Rcpp::stop("incompatible dimensions or nugget in chol_nested_stats_inplace");

  for (int i = 0; i < N; ++i) K(i, i) += nugget;
  char uplo = 'L';
  int info = 0, lda = N, n_full = N;
  F77_CALL(dpotrf)(&uplo, &n_full, REAL(K), &lda, &info FCONE);
  if (info != 0) Rcpp::stop("dpotrf failed with info=%d", info);

  double logdet_small = 0.0, logdet_full = 0.0;
  for (int i = 0; i < N; ++i) {
    const double term = 2.0 * std::log(K(i, i));
    logdet_full += term;
    if (i < n_small) logdet_small += term;
  }

  NumericVector rhs_small(n_small);
  for (int i = 0; i < n_small; ++i) rhs_small[i] = y[i];
  int one = 1, ldb_small = n_small, ns = n_small;
  F77_CALL(dpotrs)(&uplo, &ns, &one, REAL(K), &lda,
                   REAL(rhs_small), &ldb_small, &info FCONE);
  if (info != 0) Rcpp::stop("small dpotrs failed with info=%d", info);
  double quad_small = 0.0;
  for (int i = 0; i < n_small; ++i) quad_small += y[i] * rhs_small[i];

  NumericVector rhs_full = Rcpp::clone(y);
  int ldb_full = N;
  F77_CALL(dpotrs)(&uplo, &n_full, &one, REAL(K), &lda,
                   REAL(rhs_full), &ldb_full, &info FCONE);
  if (info != 0) Rcpp::stop("full dpotrs failed with info=%d", info);
  double quad_full = 0.0;
  for (int i = 0; i < N; ++i) quad_full += y[i] * rhs_full[i];

  return List::create(
    Rcpp::Named("logdet_small") = logdet_small,
    Rcpp::Named("quad_small") = quad_small,
    Rcpp::Named("logdet_full") = logdet_full,
    Rcpp::Named("quad_full") = quad_full);
}

// Factor a single covariance K + nugget I in place using its lower triangle.
// [[Rcpp::export]]
List chol_stats_lower_inplace(NumericMatrix K, NumericVector y,
                              double nugget) {
  const int N = K.nrow();
  if (K.ncol() != N || y.size() != N ||
      !R_finite(nugget) || nugget <= 0)
    Rcpp::stop("incompatible dimensions or nugget in chol_stats_lower_inplace");
  for (int i = 0; i < N; ++i) K(i, i) += nugget;
  char uplo = 'L';
  int info = 0, lda = N, n = N;
  F77_CALL(dpotrf)(&uplo, &n, REAL(K), &lda, &info FCONE);
  if (info != 0) Rcpp::stop("dpotrf failed with info=%d", info);
  double logdet = 0.0;
  for (int i = 0; i < N; ++i) logdet += 2.0 * std::log(K(i, i));
  NumericVector rhs = Rcpp::clone(y);
  int one = 1, ldb = N;
  F77_CALL(dpotrs)(&uplo, &n, &one, REAL(K), &lda,
                   REAL(rhs), &ldb, &info FCONE);
  if (info != 0) Rcpp::stop("dpotrs failed with info=%d", info);
  double quad = 0.0;
  for (int i = 0; i < N; ++i) quad += y[i] * rhs[i];
  return List::create(Rcpp::Named("logdet") = logdet,
                      Rcpp::Named("quad") = quad);
}

// Build I + Phi' D^-1 Phi in feature space, one column block at a time.
// C already includes the 1/sqrt(L) scaling and Phi_l = C_l (x) W_l.
// [[Rcpp::export]]
NumericMatrix build_feature_system_highp(const NumericMatrix& C,
                                         const NumericMatrix& W,
                                         const NumericVector& d,
                                         int block_size = 256) {
  const int n = C.nrow(), L = C.ncol(), P = W.ncol();
  if (W.nrow() != L || d.size() != P || block_size < 1)
    Rcpp::stop("incompatible dimensions in build_feature_system_highp");
  Map<const MatrixXd> Cm(REAL(C), n, L);
  Map<const MatrixXd> Wm(REAL(W), L, P);
  MatrixXd Wd = Wm;
  for (int j = 0; j < P; ++j) {
    if (!R_finite(d[j]) || d[j] <= 0) Rcpp::stop("d must be positive");
    Wd.col(j) /= std::sqrt(d[j]);
  }
  NumericMatrix B(L, L);
  Map<MatrixXd> Bm(REAL(B), L, L);
  for (int j0 = 0; j0 < L; j0 += block_size) {
    const int b = std::min(block_size, L - j0);
    MatrixXd G = Cm.transpose() * Cm.middleCols(j0, b);
    MatrixXd H = Wd * Wd.middleRows(j0, b).transpose();
    Bm.middleCols(j0, b) = G.cwiseProduct(H);
  }
  Bm.diagonal().array() += 1.0;
  return B;
}

// Factor a feature-space system in place and return log(det(B)) and q'B^-1q.
// [[Rcpp::export]]
List chol_feature_stats_inplace(NumericMatrix B, NumericVector q) {
  const int N = B.nrow();
  if (B.ncol() != N || q.size() != N)
    Rcpp::stop("incompatible dimensions in chol_feature_stats_inplace");
  char uplo = 'U';
  int info = 0, lda = N, n = N;
  F77_CALL(dpotrf)(&uplo, &n, REAL(B), &lda, &info FCONE);
  if (info != 0) Rcpp::stop("feature dpotrf failed with info=%d", info);
  double logdet = 0.0;
  for (int i = 0; i < N; ++i) logdet += 2.0 * std::log(B(i, i));
  NumericVector rhs = Rcpp::clone(q);
  int one = 1, ldb = N;
  F77_CALL(dpotrs)(&uplo, &n, &one, REAL(B), &lda,
                   REAL(rhs), &ldb, &info FCONE);
  if (info != 0) Rcpp::stop("feature dpotrs failed with info=%d", info);
  double quad = 0.0;
  for (int i = 0; i < N; ++i) quad += q[i] * rhs[i];
  return List::create(Rcpp::Named("logdet") = logdet,
                      Rcpp::Named("quad") = quad);
}
