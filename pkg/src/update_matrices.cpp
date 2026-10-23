// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::MappedSparseMatrix;

// Per-feature matrices for the Gibbs sampler:
//
//   Q_W[l]              = diag(M[l, ]) sigma^-1 diag(M[l, ]) = Sigma_{W_l}^-1,
//   diagonal_entries[l] = A_l^T D^-1 A_l,
//
// with M[l, j] = sqrt(g(X_l) / (2 f_j(X_l))), built in R.
//
// Every row of A_l has one nonzero (an observation belongs to one process), so
// A_l^T D^-1 A_l is diagonal. Which of `col_ss` / `col_ss_mat` is non-empty
// selects how it is formed:
//
//  * col_ss (length L), colocated: ||c_l||^2 diag(1/d_p).
//  * col_ss_mat (P x L), misaligned: col_ss_mat[p, l] / tau_p, where
//    col_ss_mat[p, l] = sum_{i observed for p} c_{il}^2.
//  * neither: the general sparse product from A_l_list, the reference used
//    only under options(SSGGP.force_general = TRUE).
//
// `d_obs` is the per-observation nugget and `tau` the per-process one.

// [[Rcpp::export]]
Rcpp::List update_matrices_cpp(
    const MatrixXd&                    M,           // L x P
    const MappedSparseMatrix<double>&  sigma_inv,   // P x P
    const VectorXd&                    d_obs,       // length n_obs
    const Rcpp::List&                  A_l_list,    // length L, each n_obs x P
    const Rcpp::NumericVector&         col_ss =
      Rcpp::NumericVector::create(),                // length L, or empty
    const Rcpp::NumericVector&         col_ss_mat =
      Rcpp::NumericVector::create(),                // P x L (col-major), or empty
    const Rcpp::NumericVector&         tau =
      Rcpp::NumericVector::create()                 // length P, or empty
) {
  const int L = M.rows();
  const int P = M.cols();
  const bool masked    = (col_ss_mat.size() > 0);
  const bool colocated = !masked && (col_ss.size() > 0);

  const VectorXd inv_sqrt_d = d_obs.array().sqrt().inverse();
  Rcpp::List Q_W(L), diagonal_entries(L);

  for (int l = 0; l < L; ++l) {
    const VectorXd m = M.row(l).transpose();

    SparseMatrix<double> QW = sigma_inv;
    for (int k = 0; k < QW.outerSize(); ++k) {
      for (SparseMatrix<double>::InnerIterator it(QW, k); it; ++it) {
        it.valueRef() *= m(it.row()) * m(it.col());
      }
    }
    QW = 0.5 * (QW + SparseMatrix<double>(QW.transpose()));
    QW.makeCompressed();
    Q_W[l] = QW;

    SparseMatrix<double> tAA(P, P);
    if (masked || colocated) {
      tAA.reserve(Eigen::VectorXi::Constant(P, 1));
      const R_xlen_t off = static_cast<R_xlen_t>(l) * P;
      for (int p = 0; p < P; ++p) {
        // colocated: d_obs is process-fastest, so its first P entries are tau
        tAA.insert(p, p) = masked ? col_ss_mat[off + p] / tau[p]
                                  : col_ss[l] / d_obs(p);
      }
    } else {
      SparseMatrix<double> A_l = Rcpp::as<MappedSparseMatrix<double> >(A_l_list[l]);
      for (int k = 0; k < A_l.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(A_l, k); it; ++it) {
          it.valueRef() *= inv_sqrt_d(it.row());
        }
      }
      tAA = A_l.transpose() * A_l;
      tAA = 0.5 * (tAA + SparseMatrix<double>(tAA.transpose()));
    }
    tAA.makeCompressed();
    diagonal_entries[l] = tAA;
  }

  return Rcpp::List::create(
    Rcpp::Named("Q_W")              = Q_W,
    Rcpp::Named("diagonal_entries") = diagonal_entries
  );
}
