// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::MappedSparseMatrix;

// Marginal log-likelihood of Y after integrating out the latent W, for a fixed
// draw of random Fourier features.
//
//   Y ~ N(0, Sigma_L + D),  Sigma_L = A Sigma_W A^T,
//   Sigma_W = block-diagonal(Sigma_{W_l}),  Sigma_{W_l}^{-1} = D_l sigma^-1 D_l
//
// with D_l = diag(M[l, ]). Evaluated by the Woodbury identity and the matrix
// determinant lemma, so the cost is driven by an (L p) x (L p) sparse Cholesky
// rather than an (n p) x (n p) dense one:
//
//   log det(Sigma_L + D) = log det(Sigma_W^-1 + A^T D^-1 A)
//                          - log det(Sigma_W^-1) + log det(D)
//   Y^T (Sigma_L + D)^-1 Y = Y^T D^-1 Y
//                            - u^T (Sigma_W^-1 + A^T D^-1 A)^-1 u,  u = A^T D^-1 Y
//
// M(l, j) = sqrt( g(X_l) / (2 f_j(X_l)) ), as in update_matrices_cpp. A_obs
// holds the observed rows only.

// [[Rcpp::export]]
double ll_cpp(
    const MatrixXd&                   M,          // L x P
    const MappedSparseMatrix<double>& sigma_inv,  // P x P
    const MappedSparseMatrix<double>& A_obs,      // n_obs x (L * P)
    const VectorXd&                   Y_obs,      // n_obs
    const VectorXd&                   d_obs       // n_obs, strictly positive
) {
  const int L = M.rows();
  const int P = M.cols();
  const int n_obs = Y_obs.size();

  // sigma^-1 = Lambda Lambda^T
  Eigen::SimplicialLLT<SparseMatrix<double>, Eigen::Lower,
                       Eigen::NaturalOrdering<int> > chol_sigma_inv;
  chol_sigma_inv.compute(sigma_inv);
  if (chol_sigma_inv.info() != Eigen::Success) {
    Rcpp::stop("Cholesky of sigma_inv failed; is it positive definite?");
  }
  const SparseMatrix<double> Lambda = chol_sigma_inv.matrixL();

  // L_mat = block-diagonal( diag(M[l, ]) Lambda ), so that
  // L_mat L_mat^T = block-diagonal( Sigma_{W_l}^{-1} ) = Sigma_W^{-1}
  SparseMatrix<double> L_mat(L * P, L * P);
  {
    std::vector<Eigen::Triplet<double> > trip;
    trip.reserve(static_cast<size_t>(Lambda.nonZeros()) * L);
    for (int l = 0; l < L; ++l) {
      const int off = l * P;
      for (int k = 0; k < Lambda.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(Lambda, k); it; ++it) {
          trip.emplace_back(it.row() + off, it.col() + off,
                            it.value() * M(l, it.row()));
        }
      }
    }
    L_mat.setFromTriplets(trip.begin(), trip.end());
    L_mat.makeCompressed();
  }

  const VectorXd dinv = d_obs.cwiseInverse();

  SparseMatrix<double> Sinv_ATDinvA =
    SparseMatrix<double>(L_mat * L_mat.transpose()) +
    SparseMatrix<double>(A_obs.transpose() * dinv.asDiagonal() * A_obs);

  Eigen::SimplicialLLT<SparseMatrix<double> > chol_S(Sinv_ATDinvA);
  if (chol_S.info() != Eigen::Success) {
    Rcpp::stop("Cholesky of Sigma_W^-1 + A^T D^-1 A failed");
  }
  const SparseMatrix<double> L_S = chol_S.matrixL();

  const VectorXd u = A_obs.transpose() * (dinv.array() * Y_obs.array()).matrix();
  const double quad = u.dot(chol_S.solve(u));

  return -0.5 * n_obs * std::log(2.0 * M_PI)
         - L_S.diagonal().array().log().sum()
         + L_mat.diagonal().array().log().sum()
         - 0.5 * d_obs.array().log().sum()
         - 0.5 * (Y_obs.array().square() * dinv.array()).sum()
         + 0.5 * quad;
}
