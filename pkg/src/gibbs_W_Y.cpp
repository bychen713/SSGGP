// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "utils.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MappedSparseMatrix;

using Sp   = Eigen::SparseMatrix<double>;
using Chol = Eigen::SimplicialLLT<Sp, Eigen::Lower, Eigen::NaturalOrdering<int>>;

// Gibbs sampler for W given the training data, recording the latent prediction
// A_test W after every sweep. The sweep is the one in gibbs_W.cpp; see there
// for the dense identities. Test locations are always fully observed, so
// A_test = c_test (x) I_P.
//
// W_vals_sample is a two-column scratch buffer (current and previous state);
// only Y_pred is returned.
// [[Rcpp::export]]
Rcpp::List gibbs_W_Y_cpp(
  Eigen::Map<MatrixXd> W_vals_sample,          // (L * P) x 2
  const VectorXd& Y_reordered,                 // length nP
  const Rcpp::List& updated_list,              // Q_W, diagonal_entries
  const VectorXd& D_vector,                    // per-observation nugget
  const int MC,
  const int L,
  const MappedSparseMatrix<double>& A_obs,     // general branch only
  const MappedSparseMatrix<double>& A_test,    // general branch only
  const Rcpp::List& A_l_list,                  // general branch only
  const MatrixXd& c_obs,                       // n_loc x L; 0x0 = general
  const MatrixXd& c_test,                      // n_test_loc x L
  const MatrixXd& obs_mask,                    // P x n_loc; 0x0 = colocated
  const MatrixXd& col_ss_mat,                  // P x L masked column sums
  const VectorXd& tau                          // length P, per-process nugget
) {
  Rcpp::RNGScope scope;

  const int WP_rows = W_vals_sample.rows();
  if (W_vals_sample.cols() != 2) Rcpp::stop("W_vals_sample must have 2 columns");
  const int P = WP_rows / L;
  const int n = Y_reordered.size();

  const bool dense  = (c_obs.size() > 0);
  const bool masked = dense && (obs_mask.size() > 0);
  const int n_test  = dense ? static_cast<int>(c_test.rows()) * P
                            : static_cast<int>(A_test.rows());

  Rcpp::List Q_W_list          = updated_list["Q_W"];
  Rcpp::List diag_entries_list = updated_list["diagonal_entries"];

  std::vector<Chol> chol_QW_add(L);
  for (int l = 0; l < L; ++l) {
    const Sp QW_add = Rcpp::as<Sp>(Q_W_list[l]) + Rcpp::as<Sp>(diag_entries_list[l]);
    chol_QW_add[l].compute(QW_add);
    if (chol_QW_add[l].info() != Eigen::Success) {
      Rcpp::stop("Cholesky failed for Q_W_added at block " + std::to_string(l + 1));
    }
  }

  MatrixXd Y_pred(n_test, MC);

  VectorXd AW;
  const int n_loc = dense ? n / P : 0;
  MatrixXd resid;
  VectorXd tau_inv, css;
  if (dense) {
    tau_inv = tau.cwiseInverse();
    Eigen::Map<const MatrixXd> Wmat0(W_vals_sample.data(), P, L);
    resid = Eigen::Map<const MatrixXd>(Y_reordered.data(), P, n_loc)
          - Wmat0 * c_obs.transpose();
    if (masked) {
      resid.array() *= obs_mask.array();
    } else {
      css.resize(L);
      for (int l = 0; l < L; ++l) css[l] = c_obs.col(l).squaredNorm();
    }
  } else {
    AW = A_obs * W_vals_sample.col(0);
  }

  for (int mc = 0; mc < MC; ++mc) {
    // read the state written at mc - 1 and write the other column; at mc == 0
    // both are column 0, an in-place update of the initial state
    const int idx_prev  = (mc == 0) ? 0 : ((mc - 1) % 2);
    const int idx_store = mc % 2;

    for (int l = 0; l < L; ++l) {
      const int row_start = l * P;
      const VectorXd W_l_old = W_vals_sample.block(row_start, idx_prev, P, 1);

      // t(A_l) D^-1 (Y - AW + A_l W_l_old)
      VectorXd info(P);
      if (dense) {
        info.noalias() = resid * c_obs.col(l);
        if (masked) {
          info.array() += W_l_old.array() * col_ss_mat.col(l).array();
        } else {
          info += W_l_old * css[l];
        }
        info.array() *= tau_inv.array();
      } else {
        MappedSparseMatrix<double> A_l =
          Rcpp::as<MappedSparseMatrix<double>>(A_l_list[l]);
        VectorXd v = Y_reordered - AW + A_l * W_l_old;
        v.array() /= D_vector.array();
        info = A_l.transpose() * v;
      }

      const VectorXd mean    = chol_QW_add[l].solve(info);
      const VectorXd W_l_new = mean + chol_QW_add[l].matrixU().solve(rnorm_vec(P));
      W_vals_sample.block(row_start, idx_store, P, 1) = W_l_new;

      const VectorXd delta_W = W_l_new - W_l_old;
      if (dense) {
        if (masked) {
          for (int i = 0; i < n_loc; ++i) {
            resid.col(i).array() -=
              obs_mask.col(i).array() * delta_W.array() * c_obs(i, l);
          }
        } else {
          resid.noalias() -= delta_W * c_obs.col(l).transpose();
        }
      } else {
        MappedSparseMatrix<double> A_l =
          Rcpp::as<MappedSparseMatrix<double>>(A_l_list[l]);
        AW.noalias() += A_l * delta_W;
      }
    }

    // A_test W after this sweep; dense: vec(Wmat c_test^T)
    if (dense) {
      Eigen::Map<const MatrixXd> Wmat(
        W_vals_sample.data() + static_cast<std::size_t>(idx_store) * WP_rows, P, L);
      const MatrixXd Mtest = Wmat * c_test.transpose();
      Y_pred.col(mc) = Eigen::Map<const VectorXd>(Mtest.data(), n_test);
    } else {
      Y_pred.col(mc) = A_test * W_vals_sample.col(idx_store);
    }
  }

  return Rcpp::List::create(Rcpp::Named("Y_pred") = Y_pred);
}
