// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "utils.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MappedSparseMatrix;

using Sp   = Eigen::SparseMatrix<double>;
using Chol = Eigen::SimplicialLLT<Sp, Eigen::Lower, Eigen::NaturalOrdering<int>>;

// Gibbs sampler for the latent W, one feature block at a time:
//
//   W_l | W_-l, Y ~ N( Q_add^-1 A_l^T D^-1 (Y - A_-l W_-l),  Q_add^-1 ),
//   Q_add = Sigma_{W_l}^-1 + A_l^T D^-1 A_l.
//
// Column 0 of W_vals_sample is the starting state; column mc holds the state
// after sweep mc + 1 (the start is overwritten by the first sweep).
//
// Dense branches (c_vals supplied). With the residual Y - AW viewed as a
// P x n_loc matrix (Y is process-fastest, so this is a reshape),
//
//   A_l^T D^-1 (Y - AW + A_l W_l_old) = [resid c_l + W_l_old * C_l] / tau,
//
// where C_l = ||c_l||^2 when colocated and the masked column sums
// C_{pl} = sum_{i in O_p} c_{il}^2 when misaligned. Unobserved entries of the
// residual are held at zero by the mask, which is exactly the row selection
// the sparse A_l performs. The residual update is a rank-1 subtraction.
//
// General branch (c_vals empty): the sparse A_obs / A_l_list reference, used
// only under options(SSGGP.force_general = TRUE) by the parity tests.
// [[Rcpp::export]]
Rcpp::List gibbs_W_cpp(
  MatrixXd W_vals_sample,                     // (L * P) x MC
  const Rcpp::List& updated_list,             // Q_W, diagonal_entries
  const VectorXd& D_vector,                   // per-observation nugget
  const int MC,
  const int L,
  const VectorXd& Y_reordered,                // length n_obs (or n_loc * P)
  const MappedSparseMatrix<double>& A_obs,    // general branch only
  const Rcpp::List& A_l_list,                 // general branch only
  const MatrixXd& c_vals,                     // n_loc x L; 0x0 = general
  const MatrixXd& obs_mask,                   // P x n_loc; 0x0 = colocated
  const MatrixXd& col_ss_mat,                 // P x L masked column sums
  const VectorXd& tau                         // length P, per-process nugget
) {
  Rcpp::RNGScope scope;

  const int WP_rows = W_vals_sample.rows();
  const int P = WP_rows / L;
  const int n = Y_reordered.size();

  const bool dense  = (c_vals.size() > 0);
  const bool masked = dense && (obs_mask.size() > 0);

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

  VectorXd AW;
  const int n_loc = dense ? n / P : 0;
  MatrixXd resid;
  VectorXd tau_inv, css;
  if (dense) {
    tau_inv = tau.cwiseInverse();
    Eigen::Map<const MatrixXd> Wmat0(W_vals_sample.data(), P, L);
    resid = Eigen::Map<const MatrixXd>(Y_reordered.data(), P, n_loc)
          - Wmat0 * c_vals.transpose();
    if (masked) {
      resid.array() *= obs_mask.array();
    } else {
      css.resize(L);
      for (int l = 0; l < L; ++l) css[l] = c_vals.col(l).squaredNorm();
    }
  } else {
    AW = A_obs * W_vals_sample.col(0);
  }

  for (int mc = 0; mc < MC; ++mc) {
    const int idx_prev = (mc == 0) ? 0 : (mc - 1);
    for (int l = 0; l < L; ++l) {
      const int row_start = l * P;
      const VectorXd W_l_old = W_vals_sample.block(row_start, idx_prev, P, 1);

      // A_l^T D^-1 (Y - AW + A_l W_l_old)
      VectorXd info(P);
      if (dense) {
        info.noalias() = resid * c_vals.col(l);
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

      // Q_add = L L^T and matrixU() is L^T, so U^-1 z ~ N(0, Q_add^-1)
      const VectorXd mean    = chol_QW_add[l].solve(info);
      const VectorXd W_l_new = mean + chol_QW_add[l].matrixU().solve(rnorm_vec(P));
      W_vals_sample.block(row_start, mc, P, 1) = W_l_new;

      const VectorXd delta_W = W_l_new - W_l_old;
      if (dense) {
        if (masked) {
          for (int i = 0; i < n_loc; ++i) {
            resid.col(i).array() -=
              obs_mask.col(i).array() * delta_W.array() * c_vals(i, l);
          }
        } else {
          resid.noalias() -= delta_W * c_vals.col(l).transpose();
        }
      } else {
        MappedSparseMatrix<double> A_l =
          Rcpp::as<MappedSparseMatrix<double>>(A_l_list[l]);
        AW.noalias() += A_l * delta_W;
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("W_vals_sample") = W_vals_sample);
}
