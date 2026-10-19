
#ifndef DAGGP
#define DAGGP

#include <RcppArmadillo.h>
#include <Eigen/SparseCore>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <string>

#include "omp_import.h"
#include "covariance.h"

using namespace std;

// DagGP no longer stores the sparse H or Ci as members.  Vecchia factors
// live as the row-wise (hrows / ax / sqrtR / h) representation, which is
// what the operator interfaces (H_times_A / Ht_times_A / H_solve_A /
// Ht_solve_A) consume.  Callers that genuinely need the sparse matrix
// (R-exported daggp_build, the latent_model=2 single-site sampler,
// PPCG sparse solves, etc.) build one on demand via make_H() / make_Ci().
class DagGP {
public:
  int nr;
  arma::mat coords;
  arma::vec theta;

  // Optional per-location nugget added to the diagonal of the covariance
  // blocks consumed by compute_comps.  Empty (default) means "no nugget".
  arma::vec nugget;

  double precision_logdeterminant;
  double logdens(const arma::vec& x);
  void update_theta(const arma::vec& newtheta);

  arma::field<arma::uvec> dag;

  // compute and store markov blanket (= parents ∪ children ∪ co-parents
  // of children); built directly from `dag` without going through H/Ci.
  arma::field<arma::uvec> mblanket;
  // coloring (for the latent_model=2 single-site sampler)
  arma::field<arma::uvec> colors;

  // Row-wise H representation: nonzeros of H row i are at columns ax(i)
  // with values hrows(i).  ax(i) = [i, dag(i)...]; the diagonal entry
  // 1/sqrtR(i) is at position 0 of hrows(i).
  arma::field<arma::uvec> ax;
  arma::field<arma::vec> hrows;
  arma::vec sqrtR;
  arma::field<arma::vec> h;

  // Eigen mirror of H (and its transpose) for the four hot-path operators
  // H_times_A / Ht_times_A / H_solve_A / Ht_solve_A.  Column-major sparse so
  // the triangular solves (forward on Lower H, back on Upper H^T) take the
  // column-traversal path Eigen specialises for.  Rebuilt at the end of
  // every compute_comps() — cost is O(nnz), amortised against the Vecchia
  // factorisation that just ran in the same call.
  Eigen::SparseMatrix<double> H_eigen;    // n x n, lower triangular
  Eigen::SparseMatrix<double> Ht_eigen;   // n x n, upper triangular (= H_eigen^T)
  void build_H_eigen();

  void compute_comps();

  // Operator-style: all four dispatch to the Eigen mirror built at the end
  // of compute_comps().  The use_spmat parameter is retained for ABI
  // compatibility with existing call sites but is ignored — the Eigen path
  // is always taken.
  arma::mat H_times_A(const arma::mat& A, bool use_spmat=false) const;
  arma::mat Ht_times_A(const arma::mat& A, bool use_spmat=false) const;
  arma::mat H_solve_A(const arma::mat& A, bool use_spmat=false) const;
  arma::mat Ht_solve_A(const arma::mat& A, bool use_spmat=false) const;

  // Build the lower-triangular Vecchia precision Cholesky H on demand,
  // from sqrtR/h/dag/cache_map.  O(n·m) cost.  Use sparingly — only in
  // sites that need structural sparse access (R export, the latent_model=2
  // single-site sampler, PPCG sparse solves).
  arma::sp_mat make_H() const;

  // diag(H^T H) computed in O(n·m) without materialising H.  Used by
  // Jacobi preconditioners that need the per-column squared norm of H.
  arma::vec H_col_squared_norms() const;

  // info about covariance model:
  int matern; // 0: pexp; 1: matern; 2: wave
  double * bessel_ws;

  DagGP(){};

  // dag_opts:  0 = non-gridded (no cache amortisation, each node computes
  //                its own parent block);
  //           -1 = gridded coords + translation-invariant kernel.
  DagGP(
    const arma::mat& coords_in,
    const arma::vec& theta_in,
    const arma::field<arma::uvec>& custom_dag,
    int dag_opts=0,
    int covariance_matern=1,
    int num_threads_in=1);

  // utils
  arma::uvec oneuv;
  int n_threads;

  // Gridded amortisation: when true, compute_comps reuses Pinv / h / CPC
  // across nodes that share a parent-offset pattern.
  bool gridded;
  arma::umat cache_map;
  arma::field<arma::uvec> dag_cache;
  arma::uvec child_cache;
  void build_grid_exemplars();
  void color_from_mblanket();
  // Compute the Markov blanket from the DAG (parents + children +
  // co-parents-of-children) without materialising H or Ci.
  void build_mblanket_from_dag();

  arma::mat Corr_export(const arma::mat& cx, const arma::uvec& ix, const arma::uvec& jx, int matern, bool same);
};


static inline std::vector<double> sorted_unique_col(const arma::vec& x) {
  arma::vec u = arma::unique(x);
  std::vector<double> v(u.begin(), u.end());
  std::sort(v.begin(), v.end());
  return v;
}

static inline arma::uvec indexify_col_exact(const arma::vec& x,
                                            const std::vector<double>& uniq_sorted) {
  arma::uvec out(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    double val = x[i];
    auto it = std::lower_bound(uniq_sorted.begin(), uniq_sorted.end(), val);
    if (it == uniq_sorted.end() || *it != val) {
      Rcpp::stop("Grid indexing failed: non-gridded or jittered coordinate.");
    }
    out[i] = static_cast<arma::uword>(std::distance(uniq_sorted.begin(), it));
  }
  return out;
}

static inline std::string key_from_offsets(const arma::Mat<arma::uword>& idx,
                                           arma::uword node,
                                           const arma::uvec& parents) {
  if (parents.n_elem == 0) return std::string("[]");
  std::string s; s.reserve(parents.n_elem * (idx.n_cols * 4 + 1));
  for (arma::uword t = 0; t < parents.n_elem; ++t) {
    arma::uword p = parents[t];
    for (arma::uword k = 0; k < idx.n_cols; ++k) {
      long long off = (long long)idx(p, k) - (long long)idx(node, k);
      s += std::to_string(off);
      if (k + 1 < idx.n_cols) s.push_back(',');
    }
    if (t + 1 < parents.n_elem) s.push_back('|');
  }
  return s;
}


#endif
