#include "omp_import.h"
#include "sparse_solvers.h"
#include "daggp.h"
#include "ramadapt.h"

using namespace std;

int time_count(std::chrono::steady_clock::time_point tstart);

class SpIOX {
public:
  // Y = XB + Z where Z is multivariate q,q IOC-GP
  // S^T*S = Sigma = Q^-1 = (Lambda*Lambda^T + Delta)^-1. Lambda, Q are sparse
  
  // -------------- data
  
  // matrix of outcomes dim n, q 
  arma::mat Y;
  // matrix of predictors dim n, p
  arma::mat X;
  
  // metadata
  unsigned int n, q, p;
  int intercept; // location of intercept in X, if at all
  int num_threads;
  //double spatial_sparsity;
  
  // -------------- model parameters
  // objects that depend on B
  arma::mat B;
  arma::mat YXB;
  arma::mat B_Var; // prior variance on B, element by element
  
  // objects that depend on Sigma
  arma::mat S, Si, Sigma, Q; // S^T * S = Sigma = Q^-1 = (Lambda*Lambda^T + Delta)^-1 = (Si * Si^T)^-1
  
  // future:
  //arma::mat A, Aplus, AplusT; 
  
  // objects that depend on theta RadGP for spatial dependence
  std::vector<DagGP> daggps, daggps_alt;
  arma::mat theta; // each column is one alternative value for theta
  bool gridded;
  
  // objects that depend on W
  arma::mat V; 
  
  // -------------- utilities
  int matern;
  void update_B(); // 
  void update_Sigma_iwishart();
  
  void compute_V(); 
  bool upd_theta_metrop();
  arma::uvec upd_theta_metrop_conditional(); // returns uvec with changes to thetaj
  void init_theta_adapt();
  
  // data with with misalignment
  arma::field<arma::uvec> avail_by_outcome;
  arma::umat missing_mat;
  arma::uvec rows_with_missing;
  arma::uvec rows_some;
  int n_some; // n where at least one is nonmissing
  bool Y_needs_filling;
  arma::uvec Y_na_indices;
  void manage_missing_data(){
    // managing misalignment, i.e. not all outcomes observed at all locations
    // indices of non-NAs 
    avail_by_outcome = arma::field<arma::uvec>(q);
    missing_mat = arma::zeros<arma::umat>(n, q);
    arma::uvec row_miss_01 = arma::zeros<arma::uvec>(n);
    arma::uvec row_complete_miss_01 = arma::zeros<arma::uvec>(n);
    Y_na_indices = arma::find_nonfinite(Y);
    for(int j=0; j<q; j++){
      avail_by_outcome(j) = arma::find_finite(Y.col(j));
      for(int i=0; i<n; i++){
        missing_mat(i,j) = !arma::is_finite(Y(i,j));
      }
    }
    Y.elem(arma::find(missing_mat)).zeros(); // set missing to zero
    
    for(int i=0; i<n; i++){
      row_miss_01(i) = arma::any(missing_mat.row(i) == 1);
      row_complete_miss_01(i) = arma::all(missing_mat.row(i) == 1); 
    }
    rows_with_missing = arma::find(row_miss_01 == 1);
    rows_some = arma::find(row_complete_miss_01 == 0);
    n_some = rows_some.n_elem;
    if(n_some < n){
      if(latent_model != 1){
        Rcpp::stop("Fully missing rows in Y detected. Use latent_model=1 or exclude from data and rerun.\n");
      }
    }
    //if(Y.has_nonfinite() & (latent_model==1)){
    //  Rcpp::stop("nq block not implemented for misaligned data.\n");
    //}
    if(arma::accu(missing_mat)>0){
      Y_needs_filling = true;
    } else {
      Y_needs_filling = false;
    }
  }
  void sample_Y_misaligned(const arma::uvec& theta_changed);
  
  // Preconditioner choices for the joint BW PCG in gibbs_BW_block.
  // Only two real PCs survive — the adaptive PROBE selects between them.
  //
  //   JACOBI    : diagonal of the joint operator (precision form).
  //               O(nq + pq) per apply; nearly free; tends to be a weak PC.
  //
  //   POSTERIOR : block-diagonal-on-(B,W) PC with cross-outcome Σ-mix on
  //               the W half.  B half: per-outcome exact dense Cholesky of
  //                 M_B,j = diag(1/B_Var.col(j)) + X^T·diag(invD.col(j))·X
  //               W half: blkdiag(H_A,j^T)·(R_corr ⊗ I)·blkdiag(H_A,j)
  //                 with R_corr = D_σ^{-1}·Σ·D_σ^{-1}, the correlation matrix
  //                 of Σ.  H_A,j is the per-outcome Vecchia precision factor
  //                 of A_j = Q_jj·H_j^T H_j + diag(invD.col(j)), built once
  //                 per chain via local regression on (m+1)x(m+1) sub-blocks
  //                 of A_j and mirrored as col-major sparse Eigen for the
  //                 fast double-mult apply.
  //               R_corr chosen so the PC equals A_W^{-1} exactly in the
  //               D→0 limit (modulo Vecchia); reduces to a per-outcome block
  //               PC at Σ diagonal.  Extra cost vs the diagonal-Σ version:
  //               one n×q · q×q dense multiply per CG iter (cheap for q small).
  //
  //   PROBE     : adaptive default.  Runs POSTERIOR for the first 5 sweeps,
  //               then JACOBI for the next 5 with maxit capped at
  //               2·max(posterior CG iters).  JACOBI wins iff it converges
  //               (within the cap) in fewer iters on average; otherwise
  //               POSTERIOR is locked in.  POSTERIOR is the safe first
  //               choice — without its iter count, JACOBI could waste many
  //               sweeps before we'd know it's failing.
  enum PrecondChoice {
    PRECOND_PROBE     = 0,
    PRECOND_JACOBI    = 1,
    PRECOND_POSTERIOR = 2
  };
  PrecondChoice precond_choice = PRECOND_PROBE;

  // Probe state.  Two phases of 5 sweeps each:
  //   probe_count ∈ [0,5)  : run POSTERIOR, accumulate post_iter_sum and
  //                          track post_iter_max.
  //   probe_count ∈ [5,10) : run JACOBI with maxit cap = 2·post_iter_max,
  //                          accumulate jac_iter_sum and jac_converged_all.
  //   probe_count == 10    : decide and lock precond_choice.
  // After lock-in, the probe machinery is dormant (precond_choice ≠ PROBE).
  int           probe_count        = 0;
  static constexpr int probe_per_pc = 5;
  int           post_iter_sum      = 0;
  int           post_iter_max      = 0;
  int           jac_iter_sum       = 0;
  int           jac_iter_cap       = 0;   // = 2 * post_iter_max; set after POSTERIOR phase
  bool          jac_converged_all  = true; // false if any sweep hit the cap without converging

  // Per-outcome state for the POSTERIOR PC: Vecchia precision factor of
  // A_j = Q_jj·H_j^T H_j + diag(invD.col(j)), built once via local regression
  // on (m+1)x(m+1) sub-blocks of A_j.  Mirrored as col-major sparse Eigen
  // for the fast double-mult apply.  Idempotent build via vaprec_n_builds.
  // (Naming kept as "vaprec" for historical continuity with the build code.)
  int vaprec_n_builds = 0;
  std::vector<arma::vec>              vaprec_sqrtR;     // size q, each length n
  std::vector<arma::field<arma::vec>> vaprec_h;         // size q, each length n
  arma::field<arma::umat>             vaprec_children;  // length n, each rows = [k, t_of_i_in_parents_of_k]
  std::vector<Eigen::SparseMatrix<double>> vaprec_H_eigen;   // size q (col-major lower)
  std::vector<Eigen::SparseMatrix<double>> vaprec_Ht_eigen;  // size q (col-major upper, = transpose)
  // Build the Vecchia factors of A_j (no-op if already built).  Called once
  // on the first POSTERIOR apply.
  void build_vaprec_factors();

  // Telemetry: number of CG iterations used in the most-recent gibbs_BW_block
  // call, plus an integer code for which preconditioner ran
  // (0 unset / 1 jacobi / 2 posterior).  Read by the outer MCMC driver
  // and surfaced back to R.
  int last_cg_iter      = 0;
  int last_precond_used = 0;
  
  // latent model 
  int latent_model; // 0: response, 1: block, 2: row seq, 3: col seq
  arma::mat W;
  void w_sequential_singlesite(const arma::uvec& theta_changed);

  // Per-outcome sequential sampler for W (latent_model = 3).  Uses a fixed
  // sparse-precision-form CG with Jacobi PC (pcg_diag_solve) — the precond
  // enum is accepted for ABI compatibility but ignored.
  void gibbs_w_sequential_byoutcome(int& cg_iter, PrecondChoice precond);

  // β-conditional W update (ASIS reparameterisation) used by the non-block
  // latent samplers (latent_model = 2, 3) when sample_Beta is on.  Carries
  // its own inline Jacobi PC — doesn't touch the PrecondChoice enum.
  void update_BW_asis(int& cg_iter, arma::mat& B, arma::mat& W, bool sampling);

  // Joint BW PCG sampler — the only block sampler that survives the cleanup.
  // PC dispatched on `precond`: PRECOND_JACOBI (diagonal of the joint
  // precision) or PRECOND_POSTERIOR (block-diagonal-on-(B,W) with cross-
  // outcome Σ-mix on the W half via VAPREC factors).  `sampling` toggles
  // between MCMC sampling mode (RHS includes Bhattacharya noise terms) and
  // posterior-mean mode (deterministic).  `cg_maxit_override` > 0 caps the
  // PCG iter count (used by the probe phase to budget JACOBI against
  // POSTERIOR's measured iter count); 0 means use the default cap (n).
  void gibbs_BW_block(int& cg_iter, PrecondChoice precond, bool sampling=true,
                      int cg_maxit_override=0);
    
  arma::vec Ddiag;
  
  // centering of W and move to intercept (if there is one)
  void W_centering();
  // update running means
  void update_running_means(arma::mat&, const arma::mat&, bool pr=false); 
  
  // utilities for gibbs
  void update_BWSigma_px();
  void update_Sigma_gibbs();
  void update_Ddiag_gibbs();
  
  // whitened Y and X (that is, applying the same operation that makes W white noise)
  arma::mat Ytilde;
  arma::cube HX;
  
  // utilities for vi
  int vi_min_iter; // burn-in for vi in latent models
  int vi_it; // internal iteration counter
  //void update_B_vi();
  void update_Sigma_vi();
  void update_Ddiag_vi();
  bool vi;
  //arma::mat B_post_cov;
  arma::mat VTV, VTV_ma; // ma for moving average
  arma::mat ETE, ETE_ma;  
  arma::mat E_B; 
  arma::mat W_RB;
  arma::mat E_W;
  
  void vi_Beta_UQ(); // for computing Beta_UQ
  arma::vec delta_t, beta_running_mean;
  arma::mat Beta_UQ;
  arma::vec Ddiag_UQ;
  arma::mat Sigma_UQ;
  
  
  // utility for latent model and misaligned response model
  arma::field<arma::mat> Rw_no_Q;
  arma::field<arma::mat> Pblanket_no_Q;
  void cache_blanket_comps(const arma::uvec& theta_changed);
  
  // which theta updates are we doing
  bool phi_sampling, sigmasq_sampling, nu_sampling, alpha_sampling;
  
  // adaptive metropolis to update theta atoms
  int theta_mcmc_counter;
  arma::uvec which_theta_elem;
  arma::mat theta_unif_bounds;
  //arma::mat theta_metrop_sd;
  RAMAdapt theta_adapt;
  bool theta_adapt_active;
  // adaptive metropolis (conditional update) to update theta atoms
  // assume shared covariance functions and unknown parameters across variables
  arma::mat c_theta_unif_bounds;
  std::vector<RAMAdapt> c_theta_adapt;
  // --------
  
  
  // -------------- run 1 gibbs iteration based on current values
  void response_gibbs(int it, int sample_sigma, bool sample_beta, bool update_theta, bool sample_tausq=false);
  void latent_gibbs(int it, int sample_sigma, bool sample_beta, bool update_theta, bool sample_tausq=false);
  void response_vi();
  void latent_vi();
  double latent_fit_eval();
  
  std::chrono::steady_clock::time_point tstart;
  arma::vec timings;
  
  // -------------- constructor for building C (dense! beware)
  SpIOX(const arma::mat& _coords,
        const arma::field<arma::uvec>& custom_dag,
        int dag_opts,
        const arma::mat& daggp_theta, 
        const arma::mat& Sigma,
        int cov_model_matern,
        int num_threads_in){
    
    theta = daggp_theta;
    
    q = theta.n_cols;
    n = _coords.n_rows;
    
    daggps = std::vector<DagGP>(q);

    matern = cov_model_matern;

    // make n_threads depend on whether data are gridded, since behavior is opposite
    gridded = dag_opts==-1;
    int daggp_n_threads = gridded ? 1 : num_threads;
    for(unsigned int i=0; i<q; i++){
      daggps[i] = DagGP(_coords, theta.col(i), custom_dag,
                        dag_opts,
                        matern,
                        daggp_n_threads);
    }
  }
  
  // -------------- constructor
  SpIOX(const arma::mat& _Y,
        const arma::mat& _X,
        const arma::mat& _coords,
        const arma::field<arma::uvec>& custom_dag,
        int dag_opts,
        int latent_model_choice,
        const arma::mat& Beta_start,
        const arma::mat& W_start,
        const arma::mat& Sigma_start,
        const arma::mat& daggp_theta,
        const arma::uvec& update_theta_which,
        const arma::vec& tausq_start,
        int cov_model_matern,
        int num_threads_in,
        int _vi_min_iter = 100)
  {
    num_threads = num_threads_in;

    Y = _Y;
    X = _X;

    n = Y.n_rows;
    q = Y.n_cols;
    p = X.n_cols;
    
    // intercept? 
    intercept = -1;
    for(int j=0; j<p; j++){
      if(arma::all(X.col(j) == 1.0)){
        intercept = j;
        break;
      }
    }
    
    latent_model = latent_model_choice; // latent model?

    manage_missing_data();
    
    if(latent_model>0){
      W = W_start;
      W(arma::find_nonfinite(W)).fill(0);
      E_W = W;
      W_RB = W;
      Ddiag = tausq_start;
    }
    
    
    // mcvi params
    vi = _vi_min_iter > 0;
    vi_min_iter = _vi_min_iter;
    vi_it = 0;
    
    VTV = arma::zeros(q, q);
    ETE = arma::zeros(q, q);
    delta_t = arma::zeros(p*q);
    beta_running_mean = arma::zeros(p*q);
    Beta_UQ = arma::zeros(p*q, p*q);
    
    VTV_ma = arma::zeros(q, q);
    ETE_ma = arma::zeros(q, q);
    
    Ytilde = arma::zeros(n, q);
    HX = arma::zeros(n, p, q);
    
    // continue with other init
    B = Beta_start;
    E_B = B;
    YXB = Y - X * B;
    
    B_Var = 1000 * arma::ones(arma::size(B));
    
    theta = daggp_theta;
    daggps = std::vector<DagGP>(q);

    // if multiple nu options, interpret as wanting to sample smoothness for matern
    // otherwise, power exponential with fixed exponent.
    phi_sampling = update_theta_which(0) == 1;
    sigmasq_sampling = update_theta_which(1) == 1;
    nu_sampling = update_theta_which(2) == 1;
    alpha_sampling = update_theta_which(3) == 1;

    matern = cov_model_matern;

    // make n_threads depend on whether data are gridded, since behavior is opposite
    gridded = dag_opts==-1;
    int daggp_n_threads = gridded ? 1 : num_threads;
    for(unsigned int i=0; i<q; i++){
      daggps[i] = DagGP(_coords, theta.col(i), custom_dag,
                               dag_opts,
                               matern,
                               daggp_n_threads);
    }
    daggps_alt = daggps;

    init_theta_adapt();
    
    S = arma::chol(Sigma_start, "upper");
    Si = arma::inv(arma::trimatu(S));
    Sigma = S.t() * S;
    Q = Si * Si.t();
    
    // future
    //A = S.t();
    //Aplus = arma::pinv(A);
    //AplusT = Aplus.t();
    
    compute_V();
    
    int nfill = latent_model == 2 ? n : rows_with_missing.n_elem;
    Rw_no_Q = arma::field<arma::mat> (nfill);
    Pblanket_no_Q = arma::field<arma::mat> (nfill);
    
    arma::uvec updater = arma::ones<arma::uvec>(q);
    if(latent_model | Y_needs_filling){
      // first time making markov blanket cache
      cache_blanket_comps(updater);
    }
    
    timings = arma::zeros(10);
  }
};
