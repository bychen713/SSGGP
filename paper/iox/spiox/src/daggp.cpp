#include "daggp.h"
#include <thread>
#include <Eigen/Core>

DagGP::DagGP(
  const arma::mat& coords_in,
  const arma::vec& theta_in,
  const arma::field<arma::uvec>& custom_dag,
  int dag_opts,
  int covariance_matern,
  int num_threads_in){

  coords = coords_in;
  theta = theta_in;
  nr = coords.n_rows;

  dag = custom_dag;
  
  oneuv = arma::ones<arma::uvec>(1);
  
  matern = covariance_matern; // 0 pexp or 1 matern or 2 wave
  //Rcpp::Rcout << "DagGP covariance choice " << matern << endl;
  //thread safe stuff
  n_threads = num_threads_in;
  
  int bessel_ws_inc = MAT_NU_MAX;//see bessel_k.c for working space needs
  // bessel_ws is indexed by omp_get_thread_num() at call-time.  Methods on
  // this DagGP may be invoked from any outer parallel context — including
  // spiox-level omp parallel-for loops that spawn more threads than the
  // current omp_get_max_threads() setting (e.g. when a previous DagGP was
  // built with n_threads=1, which sets the OMP default to 1 globally).
  // Size the buffer for the hardware-concurrency upper bound so it covers
  // any reasonable outer fan-out and the locally-requested n_threads_in.
#ifdef _OPENMP
  int max_threads = std::max({omp_get_num_procs(),
                              omp_get_max_threads(),
                              num_threads_in});
#else
  int max_threads = std::max(1, num_threads_in);
#endif
  bessel_ws = (double *) R_alloc(max_threads * bessel_ws_inc, sizeof(double));

  //
  mblanket = arma::field<arma::uvec>(nr);
  hrows = arma::field<arma::vec>(nr);
  ax = arma::field<arma::uvec>(nr);
  
  if (dag_opts == -1) {
    // Gridded coords + translation-invariant kernel: amortise per-group
    // Pinv / h / CPC across nodes that share a relative-offset pattern.
    gridded = true;
    cache_map = arma::zeros<arma::umat>(nr, 2);
    build_grid_exemplars();
  } else if (dag_opts == 0) {
    // Non-gridded: every node computes its own block.  dag_cache is just
    // `dag` and cache_map is identity, so the per-group loop in
    // compute_comps still works without special-casing the non-gridded path.
    gridded = false;
    dag_cache = dag;
    cache_map.set_size(nr, 2);
    for (unsigned int i = 0; i < nr; i++) {
      cache_map(i, 0) = i;
      cache_map(i, 1) = i;
    }
  } else {
    Rcpp::stop("DagGP: dag_opts must be 0 (non-gridded) or -1 (gridded).");
  }
  
  
  compute_comps();
  build_mblanket_from_dag();
  color_from_mblanket();
}



void DagGP::build_grid_exemplars() {
  // this function computes exemplar sets (parents, i) for gridded data
  // and assuming isotropic covariance
  // in this case the relationships between locations replicates across the domain
  // so we can massively reduce the number of matrix operations
  // this function is best paired with a neighbor search on an ordering that
  // makes it more likely to find many ties.
  
  const arma::uword N = coords.n_rows;
  const arma::uword d = coords.n_cols;
  if (N != dag.n_elem) Rcpp::stop("coords / dag size mismatch");
  
  arma::field<arma::uvec> dagplus = dag;
  if(gridded){
    // gridded coords, so cache the child too
    for(unsigned int i=0; i<nr; i++){
      dagplus(i) = arma::join_vert(dag(i), oneuv*i);
    }
  }
  
  
  // 1) grid index per dimension (0-based)
  arma::Mat<arma::uword> idx(N, d);
  for (arma::uword k = 0; k < d; ++k) {
    auto levels = sorted_unique_col(coords.col(k));
    idx.col(k) = indexify_col_exact(coords.col(k), levels);
  }
  
  // 2) hash by ordered offset pattern
  std::unordered_map<std::string, int> key2gid;
  key2gid.reserve(N * 2);
  
  cache_map.set_size(N, 2);
  std::vector< arma::uvec > exemplar_sets; exemplar_sets.reserve(N);
  
  for (arma::uword i = 0; i < N; ++i) {
    const arma::uvec& Pi = dagplus(i);           // 0-based parents for node i
    std::string key = key_from_offsets(idx, i, Pi);
    auto it = key2gid.find(key);
    int gid;
    if (it == key2gid.end()) {
      gid = (int)exemplar_sets.size();
      key2gid.emplace(std::move(key), gid);
      exemplar_sets.push_back(Pi);            // store this node's parents as exemplar
    } else {
      gid = it->second;
    }
    cache_map(i, 0) = i;                      // node id (0-based)
    cache_map(i, 1) = gid;                    // exemplar id (0-based)
  }
  
  // 3) materialize dag_cache from exemplars
  dag_cache.set_size(exemplar_sets.size());
  for (size_t g = 0; g < exemplar_sets.size(); ++g) {
    dag_cache(g) = exemplar_sets[g];
  }
  
  const arma::uword G = exemplar_sets.size();
  dag_cache.set_size(G);
  child_cache = arma::zeros<arma::uvec>(G);  // arma::uvec
  
  for (arma::uword g = 0; g < G; ++g) {
    if(gridded){
      const arma::uvec& v = exemplar_sets[g];      // expects [parents..., child]
      const arma::uword n = v.n_elem;
      
      if (n == 0) {
        dag_cache(g)   = arma::uvec();             // empty parents
        child_cache(g) = arma::uword(-1);          // sentinel (or choose 0)
        continue;
      }
      
      // parents = first n-1 elements (copy); handles n==1 by making empty uvec
      if (n > 1) dag_cache(g) = v.head(n - 1);
      else       dag_cache(g) = arma::uvec();      // no parents
      // child = last element
      child_cache(g) = v[n - 1];                   // or v.back()
    } else {
      dag_cache(g) = exemplar_sets[g];
    }
  }
}

void DagGP::color_from_mblanket(){
  
  arma::uword n_nodes = mblanket.n_elem;
  
  // Stores the color assigned to each node (-1 = unassigned)
  // This allows O(1) lookup to see a neighbors color.
  std::vector<int> node_colors(n_nodes, -1);
  
  // Intermediate storage: groups[c] is the list of nodes with color c
  std::vector<std::vector<arma::uword>> groups;
  
  for (arma::uword i = 0; i < n_nodes; ++i) {
    
    // 1. Identify which colors are forbidden for this node
    // We assume the number of colors wont exceed the current max + 1
    std::vector<bool> forbidden(groups.size() + 1, false);
    
    const arma::uvec& neighbors = mblanket(i);
    for (arma::uword neighbor_idx : neighbors) {
      // Safety check for index and check if neighbor is already colored
      if (neighbor_idx < n_nodes && node_colors[neighbor_idx] != -1) {
        int neighbor_c = node_colors[neighbor_idx];
        
        // If neighbors color is within our current range, mark it
        if (neighbor_c < forbidden.size()) {
          forbidden[neighbor_c] = true;
        }
      }
    }
    
    // 2. Pick the first color that is NOT forbidden
    int chosen_color = 0;
    while (chosen_color < forbidden.size() && forbidden[chosen_color]) {
      chosen_color++;
    }
    
    // 3. Assign the color
    node_colors[i] = chosen_color;
    
    // 4. Add to the corresponding group
    // If we picked a new color we havent used before, resize the groups vector
    if (chosen_color >= groups.size()) {
      groups.resize(chosen_color + 1);
    }
    groups[chosen_color].push_back(i);
  }
  
  // 5. Convert std::vector to arma::field
  colors = arma::field<arma::uvec>(groups.size());
  for (size_t c = 0; c < groups.size(); ++c) {
    colors(c) = arma::uvec(groups[c]);
  }
  
}

void DagGP::compute_comps(){
  // Compute the row-wise Vecchia factors (sqrtR, h, hrows, ax) plus the
  // precision log-determinant.  H is never assembled here — callers that
  // need a sparse matrix build one on demand via make_H().
  sqrtR = arma::zeros(nr);
  h = arma::field<arma::vec>(nr);
  arma::vec logdetvec(nr);
  
  arma::uvec errors = arma::zeros<arma::uvec>(nr);
  
  std::chrono::steady_clock::time_point tstart;
  std::chrono::steady_clock::time_point tend;
  int timed;
  
  tstart = std::chrono::steady_clock::now();
  const int G = static_cast<int>(dag_cache.n_elem);
  arma::field<arma::mat> Pinv_cache(G);
  arma::field<arma::mat> h_cache(G);
  arma::vec CPC_cache(G);

  bool cache_error = false;
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) 
#endif
  for (int g = 0; g < G; ++g) {
    try {
      const arma::uvec& pxg = dag_cache(g);          // 0-based indices
      
      if (pxg.n_elem == 0u) { Pinv_cache(g).reset(); continue; }
      
      arma::mat Pg  = Correlationf(coords, pxg, pxg, theta, bessel_ws, matern, /*same=*/true);
      if (nugget.n_elem > 0) Pg.diag() += nugget(pxg);
      Pinv_cache(g) = arma::inv_sympd(Pg);
      
      if(gridded){
        const arma::uvec& ixg = oneuv*child_cache(g);
        arma::mat CPt = Correlationf(coords, pxg, ixg, theta, bessel_ws, matern, false);
        h_cache(g) = Pinv_cache(g) * CPt;
        CPC_cache(g) = arma::conv_to<double>::from(CPt.t() * h_cache(g));
      }
    } catch (...) {
      cache_error = true;
    }
  }
  tend = std::chrono::steady_clock::now();
  timed = std::chrono::duration_cast<std::chrono::microseconds>(tend - tstart).count();
  
  if(cache_error){
    Rcpp::stop("Failed to compute cached elements.\n");
  }
  
  //Rcpp::Rcout << "First loop: " << timed << endl; 
  
  tstart = std::chrono::steady_clock::now();
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads)
#endif
  for(int i=0; i<nr; i++){
    try {
      arma::uvec ix = oneuv * i;                 
      const int gid = static_cast<int>(cache_map(i,1));   
      
      arma::uvec px;
      if(gridded){
        // gridded, look up original dag
        // because dag_cache has the exemplar parents
        px = dag(i);
      } else {
        px = dag_cache(gid);   
      }   
      
      ax(i) = arma::join_vert(ix, px);
      
      arma::mat CC  = Correlationf(coords, ix, ix, theta, bessel_ws, matern, true);
      if (nugget.n_elem > 0) CC(0, 0) += nugget(i);

      arma::mat CPC;
      arma::vec hi;
      
      if (px.n_elem == 0u) {
        hi.set_size(0);
        sqrtR(i) = std::sqrt( arma::conv_to<double>::from( CC ) );
      } else {
        //const arma::mat& Pinv = Pinv_cache(gid);   
        // gridded data
        if(gridded){
          hi = h_cache(gid);   
          CPC = CPC_cache(gid);
        } else {
          arma::mat CPt = Correlationf(coords, px, ix, theta, bessel_ws, matern, false);
          hi = Pinv_cache(gid) * CPt;             
          CPC = CPt.t() * hi;
        }
        sqrtR(i) = std::sqrt( arma::conv_to<double>::from( CC - CPC ) );
      }
      
      h(i) = hi;
      
      arma::vec mhR = (px.n_elem ? -hi / sqrtR(i) : arma::vec());
      arma::vec rowoneR = arma::ones<arma::vec>(1) / sqrtR(i);
      
      hrows(i) = (px.n_elem ? arma::join_vert(rowoneR, mhR) : rowoneR);
      logdetvec(i) = -std::log(sqrtR(i));
    } catch (...) {
      errors(i) = 1;
    }
  }
  tend = std::chrono::steady_clock::now();
  
  timed = std::chrono::duration_cast<std::chrono::microseconds>(tend - tstart).count();
  
  //Rcpp::Rcout << "Second loop: " << timed << endl; 
  
  if(arma::any(errors)){
    Rcpp::stop("Failure in building sparse DAG GP. Check coordinates/values of theta.");
  }

  precision_logdeterminant = 2 * arma::accu(logdetvec);

  // Build the Eigen mirror of H once the row-wise representation is final.
  // All four operators below dispatch to this; the arma row-wise rep stays
  // around for make_H() / H_col_squared_norms() and is also the source of
  // truth here.
  build_H_eigen();
}

void DagGP::build_H_eigen() {
  // Assemble H_eigen (col-major, lower triangular) from the row-wise rep
  // (ax / hrows / sqrtR).  Triplets are (row, col, value); setFromTriplets
  // sorts and packs them into CSC.  Cost ~ O(nnz log nnz) — negligible
  // against the surrounding Vecchia factorisation.
  std::vector<Eigen::Triplet<double>> trips;
  arma::uword nnz_guess = nr;
  for (int i = 0; i < nr; ++i) nnz_guess += ax(i).n_elem;
  trips.reserve(nnz_guess);
  for (int i = 0; i < nr; ++i) {
    const arma::vec&  vals = hrows(i);
    const arma::uvec& cols = ax(i);
    for (arma::uword idx = 0; idx < vals.n_elem; ++idx) {
      trips.emplace_back((int)i, (int)cols(idx), vals(idx));
    }
  }
  H_eigen.resize(nr, nr);
  H_eigen.setFromTriplets(trips.begin(), trips.end());
  H_eigen.makeCompressed();
  // Explicit upper-triangular transpose for the H^T paths.  Eigen could
  // form .transpose() lazily, but materialising once gives the back-solve
  // a native column-major Upper view with no per-call wrapper overhead.
  Ht_eigen = H_eigen.transpose();
  Ht_eigen.makeCompressed();
}


double DagGP::logdens(const arma::vec& x){
  // x^T H^T H x = ||H x||^2 ; use the operator path (no sparse H).
  arma::vec Hx = H_times_A(x);
  double loggausscore = arma::dot(Hx, Hx);
  return 0.5 * ( precision_logdeterminant - loggausscore );
}

void DagGP::update_theta(const arma::vec& newtheta){
  theta = newtheta;
  compute_comps();
}

arma::vec DagGP::H_col_squared_norms() const {
  // diag(H^T H)(c) = sum_k H(k, c)^2.  Row k of H stores values hrows(k)
  // at columns ax(k); scatter the squared values into the result.
  arma::vec result = arma::zeros(nr);
  for (int k = 0; k < nr; ++k) {
    const arma::vec&  vals = hrows(k);
    const arma::uvec& cols = ax(k);
    for (arma::uword idx = 0; idx < vals.n_elem; ++idx) {
      result(cols(idx)) += vals(idx) * vals(idx);
    }
  }
  return result;
}

arma::sp_mat DagGP::make_H() const {
  // Materialise H from the row-wise representation.  Cost ~ O(n·m).
  unsigned int H_n_elem = nr;
  for (int i = 0; i < nr; i++) {
    if (gridded) {
      H_n_elem += dag(i).n_elem;
    } else {
      int u = cache_map(i, 1);
      H_n_elem += dag_cache(u).n_elem;
    }
  }
  arma::umat H_locs(2, H_n_elem);
  arma::vec  H_values(H_n_elem);
  unsigned int ix = 0;
  for (int i = 0; i < nr; i++) {
    H_locs(0, ix) = i;
    H_locs(1, ix) = i;
    H_values(ix)  = 1.0 / sqrtR(i);
    ix++;
    arma::uvec daghere;
    if (gridded) {
      daghere = dag(i);
    } else {
      daghere = dag_cache(cache_map(i, 1));
    }
    for (unsigned int j = 0; j < daghere.n_elem; j++) {
      H_locs(0, ix) = i;
      H_locs(1, ix) = daghere(j);
      H_values(ix)  = -h(i)(j) / sqrtR(i);
      ix++;
    }
  }
  return arma::sp_mat(H_locs, H_values, nr, nr);
}

void DagGP::build_mblanket_from_dag() {
  // Markov blanket of i = parents(i) ∪ children(i) ∪ co-parents-of-children(i).
  // Computed directly from the DAG without going through Ci = H^T H.
  std::vector< std::vector<int> > children(nr);
  for (int k = 0; k < nr; ++k) {
    const arma::uvec& pk = dag(k);
    for (unsigned int t = 0; t < pk.n_elem; ++t) {
      children[(int)pk(t)].push_back(k);
    }
  }

  for (int i = 0; i < nr; ++i) {
    std::vector<int> blanket;
    const arma::uvec& pi = dag(i);
    blanket.reserve(pi.n_elem + 4 * children[i].size());
    for (unsigned int t = 0; t < pi.n_elem; ++t) blanket.push_back((int)pi(t));
    for (int c : children[i]) {
      blanket.push_back(c);
      const arma::uvec& pc = dag(c);
      for (unsigned int t = 0; t < pc.n_elem; ++t) {
        int p = (int)pc(t);
        if (p != i) blanket.push_back(p);
      }
    }
    std::sort(blanket.begin(), blanket.end());
    blanket.erase(std::unique(blanket.begin(), blanket.end()), blanket.end());

    mblanket(i) = arma::uvec(blanket.size());
    for (size_t t = 0; t < blanket.size(); ++t) {
      mblanket(i)(t) = (arma::uword) blanket[t];
    }
  }
}

// ---------------------------------------------------------------------------
// The four sparse operators.  All dispatch to the Eigen mirror H_eigen
// (col-major lower triangular) and Ht_eigen (col-major upper triangular,
// stored as the explicit transpose).  Inputs/outputs are arma::mat with
// column-major contiguous storage; we wrap them in zero-copy Eigen::Map
// views.  Eigen expression templates handle SIMD and unrolling, and the
// triangular solves dispatch to its column-walking sparse forward/back
// substitution — typically several times faster than the row-wise hand-rolled
// loop on realistic Vecchia DAGs.

arma::mat DagGP::H_solve_A(const arma::mat& A, bool /*use_spmat*/) const {
  arma::mat Y = A;
  if (Y.n_rows == 0u || Y.n_cols == 0u) return Y;
  Eigen::Map<Eigen::MatrixXd> Ye(Y.memptr(), Y.n_rows, Y.n_cols);
  H_eigen.triangularView<Eigen::Lower>().solveInPlace(Ye);
  return Y;
}

arma::mat DagGP::Ht_solve_A(const arma::mat& A, bool /*use_spmat*/) const {
  arma::mat Y = A;
  if (Y.n_rows == 0u || Y.n_cols == 0u) return Y;
  Eigen::Map<Eigen::MatrixXd> Ye(Y.memptr(), Y.n_rows, Y.n_cols);
  Ht_eigen.triangularView<Eigen::Upper>().solveInPlace(Ye);
  return Y;
}

arma::mat DagGP::H_times_A(const arma::mat& A, bool /*use_spmat*/) const {
  arma::mat result(nr, A.n_cols);
  if (A.n_cols == 0u) return result;
  Eigen::Map<const Eigen::MatrixXd> Ae(A.memptr(), A.n_rows, A.n_cols);
  Eigen::Map<Eigen::MatrixXd>       Re(result.memptr(), nr, A.n_cols);
  Re.noalias() = H_eigen * Ae;
  return result;
}

arma::mat DagGP::Ht_times_A(const arma::mat& A, bool /*use_spmat*/) const {
  arma::mat result(nr, A.n_cols);
  if (A.n_cols == 0u) return result;
  Eigen::Map<const Eigen::MatrixXd> Ae(A.memptr(), A.n_rows, A.n_cols);
  Eigen::Map<Eigen::MatrixXd>       Re(result.memptr(), nr, A.n_cols);
  Re.noalias() = Ht_eigen * Ae;
  return result;
}


arma::mat DagGP::Corr_export(const arma::mat& these_coords, 
                             const arma::uvec& ix, const arma::uvec& jx, 
                             int matern, bool same){
  return Correlationf(these_coords, ix, jx, theta, bessel_ws, matern, same);
}

//[[Rcpp::export]]
Rcpp::List daggp_build(const arma::mat& coords, const arma::field<arma::uvec>& dag,
                       double phi, double sigmasq, double nu, double tausq,
                       int matern=1, int num_threads=1, bool gridded=false){

  arma::vec theta(4);
  theta(0) = phi;
  theta(1) = sigmasq;
  theta(2) = nu;
  theta(3) = tausq;

  // dag_opts encodes the gridded flag for DagGP: -1 = gridded (cache by
  // relative-offset pattern), 0 = non-gridded (each node its own block).
  const int dag_opts = gridded ? -1 : 0;
  DagGP adag(coords, theta, dag, dag_opts, matern, num_threads);

  // Materialise H and Ci on demand for the R-side export.
  arma::sp_mat H  = adag.make_H();
  arma::sp_mat Ci = H.t() * H;

  return Rcpp::List::create(
    Rcpp::Named("H") = H,
    Rcpp::Named("Cinv") = Ci,
    Rcpp::Named("Cinv_logdet") = adag.precision_logdeterminant,
    Rcpp::Named("dag") = adag.dag,
    Rcpp::Named("dag_cache") = adag.dag_cache,
    Rcpp::Named("cache_map") = adag.cache_map,
    Rcpp::Named("markov_blankets") = adag.mblanket
  );
}



/*

using EigenFromArma = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>;
 
arma::sp_mat sparse_convert(const Eigen::SparseMatrix<double>& XE){
 unsigned int num_cols = XE.cols();
 unsigned int num_rows = XE.rows();
 
 std::vector<unsigned int> rowind_vect(
 XE.innerIndexPtr(), XE.innerIndexPtr() + XE.nonZeros());
 std::vector<unsigned int> colptr_vect(
 XE.outerIndexPtr(), XE.outerIndexPtr() + XE.outerSize() + 1);
 std::vector<double> values_vect(
 XE.valuePtr(), XE.valuePtr() + XE.nonZeros());
 
 arma::vec values(values_vect.data(), values_vect.size(), false);
 arma::uvec rowind(rowind_vect.data(), rowind_vect.size(), false);
 arma::uvec colptr(colptr_vect.data(), colptr_vect.size(), false);
 
 arma::sp_mat result(rowind, colptr, values, num_rows, num_cols);
 
 return result; 
} */

