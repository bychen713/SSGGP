#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <cmath>


// --- Helper: Generate Stencil with Pre-computed Linear Offsets ---
struct Stencil {
  arma::imat offsets;       // (d x K) matrix of coordinate offsets
  arma::vec linear_offsets; // (K) vector of linear index offsets
};

Stencil get_precomputed_stencil(int d, int m, const arma::uvec& strides) {
  // 1. Determine Radius
  int r = 1;
  while (std::pow(2 * r + 1, d) <= 3 * m + 1) r++; 
  if (r > 4 && d > 3) r = 4; // Safety cap
  
  int side = 2 * r + 1;
  int n_points = (int)std::pow(side, d);
  
  arma::imat offsets(d, n_points);
  
  // 2. Generate Offsets
  for(int i = 0; i < n_points; ++i) {
    int temp = i;
    for(int k = 0; k < d; ++k) {
      offsets(k, i) = (temp % side) - r;
      temp /= side;
    }
  }
  
  // 3. Sort by Distance
  arma::rowvec dists = arma::conv_to<arma::rowvec>::from(arma::sum(arma::square(offsets), 0));
  arma::uvec sorted_indices = arma::stable_sort_index(dists.t());
  
  // Remove self (index 0)
  if (n_points <= 1) return {arma::imat(), arma::vec()};
  arma::uvec keep = sorted_indices.subvec(1, n_points - 1);
  
  arma::imat sorted_offsets = offsets.cols(keep);
  
  // 4. Pre-compute Linear Offsets (dot product with strides)
  // We use double to prevent overflow during calculation, then cast to long long later
  arma::vec lin_offs(sorted_offsets.n_cols);
  for(int k=0; k<sorted_offsets.n_cols; ++k) {
    double off_val = 0;
    for(int dim=0; dim<d; ++dim) {
      off_val += (double)sorted_offsets(dim, k) * (double)strides(dim);
    }
    lin_offs(k) = off_val;
  }
  
  return {sorted_offsets, lin_offs};
}

// [[Rcpp::export]]
arma::field<arma::uvec> dag_for_gridded_cols(const arma::mat& coords, int m = 20) {
  if (coords.is_empty()) Rcpp::stop("coords is empty");
  
  int N = coords.n_rows;
  int d = coords.n_cols;
  
  // --- 1. Grid Discretization ---
  arma::umat grid_pts(N, d);
  arma::uvec dims(d);
  
  for(int k = 0; k < d; ++k) {
    arma::vec uniq = arma::unique(coords.col(k));
    dims(k) = uniq.n_elem;
    auto begin = uniq.begin();
    auto end = uniq.end();
    for(int i = 0; i < N; ++i) {
      auto it = std::lower_bound(begin, end, coords(i, k));
      grid_pts(i, k) = (int)(it - begin);
    }
  }
  
  // --- 2. Strides & Linear Indices ---
  arma::uvec strides(d);
  strides(0) = 1;
  for(int k = 1; k < d; ++k) strides(k) = strides(k - 1) * dims(k - 1);
  
  // Calculate max possible linear index to decide lookup strategy
  double max_lin_idx_d = 0;
  for(int k=0; k<d; ++k) max_lin_idx_d += (double)(dims(k)-1) * (double)strides(k);
  
  // --- 3. Lookup Table Strategy ---
  // If grid space is reasonable (< 200 million points), use direct vector lookup.
  // Otherwise fall back to hash map.
  bool use_direct_lookup = (max_lin_idx_d < 200000000.0); 
  
  std::vector<int> direct_lookup;
  std::unordered_map<int, int> map_lookup;
  
  arma::uvec lin_indices(N);
  
  if (use_direct_lookup) {
    int max_idx = (int)max_lin_idx_d;
    direct_lookup.resize(max_idx + 1, N + 1); // Init with sentinel
  } else {
    map_lookup.reserve(N * 2);
  }
  
  for(int i = 0; i < N; ++i) {
    double lin = 0; 
    for(int k=0; k<d; ++k) lin += (double)grid_pts(i, k) * (double)strides(k);
    lin_indices(i) = (int)lin;
  }
  
  // Sort to establish DAG order
  arma::uvec order = arma::sort_index(lin_indices);
  
  // Fill Lookup Table
  for(int i = 0; i < N; ++i) {
    int lin = lin_indices(order(i));
    if (use_direct_lookup) {
      direct_lookup[lin] = i; // Store 'rank' in sorted order
    } else {
      map_lookup[lin] = i;
    }
  }
  
  // --- 4. Process Stencil ---
  Stencil sten = get_precomputed_stencil(d, m, strides);
  arma::field<arma::uvec> parents(N);
  
  // Reuse vector to avoid allocation inside loop
  std::vector<int> found_buffer; 
  found_buffer.reserve(m + 5);
  
  for(int i = 0; i < N; ++i) {
    int original_idx = order(i);
    int current_lin = lin_indices(original_idx);
    
    // Direct pointer to current grid row for fast access
    arma::umat grid_pts_orig = grid_pts.row(original_idx);
    const arma::uword* current_coords = grid_pts_orig.memptr();
    
    found_buffer.clear();
    
    for(int k = 0; k < sten.offsets.n_cols; ++k) {
      if (found_buffer.size() >= (size_t)m) break;
      
      // 1. Fast Linear Offset Check
      // Calculate potential linear index using pre-computed offset
      double potential_lin_d = (double)current_lin + sten.linear_offsets(k);
      
      // Basic range check
      if (potential_lin_d < 0 || potential_lin_d > max_lin_idx_d) continue;
      
      int neighbor_lin = (int)potential_lin_d;
      
      // 2. Lookup Check (Does a node exist there?)
      int neighbor_rank = N + 1;
      
      if (use_direct_lookup) {
        neighbor_rank = direct_lookup[neighbor_lin];
        if (neighbor_rank > N) continue; // Empty slot
      } else {
        auto it = map_lookup.find(neighbor_lin);
        if (it == map_lookup.end()) continue;
        neighbor_rank = it->second;
      }
      
      // 3. DAG Check (Must be "past")
      if (neighbor_rank >= i) continue;
      
      // 4. Exact Coordinate Boundary Check (Expensive, so done last)
      // We only do this if we found a candidate in the "past".
      // This prevents "wrapping" around grid edges.
      bool valid_coords = true;
      for(int dim = 0; dim < d; ++dim) {
        int off = sten.offsets(dim, k);
        if (off == 0) continue;
        
        long long nc = (long long)current_coords[dim] + off;
        if (nc < 0 || nc >= (long long)dims(dim)) {
          valid_coords = false; 
          break;
        }
      }
      
      if (valid_coords) {
        found_buffer.push_back(order(neighbor_rank));
      }
    }
    parents(original_idx) = arma::conv_to<arma::uvec>::from(found_buffer);
  }
  
  return parents;
}