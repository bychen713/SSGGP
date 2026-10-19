#include "sparse_solvers.h"
#include "daggp.h" 


arma::vec pcg_diag_solve(
    const arma::sp_mat& A,
    const arma::vec& b,
    const arma::vec& x0,
    const double tol,
    const int maxit,
    const double diag_floor,
    int num_threads
){
  const arma::uword n = A.n_rows;
  
  // CSC arrays -- reused as CSR because A is symmetric
  A.sync();
  const double*      vals = A.values;
  const arma::uword* rows = A.row_indices;
  const arma::uword* cols = A.col_ptrs;
  
  // Parallel symmetric SpMV:  y = A x
  auto spmv = [&](const arma::vec& xin, arma::vec& yout){
    const double* xp = xin.memptr();
    double*       yp = yout.memptr();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(num_threads)
#endif
    for(arma::uword i = 0; i < n; ++i){
      double yi = 0.0;
      const arma::uword ks = cols[i];
      const arma::uword ke = cols[i+1];
      for(arma::uword k = ks; k < ke; ++k){
        yi += vals[k] * xp[ rows[k] ];
      }
      yp[i] = yi;
    }
  };
  
  arma::vec x = (x0.n_elem == n) ? x0 : arma::zeros(n);
  arma::vec r(n), z(n), p(n), Ap(n);
  
  spmv(x, Ap);          // Ap = A x0
  r = b - Ap;
  
  const double bnorm = std::max(arma::norm(b, 2), 1e-300);
  
  arma::vec Mdiag(A.diag());
  for(arma::uword i=0; i<n; i++) if(!(Mdiag(i) > diag_floor)) Mdiag(i) = diag_floor;
  
  z = r / Mdiag;
  p = z;
  double rz_old = arma::dot(r, z);
  
  int k = 0;
  for(; k < maxit; k++){
    spmv(p, Ap);
    
    const double denom = arma::dot(p, Ap);
    if(!(denom > 0.0)) break;
    
    const double alpha = rz_old / denom;
    
    // Fused update: x += alpha p ; r -= alpha Ap ; and compute ||r|| on the fly
    double rr = 0.0;
    double*       xp = x.memptr();
    double*       rp = r.memptr();
    const double* pp = p.memptr();
    const double* ap = Ap.memptr();
#ifdef _OPENMP
#pragma omp parallel for reduction(+:rr) schedule(static) num_threads(num_threads)
#endif
    for(arma::uword i = 0; i < n; ++i){
      xp[i] += alpha * pp[i];
      rp[i] -= alpha * ap[i];
      rr    += rp[i] * rp[i];
    }
    
    if(std::sqrt(rr) / bnorm <= tol){ k++; break; }
    
    z = r / Mdiag;
    const double rz_new = arma::dot(r, z);
    const double beta = rz_new / rz_old;
    p = z + beta * p;
    rz_old = rz_new;
  }
  
  // (optional) 
  //Rcpp::Rcout << "[pcg] iters=" << k << "\n";
  return x;
}

