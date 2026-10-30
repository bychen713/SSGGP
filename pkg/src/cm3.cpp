// Fused CM3 objective, and the graphical-lasso statistic contraction.
//
// CM3 maximizes, for one coordinate j at a time,
//
//   Q_j(iota_j) = - sum_l log s_lj
//                 - (1/4) sum_l g_l [ Q_jj u_l,jj / s_lj^2 + 2 b_lj / s_lj ]
//
// with s_lj = sqrt(f_j(X_l)) and b_lj = sum_{k != j} Q_jk u_l,jk / s_lk. The
// optimizer needs only the sums, so they are accumulated in one pass without
// materializing the density. Sums use long double, as R's sum() does.
//
// A trial point is rejected (value NA, which the R side maps to a large
// constant) when any f is non-finite or non-positive, or the objective is not
// finite -- the same rule as the R evaluator in cm3_sweep().
//
// Threading: only the CH density map is parallel. The reduction stays serial
// and in order, so the result does not depend on the thread count.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <gsl/gsl_sf_psi.h>
#include "ch_kernel.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;

// exp(x) underflows to zero below the first and overflows above the second
static const double LOG_UNDERFLOW = -745.1332191019411;
static const double LOG_OVERFLOW  =  709.782712893384;

enum { FAM_MATERN = 0, FAM_SE = 1, FAM_CH = 2 };

// Clamp to the processors this process may run on (the affinity mask, which is
// what Slurm's --cpus-per-task sets), not to OMP_NUM_THREADS.
static int clamp_threads(int nthread) {
#ifdef _OPENMP
  const int avail = omp_get_num_procs();
  return (avail > 0 && nthread > avail) ? avail : nthread;
#else
  (void) nthread;
  return 1;
#endif
}

// The CH density pref * U(A, B, sc r^2) at every frequency.
static std::vector<double> ch_density_map(const NumericVector& r2, double pref,
                                          double A, double B, double sc,
                                          int nthread) {
  const R_xlen_t L = r2.size();
  std::vector<double> f(static_cast<std::size_t>(L));
  double* fp = f.data();
  const double* r2p = &r2[0];
  const int nt = clamp_threads(nthread);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if(nt > 1)
#endif
  for (R_xlen_t l = 0; l < L; ++l) {
    // sqrt then square, as ch_spectral_vec's sqrt(rowSums(X^2)) does
    const double r = std::sqrt(r2p[l]);
    fp[l] = pref * ssggp_U_safe(A, B, sc * r * r);
  }
  return f;
}

//' Fused CM3 objective for the built-in spectral families
//'
//' @param family 0 Matern (a, nu), 1 squared exponential (a), 2 CH (beta, nu, alpha).
//' @param par The family's parameters, on the natural scale, in `par_names`
//'   order.
//' @param r2 Squared frequency norms, `rowSums(Xmat^2)`; constant for a fit.
//' @param gvals,u_jj,b,Qjj The coordinate's fixed pieces of `Q_j`.
//' @param d Spatial dimension.
//' @param want_grad Also return the gradient with respect to `log(par)`.
//'   Unavailable for CH.
//' @param nthread Threads for the CH density map.
//' @param log_cap Ceiling on `log R`, `R = g/f`. Negative (the default) selects
//'   the plain form. Otherwise `gvals` is read as `log g`, `b` as the
//'   ratio-scaled neighbor sum `sum_k Q_jk u_l,jk sqrt(R_lk)`, and the
//'   objective is evaluated as
//'   `Q_j = 1/2 sum log R - 1/4 sum (Qjj u R + 2 sqrt(R) b)`,
//'   which never forms `1/f`. The dropped `-1/2 sum log g` is constant in
//'   `iota`.
//' @return A list with `value`, the NEGATED `Q_j` that `optim` minimizes
//'   (`NA_real_` for a rejected point), `grad` when asked, and `capped` under
//'   the ceiling.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cm3_objective_cpp(int family, const NumericVector& par,
                             const NumericVector& r2, const NumericVector& gvals,
                             const NumericVector& u_jj, const NumericVector& b,
                             double Qjj, int d, bool want_grad,
                             int nthread = 1, double log_cap = -1.0) {
  const R_xlen_t L = r2.size();
  const double hd = 0.5 * d;
  const Rcpp::List bad = Rcpp::List::create(Rcpp::Named("value") = NA_REAL);

  // ---- the capped, ratio-only form ----------------------------------------
  if (log_cap > 0.0) {
    if (want_grad) Rcpp::stop("no analytic gradient under the ratio ceiling");
    double kc = 0.0, kp = 0.0, ka2 = 0.0, kinv = 0.0;   // matern / squared exponential
    std::vector<double> fch;                            // CH densities
    if (family == FAM_MATERN) {
      const double a = par[0], nu = par[1];
      if (!(a > 0.0) || !(nu > 0.0)) return bad;
      kc = 2.0 * nu * std::log(a) + std::lgamma(nu + hd) - std::lgamma(nu)
           - hd * std::log(M_PI);
      kp = nu + hd; ka2 = a * a;
    } else if (family == FAM_SE) {
      const double a = par[0];
      if (!(a > 0.0)) return bad;
      kc = -d * std::log(a) - hd * std::log(2.0 * M_PI);
      kinv = 1.0 / (2.0 * a * a);
    } else {
      ssggp_gsl_handler_off();
      const double beta = par[0], nu = par[1], alpha = par[2];
      if (!(beta > 0.0) || !(nu > 0.0) || !(alpha > 0.0)) return bad;
      const double pref = ssggp_ch_pref(beta, nu, alpha, d, 1.0);
      if (!R_finite(pref)) return bad;
      fch = ch_density_map(r2, pref, nu + hd, 1.0 - alpha + hd,
                           0.5 * beta * beta, nthread);
    }

    long double s_lr = 0.0L, s_q = 0.0L;
    R_xlen_t n_capped = 0;
    for (R_xlen_t l = 0; l < L; ++l) {
      double lf;
      if (family == FAM_MATERN) {
        lf = kc - kp * std::log(ka2 + r2[l]);
      } else if (family == FAM_SE) {
        lf = kc - r2[l] * kinv;
      } else {
        const double f = fch[l];
        if (ISNAN(f) || !R_finite(f) || !(f > 0.0)) return bad;
        lf = std::log(f);
      }
      if (ISNAN(lf)) return bad;
      double lr = gvals[l] - lf;           // gvals holds log g here
      if (lr > log_cap) { lr = log_cap; ++n_capped; }
      const double R = std::exp(lr);
      s_lr += (long double) lr;
      s_q  += (long double) (Qjj * u_jj[l] * R + 2.0 * std::sqrt(R) * b[l]);
    }
    const double value = -(0.5 * (double) s_lr - 0.25 * (double) s_q);
    if (!R_finite(value)) return bad;
    return Rcpp::List::create(Rcpp::Named("value")  = value,
                              Rcpp::Named("capped") = (double) n_capped);
  }

  long double s_quad = 0.0L;   // sum_l g_l (Qjj u_l / f_l + 2 b_l / sqrt(f_l))
  long double s_lf   = 0.0L;   // sum_l log f_l
  NumericVector grad;

  if (family == FAM_MATERN || family == FAM_SE) {
    const bool matern = (family == FAM_MATERN);
    const double a  = par[0];
    const double nu = matern ? par[1] : 0.0;
    if (!(a > 0.0) || (matern && !(nu > 0.0))) return bad;

    // f = C * t^-p  with t = a^2 + r^2, p = nu + d/2      (Matern)
    // f = C * exp(-r^2 / (2 a^2))                          (squared exponential)
    const double p    = matern ? nu + hd : 0.0;
    const double logC = matern
      ? 2.0 * nu * std::log(a) + std::lgamma(p) - std::lgamma(nu) - hd * std::log(M_PI)
      : -d * std::log(a) - hd * std::log(2.0 * M_PI);
    const double inv2a2 = matern ? 0.0 : 1.0 / (2.0 * a * a);
    if (!R_finite(logC)) return bad;

    // d log f / d a and d log f / d nu: frequency-free part, then per feature
    const double dA0 = want_grad ? (matern ? 2.0 * nu / a : -((double) d) / a) : 0.0;
    const double dN0 = (want_grad && matern)
      ? 2.0 * std::log(a) + gsl_sf_psi(p) - gsl_sf_psi(nu) : 0.0;
    const double a3 = a * a * a;

    long double sDa = 0.0L, sDn = 0.0L, qDa = 0.0L, qDn = 0.0L;
    for (R_xlen_t l = 0; l < L; ++l) {
      const double t  = matern ? a * a + r2[l] : 0.0;
      const double w  = matern ? std::log(t) : 0.0;
      const double lf = matern ? logC - p * w : logC - r2[l] * inv2a2;
      if (ISNAN(lf) || lf <= LOG_UNDERFLOW || lf >= LOG_OVERFLOW) return bad;
      const double inv_sf = std::exp(-0.5 * lf);
      const double inv_f  = inv_sf * inv_sf;
      // half-weighted form for the gradient; the objective doubles the b term
      const double h = Qjj * u_jj[l] * inv_f + b[l] * inv_sf;
      s_lf   += (long double) lf;
      s_quad += (long double) (gvals[l] * (h + b[l] * inv_sf));
      if (want_grad) {
        const double Da = matern ? dA0 - 2.0 * a * p / t : dA0 + r2[l] / a3;
        sDa += (long double) Da;
        qDa += (long double) (Da * gvals[l] * h);
        if (matern) {
          const double Dn = dN0 - w;
          sDn += (long double) Dn;
          qDn += (long double) (Dn * gvals[l] * h);
        }
      }
    }

    if (want_grad) {
      // dQ/d(log par) = par * dQ/d(par); the objective returned is -Q_j
      grad = NumericVector(matern ? 2 : 1);
      grad[0] = a * (0.5 * (double) sDa - 0.25 * (double) qDa);
      if (matern) grad[1] = nu * (0.5 * (double) sDn - 0.25 * (double) qDn);
    }
  } else if (family == FAM_CH) {
    if (want_grad) Rcpp::stop("no analytic gradient for the CH family");
    ssggp_gsl_handler_off();
    const double beta = par[0], nu = par[1], alpha = par[2];
    if (!(beta > 0.0) || !(nu > 0.0) || !(alpha > 0.0)) return bad;
    const double pref = ssggp_ch_pref(beta, nu, alpha, d, 1.0);
    if (!R_finite(pref)) return bad;
    const std::vector<double> fch = ch_density_map(
      r2, pref, nu + hd, 1.0 - alpha + hd, 0.5 * beta * beta, nthread);

    for (R_xlen_t l = 0; l < L; ++l) {
      const double f = fch[l];
      if (ISNAN(f) || !R_finite(f) || !(f > 0.0)) return bad;
      s_lf += (long double) std::log(f);
      const double inv_sf = 1.0 / std::sqrt(f);
      s_quad += (long double) (gvals[l] * (Qjj * u_jj[l] / f + 2.0 * b[l] * inv_sf));
    }
  } else {
    Rcpp::stop("unknown family id");
  }

  const double value = -(-0.5 * (double) s_lf - 0.25 * (double) s_quad);
  if (!R_finite(value)) return bad;
  if (want_grad) {
    for (R_xlen_t k = 0; k < grad.size(); ++k) {
      if (!R_finite(grad[k])) return bad;
    }
    return Rcpp::List::create(Rcpp::Named("value") = value,
                              Rcpp::Named("grad")  = grad);
  }
  return Rcpp::List::create(Rcpp::Named("value") = value);
}

//' The graphical LASSO statistic, from the per-feature moments
//'
//' `S[a, b] = sum_l M[l, a] M[l, b] u[l, k] / L`, accumulated in long double.
//'
//' @param u `L` by `P(P+1)/2` packed per-feature moments.
//' @param M `L` by `P` importance scaling.
//' @param ti,tj 1-based row and column of each packed column.
//' @return A `P` by `P` symmetric matrix.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix suffstat_to_s_cpp(const NumericMatrix& u, const NumericMatrix& M,
                                const Rcpp::IntegerVector& ti,
                                const Rcpp::IntegerVector& tj) {
  const R_xlen_t L = M.nrow();
  const int P = M.ncol();
  const R_xlen_t ntri = ti.size();
  NumericMatrix s(P, P);
  for (R_xlen_t k = 0; k < ntri; ++k) {
    const int A = ti[k] - 1, B = tj[k] - 1;
    const double* Ma = &M(0, A);
    const double* Mb = &M(0, B);
    const double* uk = &u(0, k);
    long double acc = 0.0L;
    for (R_xlen_t l = 0; l < L; ++l) acc += (long double) (Ma[l] * Mb[l] * uk[l]);
    const double v = ((double) acc) / (double) L;
    s(A, B) = v;
    s(B, A) = v;
  }
  return s;
}

//' Threads available to the compiled kernels
//'
//' @return `omp_get_max_threads()`, or `1` when built without OpenMP.
//' @keywords internal
// [[Rcpp::export]]
int ssggp_max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

//' Processors this process may run on
//'
//' The ceiling the CH density map clamps its thread count to.
//'
//' @return A single integer.
//' @keywords internal
// [[Rcpp::export]]
int ssggp_num_procs() {
#ifdef _OPENMP
  return omp_get_num_procs();
#else
  return 1;
#endif
}

//' Whether Eigen would parallelize its own dense products
//'
//' Must be `FALSE`: Makevars.in defines `EIGEN_DONT_PARALLELIZE` so that the
//' Gibbs kernels' dense products do not depend on the thread count.
//'
//' @return `TRUE` if Eigen would thread its dense products.
//' @keywords internal
// [[Rcpp::export]]
bool ssggp_eigen_parallel() {
#if defined(_OPENMP) && !defined(EIGEN_DONT_PARALLELIZE)
  return true;
#else
  return false;
#endif
}
