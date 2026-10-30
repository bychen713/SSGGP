#ifndef SSGGP_CH_KERNEL_H
#define SSGGP_CH_KERNEL_H

// The confluent hypergeometric density, factored out so ch_spectral_cpp (which
// returns the whole vector) and the fused CM3 objective (which returns only
// sums) evaluate the SAME arithmetic. Two copies of a prefactor would be two
// densities that agree until someone edits one of them.

#include <Rcpp.h>
#include <cmath>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>

// GSL's default error handler calls abort(), which would take the R session
// down the first time the optimizer probes an awkward parameter.
inline void ssggp_gsl_handler_off() {
  static bool done = false;
  if (!done) { gsl_set_error_handler_off(); done = true; }
}

// U(a, b, z), returning NA rather than aborting when GSL cannot evaluate it.
inline double ssggp_U_safe(double a, double b, double z) {
  if (z <= 0.0) {
    // As z -> 0: finite for b < 1, log divergence at b = 1, ~ z^(1 - b) for b > 1
    if (b < 1.0) return gsl_sf_gamma(1.0 - b) / gsl_sf_gamma(a + 1.0 - b);
    return R_PosInf;
  }
  gsl_sf_result r;
  if (gsl_sf_hyperg_U_e(a, b, z, &r)) return NA_REAL;
  return r.val;
}

// The parameter-only part of f: everything outside U.
inline double ssggp_ch_pref(double beta, double nu, double alpha, int d,
                            double variance) {
  const double hd = 0.5 * d;
  // B(alpha, nu) = Gamma(alpha) Gamma(nu) / Gamma(alpha + nu)
  const double Beta = gsl_sf_gamma(alpha) * gsl_sf_gamma(nu) /
                      gsl_sf_gamma(alpha + nu);
  return variance * gsl_sf_gamma(nu + hd) * std::pow(beta, (double) d) /
         (std::pow(2.0 * M_PI, hd) * Beta);
}

#endif
