// Confluent Hypergeometric (CH) covariance and spectral density.
//
//   C(h) = sigma * Gamma(nu+alpha)/Gamma(nu) * U(alpha, 1-nu, |h|^2/(2 beta^2))
//   f(x) = sigma * Gamma(nu+d/2) beta^d / ((2 pi)^{d/2} B(alpha,nu))
//                * U(nu+d/2, 1-alpha+d/2, beta^2 |x|^2 / 2)
//
// where U is Tricomi's confluent hypergeometric function of the second kind.
// The (2 pi)^{-d/2} factor makes the spectral density integrate to one.
//
// C(0) = sigma exactly, since U(a,b,0) = Gamma(1-b)/Gamma(a+1-b) for b < 1 and
// with b = 1-nu that cancels the prefactor.
//
// NOTE alpha > d/2. Below that the spectral density is not finite at the
// origin, so the model is not defined. The R side enforces
// the bound; this file only refuses to return nonsense.

#include <Rcpp.h>
#include "ch_kernel.h"

//' Tricomi's U on a vector
//'
//' Exposed for testing the covariance and spectral density against each other
//' and against known values; not part of the modeling interface.
//'
//' @param a,b Parameters of U.
//' @param z Non-negative numeric vector.
//' @return U(a, b, z), elementwise.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector ch_U(double a, double b, Rcpp::NumericVector z) {
  ssggp_gsl_handler_off();
  const R_xlen_t n = z.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) out[i] = ssggp_U_safe(a, b, z[i]);
  return out;
}

//' CH covariance matrix from a distance matrix
//'
//' Symmetric, so only the upper triangle is evaluated: n(n+1)/2 evaluations of
//' U rather than n^2.
//'
//' @param D Symmetric distance matrix.
//' @param beta Range.
//' @param nu Smoothness.
//' @param alpha Tail decay; must exceed d/2 for the model to be defined.
//' @param variance Marginal variance.
//' @return The covariance matrix.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix ch_cov_cpp(Rcpp::NumericMatrix D, double beta, double nu,
                               double alpha, double variance) {
  ssggp_gsl_handler_off();
  const int n = D.nrow();
  if (D.ncol() != n) Rcpp::stop("D must be square");
  Rcpp::NumericMatrix C(n, n);

  const double pref  = variance * gsl_sf_gamma(nu + alpha) / gsl_sf_gamma(nu);
  const double scale = 1.0 / (2.0 * beta * beta);
  const double b     = 1.0 - nu;

  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      const double h = D(i, j);
      // h == 0 handled exactly rather than through the Gamma-ratio limit
      const double v = (h == 0.0) ? variance
                                  : pref * ssggp_U_safe(alpha, b, h * h * scale);
      C(i, j) = v;
      C(j, i) = v;
    }
  }
  return C;
}

//' CH spectral density on a vector of radii
//'
//' Isotropic, so it depends on the frequency only through its length.
//'
//' @param r Non-negative radii, `|x|`.
//' @param beta,nu,alpha Range, smoothness and tail decay.
//' @param d Spatial dimension.
//' @param variance Marginal variance; leave at one for a unit-integral density.
//' @return The spectral density at each radius.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector ch_spectral_cpp(Rcpp::NumericVector r, double beta,
                                    double nu, double alpha, int d,
                                    double variance = 1.0) {
  ssggp_gsl_handler_off();
  const double hd = 0.5 * d;
  const double pref = ssggp_ch_pref(beta, nu, alpha, d, variance);
  const double a = nu + hd;
  const double b = 1.0 - alpha + hd;
  const double s = 0.5 * beta * beta;

  const R_xlen_t n = r.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) out[i] = pref * ssggp_U_safe(a, b, s * r[i] * r[i]);
  return out;
}
