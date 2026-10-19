#ifndef UTILS_HPP
#define UTILS_HPP

#include <RcppEigen.h>

inline Eigen::VectorXd rnorm_vec(int n) {
  Eigen::VectorXd z(n);
  for (int i = 0; i < n; ++i) {
    z[i] = R::rnorm(0.0, 1.0);
  }
  return z;
}

#endif
