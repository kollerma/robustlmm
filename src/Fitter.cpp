#include "Fitter.h"

#define DEBUG(STRING)                                          \
// Rcpp::Rcout << STRING << std::endl;

using namespace Rcpp;

namespace Rcpp {
  template <typename T> SEXP wrapFit(const Fit<T>& fit) {
    return wrap(List::create(Named("value", fit.getValue()),
                             Named("convergence", fit.getConvergenceStatus()),
                             Named("iter", fit.getNumberOfOperations()),
                             Named("message", fit.getMessage())));
  }

  template <> SEXP wrap(const Fit<VectorXd>& fit) {
    return wrapFit(fit);
  }

  template <> SEXP wrap(const Fit<MatrixXd>& fit) {
    return wrapFit(fit);
  }

  template <> SEXP wrap(const Fit<double>& fit) {
    return wrapFit(fit);
  }
}
