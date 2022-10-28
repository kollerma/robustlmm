#include "misc.h"

#define DEBUG(STRING)                                          \
// Rcpp::Rcout << STRING << std::endl;

void warn(std::string msg) {
  static Rcpp::Function *fun = NULL;
  if (fun == NULL)
     fun = new Rcpp::Function("warning");
  (*fun)(msg);
}

bool isDifferent(const double v1, const double v2, const double relativeTolerance) {
  double diff = v1 - v2;
  double scale = (std::abs(v1) + std::abs(v2)) / 2.;
  return std::abs(diff) > std::max(relativeTolerance * scale, relativeTolerance);
}
