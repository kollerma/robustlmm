#include "misc.h"

void warn(std::string msg) {
  static Rcpp::Function *fun = NULL;
  if (fun == NULL) 
    fun = new Rcpp::Function("warning");
  (*fun)(msg);
}
