#if !defined  ROBUSTLMM_MISC_H__
#define  ROBUSTLMM_MISC_H__

#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <sstream>
#include <string>
#include <iostream>

template <typename T>
std::string to_string(T value)
{
  std::ostringstream os ;
  os << value ;
  return os.str() ;
}

#endif