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

void warn(std::string msg);

template <typename T1>
inline void warn(const char* fmt, const T1& arg1) {
  // Rcpp::Rcout << "Formatting message" << std::endl;
  // This will cause a stackoverflow if -fno-elide-constructors is set
  std::string msg = tfm::format(fmt, arg1);
  // Rcpp::Rcout << "Issuing warning" << std::endl;
  warn(msg);
}

template <typename T1, typename T2>
inline void warn(const char* fmt, const T1& arg1, const T2& arg2) {
  // Rcpp::Rcout << "Formatting message" << std::endl;
  // This will cause a stackoverflow if -fno-elide-constructors is set
  std::string msg = tfm::format(fmt, arg1, arg2);
  // Rcpp::Rcout << "Issuing warning" << std::endl;
  warn(msg);
}

#endif
