#if !defined  ROBUSTLMM_MISC_H__
#define  ROBUSTLMM_MISC_H__

#include <memory>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <Rcpp.h>
#include <sstream>
#include <string>
#include <iostream>

/*
#ifdef NDEBUG
#define R_ASSERT(EX)
#else
*/
#define R_ASSERT(EX)                                            \
  if (EX) {} else                                              \
  throw Rcpp::exception(tfm::format("Assertion '%s' failed at %s, line %i", \
    #EX, __FILE__, __LINE__).c_str())
// #endif

template <typename T> int sgn(T val);

#define REP_DO(TYPE_)                                          \
TYPE_ result(rep.size());                                      \
for (unsigned i = 0; i < rep.size(); i++) {                         \
  result[i] = input[rep[i] - 1];                               \
}                                                              \
return result;                                                 \

#define REP_DO_P(TYPE_)                                        \
TYPE_ result(rep.size());                                      \
for (unsigned i = 0; i < rep.size(); i++) {                         \
  result[i] = input(rep[i] - 1);                               \
}                                                              \
return result;

namespace std {
template<typename T>
std::string to_string(const T &n) {
  std::ostringstream s;
  s << n;
  return s.str();
}

template<typename T>
std::string to_string(const vector<T> &v) {
  std::stringstream s;
  for(size_t i = 0; i < v.size(); ++i)
  {
    if(i != 0)
      s << ", ";
    s << v[i];
  }
  return s.str();
}
} // namespace std

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

bool isDifferent(const double v1, const double v2, const double relativeTolerance);

#endif
