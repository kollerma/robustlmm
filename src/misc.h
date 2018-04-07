#if !defined  ROBUSTLMM_MISC_H__
#define  ROBUSTLMM_MISC_H__

#include <memory>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <Rcpp.h>
#include <RcppEigen.h>
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

using Eigen::Lower;
using Eigen::ArrayXd;
using Eigen::LLT;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

typedef Eigen::Map<MatrixXd>                      MMap;
typedef Eigen::Map<VectorXd>                      MVec;
typedef Eigen::Map<VectorXi>                      MiVec;
typedef MatrixXd::Scalar                          Scalar;
typedef MatrixXd::Index                           Index;
typedef Eigen::SparseMatrix<double>               SpMatrixd;
typedef Eigen::MappedSparseMatrix<double>         MSpMatrixd;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic>    DMatrixXd;

typedef const Eigen::Ref<const Eigen::VectorXd>   ConstVectorXd;
typedef const Eigen::Ref<const Eigen::VectorXi>   ConstVectorXi;
typedef const Eigen::Ref<const Eigen::MatrixXd>   ConstMatrixXd;
// typedef const Eigen::Ref<const DMatrixXd> ConstDMatrixXd;
typedef const SpMatrixd ConstSpMatrixd;

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

template<class T> static T repDo(const T& input, ConstVectorXi& rep) {
  REP_DO(T);
}

inline static VectorXd repDo(const std::vector<double>& input, ConstVectorXi& rep) {
  REP_DO(VectorXd);
}

inline static VectorXd repDo(const MMap& input, ConstVectorXi& rep) {
  REP_DO_P(VectorXd);
}

inline static VectorXi repDo(ConstVectorXi& input, ConstVectorXi& rep) {
  REP_DO(VectorXi);
}

inline static VectorXi repDo(const MiVec& input, ConstVectorXi& rep) {
  REP_DO_P(VectorXi);
}


template<class T1, typename T2>
int vecApply(const std::vector<T1>& mVec, T2 (Eigen::PlainObjectBase<T1>::*fptr)() const);
template<class T1, class T2>
void setBlock(const T1& target, const int colOffset, const int rowOffset, const T2& source);
SpMatrixd bdiag(const std::vector<MatrixXd>& mVec);
MatrixXd zeroDiag(ConstMatrixXd& m);
void zeroNaNs(MatrixXd& m);

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

inline DMatrixXd crossprod(const DMatrixXd& A) {
  return DMatrixXd(A.diagonal().cwiseProduct(A.diagonal()));
}

inline MatrixXd crossprod(ConstMatrixXd& A) {
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
                      .rankUpdate(A.adjoint());
}

inline SpMatrixd crossprod(ConstSpMatrixd& A) {
  return A.adjoint() * A;
}

inline MatrixXd crossprod(ConstMatrixXd& A, ConstMatrixXd& B) {
  return A.adjoint() * B;
}

inline MatrixXd crossprod(ConstMatrixXd& A, ConstSpMatrixd& B) {
  return A.adjoint() * B;
}

inline MatrixXd crossprod(ConstSpMatrixd& A, ConstMatrixXd& B) {
  return A.adjoint() * B;
}

inline SpMatrixd crossprod(ConstSpMatrixd& A, ConstSpMatrixd& B) {
  return A.adjoint() * B;
}

inline MatrixXd tcrossprod(ConstMatrixXd& A) {
  int n(A.rows());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>()
                      .rankUpdate(A);
}

inline SpMatrixd tcrossprod(ConstSpMatrixd& A) {
  return A * A.adjoint();
}

inline MatrixXd tcrossprod(ConstMatrixXd& A, ConstMatrixXd& B) {
  return A * B.adjoint();
}

inline MatrixXd tcrossprod(ConstMatrixXd& A, ConstSpMatrixd& B) {
  return A * B.adjoint();
}

inline MatrixXd tcrossprod(ConstSpMatrixd& A, ConstMatrixXd& B) {
  return A * B.adjoint();
}

inline SpMatrixd tcrossprod(ConstSpMatrixd& A, ConstSpMatrixd& B) {
  return A * B.adjoint();
}

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
bool isDifferent(const MatrixXd& m1, const MatrixXd& m2, const double relativeTolerance);

#endif
