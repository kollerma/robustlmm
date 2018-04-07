#include "misc.h"

#define DEBUG(STRING)                                          \
// Rcpp::Rcout << STRING << std::endl;

template<class T1, typename T2>
int vecApply(const std::vector<T1>& mVec, T2 (Eigen::PlainObjectBase<T1>::*fptr)() const) {
  int sum = 0;
  for (unsigned i = 0; i < mVec.size(); i++) {
    sum += (mVec[i].*fptr)();
  }
  return sum;
}

template<class T1, class T2>
void setBlock(const T1& target, const unsigned colOffset, const unsigned rowOffset, const T2& source) {
  for (unsigned i = 0; i < source.cols(); i++) {
    for (unsigned j = 0; j < source.rows(); j++) {
      DEBUG("setBlock (" << j << ", " << i << ") value " << source.coeff(j, i))
      target->insert(rowOffset + j, colOffset + i) = source.coeff(j, i);
    }
  }
}

SpMatrixd bdiag(const std::vector<MatrixXd>& mVec) {
  DEBUG("bdiag 10 with vector of length " << mVec.size())
  unsigned ncols = vecApply(mVec, &MatrixXd::cols);
  unsigned nrows = vecApply(mVec, &MatrixXd::rows);
  DEBUG("bdiag 20 ncols = " << ncols << " and nrows = " << nrows)
  SpMatrixd result(ncols, nrows);
  result.reserve(VectorXi::Constant(ncols, nrows / mVec.size()));
  DEBUG("bdiag 30")
    unsigned colIndex = 0, rowIndex = 0;
  for (unsigned i = 0; i < mVec.size(); i++) {
    DEBUG("bdiag 40 (" << i << "): colIndex = " << colIndex << " rowIndex = " << rowIndex)
    setBlock(&result, colIndex, rowIndex, mVec[i]);
    colIndex += mVec[i].cols();
    rowIndex += mVec[i].rows();
  }
  DEBUG(result)
  return result;
}

MatrixXd zeroDiag(ConstMatrixXd& m) {
  R_ASSERT(m.rows() == m.cols());
  MatrixXd out(m);
  out.diagonal().setZero();
  return out;
}

void zeroNaNs(MatrixXd& m) {
  if (!m.hasNaN()) return;
  unsigned cols = m.cols();
  for (unsigned i = 0, rows = m.rows(); i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      if (ISNAN(m.coeff(i, j))) {
        m(i, j) = 0.;
      }
    }
  }
}

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

bool isDifferent(const MatrixXd& m1, const MatrixXd& m2, const double relativeTolerance) {
  if (m1.rows() != m2.rows() || m1.cols() != m2.cols()) {
    return true;
  }

  for (int row = 0, rows = m1.rows(); row < rows; ++row) {
    for (int col = 0, cols = m1.cols(); col < cols; ++col) {
      if (isDifferent(m1.coeff(row, col), m2.coeff(row, col), relativeTolerance)) {
        return true;
      }
    }
  }

  return(false);
}

