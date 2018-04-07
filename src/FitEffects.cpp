#include "FitEffects.h"

#define DEBUG(STRING)                                          \
 Rcpp::Rcout << STRING << std::endl;

FitEffects::FitEffects(rlmerPredD* pp, rlmerResp* resp,
                       const double relativeTolerance, const int maxOperations) :
  SimpleIterativeFitter(pp->getEffects(), relativeTolerance, maxOperations),
  pp_(pp), resp_(resp), n_(pp_->n_), p_(pp_->p_), q_(pp_->q_),
  mat1_(p_+q_, n_+q_), mat2_(p_+q_, n_),
  invU_ey_(pp_->getInvU_e() * resp_->y_), W_(n_+q_)  {
  init();
}

FitEffects::FitEffects(rlmerPredDXPtr pp, rlmerRespXPtr resp,
           const double relativeTolerance, const int maxOperations) :
  SimpleIterativeFitter(pp->getEffects(), relativeTolerance, maxOperations),
  pp_(static_cast<rlmerPredD*> (R_ExternalPtrAddr(pp))),
  resp_(static_cast<rlmerResp*> (R_ExternalPtrAddr(resp))),
  n_(pp_->n_), p_(pp_->p_), q_(pp_->q_),
  mat1_(p_+q_, n_+q_), mat2_(p_+q_, n_),
  invU_ey_(pp_->getInvU_e() * resp_->y_), W_(n_+q_)  {
  init();
}

std::string FitEffects::init() {
  initMat1();
  // decomp_.cholmod().final_ll = 1; // force LL' decomposition.
  decomp_.analyzePattern(mat1_ * mat1_.adjoint());
  if (decomp_.info() != Eigen::Success) {
    return "CholeskyDecomposition.analyzePattern failed";
  }
  initMat2();
  return "";
}

void FitEffects::initMat1() {
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  const MatrixXd invU_eX = pp_->getInvU_eX();
  tripletList.reserve(invU_eX.size() +
    pp_->getInvU_btZtU_et().nonZeros() + pp_->getSqrtLambda_b().rows());
  // mat1_.block(0, 0, p_, n_) = pp_->getInvU_eX().adjoint();
  for (int i = 0, nRows = invU_eX.rows(), nCols = invU_eX.cols(); i < nCols; ++i) {
    for (int j = 0; j < nRows; ++j) {
      double value = invU_eX(j,i);
      if (value != 0.) {
        tripletList.push_back(T(i, j, value));
      }
    }
  }

  // mat1_.block(p_, 0, q_, n_) =  pp_->getInvU_btZtU_et();
  for (int k=0; k<pp_->getInvU_btZtU_et().outerSize(); ++k) {
    for (SpMatrixd::InnerIterator it(pp_->getInvU_btZtU_et(),k); it; ++it)
    {
      tripletList.push_back(T(p_ + it.row(), it.col(), it.value()));
    }
  }

  // mat1_.block(p_, n_, q_, q_) = pp_->getSqrtLambda_b().adjoint();
  const double* data = pp_->getSqrtLambda_b().diagonal().data();
  for (size_t i = 0, nRows = pp_->getSqrtLambda_b().rows(); i < nRows; ++i)
    tripletList.push_back(T(p_ + i, n_ + i, *(data + i)));

  mat1_.setFromTriplets(tripletList.begin(), tripletList.end());
}

void FitEffects::initMat2() {
  mat2_.block(0, 0, p_, n_) = pp_->getInvU_eX().adjoint();
  mat2_.block(p_, 0, q_, n_) = pp_->getInvU_btZtU_et();
}

FitEffects::~FitEffects() {}

const Fit<VectorXd>& FitEffects::fit() {
  return SimpleIterativeFitter::fit();
}

std::string FitEffects::doIteration() {
  setEffects(pp_, resp_, fit_.getValue());
  return computeAndSetNextValue();
}

std::string FitEffects::computeAndSetNextValue() {
  updateW();
  SpMatrixd mat(mat1_ * W_.cwiseSqrt().asDiagonal());
  VectorXd w_y(W_.segment(0, n_).cwiseProduct(invU_ey_));
  // FIXME need to get rid of this * transpose.
  decomp_.factorize(mat * mat.transpose());
  if (decomp_.info() != Eigen::Success) {
    return "CholeskyDecomposition.factorize failed";
  }
  VectorXd newValue(decomp_.solve(mat2_ * w_y));
  if (decomp_.info() != Eigen::Success) {
    return "CholeskyDecomposition.solve failed";
  }
  this->update(newValue);
  return "";
}

void FitEffects::updateW() {
  W_.segment(0, n_) = wgt_e(pp_, resp_);
  W_.segment(n_, q_) = wgt_b(pp_);
}

SpMatrixd FitEffects::getMat1_copy() const {
  return SpMatrixd(mat1_);
}

MatrixXd FitEffects::getMat2_copy() const {
  return MatrixXd(mat2_);
}

const VectorXd FitEffects::getInvU_ey() const {
  return invU_ey_;
}

const VectorXd FitEffects::getW() {
  updateW();
  return W_;
}

const VectorXd FitEffects::getNextValue() {
  computeAndSetNextValue();
  return fit_.getValue();
}

// RCPP_EXPOSED_CLASS(Fit);
// RCPP_EXPOSED_CLASS(FitEffects);
RCPP_MODULE(fitEffects_module) {

  class_<FitEffects>("FitEffects")
  .constructor<rlmerPredDXPtr, rlmerRespXPtr, double, double>()
  .method("mat1", &FitEffects::getMat1_copy)
  .method("mat2", &FitEffects::getMat2_copy)
  .method("invU_ey", &FitEffects::getInvU_ey)
  .method("W", &FitEffects::getW)
  .method("nextValue", &FitEffects::getNextValue)
  .method("fit", &FitEffects::fit)
  ;

}
