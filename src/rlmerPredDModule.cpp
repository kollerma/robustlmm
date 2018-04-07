#include "rlmerPredDModule.h"

using namespace Rcpp;

using std::invalid_argument;
using Eigen::LLT;
using Eigen::SimplicialLDLT;

#define DEBUG(STRING)                                          \
//  Rcpp::Rcout << STRING << std::endl;

#define DEBUG2(STRING)                                         \
//   Rcpp::Rcout << STRING << std::endl;

#define DEBUG0(STRING)                                         \
//  Rcpp::Rcout << STRING << std::endl;

namespace Rcpp {

template <> std::vector<PsiFuncXPtr> as(SEXP arg) {
  List list(as<List>(arg));
  std::vector<PsiFuncXPtr> result;
  result.reserve(list.size());
  for (int i = 0; i < list.size(); i++) {
    PsiFuncXPtr p(as<PsiFuncXPtr>(list[i]));
    result.push_back(p);
  }
  return result;
}

}

// class rlmerPredD

rlmerPredD::rlmerPredD(List args, MMap X, MSpMatrixd Zt)
  : X_(X), Zt_(Zt),
    n_(X_.rows()), p_(X_.cols()), q_(Zt_.rows()),
    maxOperations_(args["maxOperations"]),
    relativeTolerance_(args["relativeTolerance"]),
    indexMapper_(args), Lambdat_(as<MSpMatrixd>(args["Lambdat"])),
    lower_(as<MVec>(args["lower"])), v_e_(as<MVec>(args["v_e"])),
    theta_(as<MVec>(args["theta"])), beta_(as<MVec>(args["beta"])),
    b_s_(as<MVec>(args["b.s"])), sigma_(as<MVec>(args["sigma"])),
    Lind_(as<MiVec>(args["Lind"])),
    rho_e_(as<PsiFuncXPtr>(args["rho_e"])),
    rhoSigma_e_(as<PsiFuncXPtr>(args["rhoSigma_e"])),
    rho_b_(as<std::vector<PsiFuncXPtr> >(args["rho_b"])),
    rhoSigma_b_(as<std::vector<PsiFuncXPtr> >(args["rhoSigma_b"])),
    U_e_(v_e_.asDiagonal()), invU_e_(U_e_.inverse()), V_e_(crossprod(U_e_)),
    b_r_(Lambdat_.adjoint() * b_s_),
    integration_(*(new DqagIntegration())),
    calledInit_(false), setM_(false), setUnsc_(false) {
  // FIXME add dimension checks
  /*
   * stopifnot(length(theta) > 0L,
    length(Lind) > 0L,
    all(sort(unique(Lind)) == seq_along(theta)),
    length(sigma) == 1,
    sigma > 0,
    length(lower) == length(theta),
    length(v_e) == n)
   */
  setZeroB();
}

rlmerPredD::~rlmerPredD() {
  delete &integration_;
}

void rlmerPredD::initMatrices() {
  if (calledInit_) return;
  invU_eX_ = invU_e_ * X_;
  invU_eZ_ = invU_e_ * Zt_.adjoint();
  this->initRho();
  calledInit_ = true;
  this->updateMatrices();
}

void rlmerPredD::initRho() {
  D_e_ = VectorXd(n_).setOnes().asDiagonal() * rho_e_->EDpsi();
  // calculate D_b_ per blocktype
  VectorXd tmp = btapplyV(rho_b_, &calcED_re, indexMapper_.getBBlockMap());
  D_b_ = tmp.asDiagonal();
  Lambda_b_ = D_b_.inverse() * rho_e_->EDpsi();
  sqrtLambda_b_ = Lambda_b_.diagonal().cwiseSqrt().asDiagonal();
  Lambda_bD_b_ = VectorXd(q_).setOnes().asDiagonal() * rho_e_->EDpsi();
  sqrtD_e_ = D_e_.diagonal().cwiseSqrt().asDiagonal();
  MatrixXd tmp0(sqrtD_e_ * invU_eX_);
  M_XX_ = crossprod(tmp0);
  lltM_XX_ = LLT<MatrixXd>(M_XX_);
  SpMatrixd tmp1(sqrtD_e_ * invU_eZ_);
  M_XZ0_ = tmp0.adjoint() * tmp1;
  M_ZZ0_ = crossprod(tmp1);
  LLT<MatrixXd> lltM_XX(M_XX_);
  invM_XXMZZ0_ = lltM_XX.solve(M_XZ0_);
  M_ZZ0_sub_M_ZX0invM_XXMZZ0_ = -1 * M_XZ0_.adjoint() * invM_XXMZZ0_;
  M_ZZ0_sub_M_ZX0invM_XXMZZ0_ += M_ZZ0_;
  Epsi_bbt_ = btapplyM(rho_b_, &calcEpsi_bbt, indexMapper_.getInd());
  Epsi_bpsi_bt_ = btapplyM(rho_b_, &calcEpsi_bpsi_bt, indexMapper_.getInd());
  Epsi2_b_ = Epsi_bpsi_bt_.diagonal();
}

void rlmerPredD::updateMatrices() {
  this->initMatrices();
  invU_btZtU_et_ = Lambdat_ * invU_eZ_.adjoint();
}

void rlmerPredD::setZeroB() {
  indexMapper_.updateBlockTypesDropped(theta_);
}

VectorXi rlmerPredD::getZeroB() const {
  VectorXi zeroB(indexMapper_.getNumberOfRandomEffects());
  zeroB.setZero();
  for (const auto &randomEffect : indexMapper_.getRandomEffects()) {
    if (randomEffect->isBlockTypeDropped()) {
      zeroB[randomEffect->getIndex()] = 1;
    }
  }
  return zeroB;
}

void rlmerPredD::setSigma(double sigma) {
  sigma_[0] = sigma;
}

void rlmerPredD::setU(const VectorXd& u) {
  setB_s(u);
}

void rlmerPredD::setB_s(const VectorXd& b_s) {
  std::copy(b_s.data(), b_s.data() + b_s.size(),
            b_s_.data());
  b_r_ = Lambdat_.adjoint() * b_s;
}

void rlmerPredD::setB_r(const VectorXd& b_r) {
  std::copy(b_r.data(), b_r.data() + b_r.size(),
            b_r_.data());
  b_s_ = stdB(b_r);
}

VectorXd rlmerPredD::stdB(ConstVectorXd& vector) {
  VectorXd out(vector);
  if (allRandomEffectsDropped()) {
    out.setZero();
  } else {
    MatrixXd U_b = getU_b_copy();
    if (indexMapper_.isAnyBlockTypeDropped()) {
      for (const auto &randomEffect : indexMapper_.getRandomEffects()) {
        if (randomEffect->isBlockTypeDropped()) {
          U_b(randomEffect->getIndex(), randomEffect->getIndex()) = 1.;
        }
      }
    }
    U_b.triangularView<Eigen::Lower>().solveInPlace(out);
    if (indexMapper_.isAnyBlockTypeDropped()) {
      for (const auto &randomEffect : indexMapper_.getRandomEffects()) {
        if (randomEffect->isBlockTypeDropped()) {
          out[randomEffect->getIndex()] = 0.;
        }
      }
    }
  }
  return out;
}

void rlmerPredD::setBeta(const VectorXd& beta) {
  if (beta.size() != beta_.size()) {
    Rcpp::Rcout << "(" << beta.size() << "!=" <<
      beta_.size() << ")" << std::endl;
    throw invalid_argument("beta size mismatch");
  }
  std::copy_n(beta.data(), beta.size(), beta_.data());
}

void rlmerPredD::setTheta(const VectorXd& theta) {
  // follow lme4 closely here
  if (theta.size() != theta_.size()) {
    Rcpp::Rcout << "(" << theta.size() << "!=" <<
      theta_.size() << ")" << std::endl;
    throw invalid_argument("theta size mismatch");
  }
  // update theta
  std::copy_n(theta.data(), theta.size(), theta_.data());
  // update Lambdat
  const int *lipt = Lind_.data();
  double *LamX = Lambdat_.valuePtr(), *thpt = theta_.data();
  for (int i = 0; i < Lind_.size(); ++i) {
    LamX[i] = thpt[lipt[i] - 1];
  }
  // and finally:
  setZeroB();
  this->resetCaches();
  this->updateMatrices();
}

const MSpMatrixd& rlmerPredD::Lambdat() const {
  return Lambdat_;
}

M* rlmerPredD::Mobj() {
  if (!setM_) {
    M_.init(this);
    setM_ = true;
  }
  return &M_;
}

const M& rlmerPredD::getM() {
  Mobj();
  return M_;
}

SpMatrixd rlmerPredD::getU_b() const {
  return Lambdat_.adjoint();
}

SpMatrixd rlmerPredD::getU_b_copy() const {
  return SpMatrixd(getU_b());
}

MatrixXd rlmerPredD::getU_e_copy() const {
  return MatrixXd(U_e_);
}

MatrixXd rlmerPredD::getV_e_copy() const {
  return MatrixXd(V_e_);
}

const MVec& rlmerPredD::getTheta() const {
  return theta_;
}

const MVec& rlmerPredD::getBeta() const {
  return beta_;
}

const MVec& rlmerPredD::getB_s() const {
  return b_s_;
}

double rlmerPredD::getSigma() const {
  return sigma_[0];
}

const IndexMapper& rlmerPredD::getIndexMapper() const {
  return indexMapper_;
}

const VectorXi& rlmerPredD::getBBlockMap() const {
  return indexMapper_.getBBlockMap();
}

MatrixXd rlmerPredD::getM_XZ0_copy() const {
  return MatrixXd(M_XZ0_);
}

SpMatrixd rlmerPredD::getM_ZZ0_copy() const {
  return SpMatrixd(M_ZZ0_);
}

MatrixXd rlmerPredD::getM_ZZ0_sub_M_ZX0invM_XXMZZ0_copy() const {
  return MatrixXd(M_ZZ0_sub_M_ZX0invM_XXMZZ0_);
}

MatrixXd rlmerPredD::getLambda_b_copy() const {
  return MatrixXd(Lambda_b_);
}

MatrixXd rlmerPredD::getLambda_bD_b_copy() const {
  return MatrixXd(Lambda_bD_b_);
}

SpMatrixd rlmerPredD::getEpsi_bbt_copy() const {
  return SpMatrixd(Epsi_bbt_);
}

SpMatrixd rlmerPredD::getEpsi_bpsi_bt_copy() const {
  return SpMatrixd(Epsi_bpsi_bt_);
}

const VectorXd& rlmerPredD::getEpsi2_b() const {
  return Epsi2_b_;
}

const SpMatrixd& rlmerPredD::getInvU_btZtU_et() const {
  return invU_btZtU_et_;
}

SpMatrixd rlmerPredD::getInvU_btZtU_et_copy() const {
  return SpMatrixd(invU_btZtU_et_);
}

const MatrixXd& rlmerPredD::getInvU_eX() const {
  return invU_eX_;
}

MatrixXd rlmerPredD::getInvU_eX_copy() const {
  return MatrixXd(invU_eX_);
}

const DMatrixXd& rlmerPredD::getSqrtLambda_b() const {
  return sqrtLambda_b_;
}

const DMatrixXd& rlmerPredD::getLambda_b() const {
  return Lambda_b_;
}

const DMatrixXd& rlmerPredD::getInvU_e() const {
  return invU_e_;
}

VectorXd rlmerPredD::getMu() const {
  return crossprod(Zt_, b_r_) + X_* beta_;
}

VectorXd rlmerPredD::getDist_b() const {
  VectorXd db = getD_k();
  if (indexMapper_.isAnyBlockTypeNonDiagonal()) {
    for (const auto &block : indexMapper_.getBlocks()) {
      if (block->isNonDiagonalBlockType()) {
        db[block->getIndex()] = std::sqrt(db[block->getIndex()]);
      }
    }
  }
  return repDo(db, indexMapper_.getK());
}

VectorXd rlmerPredD::getD_k() const {
  VectorXd dk(indexMapper_.getNumberOfBlocks());
  dk.setZero();
  for (const auto &randomEffect : indexMapper_.getRandomEffects()) {
    if (randomEffect->isNonDiagonalBlockType()) {
      dk[randomEffect->getBlockIndex()->getIndex()] +=
        b_s_[randomEffect->getIndex()] * b_s_[randomEffect->getIndex()] /
          (sigma_[0] * sigma_[0]);
    } else {
      dk[randomEffect->getBlockIndex()->getIndex()] =
        b_s_[randomEffect->getIndex()] / sigma_[0];
    }
  }
  return dk;
}

VectorXd rlmerPredD::getEffects() const {
  VectorXd effects(p_ + q_);
  effects.block(0, 0, p_, 1) = beta_;
  effects.block(p_, 0, q_, 1) = b_s_;
  if (indexMapper_.isAnyBlockTypeDropped()) {
    for (const auto &randomEffect : indexMapper_.getRandomEffects()) {
      if (randomEffect->isBlockTypeDropped()) {
        effects(p_ + randomEffect->getIndex()) = 0.;
      }
    }
  }
  return effects;
}

const PsiFuncXPtr& rlmerPredD::getRho_e() const {
  return rho_e_;
}

const PsiFuncXPtr& rlmerPredD::getRhoSigma_e() const {
  return rhoSigma_e_;
}

const std::vector<PsiFuncXPtr>& rlmerPredD::getRho_b() const {
  return  rho_b_;
}

const std::vector<PsiFuncXPtr>& rlmerPredD::getRhoSigma_b() const {
  return rhoSigma_b_;
}

MatrixXd rlmerPredD::unsc() {
  if (setUnsc_)
    return MatrixXd(unsc_);
  if (allRandomEffectsDropped()) {
    unsc_ = (sqrtD_e_ * invU_eX_).adjoint();
    lltM_XX_.solveInPlace(unsc_);
    unsc_ = tcrossprod(unsc_) * (rho_e_->Epsi2() / rho_e_->EDpsi());
  } else {
    unsc_ = (Mobj()->BB() - crossprod(Mobj()->bB(), Lambda_bD_b_ * Mobj()->bB())) *
      (rho_e_->Epsi2() / rho_e_->EDpsi());
    unsc_ += crossprod(Epsi2_b_.cwiseSqrt().asDiagonal() *
      Lambda_b_.toDenseMatrix() * Mobj()->bB());
  }
  // Force matrix to be symmetric
  unsc_ = (unsc_ + unsc_.transpose()).eval() / 2;
  setUnsc_ = true;
  return MatrixXd(unsc_);
}

const VectorXd& rlmerPredD::getB() const {
  return b_r_;
}

void rlmerPredD::resetCaches() {
  setM_ = false;
  setUnsc_ = false;
}

template<class T, class O>
std::vector<O> rlmerPredD::btapply(const std::vector<T>& input, O (*fun)(const T&)) {
  R_ASSERT(input.size() == indexMapper_.getNumberOfBlockTypes());
  std::vector<O> result(indexMapper_.getNumberOfBlockTypes());
  for (const auto &blockType : indexMapper_.getBlockTypes()) {
     result[blockType->getIndex()] = (*fun)(input[blockType->getIndex()]);
  }
  return result;
}

template<class T, class O>
std::vector<O> rlmerPredD::btapply(const std::vector<T>& input, O (*fun)(const T&, unsigned s)) {
  R_ASSERT(input.size() == indexMapper_.getNumberOfBlockTypes());
  std::vector<O> result(indexMapper_.getNumberOfBlockTypes());
  for (const auto &blockType : indexMapper_.getBlockTypes()) {
    result[blockType->getIndex()] =
      (*fun)(input[blockType->getIndex()], blockType->dim_);
  }
  return result;
}

template<class T, class O>
std::vector<O> rlmerPredD::btapply(const std::vector<T>& input, O (*fun)(const T&, unsigned s, Integration* const integration)) {
  R_ASSERT(input.size() == indexMapper_.getNumberOfBlockTypes());
  std::vector<O> result(indexMapper_.getNumberOfBlockTypes());
  for (const auto &blockType : indexMapper_.getBlockTypes()) {
    result[blockType->getIndex()] =
      (*fun)(input[blockType->getIndex()], blockType->dim_, &integration_);
  }
  return result;
}

template<class T, class Fun>
VectorXd rlmerPredD::btapplyV(const std::vector<T>& input, Fun fun, ConstVectorXi& rep) {
  std::vector<double> result0 = btapply(input, fun);
  return repDo(result0, rep);
}

template<class T, class Fun>
SpMatrixd rlmerPredD::btapplyM(const std::vector<T>& input, Fun fun, ConstVectorXi& rep) {
  std::vector<MatrixXd> result0 = btapply(input, fun);
  std::vector<MatrixXd> result1 = repDo(result0, rep);
  SpMatrixd result2(bdiag(result1));
  return result2;
}

bool rlmerPredD::allRandomEffectsDropped() const {
  return indexMapper_.areAllBlockTypesDropped();
}

bool rlmerPredD::isBlockDropped(const BlockIndex* const block) const {
  return block->getBlockTypeIndex()->isDropped();
}

// end class rlmerPredD

// class rlmerPredD_DAS

rlmerPredD_DAS::rlmerPredD_DAS(List args, MMap X, MSpMatrixd Zt) :
  rlmerPredD(args, X, Zt),
  initTau_e_(false), setTau_e_(false), initTau_b_(false), setTau_b_(false) {}

void rlmerPredD_DAS::updateMatrices() {
  rlmerPredD::updateMatrices();
  if (allRandomEffectsDropped()) {
    A_ = invU_eX_ * lltM_XX_.solve(invU_eX_.adjoint());
    Kt_ = MatrixXd(n_, q_).setZero();
    LLT<MatrixXd> llt_D_b(D_b_);
    L_ = llt_D_b.solve(MatrixXd::Identity(q_, q_));
  } else {
    MatrixXd tmp2(tcrossprod(invU_eX_, crossprod(invU_btZtU_et_, Mobj()->bB())));
    MatrixXd tmp3(crossprod(invU_btZtU_et_, Mobj()->bb()));
    A_ = tcrossprod(invU_eX_ * Mobj()->BB(), invU_eX_) + tmp2 + tmp2.adjoint() +
      tmp3 * invU_btZtU_et_;
    Kt_ = (tcrossprod(invU_eX_, Mobj()->bB()) + tmp3) * -1.;
    L_ = Mobj()->bb() * Lambda_b_;
  }
}

void rlmerPredD_DAS::setTheta(const VectorXd& theta) {
  if (getIndexMapper().isAnyBlockTypeDropped()) {
    initTau_b_ = false;
  }
  rlmerPredD::setTheta(theta);
  this->updateMatrices();
  setTau_e_ = false;
}

MatrixXd rlmerPredD_DAS::getB() {
  return Kt_ * Lambda_b_;
}

MatrixXd rlmerPredD_DAS::getK_copy() {
  return MatrixXd(Kt_.adjoint());
}

MatrixXd rlmerPredD_DAS::getA_copy() {
  return MatrixXd(A_);
}

MatrixXd& rlmerPredD_DAS::getA() {
  return A_;
}

MatrixXd rlmerPredD_DAS::getKt_copy() {
  return MatrixXd(Kt_);
}

MatrixXd rlmerPredD_DAS::getL_copy() {
  return MatrixXd(L_);
}

SpMatrixd& rlmerPredD_DAS::getTau_b() {
  if (setTau_b_) return Tau_b_;
  Tau_b_ = SpMatrixd(this->calcTau_b());
  setTau_b_ = true;
  initTau_b_ = true;
  return Tau_b_;
}

SpMatrixd rlmerPredD_DAS::getTau_b_copy() {
  return SpMatrixd(getTau_b());
}

VectorXd rlmerPredD_DAS::s(MatrixXd m1, MatrixXd m2) {
  zeroNaNs(m1);
  zeroNaNs(m2);
  VectorXd result = rho_e_->Epsi2() * (m1.rowwise().squaredNorm());
  if (!allRandomEffectsDropped()) {
    result += m2.array().square().matrix().eval() * Epsi_bpsi_bt_.diagonal();
  }
  return result.cwiseSqrt();
}

VectorXd rlmerPredD_DAS::s_e() {
  return s(zeroDiag(A_), this->getB());
}

std::vector<MatrixXd> rlmerPredD_DAS::s_b() {
  std::vector<MatrixXd> result(getIndexMapper().getNumberOfBlocks());
  MatrixXd tmpEL(L_ * Epsi_bpsi_bt_);
  for (const auto &block : getIndexMapper().getBlocks()) {
    int blockSize = block->getBlockTypeDimension();
    MatrixXd mat(blockSize, blockSize);
    if (block->isBlockTypeDropped()) {
      result[block->getIndex()] = mat.setZero();
    } else {
      R_ASSERT(isContiguousVector(block->getRandomEffects()));
      int blockStart = block->getRandomEffects()[0]->getIndex();
      // object@rho.e@Epsi2() * tcrossprod(object@pp$K()[bidx, , drop=FALSE])
      mat = crossprod(Kt_.block(0, blockStart, Kt_.rows(), blockSize)) * rho_e_->Epsi2() +
        // tcrossprod(object@pp$L[bidx, !bidx, drop=FALSE], tmpEL[bidx, !bidx, drop=FALSE])
        tcrossprod(L_.block(blockStart, 0, blockSize, L_.cols()),
                   tmpEL.block(blockStart, 0, blockSize, tmpEL.cols())) -
        tcrossprod(L_.block(blockStart, blockStart, blockSize, blockSize),
                   tmpEL.block(blockStart, blockStart, blockSize, blockSize));
     LLT<MatrixXd> llt(mat);
     // TODO should have a fallback to sqrt(diag) if this fails as in R?!
     result[block->getIndex()] = llt.matrixU();
    }
  }
  return(result);
}

VectorXd& rlmerPredD_DAS::getTau_e() {
  if (setTau_e_) return tau_e_;
  tau_e_ = VectorXd(this->calcTau_e());
  setTau_e_ = true;
  initTau_e_ = true;
  return tau_e_;
}

VectorXd rlmerPredD_DAS::getTau_e_copy() {
  return VectorXd(getTau_e());
}

void rlmerPredD_DAS::resetCaches() {
  rlmerPredD::resetCaches();
  setTau_e_ = false;
  setTau_b_ = false;
}

void rlmerPredD_DAS::initRho() {
  rlmerPredD::initRho();
  kappa_e_ = calcKappaTau(rhoSigma_e_, 1, &integration_);
  std::vector<double> kappa_b = btapply(rhoSigma_b_, &calcKappaTau);
  kappa_b_ = VectorXd(kappa_b.size());
  for (unsigned i = 0; i < kappa_b.size(); ++i) {
    kappa_b_[i] = kappa_b[i];
  }
}

double rlmerPredD_DAS::getKappa_e() {
  return kappa_e_;
}

VectorXd& rlmerPredD_DAS::getKappa_b() {
  return kappa_b_;
}

VectorXd rlmerPredD_DAS::getKappa_b_copy() {
  return VectorXd(this->getKappa_b());
}

VectorXd rlmerPredD_DAS::getInitTau_e() {
  if (initTau_e_) return tau_e_;
  return rlmerPredD_DAS::calcTau_e();
}

SpMatrixd rlmerPredD_DAS::getInitTau_b() {
  if (initTau_b_) {
    return Tau_b_;
  }
  return rlmerPredD_DAS::calcTau_b();
}

VectorXd rlmerPredD_DAS::calcTau_e() {
  MatrixXd B(getB());
  MatrixXd Tau_e(V_e_.toDenseMatrix() - rho_e_->EDpsi() * (A_.adjoint() + A_) +
    rho_e_->Epsi2() * tcrossprod(A_) +
    B * tcrossprod(Epsi_bpsi_bt_, B));
  return Tau_e.diagonal().cwiseSqrt();
}

SpMatrixd rlmerPredD_DAS::calcTau_b() {
  MatrixXd tmp(L_ * Epsi_bbt_);
  MatrixXd Tfull(MatrixXd::Identity(q_, q_) - tmp - tmp.adjoint() +
    rho_e_->Epsi2() * crossprod(Kt_) + L_ * crossprod(Epsi_bpsi_bt_, L_));
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(q_ + std::max(((int) Lambdat_.nonZeros() - (int) q_) * 2, (int) 0));
  for (const auto &block : getIndexMapper().getBlocks()) {
    for (const RandomEffectIndex* randomEffect1 : block->getRandomEffects()) {
      int randomEffectIndex1 = randomEffect1->getIndex();
      if (block->isDiagonalBlockType()) {
        tripletList.push_back(T(randomEffectIndex1, randomEffectIndex1,
                              Tfull.coeff(randomEffectIndex1, randomEffectIndex1)));
      } else {
        for (const RandomEffectIndex* randomEffect2 : block->getRandomEffects()) {
          int randomEffectIndex2 = randomEffect2->getIndex();
          tripletList.push_back(T(randomEffectIndex1, randomEffectIndex2,
                                Tfull.coeff(randomEffectIndex1, randomEffectIndex2)));
        }
      }
    }
  }
  SpMatrixd Tau_b(q_, q_);
  Tau_b.setFromTriplets(tripletList.begin(), tripletList.end());
  return(Tau_b);
}

// end class rlmerPredD_DAS

// class M

M::M() {}

void M::init(const rlmerPredD* const obj) {
  int p(obj->p_), q(obj->q_);
  if (obj->allRandomEffectsDropped()) {
    bbinv_ = obj->Lambda_bD_b_;
    lltM_bbinv_ = LLT<MatrixXd>(bbinv_);
    XZ_.setZero(p, q);
    bB_.setZero(q, p);
  } else {
    bbinv_ = obj->Lambdat_ * (obj->M_ZZ0_sub_M_ZX0invM_XXMZZ0_) * obj->Lambdat_.adjoint() +
      obj->Lambda_bD_b_.toDenseMatrix();
    lltM_bbinv_ = LLT<MatrixXd>(bbinv_);
    XZ_ = obj->M_XZ0_ * obj->Lambdat_.adjoint();
    MatrixXd M_XZM_XXinv(obj->lltM_XX_.solve(XZ_).transpose());
    bB_ = -lltM_bbinv_.solve(M_XZM_XXinv);
  }
  BB_ = obj->lltM_XX_.solve(MatrixXd::Identity(p, p) - XZ_ * bB_);
  bb_ = lltM_bbinv_.solve(MatrixXd::Identity(q, q));
}

M::M(M* o) : XZ_(o->XZ()), bb_(o->bb()), bbinv_(o->bbinv()),
  bB_(o->bB()), BB_(o->BB()), lltM_bbinv_(o->lltM_bbinv()) {}

const MatrixXd& M::XZ() { return XZ_; }
const MatrixXd& M::bb() { return bb_; }
const MatrixXd& M::bB() { return bB_; }
MatrixXd M::Bb() { return bB_.transpose(); }
const MatrixXd& M::BB() { return BB_; }
const MatrixXd& M::bbinv() { return bbinv_; }

const LLT<MatrixXd>& M::lltM_bbinv() { return lltM_bbinv_; }

M::~M() {
  DEBUG0("Called ~M() for " << this)
}

// end class M

// class rlmerPredD_DAStau

Expectation2d* createExpectation2d(int nodes) {
  if (nodes > 0)
    return new GaussianQuadratureNormalExpectation2d(nodes);
  return new DqagNormalExpectation2d();
}

ExpectationNd* createExpectationNd(int nodes) {
  if (nodes > 0) {
    return new GaussianQuadratureNdNormalExpectation(nodes);
  }
  return new HcubatureNormalExpectation(0, 1e-5, 1e-5);
}

rlmerPredD_DAStau::rlmerPredD_DAStau(Rcpp::List args, MMap X, MSpMatrixd Zt) :
  rlmerPredD_DAS(args, X, Zt),
  expectation2d_(createExpectation2d(args["nodes"])),
  expectationNd_(createExpectationNd(args["nodes"])) {}

rlmerPredD_DAStau::~rlmerPredD_DAStau() {
  delete expectation2d_;
  delete expectationNd_;
}

VectorXd rlmerPredD_DAStau::calcTau_e() {
  VectorXd initialValue(getInitTau_e());
  return calcTau(initialValue, relativeTolerance_, maxOperations_,
                 getA().diagonal(), s_e(), getKappa_e(),
                 getRho_e(), getRhoSigma_e(), getExpectation2d());
}

SpMatrixd rlmerPredD_DAStau::calcTau_b() {
  SpMatrixd value(this->getInitTau_b());
  std::vector<MatrixXd> skbs(this->s_b());

  for (const auto &blockType : getIndexMapper().getBlockTypes()) {
    this->calcTau_b(blockType.get(), value, skbs);
  }
  return value;
}

void rlmerPredD_DAStau::calcTau_b(const BlockTypeIndex* const blockType, SpMatrixd& value,
                                  const std::vector<MatrixXd>& skbs) {
  if (blockType->isDiagonal()) {
    this->calcTau_b_diagonal(blockType, value, skbs);
  } else {
    this->calcTau_b_non_diagonal(blockType, value, skbs);
  }
}

void rlmerPredD_DAStau::calcTau_b_diagonal(const BlockTypeIndex* const blockType, SpMatrixd& value,
                                           const std::vector<MatrixXd>& skbs) {
  R_ASSERT(blockType->isDiagonal());
  if (blockType->isDropped()) {
    for (const BlockIndex* block : blockType->getBlocks()) {
      if (value.coeff(block->getIndex(), block->getIndex()) != 0.) {
        value.coeffRef(block->getIndex(), block->getIndex()) = 0.;
      }
    }
    return;
  }

  int nBlocks(blockType->getNumberOfBlocks()), i(0);
  // DEBUG("nBlocks " << nBlocks << " for L_ (" << L_.rows() << "x" << L_.cols() << ") and skbs ("
  //                  << skbs.size() << ")")
  VectorXd a(nBlocks), s(nBlocks);
  for (const BlockIndex* block : blockType->getBlocks()) {
    int bidx = block->getRandomEffects()[0]->getIndex();
    // DEBUG("Accessing bidx = " << bidx << " block index = " << block->getIndex())
    a[i] = L_(bidx, bidx);
    s[i] = skbs[block->getIndex()](0, 0);
    ++i;
  }
  R_ASSERT(isContiguousVector(blockType->getRandomEffects()));
  VectorXd valueForBlockType(value.diagonal().segment(blockType->getRandomEffects()[0]->getIndex(), nBlocks));
  valueForBlockType = calcTau(valueForBlockType, relativeTolerance_, maxOperations_,
                              a, s, getKappa_b()[blockType->getIndex()], getRho_b()[blockType->getIndex()],
                              getRhoSigma_b()[blockType->getIndex()], getExpectation2d());
  i = 0;
  for (const BlockIndex* block : blockType->getBlocks()) {
    int bidx = block->getRandomEffects()[0]->getIndex();
    double tmp = valueForBlockType[i];
    value.coeffRef(bidx, bidx) = tmp * tmp;
    ++i;
  }
}

void rlmerPredD_DAStau::calcTau_b_non_diagonal(const BlockTypeIndex* const blockType, SpMatrixd& value,
                                               const std::vector<MatrixXd>& skbs) {
  R_ASSERT(blockType->isNonDiagonal());
  R_ASSERT(isContiguousVector(blockType->getRandomEffects()));
  if (blockType->isDropped()) {
    for (const BlockIndex* block : blockType->getBlocks()) {
      for (const RandomEffectIndex* randomEffect1 : block->getRandomEffects()) {
        int index1 = randomEffect1->getIndex();
        for (const RandomEffectIndex* randomEffect2 : block->getRandomEffects()) {
          if (value.coeff(index1, randomEffect2->getIndex()) != 0.) {
            value.coeffRef(index1, randomEffect2->getIndex()) = 0.;
          }
        }
      }
    }
    return;
  }

  int status = calcTauNonDiag(blockType, value, relativeTolerance_ * 10, maxOperations_,
                              kappa_b_, L_, skbs, rho_b_[blockType->getIndex()],
                              rhoSigma_b_[blockType->getIndex()], getExpectationNd());
  if (status == 0) {
    return;
  }

  DEBUG("calcTauNonDiag failed with status " << status)
  warn("Error calcTau_b_non_diagonal for blockType %i. Setting Tau_b to DASvar version.",
       blockType->getIndex() + 1);
  SpMatrixd Tbk(rlmerPredD_DAS::calcTau_b());
  for (const BlockIndex* block : blockType->getBlocks()) {
    for (const RandomEffectIndex* randomEffect1 : block->getRandomEffects()) {
      int index1 = randomEffect1->getIndex();
      for (const RandomEffectIndex* randomEffect2 : block->getRandomEffects()) {
        value.coeffRef(index1, randomEffect2->getIndex()) =
          Tbk.coeff(index1, randomEffect2->getIndex());
      }
    }
  }
}

Expectation2d* const rlmerPredD_DAStau::getExpectation2d() {
  return expectation2d_;
}

ExpectationNd* const rlmerPredD_DAStau::getExpectationNd() {
  return expectationNd_;
}

// end class rlmerPredD_DAStau

RCPP_EXPOSED_CLASS(M)
RCPP_EXPOSED_CLASS(rlmerPredD)
RCPP_MODULE(rlmerPredD_module) {

  class_<rlmerPredD>("rlmerPredD")
  .constructor<Rcpp::List, MMap, MSpMatrixd>()
  .method("initMatrices", &rlmerPredD::initMatrices, "Initialise Matrices, this needs to be called once")
  .method("setSigma", &rlmerPredD::setSigma, "Set sigma to this value")
  .method("setU", &rlmerPredD::setU, "Set spherical random effects to this value, updating random effects accordingly")
  .method("setB_s", &rlmerPredD::setB_s, "Set spherical random effects to this value, updating random effects accordingly")
  .method("setB_r", &rlmerPredD::setB_r, "Set random effects to this value, updating spherical random effects accordingly")
  .method("setBeta", &rlmerPredD::setBeta, "Set fixed effects estimates")
  .method("setTheta", &rlmerPredD::setTheta, "Set theta and update other variables accordingly")
  .method("Lambdat", &rlmerPredD::Lambdat, "Build and return Lambdat")
  .method("M", &rlmerPredD::getM, "Build and return class M that holds helper matrices M_** (Koller 2013, Appendix C.2)")
  .method("unsc", &rlmerPredD::unsc, "the unscaled variance-covariance matrix of the fixed-effects parameters")
  .method("b", &rlmerPredD::getB, "random effects accessor")
  .method("zeroB", &rlmerPredD::getZeroB, "zeroB accessor")
  .method("U_b", &rlmerPredD::getU_b_copy, "U_b accessor")
  .method("U_e", &rlmerPredD::getU_e_copy, "U_e accessor")
  .method("V_e", &rlmerPredD::getV_e_copy, "V_e accessor")
  .method("theta", &rlmerPredD::getTheta, "theta accessor")
  .method("beta", &rlmerPredD::getBeta, "beta accessor")
  .method("b_s", &rlmerPredD::getB_s, "b_s accessor")
  .method("effects", &rlmerPredD::getEffects, "vector of beta and b_s combined")
  .method("sigma", &rlmerPredD::getSigma, "sigma accessor")
  .method("M_XZ0", &rlmerPredD::getM_XZ0_copy, "M_XZ0 accessor")
  .method("M_ZZ0", &rlmerPredD::getM_ZZ0_copy, "M_ZZ0 accessor")
  .method("M_ZZ0_sub_M_ZX0invM_XXMZZ0", &rlmerPredD::getM_ZZ0_sub_M_ZX0invM_XXMZZ0_copy, "M_ZZ0_sub_M_ZX0invM_XXMZZ0 accessor")
  .method("Lambda_b", &rlmerPredD::getLambda_b_copy, "Lambda_b accessor")
  .method("Lambda_bD_b", &rlmerPredD::getLambda_bD_b_copy, "Lambda_bD_b accessor")
  .method("Epsi_bbt", &rlmerPredD::getEpsi_bbt_copy, "Epsi_bbt accessor")
  .method("Epsi_bpsi_bt", &rlmerPredD::getEpsi_bpsi_bt_copy, "Epsi_bpsi_bt accessor")
  .method("Epsi2_b", &rlmerPredD::getEpsi2_b, "Epsi2_b accessor")
  .method("invU_btZtU_et", &rlmerPredD::getInvU_btZtU_et_copy, "invU_btZtU_et accessor")
  .method("invU_eX", &rlmerPredD::getInvU_eX_copy, "invU_eX accessor")
  .method("mu", &rlmerPredD::getMu, "compute mu")
  .method("distB", &rlmerPredD::getDist_b, "compute dist b")
  .method("bBlockMap", &rlmerPredD::getBBlockMap, "b to block map")
  .field_readonly("n", &rlmerPredD::n_)
  .field_readonly("p", &rlmerPredD::p_)
  .field_readonly("q", &rlmerPredD::q_)
  .field_readonly("X", &rlmerPredD::X_)
  .field_readonly("Zt", &rlmerPredD::Zt_)
  ;

  class_<rlmerPredD_DAS>("rlmerPredD_DAS")
  .derives<rlmerPredD>("rlmerPredD")
  .constructor<Rcpp::List, MMap, MSpMatrixd>()
  .method("s_e", &rlmerPredD_DAS::s_e, "s(obj, theta = false)")
  .method("S_b", &rlmerPredD_DAS::s_b, "S(obj)")
  .method("kappa_e", &rlmerPredD_DAS::getKappa_e, "kappa_e accessor")
  .method("kappa_b", &rlmerPredD_DAS::getKappa_b_copy, "kappa_b accessor")
  .method("Tb", &rlmerPredD_DAS::getTau_b_copy, "Tau_b accessor")
  .method("tau_e", &rlmerPredD_DAS::getTau_e_copy, "tau_e accessor")
  .method("A", &rlmerPredD_DAS::getA_copy, "A accessor")
  .method("Kt", &rlmerPredD_DAS::getKt_copy, "Kt accessor")
  .method("L", &rlmerPredD_DAS::getL_copy, "L accessor")
  ;

  class_<rlmerPredD_DAStau>("rlmerPredD_DAStau")
  .derives<rlmerPredD_DAS>("rlmerPredD_DAS")
  .constructor<Rcpp::List, MMap, MSpMatrixd>()
  ;

  class_<M>("M")
    .method("XZ", &M::XZ)
    .method("bb", &M::bb)
    .method("bB", &M::bB)
    .method("Bb", &M::Bb)
    .method("BB", &M::BB)
    .method("bbinv", &M::bbinv)
  ;

  function("calcED_re", &RcalcED_re);
  function("calcEpsi_bbt", &RcalcEpsi_bbt);
  function("calcEpsi_bpsi_bt", &RcalcEpsi_bpsi_bt);
  function("test_DqagNormalExpectation", &test_DqagNormalExpectation);
  function("test_GaussHermiteQuadrature", &test_GaussHermiteQuadrature);
  function("test_GaussHermiteNormalExpectation", &test_GaussHermiteNormalExpectation);
  function("test_DqagIntegration2d_ninfInf", &test_DqagIntegration2d_ninfInf);
  function("test_DqagIntegration2d_aInf", &test_DqagIntegration2d_aInf);
  function("test_DqagIntegration2d_ninfB", &test_DqagIntegration2d_ninfB);
  function("test_DqagIntegration2d_aB", &test_DqagIntegration2d_aB);
  function("test_DqagNormalExpectation2d_ninfInf", &test_DqagNormalExpectation2d_ninfInf);
  function("test_GaussHermiteNormalExpectation2d", &test_GaussHermiteNormalExpectation2d);
  function("test_ScalarTauParameters", &test_ScalarTauParameters);
  function("test_calcTau", &RcalcTau);
  function("test_calcTauVectorized", &RcalcTauVectorized);
  function("test_Hcubature_ninfInf", &test_Hcubature_ninfInf);
  function("test_Hcubature_ninfInf2", &test_Hcubature_ninfInf2);
  function("test_HcubatureNormalExpectation_ninfInf", &test_HcubatureNormalExpectation_ninfInf);
  function("test_Pcubature_ninfInf", &test_Pcubature_ninfInf);
  function("test_Pcubature_ninfInf", &test_Pcubature_ninfInf);
  function("test_Pcubature_aInf", &test_Pcubature_aInf);
  function("test_Pcubature_ninfB", &test_Pcubature_ninfB);
  function("test_Pcubature_aB", &test_Pcubature_aB);
  function("test_Hcubature2d_ninfInf", &test_Hcubature2d_ninfInf);
  function("test_Hcubature2d_aInf", &test_Hcubature2d_aInf);
  function("test_Hcubature2d_ninfB", &test_Hcubature2d_ninfB);
  function("test_Hcubature2d_aB", &test_Hcubature2d_aB);
  function("test_Pcubature2d_ninfInf", &test_Pcubature2d_ninfInf);
  function("test_Pcubature2d_aInf", &test_Pcubature2d_aInf);
  function("test_Pcubature2d_ninfB", &test_Pcubature2d_ninfB);
  function("test_Pcubature2d_aB", &test_Pcubature2d_aB);
  function("test_GaussianQuadratureNd_ninfInf", &test_GaussianQuadratureNd_ninfInf);
  function("test_GaussianQuadratureNdNormalExpectation_ninfInf",
           &test_GaussianQuadratureNdNormalExpectation_ninfInf);
  function("testMatrixTauParameters", &testMatrixTauParameters);

  function("wgt_e", &Rwgt_e);
  function("wgt_b", &Rwgt_b);

}
