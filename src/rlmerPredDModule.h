#if !defined  ROBUSTLMM_RLMERPREDDMODULE_H__
#define  ROBUSTLMM_RLMERPREDDMODULE_H__

#include "globals.h"
#include <math.h>

using namespace Rcpp;

class rlmerPredD;
class rlmerPredD_DAS;
class M;

class M {
private:
  MatrixXd XZ_, bb_, bbinv_, bB_, BB_;
  LLT<MatrixXd> lltM_bbinv_;

  const LLT<MatrixXd>& lltM_bbinv();
public:
  M();
  M(M* other);

  void init(const rlmerPredD* const object);

  const MatrixXd& XZ();
  const MatrixXd& bb();
  const MatrixXd& bB();
  MatrixXd Bb();
  const MatrixXd& BB();

  const MatrixXd& bbinv();

  ~M();
};

class rlmerPredD {
  friend class M;

public:
  const MMap X_;
  const MSpMatrixd Zt_;
  const unsigned n_, p_, q_, maxOperations_;
  const double relativeTolerance_;

private:
  IndexMapper indexMapper_;

protected:
  MSpMatrixd Lambdat_;
  const MVec lower_, v_e_;
  MVec theta_, beta_, b_s_, sigma_;
  const MiVec Lind_;
  const PsiFuncXPtr rho_e_, rhoSigma_e_;
  const std::vector<PsiFuncXPtr> rho_b_, rhoSigma_b_;
  MatrixXd invU_eX_, M_XX_, M_XZ0_, invM_XXMZZ0_, M_ZZ0_sub_M_ZX0invM_XXMZZ0_;
  DMatrixXd U_e_, invU_e_, V_e_, D_e_, sqrtD_e_, D_b_, Lambda_b_, Lambda_bD_b_, sqrtLambda_b_;
  SpMatrixd M_ZZ0_, invU_eZ_, invU_btZtU_et_, Epsi_bbt_, Epsi_bpsi_bt_;
  LLT<MatrixXd> lltM_XX_;
  VectorXd b_r_, Epsi2_b_;
  Integration &integration_;

private:
  MatrixXd unsc_;
  M M_;
  bool calledInit_, setM_, setUnsc_;

protected:
  M* Mobj();

public:
  rlmerPredD(Rcpp::List args, MMap X, MSpMatrixd Zt);

  void initMatrices();

  void setSigma(double sigma);
  void setU(const VectorXd& u);
  void setB_s(const VectorXd& b_s);
  void setB_r(const VectorXd& b_r);
  void setBeta(const VectorXd& beta);
  virtual void setTheta(const VectorXd& theta);

  const MSpMatrixd& Lambdat() const;
  MatrixXd unsc();

  SpMatrixd getU_b() const;
  SpMatrixd getU_b_copy() const;
  MatrixXd getU_e_copy() const;
  MatrixXd getV_e_copy() const;
  const VectorXd& getB() const;
  const M& getM();
  VectorXi getZeroB() const;
  const MVec& getTheta() const;
  const MVec& getBeta() const;
  const MVec& getB_s() const;
  double getSigma() const;
  const IndexMapper& getIndexMapper() const;
  const VectorXi& getBBlockMap() const;

  MatrixXd getM_XZ0_copy() const;
  SpMatrixd getM_ZZ0_copy() const;
  MatrixXd getM_ZZ0_sub_M_ZX0invM_XXMZZ0_copy() const;
  MatrixXd getLambda_b_copy() const;
  MatrixXd getLambda_bD_b_copy() const;
  SpMatrixd getEpsi_bbt_copy() const;
  SpMatrixd getEpsi_bpsi_bt_copy() const;
  const VectorXd& getEpsi2_b() const;
  const SpMatrixd& getInvU_btZtU_et() const;
  SpMatrixd getInvU_btZtU_et_copy() const;
  const MatrixXd& getInvU_eX() const;
  MatrixXd getInvU_eX_copy() const;
  const DMatrixXd& getSqrtLambda_b() const;
  const DMatrixXd& getLambda_b() const;
  const DMatrixXd& getInvU_e() const;

  VectorXd getMu() const;
  VectorXd getDist_b() const;
  VectorXd getD_k() const;
  VectorXd getEffects() const;

  const PsiFuncXPtr& getRho_e() const;
  const PsiFuncXPtr& getRhoSigma_e() const;
  const std::vector<PsiFuncXPtr>& getRho_b() const;
  const std::vector<PsiFuncXPtr>& getRhoSigma_b() const;

  virtual ~rlmerPredD();

  bool allRandomEffectsDropped() const;
  bool isBlockDropped(const BlockIndex* const block) const;

protected:
  void setZeroB();
  virtual void updateMatrices();
  virtual void resetCaches();
  virtual void initRho();

protected:
  VectorXd stdB(ConstVectorXd& vector);

  // Helper functions to apply a function to each block type
  template<class T, class O>
  std::vector<O> btapply(const std::vector<T>& input, O (*fun)(const T&));
  template<class T, class O>
  std::vector<O> btapply(const std::vector<T>& input, O (*fun)(const T&, unsigned s));
  template<class T, class O>
  std::vector<O> btapply(const std::vector<T>& input,
                         O (*fun)(const T&, unsigned s, Integration* const integration));
  template<class T, class Fun>
  VectorXd btapplyV(const std::vector<T>& input, Fun fun, ConstVectorXi& rep);
  template<class T, class Fun>
  SpMatrixd btapplyM(const std::vector<T>& input, Fun fun, ConstVectorXi& rep);

  // template<class T>
  // VectorXd btapply(const T& input, double *fun(const T, int s), ConstVectorXi& rep);
};

class rlmerPredD_DAS : public rlmerPredD {
  // this class corresponds to DASvar
  // overload calcTau_e and calcTau_b for other methods
protected:
  MatrixXd A_, Kt_, L_;
  double kappa_e_;
  VectorXd kappa_b_;

private:
  bool initTau_e_, setTau_e_, initTau_b_, setTau_b_;
  SpMatrixd Tau_b_;
  VectorXd tau_e_;

  VectorXd s(MatrixXd m1, MatrixXd m2);

public:
  rlmerPredD_DAS(Rcpp::List args, MMap X, MSpMatrixd Zt);

  void updateMatrices();

  void setTheta(const VectorXd& theta);

  MatrixXd getB();
  MatrixXd getK_copy();
  MatrixXd& getA();
  MatrixXd getA_copy();
  MatrixXd getKt_copy();
  MatrixXd getL_copy();
  SpMatrixd& getTau_b();
  SpMatrixd getTau_b_copy();
  VectorXd& getTau_e();
  VectorXd getTau_e_copy();
  double getKappa_e();
  VectorXd& getKappa_b();
  VectorXd getKappa_b_copy();

  VectorXd s_e();
  std::vector<MatrixXd> s_b();

protected:
  VectorXd getInitTau_e();
  SpMatrixd getInitTau_b();

private:
  virtual VectorXd calcTau_e();
  virtual SpMatrixd calcTau_b();
  void resetCaches();
  void initRho();

  friend class rlmerPredD_DAStau;
};

class rlmerPredD_DAStau : public rlmerPredD_DAS {
  friend class rlmerPredD_DAStau_gh;
  friend class rlmerPredD_DAStau_dqag;

private:
  Expectation2d* const expectation2d_;
  ExpectationNd* const expectationNd_;

public:
  rlmerPredD_DAStau(Rcpp::List args, MMap X, MSpMatrixd Zt);
  Expectation2d* const getExpectation2d();
  ExpectationNd* const getExpectationNd();
  ~rlmerPredD_DAStau();

private:
  VectorXd calcTau_e();
  SpMatrixd calcTau_b();
  void calcTau_b(const BlockTypeIndex* const blockType, SpMatrixd& value, const std::vector<MatrixXd>& skbs);
  void calcTau_b_diagonal(const BlockTypeIndex* const blockType, SpMatrixd& value, const std::vector<MatrixXd>& skbs);
  void calcTau_b_non_diagonal(const BlockTypeIndex* const blockType, SpMatrixd& value, const std::vector<MatrixXd>& skbs);
};

extern "C" SEXP _rcpp_module_boot_rlmerPredD_module();

#endif
