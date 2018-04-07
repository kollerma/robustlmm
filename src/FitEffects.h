#if !defined ROBUSTLMM_FITEFFECTS_H__
#define ROBUSTLMM_FITEFFECTS_H__

#include "globals.h"

typedef Eigen::CholmodDecomposition<SpMatrixd>    ChmDecomp;

class FitEffects : public SimpleIterativeFitter<VectorXd> {
private:
  rlmerPredD* pp_;
  rlmerResp* resp_;

  int n_, p_, q_;
  SpMatrixd mat1_;
  MatrixXd mat2_;
  VectorXd invU_ey_, W_;
  ChmDecomp decomp_;

public:
  FitEffects(rlmerPredD* pp,rlmerResp* resp,
             const double relativeTolerance, const int maxOperations);
  FitEffects(rlmerPredDXPtr pp, rlmerRespXPtr resp,
             const double relativeTolerance, const int maxOperations);
  ~FitEffects();

  const Fit<VectorXd>& fit();

  SpMatrixd getMat1_copy() const;
  MatrixXd getMat2_copy() const;
  const VectorXd getInvU_ey() const;
  const VectorXd getNextValue();
  const VectorXd getW();

private:
  std::string init();
  void initMat1();
  void initMat2();
  void updateW();
  std::string computeAndSetNextValue();

protected:
  std::string doIteration() ;
  void setParameters();
};

extern "C" SEXP _rcpp_module_boot_fitEffects_module();

#endif
