#if !defined  ROBUSTLMM_DAS_SCALE_H__
#define  ROBUSTLMM_DAS_SCALE_H__

#include "globals.h"

inline double d2(const double value) {
  return std::sqrt(value);
}

inline double Dd2(const double value) {
  return 1. / d2(value);
}

double calcED_re(const PsiFuncXPtr& rho_b, unsigned s, Integration* const integration);
MatrixXd calcEpsi_bbt(const PsiFuncXPtr& rho_b, unsigned s, Integration* const integration);
MatrixXd calcEpsi_bpsi_bt(const PsiFuncXPtr& rho_b, unsigned s, Integration* const integration);
double calcKappaTau(const PsiFuncXPtr& rho, unsigned s, Integration* const integration);

double RcalcED_re(const PsiFuncXPtr& rho_b, unsigned s);
MatrixXd RcalcEpsi_bbt(const PsiFuncXPtr& rho_b, unsigned s);
MatrixXd RcalcEpsi_bpsi_bt(const PsiFuncXPtr& rho_b, unsigned s);

class ScalarTauParameters {
  const double a_, s_, kappa_;
  const PsiFuncXPtr& rho_, rhoSigma_;
  const double *tau_;
  const int index_;

public:
  ScalarTauParameters(const double a, const double s,
                      const double kappa, const PsiFuncXPtr& rho,
                      const PsiFuncXPtr& rhoSigma, const double* const tau,
                      const unsigned index);
  ~ScalarTauParameters();

  double numerator(const double* const x, const double* const y) const;
  double denominator(const double*  const x, const double* const y) const;
  const double& getKappa() const;
  const int& getIndex() const;

private:
  double approx(const double* const x, const double* const y) const;
};

VectorXd
  test_ScalarTauParameters(const double a, const double s,
                           const double kappa, const PsiFuncXPtr& rho,
                           const PsiFuncXPtr& rhoSigma, const double tau,
                           const double x, const double y);

void calcTau_numerator(double* const x, const int n, void *const ex);
void calcTau_denominator(double* const x, const int n, void *const ex);

class ScalarTauIterativeFitter : public SimpleIterativeFitter<double> {
private:
  ScalarTauParameters parameters_;
  Expectation2d *expectation_;

public:
  ScalarTauIterativeFitter(const double initialValue, const double relativeTolerance,
                           const unsigned maxOperations, const double a, const double s,
                           const double kappa, const PsiFuncXPtr& rho,
                           const PsiFuncXPtr& rhoSigma, Expectation2d* const expectation,
                           const unsigned index);

protected:
  std::string doIteration();
};


double calcTau(const double initialValue, const double relativeTolerance,
               const unsigned maxOperations,
               const double a, const double s, const double kappa,
               const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
               Expectation2d* const expectation, const unsigned index);

double RcalcTau(const double initialValue, const double relativeTolerance,
                 const unsigned maxOperations,
                 const double a, const double s, const double kappa,
                 const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                 const unsigned nodes);

class CompareIndicesByAnotherVectorValues {
  const VectorXd *values_;
public:
  CompareIndicesByAnotherVectorValues(const VectorXd* const values);

  bool operator() (const unsigned& a, const unsigned& b) const;
};

bool isAlmostEqual(const VectorXd& a, const VectorXd& s,
                   const unsigned lastIndex, const unsigned currentIndex,
                   const double tolerance);

VectorXd calcTau(const VectorXd& initialValue, const double relativeTolerance,
                 const unsigned maxOperations,
                 const VectorXd& a, const VectorXd& s, const double kappa,
                 const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                 Expectation2d* const expectation);

VectorXd RcalcTauVectorized(const VectorXd initialValue, const double relativeTolerance,
                            const unsigned maxOperations,
                            const VectorXd a, const VectorXd s, const double kappa,
                            const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                            const unsigned nodes);

class MatrixTauParameters {
  // FIXME limit matrices to fixed size?
  const unsigned size_;
  const double skappa_;
  const MatrixXd& Lkk_, Sk_;
  const PsiFuncXPtr& rho_, rhoSigma_;
  const MatrixXd& Tbk_;
  const LLT<MatrixXd>& cholTbk_;
  bool isNotOfFullRank_;

public:
  MatrixTauParameters(const double skappa, const MatrixXd& Lkk, const MatrixXd& Sk,
                      const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                      const MatrixXd& Tbk);
  MatrixTauParameters(const double skappa, const MatrixXd& Lkk, const MatrixXd& Sk,
                      const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                      const MatrixXd& Tbk,
                      const LLT<MatrixXd>& cholTbk);
  ~MatrixTauParameters();

  double funA(const double* const u);
  MatrixXd funB(const double* const u);
  unsigned ndim() const;
  unsigned fdimB() const;

public:
  VectorXd approx(const double* const u) const;
  double crossprodSolve(const VectorXd& btilde);
  double wgtDelta(const double x) const;
  bool isNotOfFullRank() const;
};

Rcpp::List testMatrixTauParameters(const Rcpp::List & args);

// typedef int (*integrand) (unsigned ndim, const double *x, void *ex,
//             unsigned fdim, double *fval);

int calcTauNonDiagNumerator(unsigned ndim, const double *u, void *ex,
                            unsigned fdim, double *fval);

int calcTauNonDiagDenominator(unsigned ndim, const double *u, void *ex,
                              unsigned fdim, double *fval);

class MatrixTauIterativeFitter : public SimpleIterativeFitter<MatrixXd> {
private:
  MatrixXd value_;
  LLT<MatrixXd> cholTbk_;
  MatrixTauParameters parameters_;
  IntegrandNd numerator_, denominator_;
  ExpectationNd* const expectation_;
  double lastDenominator_;

public:
  MatrixTauIterativeFitter(const MatrixXd& initialValue, const double relativeTolerance,
                           const unsigned maxOperations, const double skappa_,
                           const MatrixXd& Lkk, const MatrixXd& Sk,
                           const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                           ExpectationNd* const expectation);

protected:
  std::string doIteration();
};

int calcTauNonDiag(const BlockTypeIndex* const blockType, SpMatrixd& value,
                   const double relativeTolerance, const unsigned maxOperations,
                   const VectorXd& kappa_, const MatrixXd& L,
                   const std::vector<MatrixXd>& skbs,
                   const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                   ExpectationNd* const expectation);

#endif
