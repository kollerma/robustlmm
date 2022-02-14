#include "DasScale.h"

using namespace Rcpp;

#define DEBUG(STRING)                                                  \
//    Rcpp::Rcout << STRING << std::endl;

#define SETUP_INTEGRATION(rho)                                         \
const void *exc[2] = {&rho, &s};                                       \
void **ex = const_cast<void**>(exc);

#define SETUP_INTEGRATION_PASSIN(rho)                                  \
int limit = 100;                                                       \
double epsabs = std::pow(std::numeric_limits<double>::epsilon(), .5),  \
  bound = 0., arg = 0.;                                                \
const void *exc[5] = {&rho, &s, &arg, integration, &bound};           \
void **ex = const_cast<void**>(exc);

#define PROCESS_EX_0()                                                 \
const PsiFuncXPtr **ps = static_cast<const PsiFuncXPtr **const>(ex);   \
const PsiFuncXPtr *p = ps[0];

#define PROCESS_EX()                                                   \
PROCESS_EX_0()                                                         \
const unsigned **args = static_cast<const unsigned **const>(ex);       \
const unsigned *df = args[1];

#define PROCESS_EX_PASSIN()                                            \
Integration **iargs = static_cast<Integration **>(ex);                 \
Integration *integration = iargs[3];                                   \
double **dargs = static_cast<double **>(ex);                           \
double *bound = dargs[4];                                      \

void calcED_reIntegrand1(double *x, const int n, void *const ex);
void calcED_reIntegrand2(double *x, const int n, void *const ex);
void calcEpsi_bbtIntegrand(double *x, const int n, void *const ex);
void calcEpsi_bpsi_btIntegrand(double *x, const int n, void *const ex);
void calcKappaTau_s1_numerator_integrand(double *x, const int n, void *const ex);
void calcKappaTau_s1_denominator_integrand(double *x, const int n, void *const ex);
void calcKappaTau_sgt1_integrand(double *x, const int n, void *const ex);
double calcKappaTau_sgt1_integral(double kappa, void *ex);

void calcED_reIntegrand1(double *x, const int n, void *const ex) {
  PROCESS_EX();
  for (int i = 0; i < n; i++) {
    double value = x[i];
    x[i] = (*p)->DwgtFun(d2(value)) * Dd2(value) * value * R::dchisq(value, *df, 0);
  }
  return;
}

void calcED_reIntegrand2(double *x, const int n, void *const ex) {
  PROCESS_EX();
  for (int i = 0; i < n; i++) {
    double value = x[i];
    x[i] = (*p)->wgtFun(d2(value)) * R::dchisq(value, *df, 0);
  }
  return;
}

double calcED_re(const PsiFuncXPtr& rho_b, unsigned s, Integration* const integration) {
  if (s == 1)
    return rho_b->EDpsi();

  SETUP_INTEGRATION(rho_b);
  double bound = 0.;
  double expectation = integration->aInf(calcED_reIntegrand1, ex, &bound);
  expectation /= s;
  expectation += integration->aInf(calcED_reIntegrand2, ex, &bound);
  return expectation;
}

void calcEpsi_bbtIntegrand(double *x, const int n, void *const ex) {
  PROCESS_EX();
  for (int i = 0; i < n; i++) {
    double value = x[i];
    x[i] = (*p)->wgtFun(d2(value)) * value * R::dchisq(value, *df, 0);
  }
  return;
}

MatrixXd calcEpsi_bbt(const PsiFuncXPtr& rho_b, unsigned s, Integration* const integration) {
  double expectation;
  if (s == 1) {
    expectation = rho_b->EDpsi();
  } else {
    SETUP_INTEGRATION(rho_b);
    double bound = 0.;                                                     \
    expectation = integration->aInf(calcEpsi_bbtIntegrand, ex, &bound);
    expectation /= s;
  }

  VectorXd res(s);
  res.fill(expectation);
  return res.asDiagonal();
}

void calcEpsi_bpsi_btIntegrand(double *x, const int n, void *const ex) {
  PROCESS_EX()
  for (int i = 0; i < n; i++) {
    double value = x[i];
    double wgt = (*p)->wgtFun(d2(value));
    x[i] = wgt * wgt * value * R::dchisq(value, *df, 0);
  }
  return;
}

MatrixXd calcEpsi_bpsi_bt(const PsiFuncXPtr& rho_b, unsigned s, Integration* const integration) {
  double expectation;
  if (s == 1) {
    expectation = rho_b->Epsi2();
  } else {
    SETUP_INTEGRATION(rho_b);
    double bound = 0.;                                                     \
    expectation = integration->aInf(calcEpsi_bpsi_btIntegrand, ex, &bound);
    expectation /= s;
  }

  VectorXd res;
  res.setConstant(s, expectation);
  return res.asDiagonal();
}

void calcKappaTau_s1_numerator_integrand(double *x, const int n, void *const ex) {
  PROCESS_EX_0()
  for (int i = 0; i < n; i++) {
    double value = x[i];
    x[i] = (*p)->psiFun(value) * value * stats::dnorm_0(value, 0);
  }
  return;
}

void calcKappaTau_s1_denominator_integrand(double *x, const int n, void *const ex) {
  PROCESS_EX_0()
  for (int i = 0; i < n; i++) {
    double value = x[i];
    x[i] = (*p)->wgtFun(value) * stats::dnorm_0(value, 0);
  }
  return;
}

void calcKappaTau_sgt1_integrand(double *x, const int n, void *const ex) {
  PROCESS_EX()
  const double **dargs = static_cast<const double **const>(ex);
  const double *kappa = dargs[2];
  for (int i = 0; i < n; i++) {
    double value = x[i];
    x[i] = (*p)->psiFun(value - (*df) * (*kappa)) * R::dchisq(value, *df, 0);
  }
  return;
}

double calcKappaTau_sgt1_integral(double kappa, void *ex) {
  PROCESS_EX_PASSIN()
  *dargs[2] = kappa;
  return integration->aInf(calcKappaTau_sgt1_integrand, ex, bound);
}

double calcKappaTau(const PsiFuncXPtr& rho, unsigned s, Integration* const integration) {
  double out;
  if (s > 38)
    return 1.;
  if (s == 1) {
    SETUP_INTEGRATION(rho)
    out = integration->ninfInf(calcKappaTau_s1_numerator_integrand, ex);
    out /= integration->ninfInf(calcKappaTau_s1_denominator_integrand, ex);
  } else {
    SETUP_INTEGRATION_PASSIN(rho)
    out = Rrlmm_zeroin(0., 1., calcKappaTau_sgt1_integral, ex, &epsabs, &limit);
    if (limit == -1)
      warn("calcKappaTau_sgt1_integral: zeroin failed to find a root (difference %f)", epsabs);
  }
  return out;
}

double RcalcED_re(const PsiFuncXPtr& rho_b, unsigned s) {
  DqagIntegration integration;
  return calcED_re(rho_b, s, &integration);
}

MatrixXd RcalcEpsi_bbt(const PsiFuncXPtr& rho_b, unsigned s) {
  DqagIntegration integration;
  return calcEpsi_bbt(rho_b, s, &integration);
}

MatrixXd RcalcEpsi_bpsi_bt(const PsiFuncXPtr& rho_b, unsigned s) {
  DqagIntegration integration;
  return calcEpsi_bpsi_bt(rho_b, s, &integration);
}

// class ScalarTauParameters

ScalarTauParameters::ScalarTauParameters(const double a, const double s,
                                         const double kappa, const PsiFuncXPtr& rho,
                                         const PsiFuncXPtr& rhoSigma,
                                         const double* const tau, const unsigned index) :
  a_(a), s_(s), kappa_(kappa), rho_(rho), rhoSigma_(rhoSigma),
  tau_(tau), index_(index) {}

ScalarTauParameters::~ScalarTauParameters() {}

double ScalarTauParameters::approx(const double* const x, const double* const y) const {
  // (x-a*psi(x)-s*y)/tau
  return ((*x) - a_ * rho_->psiFun(*x) - s_ * (*y)) / (*tau_);
}

double ScalarTauParameters::numerator(const double* const x, const double* const y) const {
  // rho.sigma.e@psi(t)*t*(tau*tau)
  double t = approx(x, y);
  return rhoSigma_->psiFun(t) * t * ((*tau_) * (*tau_));
}

double ScalarTauParameters::denominator(const double* const x, const double* const y) const {
  // kappa*rho.sigma.e@wgt(t)
  return rhoSigma_->wgtFun(approx(x, y));
}

const double& ScalarTauParameters::getKappa() const {
  return kappa_;
}

const int& ScalarTauParameters::getIndex() const {
  return index_;
}

VectorXd
  test_ScalarTauParameters(const double a, const double s,
                           const double kappa, const PsiFuncXPtr& rho,
                           const PsiFuncXPtr& rhoSigma, const double tau,
                           const double x, const double y) {
    ScalarTauParameters parameters(a, s, kappa, rho, rhoSigma, &tau, 1);
    VectorXd actual(2);
    actual[0] = parameters.numerator(&x, &y);
    actual[1] = parameters.denominator(&x, &y);
    return actual;
}

// end class ScalarTauParameters

void calcTau_numerator(double *x, const int n, void *const ex) {
  ScalarTauParameters** parameters = static_cast<ScalarTauParameters**>(ex);
  double **y = static_cast<double**>(ex);

  for (int i = 0; i < n; ++i) {
     x[i] = parameters[0]->numerator(&(x[i]), y[1]);
  }
}

void calcTau_denominator(double *x, const int n, void *const ex) {
  ScalarTauParameters** parameters = static_cast<ScalarTauParameters**>(ex);
  double **y = static_cast<double**>(ex);

  for (int i = 0; i < n; ++i) {
    x[i] = parameters[0]->denominator(&(x[i]), y[1]);
  }
}

// class ScalarTauIterativeFitter

ScalarTauIterativeFitter::ScalarTauIterativeFitter(
  const double initialValue, const double relativeTolerance,
  const unsigned maxOperations, const double a, const double s,
  const double kappa, const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
  Expectation2d* const expectation, const unsigned index) :
  SimpleIterativeFitter<double>(initialValue, relativeTolerance, maxOperations),
  parameters_(ScalarTauParameters(a, s, kappa, rho, rhoSigma, fit_.getValuePtr(), index)),
  expectation_(expectation) {}

std::string ScalarTauIterativeFitter::doIteration() {
  double y;
  const void *innerExc[2] = {&parameters_, &y};
  void *innerEx = const_cast<void**>(innerExc);

  double numerator = expectation_->ninfInf(calcTau_numerator, innerEx);
  if (expectation_->getErrorCode()) {
    return "Error during numerical integration.";
  }
  double denominator = expectation_->ninfInf(calcTau_denominator, innerEx) *
    parameters_.getKappa();
  if (expectation_->getErrorCode()) {
    return "Error during numerical integration.";
  }
  double newValue = std::sqrt(numerator / denominator);
  if (newValue > 1.) {
    if (fit_.getValue() != 1.)
      warn("tau[%i] > 1, setting it to 1", parameters_.getIndex() + 1);
    newValue = 1.;
  }
  this->update(newValue);
  return "";
}

// end class ScalarTauIterativeFitter

double calcTau(const double initialValue, const double relativeTolerance,
               const unsigned maxOperations,
               const double a, const double s, const double kappa,
               const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
               Expectation2d* const expectation, const unsigned index) {
  if (initialValue == 0.) {
    return 0.;
  }

  ScalarTauIterativeFitter fitter(initialValue, relativeTolerance, maxOperations,
                                  a, s, kappa, rho, rhoSigma, expectation, index);
  Fit<double> fit(fitter.fit());
  if (fit.getConvergenceStatus()) {
    warn("Convergence status: %i, %s", fit.getConvergenceStatus(), fit.getMessage().c_str());
  }
  return fit.getValue();
}

double RcalcTau(const double initialValue, const double relativeTolerance,
                const unsigned maxOperations,
                const double a, const double s, const double kappa,
                const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                const unsigned nodes) {
  if (nodes == 0) {
    DqagNormalExpectation2d expectation;
    return calcTau(initialValue, relativeTolerance, maxOperations, a, s, kappa,
                   rho, rhoSigma, &expectation, 1);
  } else {
    GaussianQuadratureNormalExpectation2d expectation(nodes);
    return calcTau(initialValue, relativeTolerance, maxOperations, a, s, kappa,
                   rho, rhoSigma, &expectation, 1);
  }
}

// class CompareIndicesByAnotherVectorValues

CompareIndicesByAnotherVectorValues::CompareIndicesByAnotherVectorValues(const VectorXd* const values) :
  values_(values) {}

bool CompareIndicesByAnotherVectorValues::operator() (const unsigned& a, const unsigned& b) const {
    return (*values_)[a] > (*values_)[b];
}

// end class CompareIndicesByAnotherVectorValues

bool isAlmostEqual(const VectorXd& a, const VectorXd& s,
                   const unsigned lastIndex, const unsigned currentIndex,
                   const double relativeTolerance) {
  if (std::abs(a[lastIndex] - a[currentIndex]) >
        relativeTolerance * std::abs(a[lastIndex]))
    return false;
  if (std::abs(s[lastIndex] - s[currentIndex]) >
        relativeTolerance * std::abs(s[lastIndex]))
    return false;
  return true;
}

VectorXd calcTau(const VectorXd& initialValue, const double relativeTolerance,
                 const unsigned maxOperations,
                 const VectorXd& a, const VectorXd& s, const double kappa,
                 const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                 Expectation2d* const expectation) {
  unsigned size = initialValue.size();
  R_ASSERT(size == a.size());
  R_ASSERT(a.size() == s.size());
  VectorXd result(size);
  result.setZero();

  std::vector<unsigned> idx(result.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in a
  std::sort(idx.begin(), idx.end(), CompareIndicesByAnotherVectorValues(&a));

  // walk and only recompute if different
  for (unsigned i = 0; i < size; ++i) {
    unsigned currentIdx = idx[i];
    if (i > 0 && isAlmostEqual(a, s, idx[i-1], currentIdx, relativeTolerance)) {
      result[currentIdx] = result[idx[i-1]];
    } else {
      result[currentIdx] =
        calcTau(initialValue[currentIdx], relativeTolerance, maxOperations,
                a[currentIdx], s[currentIdx], kappa, rho, rhoSigma, expectation,
                currentIdx);
    }
  }

  return result;
}

VectorXd RcalcTauVectorized(const VectorXd initialValue, const double relativeTolerance,
                            const unsigned maxOperations,
                            const VectorXd a, const VectorXd s, const double kappa,
                            const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                            const unsigned nodes) {
  if (nodes == 0) {
    DqagNormalExpectation2d expectation;
    return calcTau(initialValue, relativeTolerance, maxOperations, a, s, kappa,
                   rho, rhoSigma, &expectation);
  } else {
    GaussianQuadratureNormalExpectation2d expectation(nodes);
    return calcTau(initialValue, relativeTolerance, maxOperations, a, s, kappa,
                   rho, rhoSigma, &expectation);
  }
}

// class MatrixTauParameters

Rcpp::List testMatrixTauParameters(const Rcpp::List &args) {
  double skappa(Rcpp::as<double>(args["skappa"]));
  MatrixXd Lkk(Rcpp::as<MatrixXd>(args["Lkk"]));
  MatrixXd Sk(Rcpp::as<MatrixXd>(args["Sk"]));
  PsiFuncXPtr rho(Rcpp::as<PsiFuncXPtr>(args["rho"]));
  PsiFuncXPtr rhoSigma(Rcpp::as<PsiFuncXPtr>(args["rhoSigma"]));
  MatrixXd Tbk(Rcpp::as<MatrixXd>(args["Tbk"]));
  LLT<MatrixXd> cholTbk(Tbk.rows());
  cholTbk.compute(Tbk);
  std::vector<double> u(Rcpp::as<std::vector<double> >(args["u"]));

  MatrixTauParameters params(skappa, Lkk, Sk, rho, rhoSigma, Tbk, cholTbk);

  VectorXd btilde = params.approx(u.data());
  return Rcpp::List::create(Rcpp::Named("ndim") = params.ndim(),
                            Rcpp::Named("fdimB") = params.fdimB(),
                            Rcpp::Named("btilde") = btilde,
                            Rcpp::Named("crossprodSolve") = params.crossprodSolve(btilde),
                            Rcpp::Named("wgtDelta") = params.wgtDelta(u[0]),
                            Rcpp::Named("funA") = params.funA(u.data()),
                            Rcpp::Named("funB") = params.funB(u.data())
                            );
}

MatrixTauParameters::MatrixTauParameters(const double skappa,
                                         const MatrixXd& Lkk, const MatrixXd& Sk,
                                         const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                                         const MatrixXd& Tbk,
                                         const LLT<MatrixXd>& cholTbk) :
  size_(Lkk.cols()), skappa_(skappa), Lkk_(Lkk), Sk_(Sk),
  rho_(rho), rhoSigma_(rhoSigma), Tbk_(Tbk), cholTbk_(cholTbk),
  isNotOfFullRank_(false) {}

MatrixTauParameters::~MatrixTauParameters() {}

double MatrixTauParameters::funA(const double* const u) {
  return wgtDelta(crossprodSolve(approx(u)));
}

MatrixXd MatrixTauParameters::funB(const double* const u) {
  VectorXd btilde(approx(u));
  MatrixXd value(tcrossprod(btilde));
  MatrixXd result(value * rhoSigma_->wgtFun(crossprodSolve(btilde)));
  // DEBUG("  funB. btilde = " << btilde << " value = " << value << " result = " << result)
  return result;
}

unsigned MatrixTauParameters::ndim() const {
  return size_ * 2;
}

unsigned MatrixTauParameters::fdimB() const {
  return size_ * (size_ + 1) / 2;
}

VectorXd MatrixTauParameters::approx(const double* const u) const {
  VectorXd u1(size_), u2(size_);
  for (unsigned i = 0; i < size_; ++i) {
    u1(i) = u[i];
    u2(i) = u[i+size_];
  }
  // btilde <- u[1:2] - wgt(.d(u[1:2],2)) * Lkk %*% u[1:2] - crossprod(Sk, u[3:4])
  VectorXd approx(u1 - rho_->wgtFun(u1.squaredNorm()) * Lkk_ * u1 - crossprod(Sk_, u2));
  return approx;
}

double MatrixTauParameters::crossprodSolve(const VectorXd& btilde) {
  VectorXd x(cholTbk_.matrixU().solve(btilde));
  if (!btilde.isApprox(cholTbk_.matrixU()*x)) {
     DEBUG("  crossprodSolve of " << btilde.transpose() << " got x = " << x.transpose())
  }
  isNotOfFullRank_ = !btilde.isApprox(cholTbk_.matrixU()*x);
  return x.squaredNorm();
}

double MatrixTauParameters::wgtDelta(const double x) const {
  // wgtDelta <- function(u) (psi.sigma(u) - psi.sigma(u-skappa))/s
  return (rhoSigma_->psiFun(x) - rhoSigma_->psiFun(x - skappa_)) / size_;
}

bool MatrixTauParameters::isNotOfFullRank() const {
  return isNotOfFullRank_;
}

// end class MatrixTauParameters

int calcTauNonDiagNumerator(unsigned ndim, const double* const u, void *ex,
                            unsigned fdim, double *fval) {
  if (ndim * (ndim / 2 + 1) / 4 != fdim) return 1;
  MatrixTauParameters *parameters = static_cast<MatrixTauParameters*>(ex);
  MatrixXd value(parameters->funB(u));
  // DEBUG("calcTauNonDiagNumerator for " << u[0] << ", " << u[1] << ", " << u[2] << ", " << u[3] << " got " << value)
  int k = 0;
  for (unsigned i = 0, size = value.rows(); i < size; ++i) {
    for (unsigned j = i; j < size; ++j) {
      fval[k++] = value(i, j);
    }
  }
  return parameters->isNotOfFullRank();
}

int calcTauNonDiagDenominator(unsigned ndim, const double* const u, void *ex,
                              unsigned fdim, double* const fval) {
  if (fdim != 1) return 1;
  MatrixTauParameters *parameters = static_cast<MatrixTauParameters*>(ex);
  *fval = parameters->funA(u);
  // DEBUG("calcTauNonDiagDenominator for " << u[0] << ", " << u[1] << ", " << u[2] << ", " << u[3] << " got " << *fval)
  return parameters->isNotOfFullRank();
}

// class MatrixTauIterativeFitter

MatrixTauIterativeFitter::MatrixTauIterativeFitter(const MatrixXd& initialValue, const double relativeTolerance,
                                                   const unsigned maxOperations, const double skappa,
                                                   const MatrixXd& Lkk, const MatrixXd& Sk,
                                                   const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                                                   ExpectationNd* const expectation) :
  SimpleIterativeFitter<MatrixXd>(initialValue, relativeTolerance, maxOperations),
  value_(initialValue),
  cholTbk_(initialValue.rows()),
  parameters_(MatrixTauParameters(skappa, Lkk, Sk, rho, rhoSigma, value_, cholTbk_)),
  numerator_(calcTauNonDiagNumerator, parameters_.ndim(), parameters_.fdimB(),
             static_cast<void*>(&parameters_)),
  denominator_(calcTauNonDiagDenominator, parameters_.ndim(), 1,
               static_cast<void*>(&parameters_)),
  expectation_(expectation) {}

std::string MatrixTauIterativeFitter::doIteration() {
  // DEBUG(" prev iteration value_ = " << MVec(value_.data(), value_.cols()*value_.rows()).transpose())
  cholTbk_.compute(value_);
  // DEBUG(" cholTbk_ = " << cholTbk_.matrixU())
  std::vector<double> denominator = expectation_->ninfInf(&denominator_);
  if (expectation_->getErrorCode()) {
    DEBUG("Error during numerical integration (denominator, " << expectation_->getErrorCode() << ")");
    return tfm::format("Error during numerical integration (denominator, %i).",
                       expectation_->getErrorCode());
  }
  std::vector<double> numerator = expectation_->ninfInf(&numerator_);
  if (expectation_->getErrorCode()) {
    DEBUG("Error during numerical integration (numerator_, " << expectation_->getErrorCode() << ")");
    return tfm::format("Error during numerical integration (numerator, %i).",
                       expectation_->getErrorCode());
  }
  unsigned k = 0;
  for (unsigned i = 0, size = value_.rows(); i < size; ++i) {
    for (unsigned j = i; j < size; ++j) {
      if (i != j) value_(j, i) = numerator[k];
      value_(i, j) = numerator[k++];
    }
  }
  DEBUG(" will update to " << MVec(value_.data(), value_.cols()*value_.rows()).transpose() << " / " << denominator[0])
  if (std::abs(denominator[0]) >= 1e-7) {
    lastDenominator_ = denominator[0];
  }
  value_ /= lastDenominator_;
  DEBUG(" updated to " << MVec(value_.data(), value_.cols()*value_.rows()).transpose() << " hasNaN() = " << value_.hasNaN())
  this->update(value_);
  return "";
}

// end class MatrixTauIterativeFitter

int calcTauNonDiag(const BlockTypeIndex* const blockType, SpMatrixd& value,
                    const double relativeTolerance, const unsigned maxOperations,
                    const VectorXd& kappa_, const MatrixXd& L,
                    const std::vector<MatrixXd>& skbs,
                    const PsiFuncXPtr& rho, const PsiFuncXPtr& rhoSigma,
                    ExpectationNd* const expectation) {
  unsigned bidx0;
  double skappa = blockType->dim_ * kappa_[blockType->getIndex()];
  DEBUG("calcTauNonDiag: skappa = " << tfm::format("%.15f", skappa))
  MatrixXd lastLkk;
  MatrixXd lastSk;
  Fit<MatrixXd> fit;

  // sort blocks to get similar blocks one after another
  std::vector<unsigned> idx(blockType->getNumberOfBlocks());
  iota(idx.begin(), idx.end(), 0);
  VectorXd Lkk00(blockType->getNumberOfBlocks());
  for (const BlockIndex* block : blockType->getBlocks()) {
    bidx0 = block->getRandomEffects()[0]->getIndex();
    Lkk00[block->getIndex()] = L.coeff(bidx0, bidx0);
  }
  // sort indexes based on comparing values in Lkk[0,0]
  std::sort(idx.begin(), idx.end(), CompareIndicesByAnotherVectorValues(&Lkk00));

  for (const unsigned lidx : idx) {
    const BlockIndex* block(blockType->getBlocks()[lidx]);
    bidx0 = block->getRandomEffects()[0]->getIndex();
    const MatrixXd Lkk(L.block(bidx0, bidx0, blockType->dim_, blockType->dim_));
    const MatrixXd& Sk(skbs[block->getIndex()]);

    if (isDifferent(Lkk, lastLkk, relativeTolerance) ||
        isDifferent(Sk, lastSk, relativeTolerance)) {
      MatrixXd lTbk(value.block(bidx0, bidx0, blockType->dim_, blockType->dim_));
      lastLkk = Lkk;
      lastSk = Sk;
      DEBUG(" Starting from " << MVec(lTbk.data(), lTbk.cols()*lTbk.rows()).transpose())
      MatrixTauIterativeFitter fitter(lTbk, relativeTolerance, maxOperations, skappa,
                                      Lkk, Sk, rho, rhoSigma, expectation);
      fit = fitter.fit();
      if (fit.getConvergenceStatus()) {
        warn("Convergence status: %i, %s", fit.getConvergenceStatus(), fit.getMessage().c_str());
        return 1;
      }
    }
    for (const RandomEffectIndex* randomEffect1 : block->getRandomEffects()) {
      int index1 = randomEffect1->getIndex();
      int indexInBlock1 = index1 - bidx0;
      for (const RandomEffectIndex* randomEffect2 : block->getRandomEffects()) {
        int index2 = randomEffect2->getIndex();
        int indexInBlock2 = index2 - bidx0;
        double newValue = fit.getValue().coeff(indexInBlock1, indexInBlock2);
        value.coeffRef(index1, index2) = newValue;
      }
    }
  }

  return 0;
}

/*
 typedef Eigen::Triplet<double> T;
 std::vector<T> tripletList;
 tripletList.reserve(q_ + std::max((Lambdat_.nonZeros() - q_) * 2, 0));
 for (int k = 0, size = blockBMap_.size(); k < size; ++k) {
 const VectorXi bidx = blockBMap_[k];
 for (int j = 0, blockSize = bidx.size(); j < blockSize; ++j) {
 for (int i = 0; i < blockSize; ++i) {
 tripletList.push_back(T(bidx[j]-1,bidx[i]-1,Tfull.coeff(bidx[j]-1, bidx[i]-1)));
 }
 }
 }
 SpMatrixd Tau_b(q_, q_);
 Tau_b.setFromTriplets(tripletList.begin(), tripletList.end());
 */
