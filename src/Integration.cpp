#include "Integration.h"
#include <exp_cubature.h>
#include "fastGHQuad.h"

#define DEBUG(STRING)                                          \
// Rcpp::Rcout << STRING << std::endl;

// class IntegrFnEx

IntegrFnEx::IntegrFnEx(integr_fn *f, void *ex) : f_(f), ex_(ex) {}

IntegrFnEx::~IntegrFnEx() {}

// end class IntegrFnEx

// class Integration
Integration::Integration() {}

Integration::~Integration() {}

// end class Integration

// class Expectation

Expectation::Expectation() : Integration() {}

// end class NormalExpectation

// class NormalExpectation

NormalExpectation::NormalExpectation() : Expectation() {}

// end class NormalExpectation

// class DqagIntegration
DqagIntegration::DqagIntegration() : Integration(), neval_(0),
  ier_(0), limit_(100), lenw_(4*limit_), last_(0),
  epsabs_(std::pow(std::numeric_limits<double>::epsilon(), .5)),
  epsrel_(epsabs_), result_(0), noBound_(NA_REAL), abserr_(0),
  iwork_(R_Calloc(limit_, int)), work_(R_Calloc(lenw_, double)) {}

DqagIntegration::~DqagIntegration() {
  R_Free(iwork_); R_Free(work_);
}

double DqagIntegration::ninfInf(integr_fn *f, void *ex) {
  dqagi(this->wrap(f, ex), &noBound_, 2);
  return result_;
}

double DqagIntegration::aInf(integr_fn *f, void *ex, double* bound) {
  dqagi(this->wrap(f, ex), bound, 1);
  return result_;
}

double DqagIntegration::ninfB(integr_fn *f, void *ex, double* bound) {
  dqagi(this->wrap(f, ex), bound, -1);
  return result_;
}

double DqagIntegration::aB(integr_fn *f, void *ex, double* a, double* b) {
  dqags(this->wrap(f, ex), a, b);
  return result_;
}

int DqagIntegration::getNeval() {
  return neval_;
}

int DqagIntegration::getIer() {
  return ier_;
}

int DqagIntegration::getErrorCode() {
  return getIer();
}

int DqagIntegration::getLast() {
  return last_;
}

double DqagIntegration::getAbserr() {
  return abserr_;
}

IntegrFnEx DqagIntegration::wrap(integr_fn *f, void *ex) {
  return IntegrFnEx(f, ex);
}

void DqagIntegration::dqagi(IntegrFnEx fnex, double* bound, int inf) {
  Rdqagi(fnex.f_, fnex.ex_, bound, &inf, &epsabs_, &epsrel_, &result_,
         &abserr_, &neval_, &ier_, &limit_, &lenw_, &last_,
         iwork_, work_);
  checkIer();
}

void DqagIntegration::dqags(IntegrFnEx fnex, double* a, double* b) {
  Rdqags(fnex.f_, fnex.ex_, a, b, &epsabs_, &epsrel_, &result_,
         &abserr_, &neval_, &ier_, &limit_, &lenw_, &last_,
         iwork_, work_);
  checkIer();
}

void DqagIntegration::checkIer() {
  if(ier_ > 0 && ier_!=5)
    warn("integration flag %i", ier_);
}

/*
void DqagIntegration::reset() {
  neval_ = 0;
  ier_ = 0;
  last_ = 0;
  result_ = 0.;
  abserr_ = 0.;
}
*/

// end class DqagIntegration

// class DqagNormalExpectation

DqagNormalExpectation::DqagNormalExpectation() :
  DqagIntegration(), NormalExpectation() {}

IntegrFnEx DqagNormalExpectation::wrap(integr_fn *f, void *ex) {
  IntegrFnEx *original = new IntegrFnEx(f, ex);
  void *ex2 = original;
  return IntegrFnEx(dqagNormalExpectationWrapper, ex2);
}

void dqagNormalExpectationWrapper(double *x, const int n, void *const ex) {
  double* xcpy = new double[n];
  memcpy(xcpy, x, n * sizeof(double));

  IntegrFnEx *fnex = static_cast<IntegrFnEx*>(ex);
  fnex->f_(x, n, fnex->ex_);

  for (int i = 0; i < n; i++) {
    x[i] *= Rcpp::stats::dnorm_0(xcpy[i], 0);
  }
  delete[] xcpy;
}

// end class DqagNormalExpectation

// class GaussianQuadrature

GaussianQuadrature::GaussianQuadrature() : x_(20), w_(20) {
  init(gaussianQuadraturePostInit);
}

GaussianQuadrature::GaussianQuadrature(int n) : x_(n), w_(n) {
  init(gaussianQuadraturePostInit);
}

GaussianQuadrature::GaussianQuadrature(int n, void postInitFun(std::vector<double>& x, std::vector<double>& w)) :
  x_(n), w_(n) {
  init(postInitFun);
}

GaussianQuadrature::~GaussianQuadrature() {}

double GaussianQuadrature::ninfInf(integr_fn *f, void *ex) {
  IntegrFnEx fnex(this->wrap(f, ex));
  std::vector<double> x(x_);
  fnex.f_(x.data(), x.size(), fnex.ex_);

  double result = 0.;
  for (int i = 0, n = x.size(); i < n; ++i) {
    result += x[i] * w_[i];
  }
  return result;
}

double GaussianQuadrature::aInf(integr_fn *f, void *ex, double* bound) {
  throw std::logic_error("not implemented yet");
}

double GaussianQuadrature::ninfB(integr_fn *f, void *ex, double* bound) {
  throw std::logic_error("not implemented yet");
}

double GaussianQuadrature::aB(integr_fn *f, void *ex, double* a, double* b) {
  throw std::logic_error("not implemented yet");
}

int GaussianQuadrature::getNeval() {
  return x_.size();
}

int GaussianQuadrature::getErrorCode() {
  return 0;
}

double GaussianQuadrature::getAbserr() {
  throw std::logic_error("not implemented");
}

void GaussianQuadrature::init(void postInitFun(std::vector<double>& x, std::vector<double>& w)) {
  gaussHermiteDataGolubWelsch(x_.size(), &x_, &w_);
  postInitFun(x_, w_);
  DEBUG("x = " << to_string(x_))
  DEBUG("w = " << to_string(w_))
}

void gaussianQuadraturePostInit(std::vector<double>& x, std::vector<double>& w) {
  double xi;
  for (int i = 0, n = x.size(); i < n; ++i) {
    xi = x[i];
    w[i] *= exp(xi*xi);
  }
}

IntegrFnEx GaussianQuadrature::wrap(integr_fn *f, void *ex) {
  return IntegrFnEx(f, ex);
}

// end class GaussianQuadrature

// class GaussianQuadratureNormalExpectation

GaussianQuadratureNormalExpectation::GaussianQuadratureNormalExpectation() :
  GaussianQuadrature(20, gaussianQuadratureNormalExpectationPostInit) {}

GaussianQuadratureNormalExpectation::GaussianQuadratureNormalExpectation(int n) :
  GaussianQuadrature(n, gaussianQuadratureNormalExpectationPostInit) {}

IntegrFnEx GaussianQuadratureNormalExpectation::wrap(integr_fn *f, void *ex) {
  return IntegrFnEx(f, ex);
}

void gaussianQuadratureNormalExpectationPostInit(std::vector<double>& x, std::vector<double>& w) {
  for (int i = 0, n = x.size(); i < n; ++i) {
    x[i] *= 1.41421356237309514547462185873883; // sqrt(2)
    w[i] *= 0.56418958354775627928034964497783; // 1/sqrt(pi)
  }
}

// end class GaussianQuadratureNormalExpectation

void integrandRfun(double *x, const int n, void *const ex) {
  Rcpp::Function **func = static_cast<Rcpp::Function**>(ex);

  Rcpp::NumericVector xin(n);
  std::copy(x, x + n, xin.begin());
  Rcpp::NumericVector xout = (**func)(xin);
  std::copy(xout.begin(), xout.end(), x);

  return;
}

double test_DqagNormalExpectation(Rcpp::Function func) {
  DqagNormalExpectation expectation;
  const void *exc[1] = {&func};
  void *ex = const_cast<void**>(exc);
  double result = expectation.ninfInf(integrandRfun, ex);
  return result;
}

double test_GaussHermiteQuadrature(int nNodes, Rcpp::Function func) {
  GaussianQuadrature quad(nNodes);
  const void *exc[1] = {&func};
  void *ex = const_cast<void**>(exc);
  return quad.ninfInf(integrandRfun, ex);
}

double test_GaussHermiteNormalExpectation(int nNodes, Rcpp::Function func) {
  GaussianQuadratureNormalExpectation quad(nNodes);
  const void *exc[1] = {&func};
  void *ex = const_cast<void**>(exc);
  return quad.ninfInf(integrandRfun, ex);
}

// 2d integration

// class IntegrFnExInner

IntegrFnExInner::IntegrFnExInner(Integration *integration, integr_fn *f, void *ex) :
  integration_(integration), f_(f), ex_(ex), a_(NULL), b_(NULL) {}

IntegrFnExInner::IntegrFnExInner(Integration *integration, integr_fn *f, void *ex, double *a) :
  integration_(integration), f_(f), ex_(ex), a_(a), b_(NULL) {}

IntegrFnExInner::IntegrFnExInner(Integration *integration, integr_fn *f, void *ex, double *a, double *b) :
  integration_(integration), f_(f), ex_(ex), a_(a), b_(b) {}

IntegrFnExInner::~IntegrFnExInner() {}

// end class IntegrFnExInner

void integrand2d_ninfInf(double *x, const int n, void *const ex) {
  IntegrFnExInner **inner = static_cast<IntegrFnExInner**>(ex);
  double **innerEx = static_cast<double**>(inner[0]->ex_);

  for (int i = 0; i < n; ++i) {
    *innerEx[1] = x[i];
    x[i] = inner[0]->integration_->ninfInf(inner[0]->f_, inner[0]->ex_);
  }
}

void integrand2d_aInf(double *x, const int n, void *const ex) {
  IntegrFnExInner **inner = static_cast<IntegrFnExInner**>(ex);
  double **innerEx = static_cast<double**>(inner[0]->ex_);
  double *a = inner[0]->a_;

  for (int i = 0; i < n; ++i) {
    *innerEx[1] = x[i];
    x[i] = inner[0]->integration_->aInf(inner[0]->f_, inner[0]->ex_, a);
  }
}

void integrand2d_ninfB(double *x, const int n, void *const ex) {
  IntegrFnExInner **inner = static_cast<IntegrFnExInner**>(ex);
  double **innerEx = static_cast<double**>(inner[0]->ex_);
  double *b = inner[0]->b_;

  for (int i = 0; i < n; ++i) {
    *innerEx[1] = x[i];
    x[i] = inner[0]->integration_->ninfB(inner[0]->f_, inner[0]->ex_, b);
  }
}

void integrand2d_aB(double *x, const int n, void *const ex) {
  IntegrFnExInner **inner = static_cast<IntegrFnExInner**>(ex);
  double **innerEx = static_cast<double**>(inner[0]->ex_);
  double *a = inner[0]->a_;
  double *b = inner[0]->b_;

  for (int i = 0; i < n; ++i) {
    *innerEx[1] = x[i];
    x[i] = inner[0]->integration_->aB(inner[0]->f_, inner[0]->ex_, a, b);
  }
}

// class Integration2d

Integration2d::Integration2d() {}

// end class Integration2d

// class SimpleIntegration2d

double SimpleIntegration2d::ninfInf(integr_fn *f, void *ex) {
  IntegrFnExInner inner(inner_, f, ex);
  const void *outerExc[1] = {&inner};
  void *outerEx = const_cast<void**>(outerExc);
  return outer_->ninfInf(integrand2d_ninfInf, outerEx);
}

double SimpleIntegration2d::aInf(integr_fn *f, void *ex, double* bound) {
  IntegrFnExInner inner(inner_, f, ex, bound);
  const void *outerExc[1] = {&inner};
  void *outerEx = const_cast<void**>(outerExc);
  return outer_->aInf(integrand2d_aInf, outerEx, bound + 1);
}

double SimpleIntegration2d::ninfB(integr_fn *f, void *ex, double* bound) {
  IntegrFnExInner inner(inner_, f, ex, NULL, bound);
  const void *outerExc[1] = {&inner};
  void *outerEx = const_cast<void**>(outerExc);
  return outer_->ninfB(integrand2d_ninfB, outerEx, bound + 1);
}

double SimpleIntegration2d::aB(integr_fn *f, void *ex, double* a, double* b) {
  IntegrFnExInner inner(inner_, f, ex, a, b);
  const void *outerExc[1] = {&inner};
  void *outerEx = const_cast<void**>(outerExc);
  return outer_->aB(integrand2d_aB, outerEx, a + 1, b + 1);
}

int SimpleIntegration2d::getNeval() {
  return outer_->getNeval();
}

int SimpleIntegration2d::getErrorCode() {
  return outer_->getErrorCode();
}

double SimpleIntegration2d::getAbserr() {
  return outer_->getAbserr();
}

SimpleIntegration2d::SimpleIntegration2d(Integration *outer, Integration *inner) :
  outer_(outer), inner_(inner) {}

SimpleIntegration2d::~SimpleIntegration2d() {}

IntegrFnEx SimpleIntegration2d::wrap(integr_fn *f, void *ex) {
  return IntegrFnEx(f, ex);
}

void SimpleIntegration2d::deleteOuterInner() {
  delete(outer_);
  delete(inner_);
}

// end class SimpleIntegration2d

// class Expectation2d

Expectation2d::Expectation2d() {}

// end class Expectation2d

// class class SimpleExpectation2d

SimpleExpectation2d::SimpleExpectation2d(Expectation *outer, Expectation *inner) :
  SimpleIntegration2d(outer, inner) {}

// end class SimpleExpectation2d

// class NormalExpectation2d

NormalExpectation2d::NormalExpectation2d() {}

// end class NormalExpectation2d

// class SimpleNormalExpectation2d

SimpleNormalExpectation2d::SimpleNormalExpectation2d(NormalExpectation *outer, NormalExpectation *inner) :
  SimpleExpectation2d(outer, inner) {}

// end class SimpleNormalExpectation2d

// class DqagNormalExpectation2d

DqagNormalExpectation2d::DqagNormalExpectation2d() :
  SimpleNormalExpectation2d(new DqagNormalExpectation(), new DqagNormalExpectation()) {}

DqagNormalExpectation2d::~DqagNormalExpectation2d() {
  this->deleteOuterInner();
}

// end class DqagNormalExpectation2d

// class GaussianQuadratureNormalExpectation2d

GaussianQuadratureNormalExpectation2d::GaussianQuadratureNormalExpectation2d() :
  SimpleNormalExpectation2d(new GaussianQuadratureNormalExpectation(),
                      new GaussianQuadratureNormalExpectation()) {}

GaussianQuadratureNormalExpectation2d::GaussianQuadratureNormalExpectation2d(int n) :
  SimpleNormalExpectation2d(new GaussianQuadratureNormalExpectation(n),
                      new GaussianQuadratureNormalExpectation(n)) {}

GaussianQuadratureNormalExpectation2d::~GaussianQuadratureNormalExpectation2d() {
  this->deleteOuterInner();
}

// end class GaussianQuadratureNormalExpectation2d

void integrand2dRfun(double *x, const int n, void *const ex) {
  Rcpp::Function **func = static_cast<Rcpp::Function**>(ex);
  double **y = static_cast<double**>(ex);

  Rcpp::NumericVector xin(n);
  std::copy(x, x + n, xin.begin());
  Rcpp::NumericVector xout = (**func)(xin, Rcpp::wrap(*y[1]));
  std::copy(xout.begin(), xout.end(), x);

  return;
}

double test_DqagIntegration2d_ninfInf(Rcpp::Function func) {
  DqagIntegration outer;
  DqagIntegration innerIntegration;

  SimpleIntegration2d integration(&outer, &innerIntegration);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return integration.ninfInf(integrand2dRfun, innerEx);
}

double test_DqagIntegration2d_aInf(Rcpp::Function func, Rcpp::NumericVector bound) {
  DqagIntegration outer;
  DqagIntegration innerIntegration;

  SimpleIntegration2d integration(&outer, &innerIntegration);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return integration.aInf(integrand2dRfun, innerEx, bound.begin());
}

double test_DqagIntegration2d_ninfB(Rcpp::Function func, Rcpp::NumericVector bound) {
  DqagIntegration outer;
  DqagIntegration innerIntegration;

  SimpleIntegration2d integration(&outer, &innerIntegration);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return integration.ninfB(integrand2dRfun, innerEx, bound.begin());
}

double test_DqagIntegration2d_aB(Rcpp::Function func, Rcpp::NumericVector a, Rcpp::NumericVector b) {
  DqagIntegration outer;
  DqagIntegration innerIntegration;

  SimpleIntegration2d integration(&outer, &innerIntegration);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return integration.aB(integrand2dRfun, innerEx, a.begin(), b.begin());
}

double test_DqagNormalExpectation2d_ninfInf(Rcpp::Function func) {
  DqagNormalExpectation2d  expectation;

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return expectation.ninfInf(integrand2dRfun, innerEx);
}

double test_GaussHermiteNormalExpectation2d(Rcpp::Function func, int nodes) {
  GaussianQuadratureNormalExpectation2d expectation(nodes);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return expectation.ninfInf(integrand2dRfun, innerEx);
}

// Multi dimensional integration

// class IntegrandNd

IntegrandNd::IntegrandNd(integrand f, const int ndim, const int fdim, void *ex) :
  f_(f), ndim_(ndim), fdim_(fdim), ex_(ex) {}

IntegrandNd::~IntegrandNd() {}

// end class Integrand

// wrapper functions for cubature

int infiniteIntegralWrapper(unsigned ndim, const double *t, void *ex,
                            unsigned fdim, double *fval) {
  IntegrandNd *wrappedIntegrand = static_cast<IntegrandNd *>(ex);
  // compute f(t / (1 - t^2)) * (1  + t^2) / ((1 - t^2)^2)
  double tmp;
  double* x = new double[ndim];
  for (unsigned i = 0; i < ndim; ++i) {
    tmp = t[i];
    if (tmp > 1. || tmp < -1.) return 1;
    x[i] = tmp / (1. - tmp*tmp);
  }
  int status = wrappedIntegrand->f_(ndim, x, wrappedIntegrand->ex_, fdim, fval);
  for (unsigned i = 0; i < ndim; ++i) {
    tmp = t[i];
    tmp *= tmp;
    tmp = (1 + tmp) / ((1. - tmp) * (1. - tmp));
    for (unsigned j = 0; j < fdim; ++j) {
      fval[j] *= tmp;
    }
  }
  delete[] x;
  return status;
}

int semiInfiniteIntegralWrapper_aInf(unsigned ndim, const double *t, void *ex,
                                unsigned fdim, double *fval) {
  IntegrandNd **wrappedIntegrand = static_cast<IntegrandNd **>(ex);
  double **bound = static_cast<double**>(ex) + 1;
  // compute f(a + t / (1 - t)) / ((1 - t)^2)
  double tmp;
  double* x = new double[ndim];
  for (unsigned i = 0; i < ndim; ++i) {
    tmp = t[i];
    if (tmp > 1. || tmp < 0.) return 1;
    x[i] = (*bound)[i] + tmp / (1. - tmp);
  }
  int status = (*wrappedIntegrand)->f_(ndim, x, (*wrappedIntegrand)->ex_, fdim, fval);
  for (unsigned i = 0; i < ndim; ++i) {
    tmp = t[i];
    tmp = 1. - tmp;
    tmp *= tmp;
    for (unsigned j = 0; j < fdim; ++j) {
      fval[j] /= tmp;
    }
  }
  delete[] x;
  return status;
}

int semiInfiniteIntegralWrapper_ninfB(unsigned ndim, const double *t, void *ex,
                                     unsigned fdim, double *fval) {
  IntegrandNd **wrappedIntegrand = static_cast<IntegrandNd **>(ex);
  double **bound = static_cast<double**>(ex) + 1;
  // compute f(b - t / (1 - t)) / ((1 - t)^2)
  double tmp;
  double* x = new double[ndim];
  for (unsigned i = 0; i < ndim; ++i) {
    tmp = t[i];
    if (tmp > 1. || tmp < 0.) return 1;
    x[i] = (*bound)[i] - tmp / (1. - tmp);
  }
  int status = (*wrappedIntegrand)->f_(ndim, x, (*wrappedIntegrand)->ex_, fdim, fval);
  for (unsigned i = 0; i < ndim; ++i) {
    tmp = t[i];
    tmp = 1. - tmp;
    tmp *= tmp;
    for (unsigned j = 0; j < fdim; ++j) {
      fval[j] /= tmp;
    }
  }
  delete[] x;
  return status;
}

// class IntegrandNd

IntegrationNd::IntegrationNd() {}

IntegrationNd::~IntegrationNd() {}

// end class IntegrandNd

// class ExpectationNd

// empty

// end class ExpectationNd

// class NormalExpectationNd

NormalExpectationNd::NormalExpectationNd() {}

NormalExpectationNd::~NormalExpectationNd() {}

int normalExpectationNdWrapper(unsigned ndim, const double *x, void *ex,
             unsigned fdim, double *fval) {
  IntegrandNd *originalF = static_cast<IntegrandNd*>(ex);
  int status = originalF->f_(ndim, x, originalF->ex_, fdim, fval);
  for (unsigned i = 0; i < ndim; ++i) {
    double density = Rcpp::stats::dnorm_0(x[i], 0);
    for (unsigned j = 0; j < fdim; ++j) {
      fval[j] *= density;
    }
  }
  return status;
}

IntegrandNd NormalExpectationNd::wrap(IntegrandNd *f) {
  return IntegrandNd(normalExpectationNdWrapper, f->ndim_, f->fdim_, static_cast<void*>(f));
}

// end class NormalExpectationNd

// class CachedIntegrationNd

int CachedIntegrationNd::getNeval() {
  return -1;
}

int CachedIntegrationNd::getErrorCode() {
  return errorCode_;
}

std::vector<double>& CachedIntegrationNd::getAbserr() {
  return err_;
}

CachedIntegrationNd::CachedIntegrationNd() {}

CachedIntegrationNd::~CachedIntegrationNd() {}

// end class CachedIntegrationNd

// class Cubature

std::vector<double> Cubature::ninfInf(IntegrandNd *f) {
  double* xmin = new double[f->ndim_];
  double* xmax = new double[f->ndim_];
  for (int i = 0, size = f->ndim_; i < size; ++i) {
    xmin[i] = -1.;
    xmax[i] = 1.;
  }
  err_.resize(f->fdim_);
  std::vector<double> val(f->fdim_);
  IntegrandNd wrappedF = this->wrap(f);
  errorCode_ = cubature_(f->fdim_, infiniteIntegralWrapper, static_cast<void *>(&wrappedF), f->ndim_,
                         xmin, xmax, maxEval_, reqAbsError_, reqRelError_,
                         ERROR_INDIVIDUAL, val.data(), err_.data());
  delete[] xmax;
  delete[] xmin;
  return val;
}

std::vector<double> Cubature::aInf(IntegrandNd *f, double* bound) {
  double* xmin = new double[f->ndim_];
  double* xmax = new double[f->ndim_];
  for (int i = 0, size = f->ndim_; i < size; ++i) {
    xmin[i] = 0.;
    xmax[i] = 1.;
  }
  err_.resize(f->fdim_);
  std::vector<double> val(f->fdim_);
  IntegrandNd wrappedF = this->wrap(f);
  const void *exc[2] = {&wrappedF, bound};
  void *wrappedEx = const_cast<void**>(exc);
  errorCode_ = cubature_(f->fdim_, semiInfiniteIntegralWrapper_aInf, wrappedEx, f->ndim_,
                         xmin, xmax, maxEval_, reqAbsError_, reqRelError_,
                         ERROR_INDIVIDUAL, val.data(), err_.data());
  delete[] xmax;
  delete[] xmin;
  return val;
}

std::vector<double> Cubature::ninfB(IntegrandNd *f, double* bound) {
  double* xmin = new double[f->ndim_];
  double* xmax = new double[f->ndim_];
  for (int i = 0, size = f->ndim_; i < size; ++i) {
    xmin[i] = 0.;
    xmax[i] = 1.;
  }
  err_.resize(f->fdim_);
  std::vector<double> val(f->fdim_);
  IntegrandNd wrappedF = this->wrap(f);
  const void *exc[2] = {&wrappedF, bound};
  void *wrappedEx = const_cast<void**>(exc);
  errorCode_ = cubature_(f->fdim_, semiInfiniteIntegralWrapper_ninfB, wrappedEx, f->ndim_,
                         xmin, xmax, maxEval_, reqAbsError_, reqRelError_,
                         ERROR_INDIVIDUAL, val.data(), err_.data());
  delete[] xmax;
  delete[] xmin;
  return val;
}

std::vector<double> Cubature::aB(IntegrandNd *f, double* a, double* b) {
  err_.resize(f->fdim_);
  std::vector<double> val(f->fdim_);
  IntegrandNd wrappedF = this->wrap(f);
  errorCode_ = cubature_(f->fdim_, wrappedF.f_, wrappedF.ex_, f->ndim_,
                         a, b, maxEval_, reqAbsError_, reqRelError_,
                         ERROR_INDIVIDUAL, val.data(), err_.data());
  return val;
}

Cubature::Cubature(unsigned maxEval, double reqAbsError, double reqRelError,
                   cubature cubature) :
  maxEval_(maxEval), reqAbsError_(reqAbsError), reqRelError_(reqRelError),
  cubature_(cubature) {}

Cubature::~Cubature() {}

IntegrandNd Cubature::wrap(IntegrandNd *f) {
  return *f;
}

// end class Cubature

// class Hcubature

Hcubature::Hcubature(unsigned maxEval, double reqAbsError, double reqRelError) :
  Cubature(maxEval, reqAbsError, reqRelError, &hcubature) {}

Hcubature::~Hcubature() {}

// end class Hcubature

// class HcubatureNormalExpectation

HcubatureNormalExpectation::HcubatureNormalExpectation(unsigned maxEval, double reqAbsError, double reqRelError) :
  Hcubature(maxEval, reqAbsError, reqRelError) {}

HcubatureNormalExpectation::~HcubatureNormalExpectation() {}

IntegrandNd HcubatureNormalExpectation::wrap(IntegrandNd *f) {
  return NormalExpectationNd::wrap(f);
}

// end class HcubatureNormalExpectation

// class Pcubature

Pcubature::Pcubature(unsigned maxEval, double reqAbsError, double reqRelError) :
  Cubature(maxEval, reqAbsError, reqRelError, &pcubature)  {}

Pcubature::~Pcubature() {}

std::vector<double> Pcubature::ninfInf(IntegrandNd *f) {
  throw std::logic_error("infinite integration bounds not supported for pcubature");
}

std::vector<double> Pcubature::aInf(IntegrandNd *f, double* bound) {
  throw std::logic_error("infinite integration bounds not supported for pcubature");
}

std::vector<double> Pcubature::ninfB(IntegrandNd *f, double* bound) {
  throw std::logic_error("infinite integration bounds not supported for pcubature");
}

// end class Pcubature

// class PcubatureNormalExpectation

PcubatureNormalExpectation::PcubatureNormalExpectation(unsigned maxEval, double reqAbsError, double reqRelError) :
  Pcubature(maxEval, reqAbsError, reqRelError) {}

PcubatureNormalExpectation::~PcubatureNormalExpectation() {}

IntegrandNd PcubatureNormalExpectation::wrap(IntegrandNd *f) {
  return NormalExpectationNd::wrap(f);
}

// end class HcubatureNormalExpectation

// class Integrand2d

int integrand2dNdWrapper(unsigned ndim, const double *t, void *ex,
     unsigned fdim, double *fval) {
  if (ndim != 2 || fdim != 1) return 1;
  IntegrFnEx* integrFnEx = static_cast<IntegrFnEx*>(ex);
  double x = *t;
  double **inner = static_cast<double**>(integrFnEx->ex_);
  *inner[1] = t[1];
  integrFnEx->f_(&x, 1, integrFnEx->ex_);
  fval[0] = x;
  return ISNAN(x);
}

Integrand2d::Integrand2d(IntegrFnEx* integrFnEx) :
  IntegrandNd(integrand2dNdWrapper, 2, 1, static_cast<void *>(integrFnEx)) {}

Integrand2d::~Integrand2d() {}

// class IntegrationNd2d

double IntegrationNd2d::ninfInf(integr_fn *f, void *ex) {
  IntegrFnEx wrappedF = this->wrap(f, ex);
  Integrand2d integrand2d(&wrappedF);
  return integrationNd_->ninfInf(&integrand2d)[0];
}

double IntegrationNd2d::aInf(integr_fn *f, void *ex, double* bound) {
  IntegrFnEx wrappedF = this->wrap(f, ex);
  Integrand2d integrand2d(&wrappedF);
  return integrationNd_->aInf(&integrand2d, bound)[0];
}

double IntegrationNd2d::ninfB(integr_fn *f, void *ex, double* bound) {
  IntegrFnEx wrappedF = this->wrap(f, ex);
  Integrand2d integrand2d(&wrappedF);
  return integrationNd_->ninfB(&integrand2d, bound)[0];
}

double IntegrationNd2d::aB(integr_fn *f, void *ex, double* a, double* b) {
  IntegrFnEx wrappedF = this->wrap(f, ex);
  Integrand2d integrand2d(&wrappedF);
  return integrationNd_->aB(&integrand2d, a, b)[0];
}

int IntegrationNd2d::getNeval() {
  return integrationNd_->getNeval();
}

int IntegrationNd2d::getErrorCode() {
  return integrationNd_->getErrorCode();
}

double IntegrationNd2d::getAbserr() {
  return integrationNd_->getAbserr()[0];
}

IntegrationNd2d::IntegrationNd2d(IntegrationNd* integrationNd) :
  integrationNd_(integrationNd) {}

IntegrationNd2d::~IntegrationNd2d() {
  delete integrationNd_;
}

IntegrFnEx IntegrationNd2d::wrap(integr_fn *f, void *ex) {
  return IntegrFnEx(f, ex);
}

// end class IntegrationNd2d

// class Hcubature2d

Hcubature2d::Hcubature2d(unsigned maxEval, double reqAbsError, double reqRelError) :
  IntegrationNd2d(new Hcubature(maxEval, reqAbsError, reqRelError)) {}

Hcubature2d::~Hcubature2d() {}

// end class Hcubature2d

// class Pcubature2d

Pcubature2d::Pcubature2d(unsigned maxEval, double reqAbsError, double reqRelError) :
  IntegrationNd2d(new Pcubature(maxEval, reqAbsError, reqRelError)) {}

Pcubature2d::~Pcubature2d() {}

// end class Pcubature2d

int integrandNdRfun(unsigned ndim, const double *x, void *ex,
                     unsigned fdim, double *fval) {
  Rcpp::Function *func = static_cast<Rcpp::Function*>(ex);

  Rcpp::NumericVector xin(ndim);
  std::copy(x, x + ndim, xin.begin());
  Rcpp::NumericVector xout = (*func)(xin);
  std::copy(xout.begin(), xout.end(), fval);

  DEBUG("f(" << xin << ") = " << xout)

  int status = 0;
  for(const double x : xout) {
    if (Rcpp::NumericVector::is_na(x)) {
      status = 1;
      break;
    }
  }

  return status;
}

Rcpp::NumericVector test_Hcubature_ninfInf(Rcpp::Function func, int ndim, int fdim) {
  Hcubature hcubature(0, 0, 1e-6);

  IntegrandNd integrandNd(integrandNdRfun, ndim, fdim, &func);

  Rcpp::NumericVector result(fdim);
  std::vector<double> val = hcubature.ninfInf(&integrandNd);
  std::copy(val.begin(), val.end(), result.begin());
  return result;
}

Rcpp::NumericVector test_Hcubature_ninfInf2(Rcpp::Function func, int ndim, int fdim, double tol) {
  Hcubature hcubature(0, 0, tol);

  IntegrandNd integrandNd(integrandNdRfun, ndim, fdim, &func);

  Rcpp::NumericVector result(fdim);
  std::vector<double> val = hcubature.ninfInf(&integrandNd);
  std::copy(val.begin(), val.end(), result.begin());
  return result;
}

Rcpp::NumericVector test_HcubatureNormalExpectation_ninfInf(Rcpp::Function func, int ndim, int fdim) {
  HcubatureNormalExpectation hcubatureNormalExpectation(0, 0, 1e-6);

  IntegrandNd integrandNd(integrandNdRfun, ndim, fdim, &func);

  Rcpp::NumericVector result(fdim);
  std::vector<double> val = hcubatureNormalExpectation.ninfInf(&integrandNd);
  std::copy(val.begin(), val.end(), result.begin());
  return result;
}

Rcpp::NumericVector test_Pcubature_ninfInf(Rcpp::Function func, int ndim, int fdim) {
  Pcubature pcubature(0, 0, 1e-6);

  IntegrandNd integrandNd(integrandNdRfun, ndim, fdim, &func);

  Rcpp::NumericVector result(fdim);
  std::vector<double> val = pcubature.ninfInf(&integrandNd);
  std::copy(val.begin(), val.end(), result.begin());
  return result;
}

Rcpp::NumericVector test_Pcubature_aInf(Rcpp::Function func, int ndim, int fdim,
                                        Rcpp::NumericVector bound)  {
  Pcubature pcubature(0, 0, 1e-6);

  IntegrandNd integrandNd(integrandNdRfun, ndim, fdim, &func);

  Rcpp::NumericVector result(fdim);
  std::vector<double> val = pcubature.aInf(&integrandNd, bound.begin());
  std::copy(val.begin(), val.end(), result.begin());
  return result;
}

Rcpp::NumericVector test_Pcubature_ninfB(Rcpp::Function func, int ndim, int fdim,
                                         Rcpp::NumericVector bound) {
  Pcubature pcubature(0, 0, 1e-6);

  IntegrandNd integrandNd(integrandNdRfun, ndim, fdim, &func);

  Rcpp::NumericVector result(fdim);
  std::vector<double> val = pcubature.ninfB(&integrandNd, bound.begin());
  std::copy(val.begin(), val.end(), result.begin());
  return result;
}

Rcpp::NumericVector test_Pcubature_aB(Rcpp::Function func, int ndim, int fdim,
                                      Rcpp::NumericVector a, Rcpp::NumericVector b) {
  Pcubature pcubature(0, 0, 1e-6);

  IntegrandNd integrandNd(integrandNdRfun, ndim, fdim, &func);

  Rcpp::NumericVector result(fdim);
  std::vector<double> val = pcubature.aB(&integrandNd, a.begin(), b.begin());
  std::copy(val.begin(), val.end(), result.begin());
  return result;
}

double test_Hcubature2d_ninfInf(Rcpp::Function func) {
  Hcubature2d hcubature2d(0, 0, 1e-6);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return hcubature2d.ninfInf(integrand2dRfun, innerEx);
}

double test_Hcubature2d_aInf(Rcpp::Function func, Rcpp::NumericVector bound) {
  Hcubature2d hcubature2d(0, 0, 1e-6);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return hcubature2d.aInf(integrand2dRfun, innerEx, bound.begin());
}

double test_Hcubature2d_ninfB(Rcpp::Function func, Rcpp::NumericVector bound) {
  Hcubature2d hcubature2d(0, 0, 1e-6);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return hcubature2d.ninfB(integrand2dRfun, innerEx, bound.begin());
}

double test_Hcubature2d_aB(Rcpp::Function func, Rcpp::NumericVector a, Rcpp::NumericVector b) {
  Hcubature2d hcubature2d(0, 0, 1e-5);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return hcubature2d.aB(integrand2dRfun, innerEx, a.begin(), b.begin());
}

double test_Pcubature2d_ninfInf(Rcpp::Function func) {
  Pcubature2d pcubature2d(0, 0, 1e-6);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return pcubature2d.ninfInf(integrand2dRfun, innerEx);
}

double test_Pcubature2d_aInf(Rcpp::Function func, Rcpp::NumericVector bound) {
  Pcubature2d pcubature2d(0, 0, 1e-6);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return pcubature2d.aInf(integrand2dRfun, innerEx, bound.begin());
}

double test_Pcubature2d_ninfB(Rcpp::Function func, Rcpp::NumericVector bound) {
  Pcubature2d pcubature2d(0, 0, 1e-6);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return pcubature2d.ninfB(integrand2dRfun, innerEx, bound.begin());
}

double test_Pcubature2d_aB(Rcpp::Function func, Rcpp::NumericVector a, Rcpp::NumericVector b) {
  Pcubature2d pcubature2d(0, 0, 1e-5);

  double y;
  const void *innerExc[2] = {&func, &y};
  void *innerEx = const_cast<void**>(innerExc);

  return pcubature2d.aB(integrand2dRfun, innerEx, a.begin(), b.begin());
}

// class GaussianQuadratureNd

std::vector<double> GaussianQuadratureNd::ninfInf(IntegrandNd *f) {
  std::vector<double> result(f->fdim_);
  double* x = new double[f->ndim_];
  double* fval = new double[f->fdim_];
  errorCode_ = doNinfInf(f, x, fval, 1., result, f->ndim_);
  delete[] x;
  delete[] fval;
  return result;
}

int GaussianQuadratureNd::doNinfInf(IntegrandNd *f, double* x, double *fval, double carriedWgt,
                                     std::vector<double>& result, int dim) {
  int index = dim - 1;
  int ret;
  double updatedWgt;

  if (dim == 1) {
    for(int i = 0, n = x_.size(); i < n; ++i) {
      x[index] = x_[i];
      ret = f->f_(f->ndim_, x, f->ex_, f->fdim_, fval);
      if (ret) {
        return ret;
      }
      updatedWgt = carriedWgt * w_[i];
      for(int j = 0, m = f->fdim_; j < m; ++j) {
        result[j] += fval[j] * updatedWgt;
      }
    }
  } else {
    for(int i = 0, n = x_.size(); i < n; ++i) {
      x[index] = x_[i];
      updatedWgt = carriedWgt * w_[i];
      ret = doNinfInf(f, x, fval, updatedWgt, result, dim - 1);
      if (ret) {
        return ret;
      }
    }
  }

  return 0;
}

std::vector<double> GaussianQuadratureNd::aInf(IntegrandNd *f, double* bound) {
  throw new std::runtime_error("Not implemented yet");
}

std::vector<double> GaussianQuadratureNd::ninfB(IntegrandNd *f, double* bound) {
  throw new std::runtime_error("Not implemented yet");
}

std::vector<double> GaussianQuadratureNd::aB(IntegrandNd *f, double* a, double* b) {
  // TODO
  throw new std::runtime_error("Not implemented yet");
}

GaussianQuadratureNd::GaussianQuadratureNd(int n) : x_(n), w_(n) {
  init(gaussianQuadraturePostInit);
}

GaussianQuadratureNd::~GaussianQuadratureNd() {}

GaussianQuadratureNd::GaussianQuadratureNd(int n, void postInitFun(std::vector<double>& x, std::vector<double>& w)) :
  x_(n), w_(n) {
  init(postInitFun);
}

IntegrandNd GaussianQuadratureNd::wrap(IntegrandNd *f) {
  return *f;
}

void GaussianQuadratureNd::init(void postInitFun(std::vector<double>& x, std::vector<double>& w)) {
  gaussHermiteDataGolubWelsch(x_.size(), &x_, &w_);
  postInitFun(x_, w_);
}

// end class GaussianQuadratureNd

Rcpp::NumericVector test_GaussianQuadratureNd_ninfInf(Rcpp::Function func, int ndim, int fdim, int nnodes) {
  GaussianQuadratureNd gq(nnodes);
  IntegrandNd integrandNd(integrandNdRfun, ndim, fdim, &func);
  Rcpp::NumericVector result(fdim);
  std::vector<double> val = gq.ninfInf(&integrandNd);
  std::copy(val.begin(), val.end(), result.begin());

  return result;
}

// class GaussianQuadratureNdNormalExpectation

GaussianQuadratureNdNormalExpectation::GaussianQuadratureNdNormalExpectation(int n) :
  GaussianQuadratureNd(n, gaussianQuadratureNormalExpectationPostInit) {}

GaussianQuadratureNdNormalExpectation::~GaussianQuadratureNdNormalExpectation() {}

IntegrandNd GaussianQuadratureNdNormalExpectation::wrap(IntegrandNd *f) {
  return *f;
}

// end class GaussianQuadratureNdNormalExpectation

Rcpp::NumericVector test_GaussianQuadratureNdNormalExpectation_ninfInf(Rcpp::Function func, int ndim, int fdim, int nnodes) {
  GaussianQuadratureNdNormalExpectation gq(nnodes);
  IntegrandNd integrandNd(integrandNdRfun, ndim, fdim, &func);
  Rcpp::NumericVector result(fdim);
  std::vector<double> val = gq.ninfInf(&integrandNd);
  std::copy(val.begin(), val.end(), result.begin());

  return result;
}
