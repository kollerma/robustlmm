#include "Integration.h"

#define DEBUG(STRING)                                          \
Rcpp::Rcout << STRING << std::endl;

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
iwork_(Calloc(limit_, int)), work_(Calloc(lenw_, double)) {}

DqagIntegration::~DqagIntegration() {
  Free(iwork_); Free(work_);
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
    Rcpp::warning("integration flag " +  to_string(ier_)); 
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
