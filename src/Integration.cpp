#include "Integration.h"

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
