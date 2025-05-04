#include "PsiFunction.h"
#include <robustbase.h>
#include <Rinternals.h>
using namespace Rcpp;

static const double DEFAULT_K = 1.345;
static const double DEFAULT_S = 10.;

#define DEBUG(STRING)                                          \
// Rcpp::Rcout << STRING << std::endl;

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// class PsiFunction
PsiFunction::PsiFunction() {}

const std::string PsiFunction::name() const {
  return "classic (x^2/2)";
}

const std::string PsiFunction::show() const {
  return this->name() + " psi function" + this->showDefaults();
}


void PsiFunction::chgDefaults(NumericVector tuningParameters) {
  if (this->needToChgDefaults(tuningParameters)) {
    this->doChgDefaults(tuningParameters);
  }
}

bool PsiFunction::needToChgDefaults(NumericVector tuningParameters) {
  return false;
}

void PsiFunction::doChgDefaults(NumericVector tuningParameters) {}

NumericVector PsiFunction::tDefs() const {
  return NumericVector(0);
}

const std::string PsiFunction::showDefaults() const {
  return "";
}

double PsiFunction::rhoFun(const double x) {
  return x * x / 2.;
}

double PsiFunction::psiFun(const double x) {
  return x;
}

double PsiFunction::wgtFun(const double x) {
  return 1.;
}

double PsiFunction::DpsiFun(const double x) {
  return 1.;
}

double PsiFunction::DwgtFun(const double x) {
  return 0.;
}

double PsiFunction::Erho() {
  return 0.5;
}

double PsiFunction::Epsi2() {
  return 1.;
}

double PsiFunction::EDpsi() {
  return 1.;
}

PsiFunction::~PsiFunction() {}

double PsiFunction::psi2Fun(const double x) {
  double value = this->psiFun(x);
  return value * value;
}

// end class PsiFuncion

// class PsiFunctionNumIntExp
PsiFunctionNumIntExp::PsiFunctionNumIntExp() :
  PsiFunction(), integration_(*(new DqagIntegration())) {
  reset();
}

const std::string PsiFunctionNumIntExp::name() const {
  return "PsiFunction with expectations computed using numerical integration";
}

void PsiFunctionNumIntExp::doChgDefaults(NumericVector tuningParameters) {
  reset();
  PsiFunction::doChgDefaults(tuningParameters);
}

double PsiFunctionNumIntExp::Erho() {
  if (NumericVector::is_na(Erho_))
    Erho_ = computeErho();
  return Erho_;
}

double PsiFunctionNumIntExp::Epsi2() {
  if (NumericVector::is_na(Epsi2_))
    Epsi2_ = computeEpsi2();
  return Epsi2_;
}

double PsiFunctionNumIntExp::EDpsi() {
  if (NumericVector::is_na(EDpsi_))
    EDpsi_ = computeEDpsi();
  return EDpsi_;
}

PsiFunctionNumIntExp::~PsiFunctionNumIntExp() {
  delete &integration_;
}

void PsiFunctionNumIntExp::reset() {
  Erho_ = NA_REAL;
  Epsi2_ = NA_REAL;
  EDpsi_ = NA_REAL;
}

double PsiFunctionNumIntExp::computeErho() {
  DEBUG("Called computeErho()")
  return integrate(&PsiFunction::rhoFun);
}

double PsiFunctionNumIntExp::computeEpsi2() {
  return integrate(&PsiFunction::psi2Fun);
}

double PsiFunctionNumIntExp::computeEDpsi() {
  return integrate(&PsiFunction::DpsiFun);
}

double PsiFunctionNumIntExp::integrate(Fptr fptr) {
  const void *exc[2] = {this, &fptr};
  void **ex = const_cast<void**>(exc);
  return integration_.ninfInf(psiFunctionIntegrandNorm, ex);
}

// end class PsiFunctionNumIntExp

// class PsiFunctionPropII
PsiFunctionPropII::PsiFunctionPropII() :
  PsiFunctionNumIntExp(), base_(new SmoothPsi()), integration_(*(new DqagIntegration())) {
  // Rcpp::Rcout << "illegal contstructor called!!" << std::endl;
}

PsiFunctionPropII::PsiFunctionPropII(PsiFunction* base) :
  PsiFunctionNumIntExp(), base_(base), integration_(*(new DqagIntegration()))
{}

PsiFunctionPropII::~PsiFunctionPropII() {
  delete &integration_;
}

const std::string PsiFunctionPropII::name() const {
  return base_->name() + ", Proposal 2";
}

bool PsiFunctionPropII::needToChgDefaults(NumericVector x) {
  return base_->needToChgDefaults(x);
}

void PsiFunctionPropII::doChgDefaults(NumericVector x) {
  base_->doChgDefaults(x);
  PsiFunctionNumIntExp::doChgDefaults(x);
}

NumericVector PsiFunctionPropII::tDefs() const {
  return base_->tDefs();
}

double PsiFunctionPropII::rhoFun(const double x) {
  if (!R_FINITE(x))
    return x;
  return integrate(&PsiFunction::psiFun, x);
}

double PsiFunctionPropII::psiFun(const double x) {
  return base_->wgtFun(x) * base_->psiFun(x);
}

double PsiFunctionPropII::wgtFun(const double x) {
  double value = base_->wgtFun(x);
  return value * value;
}

double PsiFunctionPropII::DpsiFun(const double x) {
  return base_->wgtFun(x) * base_->DpsiFun(x) +
    base_->DwgtFun(x) * base_->psiFun(x);
}

double PsiFunctionPropII::DwgtFun(const double x) {
  return 2. * base_->wgtFun(x) * base_->DwgtFun(x);
}

const PsiFunction* PsiFunctionPropII::base() const {
  return base_;
}

const std::string PsiFunctionPropII::showDefaults() const {
  return base_->showDefaults();
}

double PsiFunctionPropII::integrate(Fptr fptr, double b) {
  const void *exc[2] = {this, &fptr};
  void **ex = const_cast<void**>(exc);
  double a = 0.;
  return integration_.aB(psiFunctionIntegrand, ex, &a, &b);
}

// end class PsiFunctionPropII

// class HuberPsi
HuberPsi::HuberPsi() : PsiFunction() {
  doChgDefaults(NumericVector(0));
}

HuberPsi::HuberPsi(NumericVector k) : PsiFunction() {
  doChgDefaults(k);
}

const std::string HuberPsi::name() const {
  return "Huber";
}

bool HuberPsi::needToChgDefaults(NumericVector k) {
  if (k.size() >= 1) {
    return k_ != k[0];
  } else {
    return k_ != DEFAULT_K;
  }
}

void HuberPsi::doChgDefaults(NumericVector k) {
  if (k.size() >= 1) {
    k_ = k[0];
  } else {
    k_ = DEFAULT_K;
  }
}

NumericVector HuberPsi::tDefs() const {
  NumericVector tDefs = NumericVector(1);
  tDefs[0] = k_;
  tDefs.names() = CharacterVector::create("k");
  return tDefs;
}

const std::string HuberPsi::showDefaults() const {
  return tfm::format(" (k = %.5g)", k_);
}

double HuberPsi::rhoFun(const double x) {
  double u = std::abs(x);
  if (u > k_) {
    return k_*(u - k_ / 2.);
  } else {
    return u * u / 2.;
  }
}

double HuberPsi::psiFun(const double x) {
  if (x < -k_) {
    return -k_;
  } else if (x > k_) {
    return k_;
  } else {
    return x;
  }
}

double HuberPsi::wgtFun(const double x) {
  if (x < -k_ || x > k_) {
    return k_ / std::abs(x);
  } else {
    return 1.;
  }
}

double HuberPsi::DpsiFun(const double x) {
  if (x < -k_ || x > k_) {
    return 0.;
  } else {
    return 1.;
  }
}

double HuberPsi::DwgtFun(const double x) {
  if (x < -k_) {
    return k_ / (x * x);
  } else if (x > k_) {
    return -k_ / (x * x);
  } else {
    return 0.;
  }
}

double HuberPsi::Erho() {
  double iP = stats::pnorm_0(k_, 0, 0);
  return .5 - iP + k_ * (stats::dnorm_0(k_, 0) - k_ * iP);
}

double HuberPsi::Epsi2() {
  if (k_ < 10.) {
    return 1. - 2.*(k_ * stats::dnorm_0(k_, 0) + (1. - k_*k_) * stats::pnorm_0(k_, 0, 0));
  } else {
    return 1.;
  }
}

double HuberPsi::EDpsi() {
  return 2. * stats::pnorm_0(k_, 1, 0) - 1.;
}

HuberPsi::~HuberPsi() {}

// end class HuberPsi

// class SmoothPsi
SmoothPsi::SmoothPsi() : PsiFunctionNumIntExp() {
  doChgDefaults(NumericVector(0));
}

SmoothPsi::SmoothPsi(NumericVector tuningParameters) : PsiFunctionNumIntExp() {
  doChgDefaults(tuningParameters);
}

const std::string SmoothPsi::name() const {
  return "smoothed Huber";
}

bool SmoothPsi::needToChgDefaults(NumericVector tuningParameters) {
  bool result = false;
  if (tuningParameters.size() >= 1) {
    result = k_ != tuningParameters[0];
  } else {
    result = k_ != DEFAULT_K;
  }
  if (result) {
    return result;
  }
  if (tuningParameters.size() >= 2) {
    result = s_ != tuningParameters[1];
  } else {
    result = s_ != DEFAULT_S;
  }
  return result;
}

void SmoothPsi::doChgDefaults(NumericVector tuningParameters) {
  PsiFunctionNumIntExp::doChgDefaults(tuningParameters);
  if (tuningParameters.size() >= 1) {
    k_ = tuningParameters[0];
  } else {
    k_ = DEFAULT_K;
  }
  if (tuningParameters.size() >= 2) {
    s_ = tuningParameters[1];
  } else {
    s_ = DEFAULT_S;
  }
  a_ = std::pow(s_, 1. / (s_ + 1.));
  c_ = k_ - std::pow(a_, -s_);
  d_ = c_ - a_;
}

NumericVector SmoothPsi::tDefs() const {
  NumericVector tDefs = NumericVector(2);
  tDefs[0] = k_;
  tDefs[1] = s_;
  tDefs.names() = CharacterVector::create("k", "s");
  return tDefs;
}

double SmoothPsi::rhoFun(const double x) {
  double ax = std::abs(x);
  if (ax <= c_) {
    return x * x / 2.;
  } else {
    return c_ * c_ / 2. + k_ * (ax - c_) -
      (std::pow(ax - d_, 1. - s_) - std::pow(a_, 1. - s_)) / (1. - s_);
  }
}

double SmoothPsi::psiFun(const double x) {
  double ax = std::abs(x);
  if (ax <= c_) {
    return x;
  } else {
    return sgn(x) * (k_ - std::pow(ax - d_, -s_));
  }
}

double SmoothPsi::DpsiFun(const double x) {
  double ax = std::abs(x);
  if (ax <= c_) {
    return 1.;
  } else {
    return s_ * std::pow(ax - d_, -s_ - 1.);
  }
}

double SmoothPsi::wgtFun(const double x) {
  double ax = std::abs(x);
  if (ax <= c_) {
    return 1.;
  } else {
    return (k_ - std::pow(ax - d_, -s_)) / ax;
  }
}

double SmoothPsi::DwgtFun(const double x) {
  double ax = std::abs(x);
  if (ax <= c_) {
    return 0.;
  } else {
    return std::pow(ax - d_, -s_ - 1.) * s_ / x -
      (k_ - std::pow(ax - d_, -s_)) / (x * ax);
  }
}

SmoothPsi::~SmoothPsi() {}

const std::string SmoothPsi::showDefaults() const {
  return tfm::format(" (k = %.5g, s = %.5g)", k_, s_);
}

// end class SmoothPsi

// class RobustbasePsi
RobustbasePsi::RobustbasePsi(NumericVector tuningParameters, int ipsi) : ipsi_(ipsi) {
  chgDefaults(tuningParameters);
}

void RobustbasePsi::chgDefaults(NumericVector tuningParameters) {
  PsiFunctionNumIntExp::chgDefaults(tuningParameters);
  initialiseTuningParametersFromDefaults();
  if (tuningParameters.hasAttribute("names")) {
    chgDefaultsUsingNamedVector(tuningParameters);
  } else {
    chgDefaultsUsingPositionInVector(tuningParameters);
  }
}

void RobustbasePsi::initialiseTuningParametersFromDefaults() {
  if (tuningParameters_ != NULL)
    return;
  const NumericVector defaults = this->getDefaults();
  tuningParameters_ = new double[defaults.size()];
  std::copy(defaults.begin(), defaults.end(), tuningParameters_);
}

void RobustbasePsi::chgDefaultsUsingNamedVector(const NumericVector &tuningParameters) {
  const NumericVector defaults = this->getDefaults();
  const std::vector<std::string> names(tuningParameters.attributeNames());
  unsigned npar = tuningParameters.size();
  R_ASSERT(names.size() == npar);
  for (unsigned i = 0; i < npar; ++i) {
    std::string name = names[i];
    if (!defaults.containsElementNamed(name.data()))
      throw std::invalid_argument("no tuning parameter for name " + name + ".");
    tuningParameters_[defaults.findName(name)] = tuningParameters[i];
  }
}

void RobustbasePsi::chgDefaultsUsingPositionInVector(const NumericVector &tuningParameters) {
  std::copy(tuningParameters.begin(), tuningParameters.end(), tuningParameters_);
}

NumericVector RobustbasePsi::tDefs() const {
  NumericVector defs(this->getDefaults());
  std::copy(tuningParameters_, tuningParameters_ + defs.size(), defs.begin());
  return defs;
}

const std::string RobustbasePsi::showDefaults() const {
  std::vector<std::string> names = this->getDefaults().attributeNames();
  std::stringstream ss;
  ss << " (";
  std::string sep = "";
  for (unsigned i = 0; i < names.size(); i++) {
    ss << sep <<  names[i] << " = " << tuningParameters_[i];
    sep = ", ";
  }
  ss << ")";
  return ss.str();
}

double RobustbasePsi::rhoFun(const double x) {
  return C_rho(x, tuningParameters_, ipsi_);
}

double RobustbasePsi::psiFun(const double x) {
  return C_psi(x, tuningParameters_, ipsi_);
}

double RobustbasePsi::DpsiFun(const double x) {
  return C_psip(x, tuningParameters_, ipsi_);
}

double RobustbasePsi::wgtFun(const double x) {
  return C_wgt(x, tuningParameters_, ipsi_);
}

double RobustbasePsi::DwgtFun(const double x) {
  if (x == 0.) return 0;
  return (DpsiFun(x) - wgtFun(x)) / x;
}

RobustbasePsi::~RobustbasePsi() {
  delete(tuningParameters_);
}

// end class RobustbasePsi

/*
 TODO: lqq, bisquare, etc, psi functions

lqqPsi <- psiFuncCached(rho = function(x, cc) Mpsi(x, cc, "lqq", -1),
  psi = function(x, cc) Mpsi(x, cc, "lqq", 0),
  Dpsi = function(x, cc) Mpsi(x, cc, "lqq", 1),
  wgt = function(x, cc) Mwgt(x, cc, "lqq"),
  Dwgt = function(x, cc)
  (Mpsi(x, cc, "lqq", 1) - Mwgt(x, cc, "lqq"))/x,
  name = "lqq",
  cc = c(-0.5, 1.5, 0.95, NA))

bisquarePsi <- psiFuncCached(rho = function(x, k) Mpsi(x, k, "biweight", -1),
  psi = function(x, k) Mpsi(x, k, "biweight", 0),
  Dpsi = function(x, k) Mpsi(x, k, "biweight", 1),
  wgt = function(x, k) (1 - (x/k)^2)^2*(abs(x) <= k),
  Dwgt = function(x, k) (-(4*(1-(x/k)^2))*x/k^2)*(abs(x) <= k),                                                                                                   name = "bisquare",
  k = 4.68)
 */

std::string name(PsiFunction* p) {
  return p->name();
}

void chgDefaults(PsiFunction* p, NumericVector x) {
  return p->chgDefaults(x);
}

NumericVector compute(PsiFunction* p, Fptr fptr, NumericVector x) {
  NumericVector result(x.size());
  for (unsigned i = 0; i < x.size(); i++) {
    result[i] = (p->*fptr)(x[i]);
  }
  return result;
}

NumericVector rho(PsiFunction* p, NumericVector x) {
  return compute(p, &PsiFunction::rhoFun, x);
}

NumericVector psi(PsiFunction* p, NumericVector x) {
  return compute(p, &PsiFunction::psiFun, x);
}

NumericVector wgt(PsiFunction* p, NumericVector x) {
  return compute(p, &PsiFunction::wgtFun, x);
}

NumericVector Dpsi(PsiFunction* p, NumericVector x) {
  return compute(p, &PsiFunction::DpsiFun, x);
}

NumericVector Dwgt(PsiFunction* p, NumericVector x) {
  return compute(p, &PsiFunction::DwgtFun, x);
}

double Erho(PsiFunction* p) {
  return p->Erho();
}

double Epsi2(PsiFunction* p) {
  return p->Epsi2();
}

double EDpsi(PsiFunction* p) {
  return p->EDpsi();
}

NumericVector tDefs(PsiFunction* p) {
  return p->tDefs();
}

void psiFunctionIntegrand(double *x, const int n, void *const ex) {
  PsiFunction **pp = static_cast<PsiFunction**>(ex);
  PsiFunction *p = pp[0];

  const Fptr **const pfptr  = static_cast<const Fptr **const>(ex);
  const Fptr *fptr = pfptr[1];

  for (int i = 0; i < n; i++) {
    double value = x[i];
    x[i] = (p->**fptr)(value);
  }

  return;
}

void psiFunctionIntegrandNorm(double *x, const int n, void *const ex) {
  PsiFunction **pp = static_cast<PsiFunction**>(ex);
  PsiFunction *p = pp[0];

  const Fptr **const pfptr  = static_cast<const Fptr **const>(ex);
  const Fptr *fptr = pfptr[1];

  for (int i = 0; i < n; i++) {
    double value = x[i];
    x[i] = (p->**fptr)(value) * stats::dnorm_0(value, 0);
  }

  return;
}

extern "C" {

  SEXP isnull(SEXP pointer) {
    return Rf_ScalarLogical(!R_ExternalPtrAddr(pointer));
  }

  SEXP deepcopy(SEXP x) {
    return(Rf_duplicate(x));
  }

}

RCPP_EXPOSED_CLASS(PsiFunction)
RCPP_MODULE(psi_function_module) {

  class_<PsiFunction>("PsiFunction")
  .constructor()
  .method("name", &name)
  .method("chgDefaults", &chgDefaults)
  .method("rho", &rho)
  .method("psi", &psi)
  .method("wgt", &wgt)
  .method("Dpsi", &Dpsi)
  .method("Dwgt", &Dwgt)
  .method("Epsi2", &Epsi2)
  .method("EDpsi", &EDpsi)
  .method("Erho", &Erho)
  .method("show", &PsiFunction::show)
  .method("tDefs", &tDefs)
  ;

  class_<HuberPsi>("HuberPsi")
    .derives<PsiFunction>("PsiFunction")
    .constructor()
    .constructor<NumericVector>()
  ;

  class_<SmoothPsi>("SmoothPsi")
    .derives<PsiFunction>("PsiFunction")
    .constructor()
    .constructor<NumericVector>()
  ;

  class_<PsiFunctionPropII>("PsiFunctionToPropIIPsiFunctionWrapper")
    .derives<PsiFunction>("PsiFunction")
    .constructor()
    .constructor<PsiFunction*>()
    .method("base", &PsiFunctionPropII::base)
  ;

  function("isnull", &isnull);
  function("deepcopy", &deepcopy);
}
