#include "globals.h"
using namespace Rcpp;

namespace std {
template<typename T>
std::string to_string(const T &n) {
  std::ostringstream s;
  s << n;
  return s.str();
}
} // namespace std

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// class PsiFunction
PsiFunction::PsiFunction() {}

const std::string PsiFunction::name() { 
  return "classic (x^2/2)"; 
}

const std::string PsiFunction::show() {
  return this->name() + " psi function" + this->showDefaults();
}

void PsiFunction::chgDefaults(NumericVector tuningParameters) {};

NumericVector PsiFunction::tDefs() {
  return NumericVector(0);
};

const std::string PsiFunction::showDefaults() {
  return "";  
}

const double PsiFunction::rhoFun(const double x) {
  return x * x / 2.;
}

const double PsiFunction::psiFun(const double x) {
  return x;
}

const double PsiFunction::wgtFun(const double x) {
  return 1.;
}

const double PsiFunction::DpsiFun(const double x) {
  return 1.;
}

const double PsiFunction::DwgtFun(const double x) {
  return 0.;
}

const double PsiFunction::Erho() {
  return 0.5;
}

const double PsiFunction::Epsi2() {
  return 1.;
}

const double PsiFunction::EDpsi() {
  return 1.;
}

PsiFunction::~PsiFunction() {}

const double PsiFunction::psi2Fun(const double x) {
  double value = this->psiFun(x);
  return value * value;
}

// end class PsiFuncion

class HuberPsi : PsiFunction {
public:
  HuberPsi() : PsiFunction() {
    chgDefaults(NumericVector(0));
  }
  
  HuberPsi(NumericVector k) : PsiFunction() {
    chgDefaults(k); 
  }
  
  const std::string name() {
    return "Huber";
  }
  
  void chgDefaults(NumericVector k) {
    if (k.size() >= 1) {
      k_ = k[0];
    } else {
      k_ = 1.345;
    }
  }
  
  NumericVector tDefs() {
    NumericVector tDefs = NumericVector(1);
    tDefs[0] = k_;
    tDefs.names() = CharacterVector::create("k");
    return tDefs;
  }
  
  const double rhoFun(const double x) {
    double u = std::abs(x);
    if (u > k_) {
      return k_*(u - k_ / 2.);
    } else {
      return u * u / 2.;
    }
  }
  
  const double psiFun(const double x) {
    if (x < -k_) {
      return -k_;
    } else if (x > k_) {
      return k_;
    } else {
      return x;
    }
  }
  
   const double wgtFun(const double x) {
     if (x < -k_ || x > k_) {
       return k_ / std::abs(x);
     } else {
       return 1.;
     }
  }
  
   const double DpsiFun(const double x) {
     if (x < -k_ || x > k_) {
       return 0.;
     } else {
       return 1.;
     }
   }
  
   const double DwgtFun(const double x) {
     if (x < -k_) {
       return k_ / (x * x);
     } else if (x > k_) {
       return -k_ / (x * x);
     } else {
       return 0.;
     }
  }
  
   const double Erho() {
     double iP = stats::pnorm_0(k_, 0, 0);
     return .5 - iP + k_ * (stats::dnorm_0(k_, 0) - k_ * iP);
  }
  
   const double Epsi2() {
     if (k_ < 10.) {
       return 1. - 2.*(k_ * stats::dnorm_0(k_, 0) + (1. - k_*k_) * stats::pnorm_0(k_, 0, 0));
     } else {
       return 1.;
     }
  }
  
   const double EDpsi() {
    return 2. * stats::pnorm_0(k_, 1, 0) - 1.;
  }
  
   ~HuberPsi() {}
  
  const std::string showDefaults() {
    char buffer[20];
    std::sprintf(buffer, " (k = %.5g)", k_);  
    return std::string(buffer);
  }
  
private:
  double k_;
};

void integrand(double *x, const int n, void *const ex) {
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

void integrandNorm(double *x, const int n, void *const ex) {
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

// class PsiFunctionNumIntExp 
PsiFunctionNumIntExp::PsiFunctionNumIntExp() : 
  PsiFunction() {
    reset();
}

const std::string PsiFunctionNumIntExp::name() {
  return "PsiFunction with expectations computed using numerical integration";
}

void PsiFunctionNumIntExp::chgDefaults(NumericVector tuningParameters) {
  reset();
  PsiFunction::chgDefaults(tuningParameters);
};

const double PsiFunctionNumIntExp::Erho() {
  if (NumericVector::is_na(Erho_))
    Erho_ = computeErho();
  return Erho_;
};

const double PsiFunctionNumIntExp::Epsi2() {
  if (NumericVector::is_na(Epsi2_))
    Epsi2_ = computeEpsi2();
  return Epsi2_;
};

const double PsiFunctionNumIntExp::EDpsi() {
  if (NumericVector::is_na(EDpsi_))
    EDpsi_ = computeEDpsi();
  return EDpsi_;
}

PsiFunctionNumIntExp::~PsiFunctionNumIntExp() {}

void PsiFunctionNumIntExp::reset() {
  Erho_ = NA_REAL;
  Epsi2_ = NA_REAL;
  EDpsi_ = NA_REAL;
}

const double PsiFunctionNumIntExp::computeErho() {
  return integrate(&PsiFunction::rhoFun);
}

const double PsiFunctionNumIntExp::computeEpsi2() {
  return integrate(&PsiFunction::psi2Fun);
}

const double PsiFunctionNumIntExp::computeEDpsi() {
  return integrate(&PsiFunction::DpsiFun);
}

double PsiFunctionNumIntExp::integrate(Fptr fptr) {
  int inf_flag = 2, neval = 0, ier = 0, limit = 100, lenw = 4*limit, last = 0;
  int* iwork = Calloc(limit, int);
  
  double epsabs = std::pow(std::numeric_limits<double>::epsilon(), .5),  
    epsrel = epsabs, result = 0, abserr = 0;
  double * work = Calloc(lenw, double);
  
  const void *exc[2] = {this, &fptr};
  void **ex = const_cast<void**>(exc);
  
  double bound = NA_REAL;
  
  Rdqagi(integrandNorm, ex, &bound, &inf_flag, &epsabs, &epsrel, &result, &abserr, &neval, &ier, 
         &limit, &lenw, &last, iwork, work);
  
  if(ier > 0 && ier!=5) warning("integration flag " +  std::to_string(ier));
  Free(iwork); Free(work);
  
  return result;
}

// end class PsiFunctionNumIntExp

class SmoothPsi : public PsiFunctionNumIntExp {
public:
  SmoothPsi() : PsiFunctionNumIntExp() {
    chgDefaults(NumericVector(0));
  }
  
  SmoothPsi(NumericVector tuningParameters) : PsiFunctionNumIntExp() {
    chgDefaults(tuningParameters);
  }

  const std::string name() {
    return "smoothed Huber";
  }
  
  void chgDefaults(NumericVector tuningParameters) {
    PsiFunctionNumIntExp::chgDefaults(tuningParameters);
    if (tuningParameters.size() >= 1) {
      k_ = tuningParameters[0];
    } else {
      k_ = 1.345;
    }
    if (tuningParameters.size() >= 2) {
      s_ = tuningParameters[1];
    } else {
      s_ = 10.;
    }
    a_ = std::pow(s_, 1. / (s_ + 1.));
    c_ = k_ - std::pow(a_, -s_);
    d_ = c_ - a_;
  }
  
  NumericVector tDefs() {
    NumericVector tDefs = NumericVector(2);
    tDefs[0] = k_;
    tDefs[1] = s_;
    tDefs.names() = CharacterVector::create("k", "s");
    return tDefs;
  }
  
  const double rhoFun(const double x) {
    double ax = std::abs(x);
    if (ax <= c_) {
      return x * x / 2.;
    } else {
      return c_ * c_ / 2. + k_ * (ax - c_) -
        (std::pow(ax - d_, 1. - s_) - std::pow(a_, 1. - s_)) / (1. - s_);
    }
  }
  
  const double psiFun(const double x) {
    double ax = std::abs(x);
    if (ax <= c_) {
      return x;
    } else {
      return sgn(x) * (k_ - std::pow(ax - d_, -s_));
    }
  }
  
  const double DpsiFun(const double x) {
    double ax = std::abs(x);
    if (ax <= c_) {
      return 1.;
    } else {
      return s_ * std::pow(ax - d_, -s_ - 1.);
    }
  }
  
  const double wgtFun(const double x) {
    double ax = std::abs(x);
    if (ax <= c_) {
      return 1.;
    } else {
      return (k_ - std::pow(ax - d_, -s_)) / ax;
    }
  }
  
  const double DwgtFun(const double x) {
    double ax = std::abs(x);
    if (ax <= c_) {
      return 0.;
    } else {
      return std::pow(ax - d_, -s_ - 1.) * s_ / x - 
        (k_ - std::pow(ax - d_, -s_)) / (x * ax);
    }
  }
  
  ~SmoothPsi() {}
  
  const std::string showDefaults() {
    char buffer[30];
    std::sprintf(buffer, " (k = %.5g, s = %.5g)", k_, s_);  
    return std::string(buffer);
  }
  
private:
  double k_;
  double s_;
  double a_;
  double c_;
  double d_;
  
};

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
  for (int i = 0; i < x.size(); i++) {
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

const double Erho(PsiFunction* p) {
  return p->Erho();
}

const double Epsi2(PsiFunction* p) {
  return p->Epsi2();
}

const double EDpsi(PsiFunction* p) {
  return p->EDpsi();
}

NumericVector tDefs(PsiFunction* p) {
  return p->tDefs();
}

class PsiFunctionPropII: public PsiFunctionNumIntExp {
public:
  PsiFunctionPropII(PsiFunction* base) : PsiFunctionNumIntExp(), base_(base) {
    init();
  }
  
  const std::string name() {
    return base_->name() + ", Proposal II";
  }
  
  void chgDefaults(NumericVector x) {
    base_->chgDefaults(x);
    PsiFunctionNumIntExp::chgDefaults(x);
  }
  
  NumericVector tDefs() {
    return base_->tDefs();
  }
  
  const double rhoFun(const double x) {
    if (!R_FINITE(x))
      return x;
    return integrate(&PsiFunction::psiFun, x);
  }
  
  const double psiFun(const double x) {
    return base_->wgtFun(x) * base_->psiFun(x);
  }
  
  const double wgtFun(const double x) {
    double value = base_->wgtFun(x);
    return value * value;
  }
  
  const double DpsiFun(const double x) {
    return base_->wgtFun(x) * base_->DpsiFun(x) + 
      base_->DwgtFun(x) * base_->psiFun(x);
  }
  
  const double DwgtFun(const double x) {
    return 2. * base_->wgtFun(x) * base_->DwgtFun(x);
  }
  
  ~PsiFunctionPropII() {
    Free(iwork_); Free(work_);
  }
  
  const std::string showDefaults() {
    return base_->showDefaults();
  }
  
private:
  PsiFunction* base_;
  int limit_;
  int lenw_;
  int* iwork_;
  double epsabs_;
  double epsrel_;
  double* work_;
  
  void init() {
    limit_ = 100;
    lenw_ = 4*limit_;
    iwork_ = Calloc(limit_, int);
    epsabs_ = std::pow(std::numeric_limits<double>::epsilon(), .5);
    epsrel_ = epsabs_;
    work_ = Calloc(lenw_, double);
  }
  
  double integrate(Fptr fptr, double b) {
    int neval = 0, ier = 0, last = 0;
    double result = 0., abserr = 0., a = 0.;
    
    const void *exc[2] = {this, &fptr};
    void **ex = const_cast<void**>(exc);
    Rdqags(integrand, ex, &a, &b, &epsabs_, &epsrel_, &result, &abserr, &neval, &ier, 
           &limit_, &lenw_, &last, iwork_, work_);
    
    if(ier > 0 && ier!=5) warning("integration flag " +  std::to_string(ier));
    
    return result;
  }
};

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
  
  class_<PsiFunctionPropII>("PsiFunction to Prop II PsiFunction wrapper")
    .derives<PsiFunction>("PsiFunction")
    .constructor<PsiFunction*>()
  ;
}


