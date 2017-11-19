#if !defined  ROBUSTLMM_PSIFUNCTION_H__
#define  ROBUSTLMM_PSIFUNCTION_H__

#include <Rcpp.h>
#include "Integration.h"

using namespace Rcpp;

class PsiFunction {
public:
  PsiFunction();
  
  virtual const std::string name() const;
  virtual const std::string show() const;
  virtual void chgDefaults(Rcpp::NumericVector tDefs);
  virtual Rcpp::NumericVector tDefs() const;
  virtual const std::string showDefaults() const;
  
  virtual double rhoFun(const double x);
  virtual double psiFun(const double x);
  virtual double wgtFun(const double x);
  virtual double DpsiFun(const double x);
  virtual double DwgtFun(const double x);
  double psi2Fun(const double x);
  
  virtual double Erho();
  virtual double Epsi2();
  virtual double EDpsi();
  
  virtual ~PsiFunction();
};

typedef double (PsiFunction::*Fptr)(const double);
typedef Rcpp::XPtr<PsiFunction> PsiFuncXPtr;

class PsiFunctionNumIntExp : public PsiFunction {
public:
  PsiFunctionNumIntExp();
  
  const std::string name() const;
  void chgDefaults(Rcpp::NumericVector tDefs);
  
  virtual double Erho();
  virtual double Epsi2();
  virtual double EDpsi();
  
  ~PsiFunctionNumIntExp();
  
private:
  double Erho_;
  double Epsi2_;
  double EDpsi_;
  Integration &integration_;
  
  void reset();
  
  double computeErho();
  double computeEpsi2();
  double computeEDpsi();
  double integrate(Fptr fptr);
};

class HuberPsi : public PsiFunction {
public:
  HuberPsi();
  HuberPsi(NumericVector k);
  
  const std::string name() const;
  void chgDefaults(NumericVector k);
  NumericVector tDefs() const;
  const std::string showDefaults() const;
  
  double rhoFun(const double x);
  double psiFun(const double x);
  double wgtFun(const double x);
  double DpsiFun(const double x);
  double DwgtFun(const double x);
  double Erho();
  double Epsi2();
  double EDpsi();
  
  ~HuberPsi();
  
private:
  double k_;
};

class SmoothPsi : public PsiFunctionNumIntExp {
public:
  SmoothPsi();
  SmoothPsi(NumericVector tuningParameters);
  const std::string name() const;
  void chgDefaults(NumericVector tuningParameters);
  NumericVector tDefs() const;
  const std::string showDefaults() const;
  
  double rhoFun(const double x);
  double psiFun(const double x);
  double DpsiFun(const double x);
  double wgtFun(const double x);
  double DwgtFun(const double x);
  
  ~SmoothPsi();
  
private:
  double k_;
  double s_;
  double a_;
  double c_;
  double d_;
};

class PsiFunctionPropII: public PsiFunctionNumIntExp {
public:
  PsiFunctionPropII();
  PsiFunctionPropII(PsiFunction* base);
  ~PsiFunctionPropII();
  
  const std::string name() const;
  void chgDefaults(NumericVector x);
  NumericVector tDefs() const;
  const std::string showDefaults() const;
  
  double rhoFun(const double x);
  double psiFun(const double x);
  double wgtFun(const double x);
  double DpsiFun(const double x);
  double DwgtFun(const double x);
  
  const PsiFunction* base() const;
  
private:
  PsiFunction* base_;
  Integration &integration_;
  
  double integrate(Fptr fptr, double b);
};

void psiFunctionIntegrand(double *x, const int n, void *const ex);
void psiFunctionIntegrandNorm(double *x, const int n, void *const ex);

std::string name(PsiFunction* p);
void chgDefaults(PsiFunction* p, NumericVector x);
NumericVector compute(PsiFunction* p, Fptr fptr, NumericVector x);
NumericVector rho(PsiFunction* p, NumericVector x);
NumericVector psi(PsiFunction* p, NumericVector x);
NumericVector wgt(PsiFunction* p, NumericVector x);
NumericVector Dpsi(PsiFunction* p, NumericVector x);
NumericVector Dwgt(PsiFunction* p, NumericVector x);
double Erho(PsiFunction* p);
double Epsi2(PsiFunction* p);
double EDpsi(PsiFunction* p);
NumericVector tDefs(PsiFunction* p);

extern "C" SEXP _rcpp_module_boot_psi_function_module();

#endif
