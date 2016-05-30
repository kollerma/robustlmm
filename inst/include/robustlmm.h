#ifndef __robustlmm_h__
#define __robustlmm_h__

#include <Rcpp.h>

class PsiFunction {
public:
  PsiFunction();
  
  virtual const std::string name();
  virtual const std::string show();
  virtual void chgDefaults(Rcpp::NumericVector tDefs);
  virtual Rcpp::NumericVector tDefs() ;
  virtual const std::string showDefaults();
  
  virtual const double rhoFun(const double x);
  virtual const double psiFun(const double x);
  virtual const double wgtFun(const double x);
  virtual const double DpsiFun(const double x);
  virtual const double DwgtFun(const double x);
  const double psi2Fun(const double x);
  
  virtual const double Erho();
  virtual const double Epsi2();
  virtual const double EDpsi();
  
  virtual ~PsiFunction();
};

typedef const double (PsiFunction::*Fptr)(const double);

class PsiFunctionNumIntExp : public PsiFunction {
public:
  PsiFunctionNumIntExp();
  
  const std::string name();
  void chgDefaults(Rcpp::NumericVector tDefs);
  
  virtual const double Erho();
  virtual const double Epsi2();
  virtual const double EDpsi();
  
  ~PsiFunctionNumIntExp();
  
private:
  double Erho_;
  double Epsi2_;
  double EDpsi_;
  
  void reset();
  
  const double computeErho();
  const double computeEpsi2();
  const double computeEDpsi();
  double integrate(Fptr fptr);
};

#endif // __robustlmm_h__