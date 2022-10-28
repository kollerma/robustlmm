#if !defined  ROBUSTLMM_INTEGRATION_H__
#define  ROBUSTLMM_INTEGRATION_H__

#include "misc.h"
#include <R_ext/Applic.h>

class IntegrFnEx {
public:
  integr_fn *f_;
  void *ex_;

  IntegrFnEx(integr_fn *f, void *ex);
  ~IntegrFnEx();

};

class Integration  {
public:

  /* typedef void integr_fn(double *x, int n, void *ex); */
  virtual double ninfInf(integr_fn *f, void *ex) = 0;
  virtual double aInf(integr_fn *f, void *ex, double* bound) = 0;
  virtual double ninfB(integr_fn *f, void *ex, double* bound) = 0;
  virtual double aB(integr_fn *f, void *ex, double* a, double* b) = 0;

  virtual int getNeval() = 0;
  virtual int getErrorCode() = 0;
  virtual double getAbserr() = 0;

  Integration();
  virtual ~Integration();

private:
  virtual IntegrFnEx wrap(integr_fn *f, void *ex) = 0;

};

class DqagIntegration : public virtual Integration {
private:
  int neval_, ier_, limit_, lenw_, last_;
  double epsabs_, epsrel_, result_, noBound_, abserr_;

  int* iwork_;
  double* work_;

public:
  DqagIntegration();
  ~DqagIntegration();

  /* typedef void integr_fn(double *x, int n, void *ex); */
  double ninfInf(integr_fn *f, void *ex);
  double aInf(integr_fn *f, void *ex, double* bound);
  double ninfB(integr_fn *f, void *ex, double* bound);
  double aB(integr_fn *f, void *ex, double* a, double* b);

  int getNeval();
  int getErrorCode();
  int getIer();
  int getLast();
  double getAbserr();

private:
  virtual IntegrFnEx wrap(integr_fn *f, void *ex);

  void dqagi(IntegrFnEx fnex, double* bound, int inf);
  void dqags(IntegrFnEx fnex, double* a, double* b);
  void checkIer();
};

#endif
