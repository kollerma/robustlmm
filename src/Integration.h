#if !defined  ROBUSTLMM_INTEGRATION_H__
#define  ROBUSTLMM_INTEGRATION_H__

#include "misc.h"
#include <R_ext/Applic.h>
#include <string.h>
#include <math.h>
#include <exp_cubature_typedefs.h>

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

class Expectation : public virtual Integration {
public:
  Expectation();

};

class NormalExpectation : public Expectation {
public:
  NormalExpectation();

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

class DqagNormalExpectation : public virtual DqagIntegration,
                              public virtual NormalExpectation {
public:
  DqagNormalExpectation();

private:
  IntegrFnEx wrap(integr_fn *f, void *ex);
};

class GaussianQuadrature : public virtual Integration {
private:
  std::vector<double> x_, w_;

public:
  GaussianQuadrature();
  GaussianQuadrature(int n);
  ~GaussianQuadrature();

  /* typedef void integr_fn(double *x, int n, void *ex); */
  double ninfInf(integr_fn *f, void *ex);
  double aInf(integr_fn *f, void *ex, double* bound);
  double ninfB(integr_fn *f, void *ex, double* bound);
  double aB(integr_fn *f, void *ex, double* a, double* b);

  int getNeval();
  int getErrorCode();
  double getAbserr();

private:
  void init(void postInitFun(std::vector<double>& x, std::vector<double>& w));
  virtual IntegrFnEx wrap(integr_fn *f, void *ex);

protected:
  GaussianQuadrature(int n, void postInitFun(std::vector<double>& x, std::vector<double>& w));
  friend class GaussianQuadratureNormalExpectation;
};

void gaussianQuadraturePostInit(std::vector<double>& x, std::vector<double>& w);

class GaussianQuadratureNormalExpectation :
  public virtual GaussianQuadrature,
  public virtual NormalExpectation {
public:
  GaussianQuadratureNormalExpectation();
  GaussianQuadratureNormalExpectation(int n);

private:
  virtual IntegrFnEx wrap(integr_fn *f, void *ex);
};

void gaussianQuadratureNormalExpectationPostInit(std::vector<double>& x, std::vector<double>& w);

void dqagNormalExpectationWrapper(double *x, const int n, void *const ex);

void integrandRfun(double *x, const int n, void *const ex);

double test_DqagNormalExpectation(Rcpp::Function func);
double test_GaussHermiteQuadrature(int nNodes, Rcpp::Function func);
double test_GaussHermiteNormalExpectation(int nNodes, Rcpp::Function func);

// 2d integration

void integrand2d_ninfInf(double *x, const int n, void *const ex);
void integrand2d_aInf(double *x, const int n, void *const ex);
void integrand2d_ninfB(double *x, const int n, void *const ex);
void integrand2d_aB(double *x, const int n, void *const ex);

class IntegrFnExInner {
public:
  Integration *integration_;
  integr_fn *f_;
  void *ex_;
  double *a_, *b_;

  IntegrFnExInner(Integration *integration, integr_fn *f, void *ex);
  IntegrFnExInner(Integration *integration, integr_fn *f, void *ex, double *a);
  IntegrFnExInner(Integration *integration, integr_fn *f, void *ex, double *a, double *b);
  ~IntegrFnExInner();
};

class Integration2d : public Integration {
public:
  Integration2d();
};

class Expectation2d : public virtual Integration2d {
public:
  Expectation2d();
};

class SimpleIntegration2d : public virtual Integration2d {
  Integration *outer_, *inner_;
public:
  double ninfInf(integr_fn *f, void *ex);
  double aInf(integr_fn *f, void *ex, double* bound);
  double ninfB(integr_fn *f, void *ex, double* bound);
  double aB(integr_fn *f, void *ex, double* a, double* b);

  int getNeval();
  int getErrorCode();
  double getAbserr();

  SimpleIntegration2d(Integration *outer, Integration *inner);
  ~SimpleIntegration2d();

  void deleteOuterInner();

protected:
  IntegrFnEx wrap(integr_fn *f, void *ex);

  friend class SimpleExpectation2d;
};

class SimpleExpectation2d : public virtual Expectation2d, public SimpleIntegration2d {
public:
  SimpleExpectation2d(Expectation *outer, Expectation *inner);
};

class NormalExpectation2d : public virtual Expectation2d {
public:
  NormalExpectation2d();
};

class SimpleNormalExpectation2d : public NormalExpectation2d, public SimpleExpectation2d {
public:
  SimpleNormalExpectation2d(NormalExpectation *outer, NormalExpectation *inner);
};

class DqagNormalExpectation2d : public SimpleNormalExpectation2d {
public:
  DqagNormalExpectation2d();
  ~DqagNormalExpectation2d();
};

class GaussianQuadratureNormalExpectation2d :
  public SimpleNormalExpectation2d {
public:
  GaussianQuadratureNormalExpectation2d();
  GaussianQuadratureNormalExpectation2d(int n);
  ~GaussianQuadratureNormalExpectation2d();
};

void integrand2dRfun(double *x, const int n, void *const ex);

double test_DqagIntegration2d_ninfInf(Rcpp::Function func);
double test_DqagIntegration2d_aInf(Rcpp::Function func, Rcpp::NumericVector bound);
double test_DqagIntegration2d_ninfB(Rcpp::Function func, Rcpp::NumericVector bound);
double test_DqagIntegration2d_aB(Rcpp::Function func, Rcpp::NumericVector a, Rcpp::NumericVector b);
double test_DqagNormalExpectation2d_ninfInf(Rcpp::Function func);
double test_GaussHermiteNormalExpectation2d(Rcpp::Function func, int nodes);

// Multi dimensional integration

class IntegrandNd {
public:
  integrand f_;
  const int ndim_;
  const int fdim_;
  void *ex_;

  IntegrandNd(integrand f, const int ndim, const int fdim, void *ex);
  ~IntegrandNd();

};

class IntegrationNd  {
public:

  virtual std::vector<double> ninfInf(IntegrandNd *f) = 0;
  virtual std::vector<double> aInf(IntegrandNd *f, double* bound) = 0;
  virtual std::vector<double> ninfB(IntegrandNd *f, double* bound) = 0;
  virtual std::vector<double> aB(IntegrandNd *f, double* a, double* b) = 0;

  virtual int getNeval() = 0;
  virtual int getErrorCode() = 0;
  virtual std::vector<double>& getAbserr() = 0;

  IntegrationNd();
  virtual ~IntegrationNd();

private:
  virtual IntegrandNd wrap(IntegrandNd *f) = 0;

};

class ExpectationNd : public virtual IntegrationNd {
protected:
  virtual IntegrandNd wrap(IntegrandNd *f) = 0;

};

class NormalExpectationNd : public ExpectationNd {
public:
  NormalExpectationNd();
 ~NormalExpectationNd();

protected:
  virtual IntegrandNd wrap(IntegrandNd *f);

  friend class HcubatureNormalExpectation;
  friend class PcubatureNormalExpectation;
};

class CachedIntegrationNd : public virtual IntegrationNd {
protected:
  int errorCode_;
  unsigned nEval_;
  std::vector<double> err_;

public:
  int getNeval();
  int getErrorCode();
  std::vector<double>& getAbserr();

  CachedIntegrationNd();
  ~CachedIntegrationNd();

  friend class Cubature;
  friend class GaussianQuadratureNd;
};

typedef int (*cubature)(
    unsigned fdim, integrand f, void *fdata,
    unsigned ndim, const double *xmin, const double *xmax,
    size_t maxEval, double reqAbsError, double reqRelError,
    error_norm norm, double *val, double *err);

class Cubature : public CachedIntegrationNd {
private:
  const unsigned maxEval_;
  const double reqAbsError_, reqRelError_;
  const cubature cubature_;

public:
  std::vector<double> ninfInf(IntegrandNd *f);
  std::vector<double> aInf(IntegrandNd *f, double* bound);
  std::vector<double> ninfB(IntegrandNd *f, double* bound);
  std::vector<double> aB(IntegrandNd *f, double* a, double* b);

  Cubature(const unsigned maxEval, const double reqAbsError,
           const double reqRelError, const cubature cubature);
  ~Cubature();

private:
  virtual IntegrandNd wrap(IntegrandNd *f);
};

class Hcubature : public Cubature {
public:
  Hcubature(unsigned maxEval, double reqAbsError, double reqRelError);
  ~Hcubature();
};

class HcubatureNormalExpectation : public Hcubature, public NormalExpectationNd {
public:
  HcubatureNormalExpectation(unsigned maxEval, double reqAbsError, double reqRelError);
  ~HcubatureNormalExpectation();
private:
  virtual IntegrandNd wrap(IntegrandNd *f);
};

class Pcubature : public Cubature {
public:
  Pcubature(unsigned maxEval, double reqAbsError, double reqRelError);
  ~Pcubature();

  std::vector<double> ninfInf(IntegrandNd *f);
  std::vector<double> aInf(IntegrandNd *f, double* bound);
  std::vector<double> ninfB(IntegrandNd *f, double* bound);
};

class PcubatureNormalExpectation : public Pcubature, public NormalExpectationNd {
public:
  PcubatureNormalExpectation(unsigned maxEval, double reqAbsError, double reqRelError);
  ~PcubatureNormalExpectation();
private:
  virtual IntegrandNd wrap(IntegrandNd *f);
};

class Integrand2d : public IntegrandNd {
public:
  Integrand2d(IntegrFnEx* integrFnEx);
  ~Integrand2d();
};

class IntegrationNd2d : public Integration2d  {
private:
  IntegrationNd* integrationNd_;

public:

  double ninfInf(integr_fn *f, void *ex);
  double aInf(integr_fn *f, void *ex, double* bound);
  double ninfB(integr_fn *f, void *ex, double* bound);
  double aB(integr_fn *f, void *ex, double* a, double* b);

  int getNeval();
  int getErrorCode();
  double getAbserr();

  IntegrationNd2d(IntegrationNd* integrationNd);
  ~IntegrationNd2d();

private:
  IntegrFnEx wrap(integr_fn *f, void *ex);

};

class Hcubature2d : public IntegrationNd2d  {
public:
  Hcubature2d(unsigned maxEval, double reqAbsError, double reqRelError);
  ~Hcubature2d();
};

class Pcubature2d : public IntegrationNd2d  {
public:
  Pcubature2d(unsigned maxEval, double reqAbsError, double reqRelError);
  ~Pcubature2d();
};

int integrandNdRfun(unsigned ndim, const double *x, void *ex,
                     unsigned fdim, double *fval);

Rcpp::NumericVector test_Hcubature_ninfInf(Rcpp::Function func, int ndim, int fdim);
Rcpp::NumericVector test_Hcubature_ninfInf2(Rcpp::Function func, int ndim, int fdim, double tol);
Rcpp::NumericVector test_HcubatureNormalExpectation_ninfInf(Rcpp::Function func, int ndim, int fdim);

Rcpp::NumericVector test_Pcubature_ninfInf(Rcpp::Function func, int ndim, int fdim);
Rcpp::NumericVector test_Pcubature_aInf(Rcpp::Function func, int ndim, int fdim, Rcpp::NumericVector bound);
Rcpp::NumericVector test_Pcubature_ninfB(Rcpp::Function func, int ndim, int fdim, Rcpp::NumericVector bound);
Rcpp::NumericVector test_Pcubature_aB(Rcpp::Function func, int ndim, int fdim, Rcpp::NumericVector a, Rcpp::NumericVector b);

double test_Hcubature2d_ninfInf(Rcpp::Function func);
double test_Hcubature2d_aInf(Rcpp::Function func, Rcpp::NumericVector bound);
double test_Hcubature2d_ninfB(Rcpp::Function func, Rcpp::NumericVector bound);
double test_Hcubature2d_aB(Rcpp::Function func, Rcpp::NumericVector a, Rcpp::NumericVector b);

double test_Pcubature2d_ninfInf(Rcpp::Function func);
double test_Pcubature2d_aInf(Rcpp::Function func, Rcpp::NumericVector bound);
double test_Pcubature2d_ninfB(Rcpp::Function func, Rcpp::NumericVector bound);
double test_Pcubature2d_aB(Rcpp::Function func, Rcpp::NumericVector a, Rcpp::NumericVector b);

class GaussianQuadratureNd : public CachedIntegrationNd  {
private:
  std::vector<double> x_, w_;

public:
  std::vector<double> ninfInf(IntegrandNd *f);
  std::vector<double> aInf(IntegrandNd *f, double* bound);
  std::vector<double> ninfB(IntegrandNd *f, double* bound);
  std::vector<double> aB(IntegrandNd *f, double* a, double* b);

  GaussianQuadratureNd(int n);
  ~GaussianQuadratureNd();

protected:
  GaussianQuadratureNd(int n, void postInitFun(std::vector<double>& x, std::vector<double>& w));

private:
  void init(void postInitFun(std::vector<double>& x, std::vector<double>& w));
  virtual IntegrandNd wrap(IntegrandNd *f);
  int doNinfInf(IntegrandNd *f, double* x, double *fval, double carriedWgt,
                std::vector<double>& result, int dim);
};

Rcpp::NumericVector test_GaussianQuadratureNd_ninfInf(Rcpp::Function func, int ndim, int fdim, int nnodes);

class GaussianQuadratureNdNormalExpectation : public GaussianQuadratureNd, public NormalExpectationNd {
public:
  GaussianQuadratureNdNormalExpectation(int n);
  ~GaussianQuadratureNdNormalExpectation();

private:
  virtual IntegrandNd wrap(IntegrandNd *f);
};

Rcpp::NumericVector test_GaussianQuadratureNdNormalExpectation_ninfInf(Rcpp::Function func, int ndim, int fdim, int nnodes);

#endif
