#if !defined  ROBUSTLMM_RLMERRESPMODULE_H__
#define  ROBUSTLMM_RLMERRESPMODULE_H__

#include "globals.h"

using namespace Rcpp;

class rlmerResp {
public:
  MVec y_, offset_, mu_, weights_, sqrtrwt_, wtres_;

  rlmerResp(MVec y, MVec weights, MVec offset, MVec mu, MVec sqrtrwt, MVec wtres);
  rlmerResp(Rcpp::XPtr<rlmerResp> toCopy, bool dummy);
  void updateMu(VectorXd mu);
};

extern "C" SEXP _rcpp_module_boot_rlmerResp_module();

#endif
