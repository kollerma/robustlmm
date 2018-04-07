#if !defined  ROBUSTLMM_HELPERS_H__
#define  ROBUSTLMM_HELPERS_H__

#include "globals.h"

class rlmerPredD;
class rlmerResp;

typedef Rcpp::XPtr<rlmerPredD>                   rlmerPredDXPtr;
typedef Rcpp::XPtr<rlmerResp>                    rlmerRespXPtr;

void setFixef(rlmerPredD *pp, rlmerResp *resp, const VectorXd& newFixef);
void setEffects(rlmerPredD *pp, rlmerResp *resp, const VectorXd& newEffects);

VectorXd wgt_e(const rlmerPredD *pp, const rlmerResp *resp);
VectorXd wgt_b(const rlmerPredD *pp);

VectorXd Rwgt_e(const rlmerPredDXPtr pp, const rlmerRespXPtr resp);
VectorXd Rwgt_b(const rlmerPredDXPtr pp);

VectorXd dist_e(const rlmerPredD *pp, const rlmerResp *resp);

#endif
