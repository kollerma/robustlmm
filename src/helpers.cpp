#include "helpers.h"

using namespace Rcpp;

#define DEBUG(STRING)                                                  \
//  Rcpp::Rcout << STRING << std::endl;

void setFixef(rlmerPredD *pp, rlmerResp *resp, const VectorXd& newFixef) {
  R_ASSERT(pp->p_ == newFixef.size());
  pp->setBeta(newFixef);
  resp->updateMu(pp->getMu());
}

void setEffects(rlmerPredD *pp, rlmerResp *resp, const VectorXd& newEffects) {
  R_ASSERT(pp->p_ + pp->q_ == newEffects.size());
  pp->setB_s(newEffects.segment(pp->p_, pp->q_));
  setFixef(pp, resp, newEffects.segment(0, pp->p_));
}

VectorXd wgt_e(const rlmerPredD *pp, const rlmerResp *resp) {
  DEBUG("Called wgt_e")
  VectorXd result(dist_e(pp, resp));
  DEBUG(" Entering for loop")
  for (int i = 0; i < result.rows(); ++i) {
    result(i) = pp->getRho_e()->wgtFun(result(i));
  }
  return result;
}

VectorXd wgt_b(const rlmerPredD *pp) {
  DEBUG("Called wgt_b")
  VectorXd result(pp->getDist_b());
  DEBUG(" Entering for loop")
  for (int i = 0; i < result.rows(); ++i) {
    result(i) = pp->getRho_b()[pp->getBBlockMap()(i) - 1]->wgtFun(result(i));
  }
  return result;
}

VectorXd Rwgt_e(const rlmerPredDXPtr pp, const rlmerRespXPtr resp) {
  return wgt_e(static_cast<rlmerPredD*> (R_ExternalPtrAddr(pp)),
               static_cast<rlmerResp*> (R_ExternalPtrAddr(resp)));
}

VectorXd Rwgt_b(const rlmerPredDXPtr pp) {
   return wgt_b(static_cast<rlmerPredD*> (R_ExternalPtrAddr(pp)));
}

VectorXd dist_e(const rlmerPredD *pp, const rlmerResp *resp) {
  return resp->wtres_ / pp->getSigma();
}
