#include "rlmerRespModule.h"

using namespace Rcpp;

rlmerResp::rlmerResp(MVec y, MVec weights, MVec offset, MVec mu, MVec sqrtrwt, MVec wtres) :
  y_(y), offset_(offset), mu_(mu), weights_(weights), sqrtrwt_(sqrtrwt), wtres_(wtres) {
  // FIXME add dimensions check
}

rlmerResp::rlmerResp(Rcpp::XPtr<rlmerResp> toCopy, bool dummy) :
  y_(toCopy->y_), offset_(toCopy->offset_),
  mu_(toCopy->mu_), weights_(toCopy->weights_),
  sqrtrwt_(toCopy->sqrtrwt_), wtres_(toCopy->wtres_)
  {}

void rlmerResp::updateMu(VectorXd mu) {
  mu_ = mu;
  wtres_ = sqrtrwt_.cwiseProduct(y_ - mu_);
}

RCPP_EXPOSED_CLASS(rlmerResp)
RCPP_MODULE(rlmerResp_module) {

  class_<rlmerResp>("rlmerResp")
  /* y, weights, offset, mu, sqrtrwt, wtres */
  .constructor<MVec, MVec, MVec, MVec, MVec, MVec>()
  .constructor<Rcpp::XPtr<rlmerResp>, bool >()
  .method("updateMu", &rlmerResp::updateMu, "update mu and weighted residuals")
  .field("y", &rlmerResp::y_)
  .field("weights", &rlmerResp::weights_)
  .field("offset", &rlmerResp::offset_)
  .field("sqrtrwt", &rlmerResp::sqrtrwt_)
  .field("wtres", &rlmerResp::wtres_)
  .field("mu", &rlmerResp::mu_)
  ;

}
