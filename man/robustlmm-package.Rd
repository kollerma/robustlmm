\name{robustlmm-package}
\alias{robustlmm}
\alias{robustlmm-package}
\docType{package}
\title{Robust linear mixed effects models}
\description{
  \code{robustlmm} provides functions for estimating linear mixed effects
  models in a robust way.

  The main workhorse is the function \code{\link{rlmer}}; it is implemented
  as direct robust analogue of the popular \code{\link{lmer}} function of
  the \code{\link{lme4}} package. The two functions have similar abilities
  and limitations. A wide range of data structures can be modeled: mixed
  effects models with hierarchical as well as complete or partially crossed
  random effects structures are possible. While the \code{\link{lmer}}
  function is optimized to handle large datasets efficiently, the
  computations employed in the \code{\link{rlmer}} function are more
  complex and for this reason also more expensive to compute. The two
  functions have the same limitations in the support of different random
  effect and residual error covariance structures. Both support only
  diagonal and unstructured random effect covariance structures.

  The \code{robustlmm} package implements most of the analysis tool chain
  as is customary in R. The usual functions such as \code{\link{summary}},
  \code{\link{coef}}, \code{\link{resid}}, etc. are provided as long as
  they are applicable for this type of models (see
  \code{\link{rlmerMod-class}} for a full list). The functions are designed
  to be as similar as possible to the ones in the \code{\link{lme4}}
  package to make switching between the two packages easy. 
  
  Details on the implementation and example analyses are provided in the
  package vignette available via \code{vignette("rlmer")} (Koller 2016).
 }
\references{
  Manuel Koller (2016).
  robustlmm: An R Package for Robust Estimation of Linear Mixed-Effects Models. 
  Journal of Statistical Software, 75(6), 1-24. doi:10.18637/jss.v075.i06
  
  Manuel Koller (2013).
  Robust estimation of linear mixed models.
  (Doctoral dissertation, Diss., Eidgenössische Technische Hochschule ETH Zürich, Nr. 20997, 2013).
}
 
