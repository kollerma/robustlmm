\name{rlmer}
\alias{rlmer}
\title{Robust linear mixed models}
\usage{
  rlmer(formula, data, REML = TRUE, ..., method = "DAStau",
    method.effects = "IRWLS", rho.e = smoothPsi,
    rho.b = smoothPsi, c.sigma.e = NULL, c.sigma.b = NULL,
    wExp.e = 2, wExp.b = 2,
    use.laplace = REML || method \%in\% c("Opt"),
    rel.tol = 1e-08, max.iter = 40 * (r + 1)^2,
    verbose = 0,
    use.grad = c(effects = TRUE, varcomp = FALSE),
    doFit = TRUE, clearCache = TRUE, init)
}
\arguments{
  \item{formula}{a two-sided linear formula object
  describing the fixed-effects part of the model, with the
  response on the left of a \code{~} operator and the
  terms, separated by \code{+} operators, on the right.
  The vertical bar character \code{"|"} separates an
  expression for a model matrix and a grouping factor.}

  \item{data}{an optional data frame containing the
  variables named in \code{formula}.  By default the
  variables are taken from the environment from which
  \code{lmer} is called.}

  \item{REML}{logical - Should the estimates be chosen to
  optimize the REML criterion (as opposed to the
  log-likelihood)?  Defaults to \code{TRUE}.}

  \item{...}{Additional parameters passed to lmer to find
  the initial estimates. See \code{\link[lme4]{lmer}}.}

  \item{method}{method to be used for estimation of theta
  and sigma, see Details.}

  \item{method.effects}{method to be used for estimation of
  the fixed and random effects.}

  \item{rho.e}{object of class psi_func, specifying the
  functions to use for robustifying the residuals part of
  the likelihood.}

  \item{rho.b}{object of class psi_func, specifying the
  functions to use for robustifying the random effects part
  of the likelihood.}

  \item{c.sigma.e}{tuning parameters for rho.e when used
  for estimating residual scale.}

  \item{c.sigma.b}{tuning parameters for rho.b when used
  for estimating variance components.}

  \item{wExp.e}{numeric specifying the exponent used for
  the robustness weights for the residuals. wExp == 0 means
  that the standard rho function is used for estimating
  theta and sigma, otherwise the weights are used with the
  specified exponent. (wExp == 2 corresponds to Proposal
  II.)}

  \item{wExp.b}{numeric specifying the exponent used for
  the robustness weights for the random effects. wExp == 0
  means that the standard rho function is used for
  estimating theta and sigma, otherwise the weights are
  used with the specified exponent. (wExp == 2 corresponds
  to Proposal II.)}

  \item{use.laplace}{whether to use the the Laplace
  approximation for calculating the deviance.}

  \item{rel.tol}{relative tolerance used as criteria in the
  fitting process.}

  \item{max.iter}{maximum number of iterations allowed.}

  \item{verbose}{verbosity of output. Ranges from 0 (none)
  to 3 (a lot of output)}

  \item{use.grad}{list of two elements, specifying whether
  to use the analytical gradients in the optimization.}

  \item{doFit}{logical scalar. When \code{doFit = FALSE}
  the model is not fit but instead a structure with the
  model matrices for the random-effects terms is returned,
  so they can be modified for special model forms. When
  \code{doFit = TRUE}, the default, the model is fit
  immediately.}

  \item{clearCache}{clear the cache before returning the
  fitted object?}

  \item{init}{optional lmerMod- or rlmerMod-object to use
  for starting values, or function producing an lmerMod
  object.}
}
\value{
  object of class rlmerMod.
}
\description{
  The method is based on the robustification of the scoring
  equations and an application of the Design Adaptive Scale
  approach.
}
\details{
  Available fitting methods for theta and sigma.e:
  \itemize{ \item DASexp: find theta using DAS-averaged
  approach, for each canditate theta, estimate beta and b
  and estimate sigma_e using DAS-averaged approach. \item
  DAStau: Analogue to the DAS-estimate in robust linear
  regression. \item Opt: optimize for (log(sigma), theta),
  for given theta and sigma_e, estimate beta and b }

  Available fitting methods for fixed and random effects:
  \itemize{ \item IRWLS: Perform iterative reweighted least
  squares using the estimating equations (eeq) \item
  Rcgmin: minimize the norm part of the logarithm of the
  density using the Rcgmin optimizer. }
}
\examples{
## dropping of VC
  system.time(rlmer(Yield ~ (1|Batch), Dyestuff2))

  ## compare classic fits of various methods,
  ## without using lmer fit as initial estimate.
  ## default method
  system.time(rfm.DAStau <- rlmer(Yield ~ (1|Batch), Dyestuff, method="DAStau",
                                  rho.e = cPsi, rho.b = cPsi, init=lmerNoFit))
  ## DASr method
  system.time(rfm.DASexp <- rlmer(Yield ~ (1|Batch), Dyestuff, method="DASexp",
                                  rho.e = cPsi, rho.b = cPsi, init=lmerNoFit))

  lf.comp <- function(rfm1, rfm2, tolerance = 1e-5) {
      ## fix calls and other slots
      rfm2@call <- rfm1@call
      rfm2@method <- rfm1@method
      all.equal(rfm1, rfm2, tolerance = tolerance)
  }

  stopifnot(lf.comp(rfm.DAStau, rfm.DASexp))
}
\seealso{
  \code{\link[lme4]{lmer}}
}

