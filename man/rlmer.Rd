\name{rlmer}
\alias{rlmer}
\title{Robust linear mixed models}
\usage{
  rlmer(formula, data, ..., method = "DAStau",
    rho.e = smoothPsi, rho.b = smoothPsi, rho.sigma.e,
    rho.sigma.b, rel.tol = 1e-08,
    max.iter = 40 * (r + 1)^2, verbose = 0, doFit = TRUE,
    init)
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

  \item{...}{Additional parameters passed to lmer to find
  the initial estimates. See \code{\link[lme4]{lmer}}.}

  \item{method}{method to be used for estimation of theta
  and sigma, see Details.}

  \item{rho.e}{object of class psi_func, specifying the
  functions to use for the huberization of the residuals.}

  \item{rho.b}{object of class psi_func or list of such
  objects (see Details), specifying the functions to use
  for the huberization of the random effects.}

  \item{rho.sigma.e}{object of class psi_func, specifying
  the functions to use for the huberization of the
  residuals when estimating the variance components.}

  \item{rho.sigma.b}{(optional) object of class psi_func or
  list of such objects, specifying the functions to use for
  the huberization of the random effects when estimating
  the variance components.}

  \item{rel.tol}{relative tolerance used as criteria in the
  fitting process.}

  \item{max.iter}{maximum number of iterations allowed.}

  \item{verbose}{verbosity of output. Ranges from 0 (none)
  to 3 (a lot of output)}

  \item{doFit}{logical scalar. When \code{doFit = FALSE}
  the model is not fit but instead a structure with the
  model matrices for the random-effects terms is returned,
  so they can be modified for special model forms. When
  \code{doFit = TRUE}, the default, the model is fit
  immediately.}

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

  Example analyses and theoretical details on the method
  are available in the vignette (\code{vignette("rlmer")}).
}
\details{
  Available fitting methods for theta and sigma.e:
  \itemize{ \item DASvar: Instead of calculating the
  expectation using the linear approximation of the
  residuals, use it to approximate the variance of the
  residuals / random effects and correct the scale estimate
  with it. \item DAStau: Analogue to the DAS-estimate in
  robust linear regression. }

  To specify separate weight functions \code{rho.b} and
  \code{rho.sigma.b} for different variance components, it
  is possible to pass a list instead of a psi_func object.
  The list entries correspond to the groups as shown by
  \code{VarCorr(.)} when applied to the model fitted with
  \code{lmer}. A set of correlated random effects count as
  just one group.
}
\examples{
## dropping of VC
  system.time(rlmer(Yield ~ (1|Batch), Dyestuff2))

  ## Default method "DAStau"
  system.time(rfm.DAStau <- rlmer(Yield ~ (1|Batch), Dyestuff))
  summary(rfm.DAStau)

  ## DASvar method
  system.time(rfm.DASvar <- rlmer(Yield ~ (1|Batch), Dyestuff, method="DASvar"))

  compare(rfm.DAStau, rfm.DAStau)
}
\seealso{
  \code{\link[lme4]{lmer}}, \code{vignette("rlmer")}
}

