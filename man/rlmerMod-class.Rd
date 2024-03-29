% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClass.R
\docType{class}
\name{rlmerMod-class}
\alias{rlmerMod-class}
\alias{coef.rlmerMod}
\alias{deviance.rlmerMod}
\alias{extractAIC.rlmerMod}
\alias{family.rlmerMod}
\alias{fitted.rlmerMod}
\alias{fixef.rlmerMod}
\alias{formula.rlmerMod}
\alias{isGLMM.rlmerMod}
\alias{isLMM.rlmerMod}
\alias{isNLMM.rlmerMod}
\alias{isREML.rlmerMod}
\alias{logLik.rlmerMod}
\alias{model.frame.rlmerMod}
\alias{model.matrix.rlmerMod}
\alias{nobs.rlmerMod}
\alias{predict.rlmerMod}
\alias{print.rlmerMod}
\alias{print.summary.rlmer}
\alias{print.VarCorr.rlmerMod}
\alias{ranef.rlmerMod}
\alias{resid.rlmerMod}
\alias{sigma.rlmerMod}
\alias{show.rlmerMod}
\alias{show,rlmerMod-method}
\alias{show.summary.rlmerMod}
\alias{summary.rlmerMod}
\alias{summary.summary.rlmerMod}
\alias{terms.rlmerMod}
\alias{update.rlmerMod}
\alias{VarCorr.rlmerMod}
\alias{VarCorr.summary.rlmerMod}
\alias{vcov.rlmerMod}
\alias{vcov.summary.rlmerMod}
\alias{weights.rlmerMod}
\title{rlmerMod Class}
\description{
Class "rlmerMod" of Robustly Fitted Mixed-Effect Models
}
\details{
A robust mixed-effects model as returned by \code{\link{rlmer}}.
}
\section{Objects from the Class}{
 Objects are created by calls to
  \code{\link{rlmer}}.
}

\section{Methods}{
 Almost all methods available from objects returned from
  \code{\link{lmer}} are also available for objects returned by
  \code{\link{rlmer}}. They usage is the same.

  It follows a list of some the methods that are exported by this package:

\itemize{
\item \code{\link{coef}}
\item \code{\link{deviance}} (disabled, see below)
\item \code{\link{extractAIC}} (disabled, see below)
\item \code{\link{family}}
\item \code{\link{fitted}}
\item \code{\link[=fixef.merMod]{fixef}}
\item \code{\link{formula}}
\item \code{\link{getInfo}}
\item \code{\link{isGLMM}}
\item \code{\link{isLMM}}
\item \code{\link{isNLMM}}
\item \code{\link{isREML}}
\item \code{\link{logLik}} (disabled, see below)
\item \code{\link{model.frame}}
\item \code{\link{model.matrix}}
\item \code{\link{nobs}}
\item \code{\link[=plot.rlmerMod]{plot}}
\item \code{\link[=predict.merMod]{predict}}
\item \code{\link[=ranef.merMod]{ranef}} (only partially implemented)
\item \code{\link[=residuals.rlmerMod]{residuals}}
\item \code{\link{sigma}}
\item \code{\link{summary}}
\item \code{\link{terms}}
\item \code{\link{update}}
\item \code{\link[=VarCorr.merMod]{VarCorr}}
\item \code{\link{vcov}}
\item \code{\link{weights}}
}
}

\section{Disabled methods}{
 A log likelihood or even a pseudo log likelihood
  is not defined for the robust estimates returned by \code{\link{rlmer}}.
  Methods that depend on the log likelihood are therefore not available. For
  this reason the methods \code{deviance}, \code{extractAIC} and
  \code{logLik} stop with an error if they are called.
}

\examples{

showClass("rlmerMod")

## convert an object of type 'lmerMod' to 'rlmerMod'
## to use the methods provided by robustlmm
fm <- lmer(Yield ~ (1|Batch), Dyestuff)
rfm <- as(fm, "rlmerMod")
compare(fm, rfm)

}
\seealso{
\code{\link{rlmer}}; corresponding class in package \code{lme4}:
  \code{\link{merMod}}
}
\keyword{classes}
