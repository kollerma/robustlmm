\docType{methods}
\name{fixef}
\alias{fixed.effects}
\alias{fixef}
\alias{fixef.rlmerMod}
\title{Extract fixed-effects estimates}
\usage{
  fixef(object, ...)
}
\arguments{
  \item{object}{any fitted model object from which fixed
  effects estimates can be extracted.}

  \item{\dots}{optional additional arguments. Currently
  none are used in any methods.}
}
\value{
  a named, numeric vector of fixed-effects estimates.
}
\description{
  Extract the fixed-effects estimates
}
\details{
  Extract the estimates of the fixed-effects parameters
  from a fitted model.
}
\examples{
## doFit=FALSE to speed up example
fixef(rlmer(Reaction ~ Days + (Days|Subject), sleepstudy, doFit=FALSE))
}
\keyword{models}

