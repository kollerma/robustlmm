\name{getInfo}
\alias{getInfo}
\alias{getInfo.lmerMod}
\alias{getInfo.mer}
\alias{getInfo.rlmerMod}
\title{getInfo of an object}
\usage{
  getInfo(object, ...)

  \method{getInfo}{lmerMod} (object, ...)

  \method{getInfo}{mer} (object, ...)

  \method{getInfo}{rlmerMod} (object, ...)
}
\arguments{
  \item{object}{object}

  \item{...}{ignored}
}
\value{
  list with estimated coefficients, estimated variance
  components, sigma, deviance and parameter configuration
  used to fit.
}
\description{
  Extract some information from objects returned from
  \code{\link{rlmer}} and \code{\link{lmer}} in the form of
  a simple list. Internally used to prepare object for
  producing a comparison chart in \code{\link{compare}}.
}
\examples{
fm <- rlmer(Yield ~ (1|Batch), Dyestuff)
str(getInfo(fm))
}

