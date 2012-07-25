\name{getInfo}
\alias{getInfo}
\alias{getInfo.lmerMod}
\alias{getInfo.rlmerMod}
\title{getInfo of an object}
\usage{
  getInfo(object, ...)

  \method{getInfo}{lmerMod} (object, ...)

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
  Prepare object for producing a comparison chart.
}

