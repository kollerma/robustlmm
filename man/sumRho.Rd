\name{sumRho.b}
\alias{sumRho.b}
\alias{sumRho.e}
\title{sumRho}
\usage{
  sumRho.b(object, sigma = object@pp$sigma, lambda = FALSE,
    ...)

  sumRho.e(object, sigma = object@pp$sigma, lambda = FALSE,
    ...)
}
\arguments{
  \item{object}{rlmerMod object}

  \item{sigma}{to use for standardization}

  \item{lambda}{apply normalization by lambda}

  \item{...}{ignored}
}
\description{
  sumRho.b: Calculate sum of rho of random effects

  sumRho.e: Calculate sum of rho of residuals
}

