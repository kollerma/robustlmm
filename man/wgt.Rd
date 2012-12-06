\name{wgt.b}
\alias{wgt.b}
\alias{wgt.e}
\title{Calculate robustness weights}
\usage{
  wgt.b(object, sigma = object@pp$sigma, center = FALSE)

  wgt.e(object, sigma = object@pp$sigma,
    use.rho.sigma = FALSE)
}
\arguments{
  \item{object}{object to use}

  \item{sigma}{scale for standardization}

  \item{center}{whether return the centered robustness
  weights, see Details.}

  \item{use.rho.sigma}{return the weights computed using
  rho.sigma.}
}
\description{
  wgt.b: Calculate the robustness weights psi(d) / d,
  standardized by sigma. The robustness weights are
  calculated with d the Mahalanobis distance. Each group of
  correlated items then gets a constant weight.  The
  robustness weights for the random effects themselves are
  different than the ones used for estimating the size of
  the covariance matrix of the random effects. Those are
  additionally centered. That way, inlier can also be
  downweighted.  If \code{center=TRUE}, then the centered
  distances are used to compute the robustness weights and
  the weight function given by rho.sigma.b is used.

  wgt.e: robustness weights of residuals
}

