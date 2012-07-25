\name{wgt.b}
\alias{wgt.b}
\alias{wgt.e}
\title{Calculate robustness weights}
\usage{
  wgt.b(object, sigma = object@pp$sigma,
    use.c.sigma = FALSE)

  wgt.e(object, sigma = object@pp$sigma,
    use.c.sigma = FALSE)
}
\arguments{
  \item{object}{object to use}

  \item{sigma}{scale for standardization}

  \item{use.c.sigma}{whether use c.sigma tuning parameters
  or the default ones.}
}
\description{
  wgt.b: Calculate the robustness weights psi(d) / d,
  standardized by sigma. The robustness weights are
  calculated with d the Mahalanobis distance. Each group of
  correlated items then gets a constant weight.

  wgt.e: robustness weights of residuals
}

