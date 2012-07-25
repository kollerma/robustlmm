\name{psiFunc2Cached}
\alias{psiFunc2Cached}
\title{psiFuncCached constructor}
\usage{
  psiFunc2Cached(rho, psi, wgt, Dpsi,
    name = NULL, ...)
}
\arguments{
  \item{rho}{rho-function}

  \item{psi}{psi-function}

  \item{wgt}{wgt-function}

  \item{Dpsi}{derivative of psi}

  \item{name}{descriptor of this function family}

  \item{...}{default values for tuning constants}
}
\value{
  psi_func2_cached-class object
}
\description{
  Create psi_func2_cached object using cached numerical
  integration for E... slots.
}
\details{
  The E... slots will not be fully functional: they just
  return the value for the current defaults and ignore
  their arguments.
}

