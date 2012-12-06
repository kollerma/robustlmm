\name{psiFuncCached}
\alias{psiFuncCached}
\title{psiFuncCached constructor}
\usage{
  psiFuncCached(rho, psi, wgt, Dwgt, Dpsi, name = NULL,
    ...)
}
\arguments{
  \item{rho}{rho-function}

  \item{psi}{psi-function}

  \item{wgt}{wgt-function}

  \item{Dwgt}{derivative of weight function}

  \item{Dpsi}{derivative of psi}

  \item{name}{descriptor of this function family}

  \item{...}{default values for tuning constants}
}
\value{
  psi_func_cached-class object
}
\description{
  Create psi_func_cached object using cached numerical
  integration for E... slots.
}
\details{
  The E... slots will not be fully functional: they just
  return the value for the current defaults and ignore
  their arguments.
}

