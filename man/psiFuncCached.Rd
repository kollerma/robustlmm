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
  Warning: the E... slots will not be fully functional:
  they just return the value for the current defaults and
  ignore their arguments.
}
\examples{
## re-define cPsi as psiFuncCached.
F0 <- function(x=1, .) rep.int(0, length(x))
F1 <- function(x=1, .) rep.int(1, length(x))
cPsi2 <- psiFuncCached(rho = function(x, .) x^2 / 2,
                       psi = function(x, .) x,
                       wgt = F1, Dwgt = F0, Dpsi = F1,
                       name = "classic (x^2/2)",
                       . = Inf ## dummy, need at least one parameter
                       )
stopifnot(all.equal(cPsi@Erho(), cPsi2@Erho()),
          all.equal(cPsi@Epsi2(), cPsi2@Epsi2()),
          all.equal(cPsi@EDpsi(), cPsi2@EDpsi()))
}
\seealso{
  \code{\link{psi_func_cached-class}}
}

