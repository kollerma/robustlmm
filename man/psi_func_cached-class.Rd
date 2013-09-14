\name{psi_func_cached-class}
\Rdversion{1.1}
\docType{class}
\alias{psi_func_cached-class}
\alias{chgDefaults,psi_func_cached-method}
\alias{psi2propII,psi_func_cached-method}

\title{Class \code{"psi_func_cached"}}
\description{
  The class \code{"psi_func_cached"} is used to store \eqn{\psi \ (psi)}{psi}
  functions for M-estimation.  In particular, an object of the class
  contains \eqn{\rho(x) \ (\code{rho})}{rho(x)}, its derivative
  \eqn{\psi(x) \ (psi)}{psi(x)}, the weight function \eqn{\psi(x)/x}, and
  first derivative of \eqn{\psi}, \code{Dpsi = } \eqn{\psi'(x)}.

  The expected values are calculated using numerical integration. For
  efficiency reasons, the results are cached. The current implementation
  only returns the expectations for the default values and ignores the
  arguments.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("psi_func_cached", ...)}
but preferably by \code{\link{psiFuncCached}(...)}.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"}, printable
      name of the psi function. }
    \item{\code{rho}:}{Object of class \code{"functionX"}, see \code{"\linkS4class{psi_func}"} }
    \item{\code{psi}:}{Object of class \code{"functionX"}, see \code{"\linkS4class{psi_func}"} }
    \item{\code{wgt}:}{Object of class \code{"functionX"}, see \code{"\linkS4class{psi_func}"} }
    \item{\code{Dwgt}:}{Object of class \code{"functionX"}, derivative of wgt}
    \item{\code{Dpsi}:}{Object of class \code{"functionX"}, see \code{"\linkS4class{psi_func}"} }
    \item{\code{tDefs}:}{Object of class \code{"numeric"}, see \code{"\linkS4class{psi_func}"} }
    \item{\code{Erho}:}{Object of class \code{"functionXal"}, see \code{"\linkS4class{psi_func}"} }
    \item{\code{Epsi2}:}{Object of class \code{"functionXal"}, see \code{"\linkS4class{psi_func}"} }
    \item{\code{EDpsi}:}{Object of class \code{"functionXal"}, see \code{"\linkS4class{psi_func}"} }
    \item{\code{name}:}{Object of class \code{"character"}, see \code{"\linkS4class{psi_func}"} }
    \item{\code{xtras}:}{Object of class \code{"list"}, see \code{"\linkS4class{psi_func}"} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{psi_func}"}, directly.
}
\section{Methods}{
  \describe{
    \item{chgDefaults}{\code{signature(object = "psi_func_cached")}:
      The method is used to change the default
      values for the tuning parameters, and returns a
      \code{\link{psi_func-class}} object, a copy of input \code{object}
      with the slot \code{tDefs} possibly changed. }
    \item{plot}{\code{signature(object = "psi_func_cached")}:
      Creates a plot of all the function slots.
    }
    \item{psi2propII}{\code{signature(object = "psi_func_cached")}:
      Converts the psi_func object into a function that
      corresponds to Proposal II, i.e., a function of the
      squared weights. The other elements of the psi_func
      object are adapted accordingly.
    }
   \item{show}{\code{signature(object = "psi_func_cached")}:
     Print function.
    }
  }
}
\author{
  Manuel Koller
}
\section{Warning}{
  The E... slots will not be fully functional: they just
  return the value for the current defaults and ignore
  their arguments.
}

\seealso{
  Inherits from \code{\linkS4class{psi_func}},
  constructor \code{\link{psiFuncCached}}.
}
\examples{
showClass("psi_func_cached")
}
\keyword{classes}
