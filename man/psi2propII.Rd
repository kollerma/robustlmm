\name{psi2propII}
\alias{psi2propII}
\alias{psi2propII,psi_func_rcpp-method}
\title{Convert to Propsal II weight function}
\usage{
\S4method{psi2propII}{psi_func_rcpp}(object, ...)
}
\arguments{
  \item{object}{PsiFunction instance to convert}

  \item{...}{optional, new default arguments passed to
  chgDefaults.}
}
\description{
  Converts the PsiFunction instance into one that
  corresponds to Proposal II, i.e., a function of the
  squared weights. The other elements of the PsiFunction
  instance are adapted accordingly.
}
\examples{
par(mfrow=c(2,1))
plot(smoothPsi)
plot(psi2propII(smoothPsi))
}
\keyword{utilities}

