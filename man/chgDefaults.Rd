\name{chgDefaults}
\alias{chgDefaults}
\alias{chgDefaults,psi_func_rcpp-method}
\title{Change default arguments}
\usage{
\S4method{chgDefaults}{psi_func_rcpp}(object, ...)
}
\arguments{
  \item{object}{PsiFunction instance to convert}
  \item{...}{arguments to change}
}
\description{
  Change the default arguments for a PsiFunction instance.
}
\note{
  Note that names of named arguments are ignored. Only the order of the arguments
  is considered when assigning new arguments.
}
\examples{
sPsi <- chgDefaults(smoothPsi, k=2)
curve(smoothPsi@psi(x), 0, 3)
curve(sPsi@psi(x), 0, 3, color="blue", add=TRUE)
}
\keyword{utilities}

