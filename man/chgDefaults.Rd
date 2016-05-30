\name{chgDefaults}
\alias{chgDefaults}
\alias{chgDefaults,Rcpp_SmoothPsi-method}
\alias{chgDefaults,Rcpp_HuberPsi-method}
\alias{chgDefaults,Rcpp_PsiFunction-method}
\alias{chgDefaults,Rcpp_PsiFunction to Prop II PsiFunction wrapper-method}
\title{Change default arguments}
\usage{
\S4method{chgDefaults}{Rcpp_SmoothPsi}(object, ...)
\S4method{chgDefaults}{Rcpp_HuberPsi}(object, ...)
\S4method{chgDefaults}{Rcpp_PsiFunction}(object, ...)
\S4method{chgDefaults}{`Rcpp_PsiFunction to Prop II PsiFunction wrapper`}(object, ...)
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
curve(smoothPsi$psi(x), 0, 3)
curve(sPsi$psi(x), 0, 3, color="blue", add=TRUE)
}
\keyword{utilities}

