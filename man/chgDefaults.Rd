\name{chgDefaults}
\alias{chgDefaults}
\title{Change default arguments}
\arguments{
  \item{...}{arguments to change}
}
\description{
  Change the default arguments for a psi_func_cached object
}
\examples{
hPsi <- chgDefaults(huberPsi, k=2)
curve(huberPsi@psi(x), 0, 3)
curve(hPsi@psi(x), 0, 3, color="blue", add=TRUE)
}
\keyword{utilities}

