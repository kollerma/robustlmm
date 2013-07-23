\docType{data}
\name{smoothPsi}
\alias{smoothPsi}
\title{smoothPsi}
\format{Formal class 'psi_func_cached' [package "test"] with 11 slots}
\usage{
  smoothPsi
}
\description{
  Smoothed Huber Function
}
\details{
  Instead of a sharp bend like the huber function, this
  function bends smoothly. The first tuning contant, k, can
  be compared to the tuning constant of the original Huber
  function. The second tuning constant, s, determines the
  smoothness of the bend.
}
\examples{
par(mfrow=c(2,1))
plot(smoothPsi)
curve(smoothPsi@psi(x, 1.345, 10), -3, 3, col="red")
curve(huberPsi@psi(x, 1.345), -3, 3, add=TRUE)
}
\keyword{datasets}

