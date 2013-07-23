\name{sigma}
\alias{sigma}
\title{Get sigma}
\usage{
  sigma(object, ...)
}
\arguments{
  \item{object}{rlmerMod object}

  \item{...}{ignored}
}
\description{
  Simple accessor for estimated residual scale sigma.
}
\examples{
fm <- rlmer(Yield ~ (1|Batch), Dyestuff)
sigma(fm)
getME(fm, "devcomp")$cmp[c("sigmaML", "sigmaREML")]
}
\seealso{
  \code{\link{getME}}
}

