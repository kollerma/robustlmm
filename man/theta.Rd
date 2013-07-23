\name{theta}
\alias{theta}
\title{Get theta}
\usage{
  theta(object)
}
\arguments{
  \item{object}{merMod object}
}
\description{
  Get theta
}
\examples{
fm <- rlmer(Yield ~ (1|Batch), Dyestuff)
stopifnot(all.equal(theta(fm), getME(fm, "theta")))
}

