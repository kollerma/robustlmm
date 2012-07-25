\name{gradient}
\alias{gradient}
\alias{gradient.rlmerMod}
\title{Gradient}
\usage{
  gradient(object, ...)

  \method{gradient}{rlmerMod} (object, ...)
}
\arguments{
  \item{object}{a lmerMod object}

  \item{...}{ignored}
}
\value{
  vector
}
\description{
  Calculates the gradient, the gradient of beta and b if a
  lmerMod object is given, the gradient of theta if a
  reModule is given.
}

