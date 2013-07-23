\name{lmerNoFit}
\alias{lmerNoFit}
\title{Construct lmerMod object but do not fit the data.}
\usage{
  lmerNoFit(..., initTheta)
}
\arguments{
  \item{...}{passed to lmer()}

  \item{initTheta}{parameter to initialize theta with
  (optional)}
}
\value{
  lmerMod object
}
\description{
  This functions replaces the deprecated option 
  \code{doFit = FALSE} in lme4.
}
\examples{
rlmer(Yield ~ (1|Batch), Dyestuff, init = lmerNoFit)
}

