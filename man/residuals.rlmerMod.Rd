\name{residuals.rlmerMod}
\alias{residuals.rlmerMod}
\title{Get residuals}
\usage{
  \method{residuals}{rlmerMod} (object,
    type = c("response", "weighted"), scaled = FALSE, ...)
}
\arguments{
  \item{object}{rlmerMod object}

  \item{type}{type of residuals}

  \item{scaled}{scale residuals by residual standard
  deviation (=scale parameter)?}

  \item{...}{ignored}
}
\description{
  The per-observation residuals are returned, i.e., the
  difference of the observation and the fitted value
  including random effects. With type one can specify
  whether the weights should be used or not.
}
\examples{
\dontrun{
  fm <- rlmer(Yield ~ (1|Batch), Dyestuff)
  stopifnot(all.equal(resid(fm, type="weighted"),
                      resid(fm) * getME(fm, "w_e")))
}
}

