\name{compare}
\alias{compare}
\alias{getInfo}
\alias{getInfo.lmerMod}
\alias{getInfo.rlmerMod}
\alias{print.xtable.comparison.table}
\alias{xtable.comparison.table}
\title{Create comparison charts for multiple fits}
\usage{
  compare(..., digits = 3, dnames = NULL,
    show.rho.functions = TRUE)

  \method{xtable}{comparison.table} (x, caption = NULL,
    label = NULL, align = NULL, digits = NULL,
    display = NULL, ...)

  \method{print}{xtable.comparison.table} (x,
    add.hlines = TRUE, latexify.namescol = TRUE,
    include.rownames = FALSE, ...)

  getInfo(object, ...)

  \method{getInfo}{lmerMod} (object, ...)

  \method{getInfo}{rlmerMod} (object, ...)
}
\arguments{
  \item{...}{objects to compare, or, for the
  \code{\link{xtable}} functions: passed to the respective
  \code{\link{xtable}} function.}

  \item{digits}{number of digits to show in output}

  \item{dnames}{names of objects given as arguments
  (optional)}

  \item{show.rho.functions}{whether to show rho functions
  in output.}

  \item{x}{object of class "comparison.table" or
  "xtable.comparison.table"}

  \item{caption}{see \code{\link{xtable}}.}

  \item{label}{see \code{\link{xtable}}.}

  \item{align}{see \code{\link{xtable}}.}

  \item{display}{see \code{\link{xtable}}.}

  \item{add.hlines}{replace empty lines in comparison table
  by hlines.  Supersedes \code{hline.after} argument of
  \code{print.xtable}.}

  \item{latexify.namescol}{replace \dQuote{sigma} and
  \dQuote{x} in the first column by latex equivalents.}

  \item{include.rownames}{include row numbers (the object
  returned by \code{xtable.comparison.table} includes names
  in the first column)}

  \item{object}{object}
}
\value{
  \code{getInfo} returns alist with estimated coefficients,
  estimated variance components, sigma, deviance and
  parameter configuration used to fit.
}
\description{
  Use \code{compare} to quickly compare the estaimated
  parameters of the fits of multiple lmerMod or rlmerMod
  objects.

  The functions \code{xtable.comparison.table} and
  \code{print.xtable.comparison.table} are wrapper
  functions for the respective \code{\link{xtable}} and
  \code{\link{print.xtable}} functions.

  The function \code{getInfo} is internally used to prepare
  object for producing a comparison chart in
  \code{compare}.
}
\examples{
fm1 <- lmer(Yield ~ (1|Batch), Dyestuff)
fm2 <- rlmer(Yield ~ (1|Batch), Dyestuff)
compare(fm1, fm2)
require(xtable)
xtable(compare(fm1, fm2))
str(getInfo(fm1))
}
\seealso{
  \code{\link{xtable}}

  \code{\link{print.xtable}}
}
\keyword{models}
\keyword{utilities}

