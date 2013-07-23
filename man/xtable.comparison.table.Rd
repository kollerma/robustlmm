\name{xtable.comparison.table}
\alias{xtable.comparison.table}
\title{Create Latex and HTML tables for output of \code{\link{compare}()}.

This function is a wrapper to \code{table} for objects of class
\code{comparison.table}. It uses the same arguments as \code{xtable()}.}
\usage{
  \method{xtable}{comparison.table} (x, caption = NULL,
    label = NULL, align = NULL, digits = NULL,
    display = NULL, ...)
}
\arguments{
  \item{x}{object of class "comparison.table".}

  \item{caption}{see \code{\link{xtable}}.}

  \item{label}{see \code{\link{xtable}}.}

  \item{align}{see \code{\link{xtable}}.}

  \item{digits}{see \code{\link{xtable}}.}

  \item{display}{see \code{\link{xtable}}.}

  \item{...}{passed to \code{\link{xtable}}.}
}
\description{
  Create Latex and HTML tables for output of
  \code{\link{compare}()}.

  This function is a wrapper to \code{table} for objects of
  class \code{comparison.table}. It uses the same arguments
  as \code{xtable()}.
}
\examples{
fm1 <- lmer(Yield ~ (1|Batch), Dyestuff)
fm2 <- rlmer(Yield ~ (1|Batch), Dyestuff)
require(xtable)
xtable(compare(fm1, fm2))
}
\seealso{
  \code{\link{xtable}},
  \code{\link{print.xtable.comparison.table}} and
  \code{\link{compare}}.
}

