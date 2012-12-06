\name{print.xtable.comparison.table}
\alias{print.xtable.comparison.table}
\title{Print Export Tables for comparison tables

Wrapper for \code{print.xtable} for objects of class
\code{xtable.comparison.table}.}
\usage{
  \method{print}{xtable.comparison.table} (x,
    add.hlines = TRUE, latexify.namescol = TRUE,
    include.rownames = FALSE, ...)
}
\arguments{
  \item{x}{An object of class
  \code{xtable.comparison.table}.}

  \item{add.hlines}{replace empty lines in comparison table
  by hlines.  Supersedes \code{hline.after} argument of
  \code{print.xtable}.}

  \item{latexify.namescol}{replace \dQuote{sigma} and
  \dQuote{x} in the first column by latex equivalents.}

  \item{include.rownames}{include row numbers (names are
  included in the first column returned by
  \code{xtable.comparison.table}).}

  \item{...}{passed to \code{\link{print.xtable}}.}
}
\description{
  Print Export Tables for comparison tables

  Wrapper for \code{print.xtable} for objects of class
  \code{xtable.comparison.table}.
}
\seealso{
  \code{\link{print.xtable}},
  \code{\link{xtable.comparison.table}}, and
  \code{\link{compare}}.
}

