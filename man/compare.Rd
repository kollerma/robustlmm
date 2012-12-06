\name{compare}
\alias{compare}
\title{Create a comparison chart for multiple fits}
\usage{
  compare(..., digits = 3, dnames = NULL,
    show.rho.functions = TRUE)
}
\arguments{
  \item{...}{objects to compare}

  \item{digits}{number of digits to show in print}

  \item{dnames}{names of objects given as arguments
  (optional)}

  \item{show.rho.functions}{whether to show rho functions
  in output.}
}
\value{
  comparison table
}
\description{
  Compare the fits of multiple lmerMod or rlmerMod objects
}

