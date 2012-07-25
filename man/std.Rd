\name{std.b}
\alias{std.b}
\alias{std.e}
\title{Standardized values}
\usage{
  std.b(object, sigma = object@pp$sigma, matrix,
    drop = TRUE, t = FALSE)

  std.e(object, sigma = object@pp$sigma, matrix,
    drop = TRUE)
}
\arguments{
  \item{object}{rlmerMod object}

  \item{sigma}{to use for standardization}

  \item{matrix}{matrix to standardize}

  \item{drop}{apply drop to result?}

  \item{t}{transpose result}
}
\description{
  std.b: Return the spherical random effects or
  "Standardize" the Matrix matrix: \eqn{\Lambda^{-1} matrix
  / \sigma}{Lambda^-1 matrix / sigma}

  std.e: Calculate the standardized residuals or
  "Standardize" the Matrix sigma: \eqn{R^{-1} matrix /
  \sigma}{R^-1 matrix / sigma}
}

