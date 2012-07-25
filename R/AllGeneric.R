##' Calculates the raw robustified log-likelihood.
##'
##' @title Raw robustified log-Likelihood
##' @param object a lmerMod object
##' @param norm.only if only the norms part should be calculated
##' @param ... parameters passed to sumRho.e() and sumRho.b()
##' @return calculated value
##' @export
robll <- function(object, ...) UseMethod("robll")

##' Calculates the matrix Q for the Laplace Approximation
##'
##' @title Matrix Q for Laplace Approximation
##' @param object a lmerMod object
##' @param determinant return determinant?
##' @param log logarithm of determinant
##' @param ... ignored
##' @return calculated value
##' @export
Q <- function(object, ...)  UseMethod("Q")

##' @rdname updateDeviance
##' @export
laplace <- function(object, ...) UseMethod("laplace")

##' Calculates the gradient, the gradient of beta and b
##' if a lmerMod object is given, the gradient of theta
##' if a reModule is given.
##'
##' @title Gradient
##' @param object a lmerMod object
##' @param ... ignored
##' @return vector
##' @export
gradient <- function(object, ...) UseMethod("gradient")

## force imports
isDiagonal <- function(...) Matrix:::isDiagonal(...)
isTriangular <- function(...) Matrix:::isTriangular(...)

#####################################################
## Accessor methods                                ##
#####################################################

##' Simple accessors for the (spherical) random effects.
##'
##' @title Random effects accessor
##' @param object merMod object
##' @param ... ignored
##' @export
b <- function(object, ...) UseMethod("b")

##' @rdname b
##' @export
u <- function(object, ...) UseMethod("u")


#####################################################
## Summary / printing methods                      ##
#####################################################

##' Prepare object for producing a comparison chart.
##'
##' @title getInfo of an object
##' @param object object
##' @param ... ignored
##' @return list with estimated coefficients, estimated
##'   variance components, sigma, deviance and
##'   parameter configuration used to fit.
##' @export
getInfo <- function(object, ...) UseMethod("getInfo")
