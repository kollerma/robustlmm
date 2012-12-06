## force imports
isDiagonal <- function(...) Matrix:::isDiagonal(...)
isTriangular <- function(...) Matrix:::isTriangular(...)

#####################################################
## Accessor methods                                ##
#####################################################

## Simple accessors for the (spherical) random effects.
##
## @title Random effects accessor
## @param object merMod object
## @param ... ignored
b <- function(object, ...) UseMethod("b")

## @rdname b
u <- function(object, ...) UseMethod("u")

##' Simple accessors for estimated residual scale sigma.
##'
##' @title Get sigma
##' @param object rlmerMod object
##' @param ... ignored
##' @export
sigma <- function(object, ...) UseMethod("sigma")

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
