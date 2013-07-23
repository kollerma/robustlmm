## force imports
isDiagonal <- function(...) Matrix:::isDiagonal(...)
isTriangular <- function(...) Matrix:::isTriangular(...)

#####################################################
## Accessor methods                                ##
#####################################################

## Simple accessor for the (spherical) random effects.
##
## @title Random effects accessor
## @param object merMod object
## @param ... ignored
b <- function(object, ...) UseMethod("b")

## @rdname b
u <- function(object, ...) UseMethod("u")

##' Simple accessor for estimated residual scale sigma.
##'
##' @title Get sigma
##' @param object rlmerMod object
##' @param ... ignored
##' @seealso \code{\link{getME}}
##' @examples
##' fm <- rlmer(Yield ~ (1|Batch), Dyestuff)
##' sigma(fm)
##' getME(fm, "devcomp")$cmp[c("sigmaML", "sigmaREML")]
##' @export
sigma <- function(object, ...) UseMethod("sigma")

#####################################################
## Summary / printing methods                      ##
#####################################################

##' Extract some information from objects returned from
##' \code{\link{rlmer}} and \code{\link{lmer}} in the
##' form of a simple list. Internally used to 
##' prepare object for producing a comparison chart in
##' \code{\link{compare}}.
##'
##' @title getInfo of an object
##' @param object object
##' @param ... ignored
##' @return list with estimated coefficients, estimated
##'   variance components, sigma, deviance and
##'   parameter configuration used to fit.
##' @examples
##' fm <- rlmer(Yield ~ (1|Batch), Dyestuff)
##' str(getInfo(fm))
##' @export
getInfo <- function(object, ...) UseMethod("getInfo")
