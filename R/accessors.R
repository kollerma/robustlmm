#######################################################
## Basic accessor functions                          ##
#######################################################

##' Get Z
##' (slot Zt in reModule might be changed to Z, so
##' don't write it explicitly)
##'
##' @title Get Z from reModule
##' @param object merMod object
##' @param t transpose or not
getZ <- function(object, t = FALSE) {
    if (t) object@pp$Zt else t(object@pp$Zt)
}

##' Get X
##'
##' @title Get X from predModule
##' @param object merMod object
##' @param t transpose or not
getX <- function(object, t = FALSE) {
    if (t) t(object@pp$X) else object@pp$X
}

##' Get rho-function used for residuals
##'
##' @title Get rho_e
##' @param object merMod object
##' @param which add "sigma" for rho.sigma.e
##' @export
rho.e <- function(object, which = "default") {
    switch(which,
           sigma=object@rho.sigma.e,
           default=,object@rho.e)
}

##' Get rho-function for used random effects
##'
##' @title Get rho_b
##' @param object merMod object
##' @param which add "sigma" for rho.sigma.e
##' @export
rho.b <- function(object, which = "default") {
    switch(which,
           sigma=object@rho.sigma.b,
           default=,object@rho.b)
}

##' Get weights exponent used for residuals
##'
##' @title Get wExp_e
##' @param object merMod object
##' @export
wExp.e <- function(object) {
    object@wExp.e
}

##' Get weights exponent used for random effects
##'
##' @title Get wExp_b
##' @param object merMod object
##' @export
wExp.b <- function(object) {
    object@wExp.b
}

##' Get theta
##'
##' @title Get theta
##' @param object merMod object
##' @export
theta <- function(object) {
    object@pp$theta
}

##' Get Lambda
##'
##' @title Get Lambda
##' @param object merMod object
Lambda <- function(object) {
    ## FIXME: which theta?
    if (class(object)[1] == "lmerMod") t(object@pp$Lambdat) else t(object@pp$Lambdat())
}

##' Get U_b
##'
##' @title Get U_b
##' @param object merMod object
U_b <- function(object) {
   if (class(object)[1] == "lmerMod") t(object@pp$Lambdat) else object@pp$U_b 
}

##' Get Lind
##'
##' @title Get Lind
##' @param object merMod object
Lind <- function(object) {
    object@pp$Lind
}

##' Get lower
##'
##' @title Get lower
##' @param object merMod object
lower <- function(object) {
    object@lower
}

##' Get indices of the random effects corresponding
##' to zero variance components
##'
##' @title Get indices of r.e. corresp. to zero v.c.
##' @param object merMod object
getZeroU <- function(object) object@pp$zeroB

##' Get various numbers of parameters / lengths of vectors.
##'
##' \itemize{
##' \item b, u: length of random effects vector
##' \item beta, coef: number of fixed effects
##' \item theta: length of vector theta
##' \item r, e: number of observations
##' }
##' 
##' @title length
##' @param x merMod object
##' @param what length is requested
len <- function(x, what) switch(what,
                                u=,
                                b=nrow(x@pp$Zt),
                                coef=,
                                beta=length(x@beta),
                                theta=length(x@theta),
                                r=,
                                e=length(x@resp$y),
                                stop("unknown length"))

.nobsLmerMod <- function(object, ...) len(object, "e")
##' @exportMethod nobs
setMethod("nobs",
          signature(object = "lmerMod"),
          .nobsLmerMod)

### Get REML (so that we are not coercing all the time)
.isREML <- function(x) {
    as.logical(x@devcomp$dims["REML"])
}

isREML.rlmerMod <- .isREML

##' The per-observation residuals are returned, i.e.,
##' the difference of the observation and the fitted value
##' including random effects. With type one can specify whether
##' the weights should be used or not. 
##' 
##' @title Get residuals
##' @param object rlmerMod object
##' @param type type of residuals
##' @param ... ignored
##' @method residuals rlmerMod
##' @S3method residuals rlmerMod
residuals.rlmerMod <- function(object, type = c("response", "weighted"), ...) {
    type <- match.arg(type)
    switch(type,
           ## FIXME: really??
           response = object@resp$wtres,
           weighted = wgt.e(object) * object@resp$wtres,
           stop("unknown type of residual"))
}

### Get sigma (so that we are not coercing all the time)
.sigma <- function(object, ...) object@pp$sigma
##' @S3method sigma rlmerMod
sigma.rlmerMod <- .sigma

### Get deviance
.deviance <- function(object, ...) object@pp$deviance
##' @S3method deviance rlmerMod
deviance.rlmerMod <- .deviance

.mu <- function(object)
{
    ## Purpose: calculate mu of respModule
    ## ----------------------------------------------------------------------
    ## Arguments: object: lmerMod object
    ## ----------------------------------------------------------------------
    ## Author: Manuel Koller, Date: 11 Apr 2011, 11:40

    ## offset will be added in updateMu
    (crossprod(object@pp$Zt, object@pp$b))@x + (object@pp$X %*% object@pp$beta)
}

### Get fixed effects
.fixef <- function(object) object@pp$beta

### Get u
b.s <- .u <- function(object, ...) object@pp$b.s
##' @rdname b
##' @method u rlmerMod
##' @S3method u rlmerMod
u.rlmerMod <- .u
##' @rdname b
##' @method u lmerMod
##' @S3method u lmerMod
u.lmerMod <- function(object, ...) object@u

### Get b
.b <- function(object, ...) object@pp$b
##' @rdname b
##' @method b rlmerMod
##' @S3method b rlmerMod
b.rlmerMod <- .b
##' @rdname b
##' @method b lmerMod
##' @S3method b lmerMod
b.lmerMod <- function(object, ...) drop(crossprod(object@pp$Lambdat, object@u))

### Get ranef
##' @S3method ranef rlmerMod
ranef.rlmerMod <- function(object)
    ranef(as(object, "lmerMod"))
