#######################################################
## Basic accessor functions                          ##
#######################################################

## Get Z
## (slot Zt in reModule might be changed to Z, so
## don't write it explicitly)
##
## @title Get Z from reModule
## @param object merMod object
## @param t transpose or not
getZ <- function(object, t = FALSE) {
    if (t) object@pp$Zt else t(object@pp$Zt)
}

## Get X
##
## @title Get X from predModule
## @param object merMod object
## @param t transpose or not
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
    ret <- switch(which,
                  sigma=object@rho.sigma.b,
                  default=,object@rho.b)
    ## backwards compatibility:
    if (inherits(ret, "psi_func"))
        ret <- rep.int(list(ret), length(object@blocks))
    ## add names
    names(ret) <- names(object@cnms)
    ret
}

##' Get theta
##'
##' @title Get theta
##' @param object merMod object
##' @export
theta <- function(object) {
    if (is(object, "rlmerMod")) {
        ## add names like lme4
        tt <- object@pp$theta
        nc <- c(unlist(mapply(function(g,e) {
            mm <- outer(e,e,paste,sep=".")
            diag(mm) <- e
            mm <- mm[lower.tri(mm,diag=TRUE)]
            paste(g,mm,sep=".")
        }, names(object@cnms),object@cnms)))
        names(tt) <- nc
        tt
    } else getME(object, "theta")
}

## Get Lambda
##
## @title Get Lambda
## @param object merMod object
Lambda <- function(object) {
    ## FIXME: which theta?
    if (class(object)[1] == "lmerMod") t(object@pp$Lambdat) else t(object@pp$Lambdat())
}

## Get U_b
##
## @title Get U_b
## @param object merMod object
U_b <- function(object) {
   if (class(object)[1] == "lmerMod") t(object@pp$Lambdat) else object@pp$U_b 
}

## Get Lind
##
## @title Get Lind
## @param object merMod object
Lind <- function(object) {
    object@pp$Lind
}

## Get lower
##
## @title Get lower
## @param object merMod object
lower <- function(object) {
    object@lower
}

## Get indices of the random effects corresponding
## to zero variance components
##
## @title Get indices of r.e. corresp. to zero v.c.
## @param object merMod object
getZeroU <- function(object) object@pp$zeroB

## Get various numbers of parameters / lengths of vectors.
##
## \itemize{
## \item b, u: length of random effects vector
## \item beta, coef: number of fixed effects
## \item theta: length of vector theta
## \item r, e: number of observations
## }
## 
## @title length
## @param x merMod object
## @param what length is requested
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
##' @importFrom stats nobs
##' @S3method nobs rlmerMod
nobs.rlmerMod <- .nobsLmerMod

### Get REML (so that we are not coercing all the time)
##' @importFrom lme4 isREML
.isREML <- function(x, ...) {
    as.logical(x@devcomp$dims["REML"])
}

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
##' @importFrom MatrixModels residuals resid
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
##' @importFrom lme4 sigma
##' @S3method sigma rlmerMod
sigma.rlmerMod <- .sigma
##' @S3method sigma mer
sigma.mer <- function(object, ...) lme4::sigma(object, ...)
##' @S3method sigma merMod
sigma.merMod <- function(object, ...) lme4:::sigma.merMod(object, ...)

### Get deviance
.deviance <- function(object, ...) NA
##' @S3method deviance rlmerMod
deviance.rlmerMod <- .deviance

.mu <- function(object)
{
    ## Purpose: calculate mu of respModule
    ## ----------------------------------------------------------------------
    ## Arguments: object: lmerMod object
    ## ----------------------------------------------------------------------
    ## Author: Manuel Koller, Date: 11 Apr 2011, 11:40

    ## FIXME: ?? offset will be added in updateMu
    drop((crossprod(object@pp$Zt, object@pp$b))@x + (object@pp$X %*% object@pp$beta))
}

### Get fixed effects
.fixef <- function(object) object@pp$beta

### Get u
b.s <- .u <- function(object, ...) object@pp$b.s
u.rlmerMod <- function(object, ...) {
    ret <- object@pp$b.s
    names(ret) <- dimnames(getZ(object))[[2]]
    ret
}

u.lmerMod <- function(object, ...) object@u

### Get b
.b <- function(object, ...) object@pp$b
b.rlmerMod <- function(object, ...) {
    ret <- object@pp$b
    names(ret) <- dimnames(getZ(object))[[2]]
    ret
}
b.lmerMod <- function(object, ...) {
    ret <- crossprod(getME(object, "Lambdat"), getME(object, "u"))
    names(ret) <- dimnames(getME(object, "Zt"))[[1]]
    ret
}

### Get coefficients
##' @importFrom stats coef
##' @S3method coef rlmerMod
coef.rlmerMod <- function(object, ...) {
    ret <- list(beta=fixef(object),b=b(object))
    class(ret) <- "coef.rlmerMod"
    ret
}

### Get ranef
##' @importFrom nlme ranef
##' @S3method ranef rlmerMod
ranef.rlmerMod <- function(object, ...) {
    ## FIXME: add postVar, drop and whichel arguments
    b <- b(object)
    ret <- lapply(object@idx, function(bidx) {
        lret <- b[bidx]
        dim(lret) <- dim(bidx)
        ## add rownames
        colnames(lret) <- names(b)[bidx[1,]]
        as.data.frame(t(lret))
    })
    ## add names and colnames
    names(ret) <- names(object@cnms)
    for (i in names(ret))
        colnames(ret[[i]]) <- object@cnms[[i]]
    class(ret) <- "ranef.rlmer"
    ret
}

## return u as list arranged like ranef
## but do not bother setting names
uArranged <- function(object, b.s = b.s(object)) {
    ret <- lapply(object@idx, function(bidx) {
        lret <- b.s[bidx]
        dim(lret) <- dim(bidx)
        t(lret)
    })
    ret
}
