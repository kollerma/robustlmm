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

## Get rho-function used for residuals
##
## @title Get rho_e
## @param object merMod object
## @param which add "sigma" for rho.sigma.e
## @export
rho.e <- function(object, which = "default") {
    switch(which,
           sigma=object@rho.sigma.e,
           default=,object@rho.e)
}

## Get rho-function for used random effects
##
## @title Get rho_b
## @param object merMod object
## @param which add "sigma" for rho.sigma.e
## @export
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
##' @param scaled scale residuals by residual standard deviation (=scale parameter)?
##' @param ... ignored
##' @examples
##' fm <- rlmer(Yield ~ (1|Batch), Dyestuff)
##' stopifnot(all.equal(resid(fm, type="weighted"),
##'                     resid(fm) * getME(fm, "w_e")))
##' @method residuals rlmerMod
##' @importFrom stats residuals resid
##' @S3method residuals rlmerMod
residuals.rlmerMod <- function(object, type = c("response", "weighted"),
                               scaled=FALSE, ...) {
    type <- match.arg(type)
    r <- switch(type,
                ## FIXME: really??
                response = object@resp$wtres,
                weighted = wgt.e(object) * object@resp$wtres,
                stop("unknown type of residual"))
    if (scaled) r <- r/sigma(object)
    r
}

### Get sigma (so that we are not coercing all the time)
.sigma <- function(object, ...) object@pp$sigma
##' @importFrom lme4 sigma
##' @S3method sigma rlmerMod
sigma.rlmerMod <- .sigma


### Get deviance
.deviance <- function(object, ...)
    stop("Deviance is not defined for rlmerMod objects")
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
    drop(crossprod(object@pp$Zt, object@pp$b.r) + (object@pp$X %*% object@pp$beta))
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
.b <- function(object, ...) object@pp$b.r
b.rlmerMod <- function(object, ...) {
    ret <- object@pp$b.r
    names(ret) <- dimnames(getZ(object))[[2]]
    ret
}
b.lmerMod <- function(object, ...) {
    ret <- crossprod(getME(object, "Lambdat"), getME(object, "u"))
    names(ret) <- dimnames(getME(object, "Zt"))[[1]]
    ret
}

### Get ranef
##' @importFrom nlme ranef
##' @S3method ranef rlmerMod
ranef.rlmerMod <- function(object, ...) {
    ## FIXME: add postVar, drop and whichel arguments
    b <- b(object)
    ret <- uArrangedNames(object, b.s = b)
    class(ret) <- "ranef.rlmerMod"
    ret
}

## return u as list arranged like ranef
uArrangedNames <- function(object, b.s = b.s(object)) {
    ret <- lapply(object@idx, function(bidx) {
        lret <- b.s[bidx]
        dim(lret) <- dim(bidx)
        ## add rownames
        colnames(lret) <- names(b.s)[bidx[1,]]
        as.data.frame(t(lret))
    })
    ## add names and colnames
    names(ret) <- names(object@cnms)
    for (i in names(ret))
        colnames(ret[[i]]) <- object@cnms[[i]]
    ret
}
## same as uArrangedNames, but do not set names
uArranged <- function(object, b.s = b.s(object)) {
    ret <- lapply(object@idx, function(bidx) {
        lret <- b.s[bidx]
        dim(lret) <- dim(bidx)
        t(lret)
    })
    ret
}

##' Extract or Get Generalize Components from a Fitted Mixed Effects Model
##'
##' Extract (or \dQuote{get}) \dQuote{components} -- in a generalized
##' sense -- from a fitted mixed-effects model, i.e. from an object
##' of class \code{"\linkS4class{rlmerMod}"} or \code{"\linkS4class{merMod}"}.
##'
##' The goal is to provide \dQuote{everything a user may want} from a fitted
##' \code{"rlmerMod"} object \emph{as far} as it is not available by methods, such
##' as \code{\link{fixef}}, \code{\link{ranef}}, \code{\link{vcov}}, etc.
##'
##' @param object a fitted mixed-effects model of class
##' \code{"\linkS4class{rlmerMod}"}, i.e. typically the result of
##' \code{\link{rlmer}()}.
##' @param name a character string specifying the name of the
##' \dQuote{component}.  Possible values are:\cr
##' \describe{
##'     \item{X}{fixed-effects model matrix}
##'     \item{Z}{random-effects model matrix}
##'     \item{Zt}{transpose of random-effects model matrix}
##'     \item{u}{conditional mode of the \dQuote{spherical} random effects variable}
##'     \item{b.s}{synonym for \dQuote{u}}
##'     \item{b}{onditional mode of the random effects variable}
##'     \item{Gp}{groups pointer vector.  A pointer to the beginning of each group
##'               of random effects corresponding to the random-effects terms.}
##'     \item{Lambda}{relative covariance factor of the random effects.}
##'     \item{U_b}{synonym for \dQuote{Lambda}}
##'     \item{Lambdat}{transpose of the relative covariance factor of the random effects.}
##'     \item{Lind}{index vector for inserting elements of \eqn{\theta}{theta} into the
##'                 nonzeros of \eqn{\Lambda}{Lambda}}
##'     \item{flist}{a list of the grouping variables (factors) involved in the random effect terms}
##'     \item{beta}{fixed-effects parameter estimates (identical to the result of \code{\link{fixef}}, but without names)}
##'     \item{theta}{random-effects parameter estimates: these are parameterized as the relative Cholesky factors of each random effect term}
##'     \item{n_rtrms}{number of random-effects terms}
##'     \item{devcomp}{a list consisting of a named numeric vector, \dQuote{cmp}, and
##'                    a named integer vector, \dQuote{dims}, describing the fitted model}
##'     \item{offset}{model offset}
##'     \item{lower}{lower bounds on model parameters (random effects parameters only)}
##'     \item{rho_e}{rho function used for the residuals}
##'     \item{rho_b}{list of rho functions used for the random effects}
##'     \item{rho_sigma_e}{rho function used for the residuals when estimating sigma}
##'     \item{rho_sigma_b}{list of rho functions used for the random effects when estimating the covariance parameters}
##'     \item{M}{list of matrices, blocks of the Henderson's equations and the matrices used for computing the linear approximations of the estimates of beta and spherical random effects.}
##'     \item{w_e}{robustness weights associated with the observations}
##'     \item{w_b}{robustness weights associated with the spherical random effects, returned in the same format as \code{\link{ranef}()}}
##'     \item{w_b_vector}{robustness weights associated with the spherical random effects, returned as one long vector}
##'     \item{w_sigma_e}{robustness weights associated with the observations when estimating sigma}
##'     \item{w_sigma_b}{robustness weights associated with the spherical random effects when estimating the covariance parameters, returned in the same format as \code{\link{ranef}()}}
##'     \item{w_sigma_b_vector}{robustness weights associated with the spherical random effects when estimating the covariance parameters, returned as one long vector}
##' }
##' @return Unspecified, as very much depending on the \code{\link{name}}.
##' @seealso \code{\link{getCall}()},
##' More standard methods for rlmerMod objects, such as \code{\link{ranef}},
##' \code{\link{fixef}}, \code{\link{vcov}}, etc.:
##' see \code{methods(class="rlmerMod")}
##' @keywords utilities
##' @usage getME(object,
##'   name = c("X", "Z", "Zt", "u", "b.s", "b", "Gp", "Lambda",
##'            "Lambdat", "U_b", "Lind", "flist", "beta", "theta",
##'            "n_rtrms", "devcomp", "offset", "lower", "rho_e",
##'            "rho_b", "rho_sigma_e", "rho_sigma_b", "M", "w_e",
##'            "w_b", "w_b_vector", "w_sigma_e", "w_sigma_b",
##'            "w_sigma_b_vector"))
##' @examples
##'
##' ## shows many methods you should consider *before* using getME():
##' methods(class = "rlmerMod")
##'
##' ## doFit = FALSE to speed up example
##' (fm1 <- rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
##'               method="DASvar", doFit=FALSE))
##' Z <- getME(fm1, "Z")
##' stopifnot(is(Z, "CsparseMatrix"),
##'           c(180,36) == dim(Z),
##' 	  all.equal(fixef(fm1), getME(fm1, "beta"),
##' 		    check.attr=FALSE, tol = 0))
##'
##' ## All that can be accessed [potentially ..]:
##' (nmME <- eval(formals(getME)$name))
##' \dontshow{
##' ## internal consistency check ensuring that all work:
##' ## "try(.)" because some are not yet implemented:
##' str(parts <- sapply(nmME, function(nm) try(getME(fm1, nm)),
##'                     simplify=FALSE))
##' }% dont..
##'
##' @export
getME <- function(object,
		  name = c("X", "Z","Zt", "u", "b.s", "b",
		  "Gp","Lambda", "Lambdat", "U_b", "Lind",
                  "flist", "beta", "theta","n_rtrms",
                  "devcomp", "offset", "lower", "rho_e",
                  "rho_b", "rho_sigma_e", "rho_sigma_b", "M",
                  "w_e", "w_b", "w_b_vector", "w_sigma_e",
                  "w_sigma_b", "w_sigma_b_vector"))
{
    if(missing(name)) stop("'name' must not be missing")
    if (is(object, "merMod"))
        return(lme4::getME(object, name))
    ## else: assume it's rlmerMod
    stopifnot(length(name <- as.character(name)) == 1,
	      is(object, "rlmerMod"))
    name <- match.arg(name)
    rsp  <- object@resp
    PR   <- object@pp
    dc   <- object@devcomp
    cmp  <- dc $ cmp
    dims <- dc $ dims
    switch(name,
	   "X" = getX(object, t=FALSE), 
	   "Z" = getZ(object, t=FALSE),
	   "Zt"= getZ(object, t=TRUE),
           "u" =,
           "b.s" = b.s(object),
           "b" = b(object),
	   "Lambda"= ,
           "U_b" = PR$ U_b,
	   "Lambdat"= PR$ Lambdat(),
           "Lind" = PR$ Lind,
           "Gp" = object@Gp,
           "flist" = object@flist,
	   "beta" = object@beta,
           "theta"= theta(object),
	   "n_rtrms" = length(object@flist), 
           "devcomp" = dc,
           "offset" = rsp$offset,
           "lower" = object@lower,
           "rho_e" = rho.e(object),
           "rho_b" = rho.b(object),
           "rho_sigma_e" = rho.e(object, "sigma"),
           "rho_sigma_b" = rho.b(object, "sigma"),
           "M" = PR$ M(),
           "w_e" = wgt.e(object),
           "w_b" = uArrangedNames(object, wgt.b(object)),
           "w_b_vector" = wgt.b(object),
           "w_sigma_e" = wgt.e(object, use.rho.sigma=TRUE),
           "w_sigma_b" = uArrangedNames(object, wgt.b(object, center=TRUE)),
           "w_sigma_b_vector" = wgt.b(object, center=TRUE),
	   "..foo.." =# placeholder!
	   stop(gettextf("'%s' is not implemented yet",
			 sprintf("getME(*, \"%s\")", name))),
	   ## otherwise
	   stop(sprintf("Mixed-Effects extraction of '%s' is not available for class \"%s\"",
			name, class(object))))
}## {getME}

##' The function \code{theta} is short for \code{getME(, "theta")}.
##'
##' @rdname getME
##' @examples
##' stopifnot(all.equal(theta(fm1), getME(fm1, "theta")))
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
