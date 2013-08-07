## This file contains functions that have been copied
## from lme4 to enable a release of this package onto CRAN.
## Copied from Revision 1785
## The authors of lme4 are:
## Douglas Bates <bates@stat.wisc.edu>,
## Martin Maechler <maechler@R-project.org> and
## Ben Bolker <bolker@mcmaster.ca>

## minimal changes: merMod -> rlmerMod
## plus some glue for mer objects, if required

##' Extract the fixed-effects estimates.
##'
##' Extract the estimates of the fixed-effects parameters from a fitted model.
##' @name fixef
##' @title Extract fixed-effects estimates
##' @aliases fixef fixed.effects fixef.rlmerMod
##' @docType methods
##' @param object any fitted model object from which fixed effects estimates can
##' be extracted.
##' @param \dots optional additional arguments. Currently none are used in any
##' methods.
##' @return a named, numeric vector of fixed-effects estimates.
##' @keywords models
##' @examples
##' ## doFit = FALSE to speed up example
##' fixef(rlmer(Reaction ~ Days + (Days|Subject), sleepstudy, doFit=FALSE))
##' @importFrom nlme fixef
##' @export fixef
##' @S3method fixef rlmerMod
##' @S3method fixef mer
##' @export
fixef.rlmerMod <- function(object, ...)
    structure(object@beta, names = dimnames(object@pp$X)[[2]])
fixef.mer <- function(object, ...) lme4::fixef(object, ...)

## Extract the conditional variance-covariance matrix of the fixed-effects
## parameters
##
## @title Extract conditional covariance matrix of fixed effects
## @param sigma numeric scalar, the residual standard error
## @param unsc matrix of class \code{"\linkS4class{dpoMatrix}"}, the
##     unscaled variance-covariance matrix
## @param nmsX character vector of column names of the model matrix
## @param correlation logical scalar, should the correlation matrix
##     also be evaluated.
## @param ... additional, optional parameters.  None are used at present.
mkVcov <- function(sigma, unsc, nmsX, correlation = TRUE, ...) {
        V <- sigma^2 * unsc
            if(is.null(rr <- tryCatch(as(V, "dpoMatrix"),
                                            error = function(e) NULL)))
                stop("Computed variance-covariance matrix is not positive definite")
            dimnames(rr) <- list(nmsX, nmsX)
            if(correlation)
                rr@factors$correlation <-
                        if(!is.na(sigma)) as(rr, "corMatrix") else rr # (is NA anyway)
            rr
    }

##' @importFrom stats vcov
##' @S3method vcov rlmerMod
vcov.rlmerMod <- function(object, correlation = TRUE, sigm = sigma(object), ...)
    mkVcov(sigm, unsc = object@pp$unsc(), nmsX = colnames(object@pp$X),
   correlation=correlation, ...)

##' @importFrom stats vcov
##' @S3method vcov summary.rlmer
vcov.summary.rlmer <- function(object, correlation = TRUE, ...) {
    if(is.null(object$vcov)) stop("logic error in summary of rlmerMod object")
    object$vcov
}

mkVarCorr <- function(sc, cnms, nc, theta, nms) {
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
    ans <- lapply(ncseq, function(i)
      {
  ## Li := \Lambda_i, the i-th block diagonal of \Lambda(\theta)
  Li <- diag(nrow = nc[i])
  Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
  rownames(Li) <- cnms[[i]]
  ## val := \Sigma_i = \sigma^2 \Lambda_i \Lambda_i', the
  val <- tcrossprod(sc * Li) # variance-covariance
  stddev <- sqrt(diag(val))
  correl <- t(val / stddev)/stddev
  diag(correl) <- 1
  attr(val, "stddev") <- stddev
  attr(val, "correlation") <- correl
  val
      })
    if(is.character(nms)) names(ans) <- nms
    attr(ans, "sc") <- sc
    ans
}

##' Extract variance and correlation components
##'
##' This function calculates the estimated variances, standard
##' deviations, and correlations between the random-effects terms in a
##' mixed-effects model, of class \code{\linkS4class{rlmerMod}}
##' (linear, generalized or nonlinear).  The within-group error
##' variance and standard deviation are also calculated.
##'
##' @name VarCorr
##' @aliases VarCorr VarCorr.rlmerMod
##' @param x a fitted model object, usually an object inheriting from class
##' \code{\linkS4class{rlmerMod}}.
##' @param sigma an optional numeric value used as a multiplier for the standard
##' deviations.  Default is \code{1}.
##' @param rdig an optional integer value specifying the number of digits used
##' to represent correlation estimates.  Default is \code{3}.
##' @return a list of matrices, one for each random effects grouping term.
##' For each grouping term, the standard deviations and correlation matrices for each grouping term
##' are stored as attributes \code{"stddev"} and \code{"correlation"}, respectively, of the
##' variance-covariance matrix, and
##' the residual standard deviation is stored as attribute \code{"sc"}
##' (for \code{glmer} fits, this attribute stores the scale parameter of the model).
##' @author This is modeled after \code{\link[nlme]{VarCorr}} from package
##' \pkg{nlme}, by Jose Pinheiro and Douglas Bates.
##' @seealso \code{\link{lmer}}, \code{\link{nlmer}}
##' @examples
##' data(Orthodont, package="nlme")
##' fm1 <- lmer(distance ~ age + (age|Subject), data = Orthodont)
##' VarCorr(fm1)
##' @keywords models
##' @importFrom nlme VarCorr
##' @export VarCorr
##' @S3method VarCorr rlmerMod
##' @S3method VarCorr mer
##' @S3method VarCorr summary.rlmer
##' @export
VarCorr.rlmerMod <- function(x, sigma, rdig)# <- 3 args from nlme
{
    ## FIXME:: add type=c("varcov","sdcorr","logs" ?)
    if (is.null(cnms <- x@cnms))
        stop("VarCorr methods require reTrms, not just reModule")
    if(missing(sigma)) # "bug": fails via default 'sigma=sigma(x)'
        sigma <- .sigma(x)
    nc <- sapply(cnms, length)  # no. of columns per term
    m <- mkVarCorr(sigma, cnms=cnms, nc=nc, theta = x@theta,
                   nms = {fl <- x@flist; names(fl)[attr(fl, "assign")]})
    attr(m,"useSc") <- as.logical(x@devcomp$dims["useSc"])
    class(m) <- "VarCorr.merMod"
    m
}
VarCorr.mer <- function(x, ...) lme4::VarCorr(x, ...)
VarCorr.summary.rlmer <- function(x, ...) x$varcor

## FIXME: should ... go to formatVC or to print ... ?
##' @S3method print VarCorr.rlmerMod
print.VarCorr.rlmerMod <- function(x,digits = max(3, getOption("digits") - 2), ...) {
    print(formatVC(x, digits = digits, useScale = attr(x,"useSc"),  ...),quote=FALSE)
}

### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
formatVC <- function(varc, digits = max(3, getOption("digits") - 2),
                          useScale) {
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"), if(useScale) list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
                   list(rep.int('', nr),
                        c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(varc, colnames)), if(useScale) "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
        maxlen <- max(reLens)
        corr <-
            do.call("rBind",
                    lapply(recorr,
                           function(x) {
                               x <- as(x, "matrix")
                               cc <- format(round(x, 3), nsmall = 3)
                               cc[!lower.tri(cc)] <- ""
                               nr <- dim(cc)[1]
                               if (nr >= maxlen) return(cc)
                               cbind(cc, matrix("", nr, maxlen-nr))
                           }))[, -maxlen, drop = FALSE]
        if (nrow(corr) < nrow(reMat))
            corr <- rbind(corr, matrix("", nrow = nrow(reMat) - nrow(corr), ncol = ncol(corr)))
        colnames(corr) <- rep.int("", ncol(corr))
        colnames(corr)[1] <- "Corr"
        cbind(reMat, corr)
    } else reMat
}

getFixedFormula <- function(form) {
    form[[3]] <- if (is.null(nb <- lme4:::nobars(form[[3]]))) 1 else nb
    form
}

##' @importFrom stats formula
##' @S3method formula rlmerMod
formula.rlmerMod <- function(x, fixed.only=FALSE, ...) {
    form <- formula(getCall(x),...)
    if (fixed.only) {
        form <- getFixedFormula(form)
    }
    form
}

##' @importFrom stats terms
##' @S3method terms rlmerMod
terms.rlmerMod <- function(x, fixed.only=TRUE, ...) {
  if (fixed.only) {
      tt <- terms.formula(formula(x,fixed.only=TRUE))
      attr(tt,"predvars") <- attr(attr(x@frame,"terms"),"predvars.fixed")
      tt
  }
  else attr(x@frame,"terms")
}

##' @importFrom stats weights
##' @S3method weights rlmerMod
weights.rlmerMod <- function(object, ...) {
  object@resp$weights
}

##' @importFrom stats fitted
##' @S3method fitted rlmerMod
fitted.rlmerMod <- function(object, ...) object@resp$mu

##' @importFrom stats model.frame
##' @S3method model.frame rlmerMod
model.frame.rlmerMod <- function(formula, ...) formula@frame

##' @importFrom stats model.matrix
##' @S3method model.matrix rlmerMod
model.matrix.rlmerMod <- function(object, ...) object@pp$X

