## This file contains functions that have been copied
## from lme4 to enable a release of this package onto CRAN.
## Copied from Revision 1785
## Updated from commit 9c7edde3f081af68f84b1e9b61be1f5cc52cfaa0
## The authors of lme4 are:
## Douglas Bates <bates@stat.wisc.edu>,
## Martin Maechler <maechler@R-project.org> and
## Ben Bolker <bolker@mcmaster.ca>
## Steven Walker

## follows the file lme4/R/lmer.R

## minimal changes: merMod -> rlmerMod
## and replaced with a stop message if not defined

## coef() method for all kinds of "mer", "*merMod", ... objects
## ------  should work with fixef() + ranef()  alone
coefRlmerMod <- function(object, ...)
{
    if (length(list(...)))
	warning('arguments named "', paste(names(list(...)), collapse = ", "),
                '" ignored')
    fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
    ref <- ranef(object)
    ## check for variables in RE but missing from FE, fill in zeros in FE accordingly
    refnames <- unlist(lapply(ref,colnames))
    nmiss <- length(missnames <- setdiff(refnames,names(fef)))
    if (nmiss >0) {
        fillvars <- setNames(data.frame(rbind(rep(0,nmiss))),missnames)
        fef <- cbind(fillvars,fef)
    }
    val <- lapply(ref, function(x)
		  fef[rep.int(1L, nrow(x)),,drop = FALSE])
    for (i in seq(a = val)) {
	refi <- ref[[i]]
	row.names(val[[i]]) <- row.names(refi)
	nmsi <- colnames(refi)
	if (!all(nmsi %in% names(fef)))
	    stop("unable to align random and fixed effects")
	for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
    }
    class(val) <- "coef.rlmerMod"
    val
} ##  {coefRlmerMod}

##' @importFrom stats coef
##' @S3method coef rlmerMod
coef.rlmerMod <- coefRlmerMod

## deviance is in accessors.R

## FIXME what about drop1

##' @importFrom stats extractAIC
##' @S3method extractAIC rlmerMod
extractAIC.rlmerMod <- function(fit, scale = 0, k = 2, ...)
   stop("AIC is not defined for rlmerMod objects")

##' @importFrom stats family
##' @S3method family rlmerMod
family.rlmerMod <- function(object, ...) gaussian()

##' @importFrom stats fitted
##' @S3method fitted rlmerMod
fitted.rlmerMod <- function(object, ...) object@resp$mu

## Extract the fixed-effects estimates.
##
## Extract the estimates of the fixed-effects parameters from a fitted model.
## @name fixef
## @title Extract fixed-effects estimates
## @aliases fixef fixed.effects fixef.rlmerMod
## @docType methods
## @param object any fitted model object from which fixed effects estimates can
## be extracted.
## @param \dots optional additional arguments. Currently none are used in any
## methods.
## @return a named, numeric vector of fixed-effects estimates.
## @keywords models
## @examples
## ## doFit = FALSE to speed up example
## fixef(rlmer(Reaction ~ Days + (Days|Subject), sleepstudy, doFit=FALSE))
##' @importFrom nlme fixef
##' @S3method fixef rlmerMod
## @export
fixef.rlmerMod <- function(object, ...)
    structure(object@beta, names = dimnames(object@pp$X)[[2]])

##' @importFrom lme4 nobars
getFixedFormula <- function(form) {
    form[[3]] <- if (is.null(nb <- nobars(form[[3]]))) 1 else nb
    form
}

##' @importFrom stats formula
##' @S3method formula rlmerMod
formula.rlmerMod <- function(x, fixed.only=FALSE, ...) {
    if (is.null(form <- attr(x@frame,"formula"))) {
        if (!grepl("rlmer$",deparse(getCall(x)[[1]])))
            stop("can't find formula stored in model frame or call")
        form <- as.formula(formula(getCall(x),...))
    }
    if (fixed.only) {
        form <- getFixedFormula(form)
    }
    form
}

##' @importFrom lme4 isREML
##' @S3method isREML rlmerMod
isREML.rlmerMod <- function(x, ...) .isREML(x, ...)

##' @importFrom lme4 isLMM
##' @S3method isLMM rlmerMod
isLMM.rlmerMod <- function(x, ...) TRUE

##' @importFrom stats logLik
##' @S3method logLik rlmerMod
logLik.rlmerMod <- function(object, REML = NULL, ...)
   stop("log-likelihood is not defined for rlmerMod objects")

##' @importFrom stats model.frame
##' @S3method model.frame rlmerMod
model.frame.rlmerMod <- function(formula, fixed.only=FALSE, ...) {
    fr <- formula@frame
    if (fixed.only) {
        ff <- formula(formula,fixed.only=TRUE)
        ## thanks to Thomas Leeper and Roman Luštrik, Stack Overflow
        vars <- rownames(attr(terms.formula(ff), "factors"))
        fr <- fr[vars]
    }
    fr
}

##' @importFrom stats model.matrix
##' @S3method model.matrix rlmerMod
model.matrix.rlmerMod <- function(object, ...) object@pp$X

## we have our own nobs.rlmerMod method

## ranef function is in accessors.R

## no refit methods

## residuals is in accessors.R

## sigma in in accessors.R

## no simulate method

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

## update is in helpers.R

## ...

.prt.call <- function(call) {
    if (!is.null(cc <- call$formula))
	cat("Formula:", deparse(cc),"\n")
    if (!is.null(cc <- call$data))
	cat("   Data:", deparse(cc), "\n")
    if (!is.null(cc <- call$subset))
	cat(" Subset:", deparse(asOneSidedFormula(cc)[[2]]),"\n")
}

.prt.VC <- function(varcor, digits, comp, ...) {
    cat("Random effects:\n")
    fVC <- if(missing(comp))
	formatVC(varcor, digits=digits)
    else
	formatVC(varcor, digits=digits, comp=comp)
    print(fVC, quote = FALSE, digits = digits, ...)
}

.prt.grps <- function(ngrps, nobs) {
    cat(sprintf("Number of obs: %d, groups: ", nobs))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
}

## .printRlmerMod is in helpers.R

## print.rlmerMod is in helpers.R

##' @exportMethod show
setMethod("show", "rlmerMod", function(object) print.rlmerMod(object))

## print.summary.rlmerMod is in helpers.R

## can import tnames if required

## getME is in helpers.R

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
##' @S3method vcov summary.rlmerMod
vcov.summary.rlmerMod <- function(object, correlation = TRUE, ...) {
    if(is.null(object$vcov)) stop("logic error in summary of rlmerMod object")
    object$vcov
}

##' @importFrom lme4 mkVarCorr

## Extract variance and correlation components
##
## This function calculates the estimated variances, standard deviations, and
## correlations between the random-effects terms in a mixed-effects model, of
## class \code{\linkS4class{merMod}} (linear, generalized or nonlinear).  The
## within-group error variance and standard deviation are also calculated.
##
## @name VarCorr
## @aliases VarCorr VarCorr.merMod
## @param x a fitted model object, usually an object inheriting from class
## \code{\linkS4class{merMod}}.
## @param sigma an optional numeric value used as a multiplier for the standard
## deviations.  Default is \code{1}.
## @param rdig an optional integer value specifying the number of digits used
## to represent correlation estimates.  Default is \code{3}.
## @return a list of matrices, one for each random effects grouping term.
## For each grouping term, the standard deviations and correlation matrices for each grouping term
## are stored as attributes \code{"stddev"} and \code{"correlation"}, respectively, of the
## variance-covariance matrix, and
## the residual standard deviation is stored as attribute \code{"sc"}
## (for \code{glmer} fits, this attribute stores the scale parameter of the model).
## @author This is modeled after \code{\link[nlme]{VarCorr}} from package
## \pkg{nlme}, by Jose Pinheiro and Douglas Bates.
## @seealso \code{\link{lmer}}, \code{\link{nlmer}}
## @examples
## data(Orthodont, package="nlme")
## fm1 <- lmer(distance ~ age + (age|Subject), data = Orthodont)
## VarCorr(fm1)
## @keywords models
##' @importFrom nlme VarCorr
## @export VarCorr
##' @method VarCorr rlmerMod
##' @export
VarCorr.rlmerMod <- function(x, sigma, rdig)# <- 3 args from nlme
{
  ## FIXME:: would like to fix nlme to add ...
  ## FIXME:: add type=c("varcov","sdcorr","logs" ?)
    if (is.null(cnms <- x@cnms))
	stop("VarCorr methods require reTrms, not just reModule")
    if(missing(sigma)) # "bug": fails via default 'sigma=sigma(x)'
	sigma <- lme4::sigma(x)  ## FIXME: do we still need lme4:: ?
    nc <- vapply(cnms, length, 1L) # no. of columns per term
    structure(mkVarCorr(sigma, cnms=cnms, nc=nc, theta = x@theta,
			nms = {fl <- x@flist; names(fl)[attr(fl, "assign")]}),
	      useSc = as.logical(x@devcomp$dims["useSc"]),
	      class = "VarCorr.rlmerMod")
}

##' @S3method VarCorr summary.rlmerMod
VarCorr.summary.rlmerMod <- function(x, ...) x$varcor

##' @S3method print VarCorr.rlmerMod
print.VarCorr.rlmerMod <- function(x, digits = max(3, getOption("digits") - 2),
		   comp = "Std.Dev.", ...)
    print(formatVC(x, digits=digits, comp=comp), quote=FALSE, ...)

## __NOT YET EXPORTED__
## "format()" the 'VarCorr' matrix of the random effects -- for
## print()ing and show()ing
##
## @title Format the 'VarCorr' Matrix of Random Effects
## @param varc a \code{\link{VarCorr}} (-like) matrix with attributes.
## @param digits the number of significant digits.
## @param comp character vector of length one or two indicating which
## columns out of "Variance" and "Std.Dev." should be shown in the
## formatted output.
## @return a character matrix of formatted VarCorr entries from \code{varc}.
formatVC <- function(varc, digits = max(3, getOption("digits") - 2),
		     comp = "Std.Dev.")
{
    c.nms <- c("Groups", "Name", "Variance", "Std.Dev.")
    avail.c <- c.nms[-(1:2)]
    if(any(is.na(mcc <- pmatch(comp, avail.c))))
	stop("Illegal 'comp': ", comp[is.na(mcc)])
    nc <- length(colnms <- c(c.nms[1:2], (use.c <- avail.c[mcc])))
    if(length(use.c) == 0)
	stop("Must *either* show variances or standard deviations")
    useScale <- attr(varc, "useSc")
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"),
		  if(useScale) list(Residual = unname(attr(varc, "sc"))))
    reLens <- vapply(reStdDev, length, 1L)
    nr <- sum(reLens)
    reMat <- array('', c(nr, nc), list(rep.int('', nr), colnms))
    reMat[1+cumsum(reLens)-reLens, "Groups"] <- names(reLens)
    reMat[,"Name"] <- c(unlist(lapply(varc, colnames)), if(useScale) "")
    if(any("Variance" == use.c))
    reMat[,"Variance"] <- format(unlist(reStdDev)^2, digits = digits)
    if(any("Std.Dev." == use.c))
    reMat[,"Std.Dev."] <- format(unlist(reStdDev),   digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x) {
			       x <- as(x, "matrix")
			       dig <- max(2, digits - 2) # use 'digits' !
			       cc <- format(round(x, dig), nsmall = dig)
			       cc[!lower.tri(cc)] <- ""
			       nr <- nrow(cc)
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }))[, -maxlen, drop = FALSE]
	if (nrow(corr) < nrow(reMat))
	    corr <- rbind(corr, matrix("", nrow(reMat) - nrow(corr), ncol(corr)))
	colnames(corr) <- c("Corr", rep.int("", max(0L, ncol(corr)-1L)))
	cbind(reMat, corr)
    } else reMat
}

## summary is in helpers.R

## plots are in plots.R

##' @importFrom stats weights
##' @S3method weights rlmerMod
weights.rlmerMod <- function(object, ...) {
  object@resp$weights
}

## no optimizer options in rlmer
