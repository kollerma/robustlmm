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

## use getS3method to copy methods from lme4
##' @importFrom utils getS3method

##' @importFrom stats coef
##' @S3method coef rlmerMod
coefMer <- getS3method("coef", "merMod")
coef.rlmerMod <- function(object, ...) {
    val <- coefMer(object, ...)
    class(val) <- "coef.rlmerMod"
    val
}

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
formula.rlmerMod <- getS3method("formula", "merMod")

##' @importFrom lme4 isREML
##' @S3method isREML rlmerMod
isREML.rlmerMod <- function(x, ...) .isREML(x, ...)

## needed for predict():
##' @importFrom lme4 isGLMM
##' @S3method isGLMM rlmerMod
isGLMM.rlmerMod <- function(x, ...) FALSE

##' @importFrom lme4 isLMM
##' @S3method isLMM rlmerMod
isLMM.rlmerMod <- function(x, ...) TRUE

##' @importFrom lme4 isNLMM
##' @S3method isNLMM rlmerMod
isNLMM.rlmerMod <- function(x, ...) FALSE

##' @importFrom stats logLik
##' @S3method logLik rlmerMod
logLik.rlmerMod <- function(object, REML = NULL, ...)
   stop("log-likelihood is not defined for rlmerMod objects")

##' @importFrom stats model.frame
##' @S3method model.frame rlmerMod
model.frame.rlmerMod <- getS3method("model.frame", "merMod")

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
terms.rlmerMod <- getS3method("terms", "merMod")

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

##' @importFrom stats vcov
##' @S3method vcov rlmerMod
vcov.rlmerMod <- getS3method("vcov", "merMod")

##' @importFrom stats vcov
##' @S3method vcov summary.rlmerMod
vcov.summary.rlmerMod <- function(object, correlation = TRUE, ...) {
    if(is.null(object$vcov)) stop("logic error in summary of rlmerMod object")
    object$vcov
}

##' @importFrom nlme VarCorr
##' @method VarCorr rlmerMod
##' @export
VarCorrMer <- getS3method("VarCorr", "merMod")
VarCorr.rlmerMod <- function(x, sigma, rdig)# <- 3 args from nlme
{
    val <- VarCorrMer(x, sigma, rdig)
    class(val) <- "VarCorr.rlmerMod"
    val
}

##' @S3method VarCorr summary.rlmerMod
VarCorr.summary.rlmerMod <- function(x, ...) x$varcor

##' @S3method print VarCorr.rlmerMod
print.VarCorr.rlmerMod <- getS3method("print", "VarCorr.merMod")

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

#######################################################
## predict method                                    ##
#######################################################

##' @importFrom stats predict
##' @S3method predict rlmerMod
predict.rlmerMod <- getS3method("predict", "merMod")
## the following is needed to get the correct getME() function:
environment(predict.rlmerMod) <- environment()
