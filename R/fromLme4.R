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
coefMer <- getS3method("coef", "merMod")
##' @export
coef.rlmerMod <- function(object, ...) {
    val <- coefMer(object, ...)
    class(val) <- "coef.rlmerMod"
    val
}

## deviance is in accessors.R

## FIXME what about drop1

##' @importFrom stats extractAIC
##' @export
extractAIC.rlmerMod <- function(fit, scale = 0, k = 2, ...)
   stop("AIC is not defined for rlmerMod objects")

##' @importFrom stats family
##' @export
family.rlmerMod <- function(object, ...) gaussian()

##' @importFrom stats fitted
##' @export
fitted.rlmerMod <- getS3method("fitted", "merMod")

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
##' @export
fixef.rlmerMod <- function(object, ...)
    structure(object@beta, names = dimnames(.X(object))[[2]])

##' @importFrom reformulas nobars
getFixedFormula <- function(form) {
    form[[3]] <- if (is.null(nb <- nobars(form[[3]]))) 1 else nb
    form
}

##' @importFrom stats formula
##' @export
formula.rlmerMod <- getS3method("formula", "merMod")

##' @importFrom lme4 isREML
##' @export
isREML.rlmerMod <- function(x, ...) .isREML(x, ...)

## needed for predict():
##' @importFrom lme4 isGLMM
##' @export
isGLMM.rlmerMod <- function(x, ...) FALSE

##' @importFrom lme4 isLMM
##' @export
isLMM.rlmerMod <- function(x, ...) TRUE

##' @importFrom lme4 isNLMM
##' @export
isNLMM.rlmerMod <- function(x, ...) FALSE

##' @importFrom stats logLik
##' @export
logLik.rlmerMod <- function(object, REML = NULL, ...)
   stop("log-likelihood is not defined for rlmerMod objects")

##' @importFrom stats model.frame
##' @export
model.frame.rlmerMod <- getS3method("model.frame", "merMod")

##' @importFrom stats model.matrix
##' @export
model.matrix.rlmerMod <- function(object, ...) .X(object)

## we have our own nobs.rlmerMod method

## ranef function is in accessors.R

## no refit methods

##' @importFrom lme4 REMLcrit
## don't add a function for now.

## residuals is in accessors.R

## sigma in in accessors.R

## no simulate method

##' @importFrom stats terms
##' @export
terms.rlmerMod <- getS3method("terms", "merMod")

## update is in helpers.R

## ...

.prt.resids <- function(resids, digits, title="Scaled residuals:", ...) {
    cat(title,"\n")
    rq <- setNames(zapsmall(quantile(resids), digits + 1L),
                   c("Min", "1Q", "Median", "3Q", "Max"))
    print(rq, digits=digits, ...)
    cat("\n")
}

.prt.call <- function(call, long=TRUE) {
    if (!is.null(cc <- call$formula))
	cat("Formula:", deparse(cc),"\n")
    if (!is.null(cc <- call$data))
	cat("   Data:", deparse(cc), "\n")
    if (!is.null(cc <- call$weights))
        cat("Weights:", deparse(cc), "\n")
    if (!is.null(cc <- call$offset))
        cat(" Offset:", deparse(cc), "\n")
    if (long && length(cc <- call$control) &&
	!identical((dc <- deparse(cc)), "lmerControl()"))
	## && !identical(eval(cc), lmerControl()))
	cat("Control:", dc, "\n")
    if (!is.null(cc <- call$subset))
	cat(" Subset:", deparse(asOneSidedFormula(cc)[[2]]),"\n")
}

.prt.VC <- function(varcor, digits, comp, formatter=format, ...) {
    cat("Random effects:\n")
    fVC <- if(missing(comp))
	formatVC(varcor, digits=digits, formatter=formatter)
    else
	formatVC(varcor, digits=digits, formatter=formatter, comp=comp)
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

## getME is in accessors.R

globalVariables("forceSymmetric", add=TRUE)

##' Variance-covariance matrix of the fixed effects of an rlmerMod fit.
##'
##' By default returns the same object as lme4's \code{vcov.merMod} (the
##' linearised model-based covariance). With \code{type = "sandwich"},
##' returns the robust cluster-sandwich \code{\link{vcov_sandwich}}.
##'
##' The cluster sandwich is exact for a single (nested) grouping factor;
##' for crossed factors it is approximate (a warning is issued, see
##' \code{\link{vcov_sandwich}}).
##'
##' \strong{Small-J caveat for the sandwich.} The G1 correction is
##' necessary but not sufficient at very small \eqn{J}: in a
##' simulation study sandwich CI coverage drops to ~0.89 at \eqn{J =
##'   8} (vs. nominal 0.95). \code{vcov_sandwich} emits a warning for
##' \eqn{J < 20}; for inference at small \eqn{J} prefer \code{type =
##'   "default"} or pair the sandwich CI with a bootstrap calibration
##' (e.g. \code{confintROB}).
##'
##' @param object An \code{rlmerMod} object.
##' @param type \code{"default"} (the lme4-inherited linearised vcov; the
##'   pre-existing behaviour) or \code{"sandwich"} (the robust
##'   cluster-sandwich).
##' @param cluster When \code{type = "sandwich"}, passed to
##'   \code{\link{vcov_sandwich}}: \code{NULL} (auto-detect for a single
##'   grouping factor), a character string naming a grouping factor of the
##'   model, or a length-\code{n} vector of cluster memberships.
##' @param correction When \code{type = "sandwich"}, \code{"G1"} (default,
##'   applies the \eqn{J/(J-1)} small-sample scaling) or \code{"none"}.
##' @param ... Additional arguments passed to the default \code{vcov}
##'   method (only used when \code{type = "default"}).
##' @return A \eqn{p \times p} covariance matrix for \eqn{\hat{\beta}}.
##' @seealso \code{\link{vcov_sandwich}}
##' @importFrom stats vcov
##' @export
vcov.rlmerMod <- function(object,
                          type = c("default", "sandwich"),
                          cluster = NULL,
                          correction = c("G1", "none"),
                          ...) {
    type <- match.arg(type)
    if (type == "default") {
        return(.vcov_merMod_default(object, ...))
    }
    correction <- match.arg(correction)
    vcov_sandwich(object, cluster = cluster, correction = correction)
}

## Cached lme4 merMod vcov method, used for type = "default" so the
## byte-for-byte pre-existing behaviour is preserved.
.vcov_merMod_default <- getS3method("vcov", "merMod")

##' @importFrom stats vcov
##' @export
vcov.summary.rlmerMod <- function(object, correlation = TRUE, ...) {
    if(is.null(object$vcov)) stop("logic error in summary of rlmerMod object")
    object$vcov
}

##' @importFrom nlme VarCorr
##' @method VarCorr rlmerMod
##' @export
VarCorrMer <- getS3method("VarCorr", "merMod")
VarCorr.rlmerMod <- function(x, ...)# <- 3 args from nlme
{
    val <- VarCorrMer(x, ...)
    class(val) <- "VarCorr.rlmerMod"
    val
}

##' @export
VarCorr.summary.rlmerMod <- function(x, ...) x$varcor

##' @export
print.VarCorr.rlmerMod <- getS3method("print", "VarCorr.merMod")

##' @export
as.data.frame.VarCorr.rlmerMod <-
    getS3method("as.data.frame", "VarCorr.merMod")

## __NOT YET EXPORTED__
## "format()" the 'VarCorr' matrix of the random effects -- for
## print()ing and show()ing
##
## @title Format the 'VarCorr' Matrix of Random Effects
## @param varc a \code{\link[lme4]{VarCorr}} (-like) matrix with attributes.
## @param digits the number of significant digits.
## @param comp character vector of length one or two indicating which
## columns out of "Variance" and "Std.Dev." should be shown in the
## formatted output.
## @param formatter the \code{\link{function}} to be used for
## formatting the standard deviations and or variances (but
## \emph{not} the correlations which (currently) are always formatted
## as "0.nnn"
## @param ... optional arguments for \code{formatter(*)} in addition
## to the first (numeric vector) and \code{digits}.
## @return a character matrix of formatted VarCorr entries from \code{varc}.
formatVC <- function(varcor, digits = max(3, getOption("digits") - 2),
                     comp = "Std.Dev.", formatter = format,
                     useScale = attr(varcor, "useSc"),
                     ...)
{
  c.nms <- c("Groups", "Name", "Variance", "Std.Dev.")
  avail.c <- c.nms[-(1:2)]
  if(anyNA(mcc <- pmatch(comp, avail.c)))
    stop("Illegal 'comp': ", comp[is.na(mcc)])
  nc <- length(colnms <- c(c.nms[1:2], (use.c <- avail.c[mcc])))
  if(length(use.c) == 0)
    stop("Must *either* show variances or standard deviations")
  reStdDev <- c(lapply(varcor, attr, "stddev"),
                if(useScale) list(Residual = unname(attr(varcor, "sc"))))
  reLens <- lengths(reStdDev)
  nr <- sum(reLens)
  reMat <- array('', c(nr, nc), list(rep.int('', nr), colnms))
  reMat[1+cumsum(reLens)-reLens, "Groups"] <- names(reLens)
  reMat[,"Name"] <- c(unlist(lapply(varcor, colnames)), if(useScale) "")
  if(any("Variance" == use.c))
    reMat[,"Variance"] <- formatter(unlist(reStdDev)^2, digits = digits, ...)
  if(any("Std.Dev." == use.c))
    reMat[,"Std.Dev."] <- formatter(unlist(reStdDev),   digits = digits, ...)
  if (any(reLens > 1)) {
    maxlen <- max(reLens)
    recorr <- lapply(varcor, attr, "correlation")
    corr <-
      do.call(rbind,
              lapply(recorr,
                     function(x) {
                       x <- as.matrix(x)
                       dig <- max(2, digits - 2) # use 'digits' !
                       ## not using formatter() for correlations
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
##' @export
weights.rlmerMod <- function(object, ...) {
  object@resp$weights
}

## no optimizer options in rlmer

#######################################################
## predict method                                    ##
#######################################################

reFormHack <- function(re.form,ReForm,REForm,REform) {
    if (!missing(ReForm)) {
        message(shQuote("re.form")," is now preferred to ",shQuote("ReForm"))
        return(ReForm)
    }
    if (!missing(REForm)) {
        message(shQuote("re.form")," is now preferred to ",shQuote("REForm"))
        return(REForm)
    }
    if (!missing(REform)) {
        message(shQuote("re.form")," is now preferred to ",shQuote("REform"))
        return(REform)
    }
    re.form
}

##' Predictions and confidence/prediction intervals for an rlmerMod fit.
##'
##' By default (\code{interval = "none"}) returns the same numeric vector
##' as the lme4-style point prediction (preserved byte-for-byte from the
##' pre-existing predict method). With \code{interval = "confidence"} or
##' \code{"prediction"}, returns a data frame with columns \code{fit},
##' \code{lwr}, \code{upr}, \code{se}.
##'
##' The fixed-effect contribution to the SE uses the robust cluster
##' sandwich \code{vcov(object, type = "sandwich")}; the random-effects
##' contribution (when REs are part of the prediction) is computed from
##' the partial influence function of \eqn{\hat{u}} (the local-shift
##' sensitivity). For \code{interval = "prediction"} the additional
##' residual variance \eqn{\hat{\sigma}^2} is added. Intervals are
##' \code{fit +/- z * se} with \eqn{z = \Phi^{-1}((1 + level)/2)}; for
##' small numbers of clusters the normal approximation may under-cover,
##' and a bootstrap CI (via \code{confint(..., method = "boot")}) is
##' preferable.
##'
##' @param object An \code{rlmerMod} object.
##' @param newdata,re.form,ReForm,REForm,REform,terms,type,allow.new.levels,na.action,...
##'   See the lme4 \code{\link[lme4]{predict.merMod}} documentation; these
##'   arguments are forwarded unchanged to the point-prediction code path.
##' @param interval One of \code{"none"} (default, returns a numeric
##'   vector for backwards compatibility), \code{"confidence"}
##'   (fit-uncertainty interval for \eqn{X\hat{\beta} + Z\hat{b}}), or
##'   \code{"prediction"} (adds residual variance \eqn{\hat{\sigma}^2}).
##' @param level Coverage level for the interval; default 0.95.
##' @return A numeric vector (when \code{interval = "none"}) or a data
##'   frame with columns \code{fit}, \code{lwr}, \code{upr}, \code{se}.
##' @seealso \code{\link[=vcov.rlmerMod]{vcov}}
##' @importFrom stats predict qnorm
##' @importFrom reformulas mkReTrms findbars
##' @export
## FIXME: the point-prediction path .predictRlmerPoint is slightly
##        modified from lme4 version 1.0-6. Newer versions have an
##        improved version.
predict.rlmerMod <- function(object, newdata=NULL,
                             re.form=NULL,
                             ReForm,
                             REForm,
                             REform,
                             terms=NULL, type=c("link","response"),
                             allow.new.levels=FALSE, na.action=na.pass,
                             interval = c("none", "confidence", "prediction"),
                             level = 0.95,
                             ...) {
    interval <- match.arg(interval)
    ## Forward to .predictRlmerPoint preserving the user's missingness
    ## pattern -- otherwise unconditionally passing type = type would
    ## trigger the "type argument ignored" warning even when the caller
    ## did not supply type. Only forward args the caller actually
    ## supplied.
    mc <- match.call(expand.dots = TRUE)
    user_args <- as.list(mc)[-1L]
    user_args[c("interval", "level")] <- NULL
    if (interval == "none")
        return(do.call(.predictRlmerPoint, user_args, envir = parent.frame()))
    pred <- do.call(.predictRlmerPoint, user_args, envir = parent.frame())
    REform_eff <- reFormHack(re.form, ReForm, REForm, REform)
    re_included <- is.null(REform_eff) ||
                   (is.language(REform_eff) && length(all.vars(REform_eff)) > 0L)
    if (is.null(newdata)) {
        X_new <- getME(object, "X")
        Z_new <- if (re_included) as.matrix(t(object@pp$Zt)) else NULL
    } else {
        RHS   <- formula(object, fixed.only = TRUE)[-2]
        Terms <- terms(object, fixed.only = TRUE)
        X_new <- model.matrix(RHS,
                              model.frame(delete.response(Terms), newdata,
                                          na.action = na.action),
                              contrasts.arg = attr(getME(object, "X"),
                                                    "contrasts"))
        Z_new <- if (re_included) .reconstruct_Z_new(object, newdata) else NULL
    }
    se <- .predictSE(object, X_new = X_new, Z_new = Z_new,
                     prediction = (interval == "prediction"))
    if (length(se) != length(pred)) {
        ## NA handling in the point path may have dropped rows; align by
        ## recycling SE = NA where pred is NA, otherwise abort.
        if (sum(!is.na(pred)) == length(se)) {
            se_full <- rep(NA_real_, length(pred))
            se_full[!is.na(pred)] <- se
            se <- se_full
        } else {
            stop("predict.rlmerMod: cannot align SE with point predictions ",
                 "(length ", length(se), " vs ", length(pred), ").")
        }
    }
    z <- stats::qnorm(0.5 + level / 2)
    data.frame(fit = as.numeric(pred),
               lwr = as.numeric(pred) - z * se,
               upr = as.numeric(pred) + z * se,
               se  = se,
               row.names = if (!is.null(names(pred))) names(pred) else NULL)
}

## The original lme4-style point-prediction body, unchanged. interval =
## "none" delegates straight to this so the byte-for-byte vector return
## is preserved.
.predictRlmerPoint <- function(object, newdata=NULL,
                               re.form=NULL,
                               ReForm,
                               REForm,
                               REform,
                               terms=NULL, type=c("link","response"),
                               allow.new.levels=FALSE, na.action=na.pass,
                               ...) {
    ## FIXME: appropriate names for result vector?
    ## FIXME: make sure behaviour is entirely well-defined for NA in grouping factors

    REform <- reFormHack(re.form,ReForm,REForm,REform)

    if (length(list(...)>0)) warning("unused arguments ignored")
    if (isLMM(object) && !missing(type)) warning("type argument ignored for linear mixed models")
    fit.na.action <- attr(object@frame,"na.action")
    type <- match.arg(type)
    if (!is.null(terms)) stop("terms functionality for predict not yet implemented")
    ## FIXME/WARNING: how do we/can we do this in an eval-safe way???
    form_orig <- formula(object)
    if (is.null(newdata) && is.null(REform)) {
        ## raw predict() call, just return fitted values (inverse-link if appropriate)
        if (isLMM(object) || isNLMM(object)) {
            pred <- fitted(object)
        } else { ## inverse-link
            pred <-  switch(type,response=object@resp$mu, ## fitted(object),
                            link=object@resp$eta)  ## fixme: getME() ?
        }
        if (!is.null(fit.na.action)) {
            pred <- napredict(fit.na.action,pred)
        }
        return(pred)

      } else { ## newdata and/or REform specified
        X_orig <- getME(object, "X")
        if (is.null(newdata)) {
            X <- X_orig
        } else {
            ## evaluate new fixed effect
            RHS <- formula(object,fixed.only=TRUE)[-2]
            Terms <- terms(object,fixed.only=TRUE)
            X <- model.matrix(RHS, mfnew <- model.frame(delete.response(Terms),
                                                        newdata, na.action=na.action),
                              contrasts.arg=attr(X_orig,"contrasts"))
        }
        pred <- drop(X %*% fixef(object))
        ## modified from predict.glm ...
        offset <- rep(0, nrow(X))
        tt <- terms(object)
        ## FIXME:: need to unname()  ?
        if (!is.null(off.num <- attr(tt, "offset"))) {
            for (i in off.num) offset <- offset + eval(attr(tt,"variables")[[i + 1]], newdata)
        }
        ## FIXME: is this redundant??
        if (!is.null(frOffset <- attr(object@frame,"offset")))
            offset <- offset + eval(frOffset, newdata)
        pred <- pred+offset
        if (is.null(REform)) {
            REform <- form_orig[-2]
        }
        ## FIXME: ??? can't apply is.na() to a 'language' object?
        ##  what's the appropriate test?
        if (is.language(REform)) {
            na.action.name <- deparse(match.call()$na.action) ## ugh
            if (!is.null(newdata) && na.action.name %in% c("na.exclude","na.omit")) {
                ## strip NAs from data for random-effects matrix construction
                if (length(nadrop <- attr(mfnew,"na.action"))>0) {
                    newdata <- newdata[-nadrop,]
                }
            }
            re <- ranef(object)
            ## ok? -- newdata used even though it was just tested for null
            if(is.null(newdata)) rfd <- object@frame else rfd <- newdata # get data for REform
            ReTrms <- mkReTrms(findbars(REform[[2]]), rfd)
            if (!allow.new.levels && any(sapply(ReTrms$flist,function(x) any(is.na(x)))))
                stop("NAs are not allowed in prediction data for grouping variables unless allow.new.levels is TRUE")
            unames <- unique(sort(names(ReTrms$cnms)))  ## FIXME: same as names(ReTrms$flist) ?
            ## convert numeric grouping variables to factors as necessary
            for (i in all.vars(REform[[2]])) {
                newdata[[i]] <- factor(newdata[[i]])
            }
            Rfacs <- setNames(lapply(unames,function(x) factor(eval(parse(text=x),envir=newdata))),
                              unames)
            new_levels <- lapply(Rfacs,function(x) levels(droplevels(factor(x))))
            ## FIXME: should this be unique(as.character(x)) instead?
            ##   (i.e., what is the proper way to protect against non-factors?)
            levelfun <- function(x,n) {
                ## find and deal with new levels
                if (any(!new_levels[[n]] %in% rownames(x))) {
                    if (!allow.new.levels) stop("new levels detected in newdata")
                    ## create an all-zero data frame corresponding to the new set of levels ...
                    newx <- as.data.frame(matrix(0,nrow=length(new_levels[[n]]),ncol=ncol(x),
                                                 dimnames=list(new_levels[[n]],names(x))))
                    ## then paste in the matching RE values from the original fit/set of levels
                    newx[rownames(x),] <- x
                    x <- newx
                }
                ## find and deal with missing old levels
                if (any(!rownames(x) %in% new_levels[[n]])) {
                    x <- x[rownames(x) %in% new_levels[[n]],,drop=FALSE]
                }
                x
            }
            ## fill in/delete levels as appropriate
            re_x <- mapply(levelfun,re,names(re),SIMPLIFY=FALSE)
            re_new <- list()
            if (any(!names(ReTrms$cnms) %in% names(re)))
                stop("grouping factors specified in REform that were not present in original model")
            ## pick out random effects values that correspond to
            ##  random effects incorporated in REform ...
            for (i in seq_along(ReTrms$cnms)) {
                rname <- names(ReTrms$cnms)[i]
                if (any(!ReTrms$cnms[[rname]] %in% names(re[[rname]])))
                    stop("random effects specified in REform that were not present in original model")
                re_new[[i]] <- re_x[[rname]][,ReTrms$cnms[[rname]]]
            }
            re_newvec <- unlist(lapply(re_new,t))  ## must TRANSPOSE RE matrices before unlisting
            if(!is.null(newdata)) pred <- pred + drop(as.matrix(re_newvec %*% ReTrms$Zt))
        } ## predictions with REform!=NA
        if (isGLMM(object) && type=="response") {
            pred <- object@resp$family$linkinv(pred)
        }
        ## fill in NAs as appropriate
        if (is.null(newdata) && !is.null(fit.na.action)) {
            pred <- napredict(fit.na.action,pred)
        } else {
            pred <- napredict(na.action,pred)
        }
        return(pred)
    }
}

forceCopy <- function(x) deepcopy(x)

# these functions pick up where findbars leaves off, in terms of sugar
# @param REtrm an element of the result of findbars
# @param REtrms the result of findbars
# @return \code{reexpr} gives a one-sided formula with the linear
# model formula for the raw model matrix. \code{grpfact} gives an
# expression with the name of the grouping factor associated with
# the raw model matrix. \code{termnms} gives a character vector with
# the names of the random effects terms.
reexpr <- function(REtrm) substitute( ~ foo, list(foo = REtrm[[2]]))
grpfact <- function(REtrm) substitute(factor(fac), list(fac = REtrm[[3]]))
termnms <- function(REtrms) vapply(REtrms, deparse1, "")

# mmList(): list of model matrices
# ------    called from getME()
mmList <- function(object, ...) UseMethod("mmList")
mmList.rlmerMod <- function(object, ...) mmList(formula(object), model.frame(object))
mmList.formula <- function(object, frame, ...) {
    bars <- findbars(object)
    mm <- setNames(lapply(bars, function(b) model.matrix(eval(reexpr(b), frame), frame)),
                   termnms(bars))
    grp <- lapply(lapply(bars, grpfact), eval, frame)
    nl <- vapply(grp, nlevels, 1L)
    if (any(diff(nl) > 0))
        mm[order(nl, decreasing=TRUE)]
    else
        mm
}
