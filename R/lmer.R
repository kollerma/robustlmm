## Define all imports for roxygen2 here

##' This functions replaces the deprecated option \code{doFit = FALSE}
##' in lme4.
##'
##' @title Construct lmerMod object but do not fit the data.
##' @param ... passed to lmer()
##' @param initTheta parameter to initialize theta with (optional)
##' @return lmerMod object
##' @importMethodsFrom Matrix diag solve determinant t crossprod tcrossprod as.vector drop rowSums rowMeans colSums colMeans chol which
##' @examples
##'   rlmer(Yield ~ (1|Batch), Dyestuff, init = lmerNoFit)
##' @export
lmerNoFit <- function(..., initTheta) {
    if (packageVersion("lme4") >= "0.99999911.0") {
        mc <- match.call()
        mc[[1]] <- quote(lmer)
        ## trick it...
        fm <- lmer(...)
        devfun <- lmer(..., devFunOnly=TRUE)
        if (!missing(initTheta)) devfun(initTheta)
        fakeOpt <- list(fval=devfun(environment(devfun)$pp$theta),
                        conv=1,
                        par=environment(devfun)$pp$theta,
                        ierr=0)
        reTrms <- list(flist=fm@flist,
                       cnms=fm@cnms,
                       Gp=fm@Gp,
                       lower=fm@lower)
        fr <- fm@frame
        ##mc <- fm@call
        lme4:::mkMerMod(environment(devfun), fakeOpt, reTrms, fr, mc)
    } else {
        ## not supported
        warning("lmerNoFit not implemented for this version of lme4, returning regular lmer fit")
        lmer(...)
    }
}

