##' The \code{lmerNoFit} function can be used to get trivial starting values.
##' This is mainly used to verify the algorithms to reproduce the fit by
##' \code{\link[lme4]{lmer}} when starting from trivial initial values.
##'
##' @rdname rlmer
##' @param initTheta parameter to initialize theta with (optional)
##' @examples
##' \dontrun{
##'   ## start from lmer's initial estimate, not its fit
##'   rlmer(Yield ~ (1|Batch), Dyestuff, init = lmerNoFit)
##' }
##' @importFrom lme4 mkMerMod
##' @export
lmerNoFit <- function(formula, data = NULL, ..., initTheta) {
    if (packageVersion("lme4") >= "0.99999911.0") {
        mc <- match.call()
        mc[[1]] <- quote(lmer)
        ## trick it...
        fm <- lmer(formula, data, ...)
        devfun <- lmer(formula, data, ..., devFunOnly=TRUE)
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
        mkMerMod(environment(devfun), fakeOpt, reTrms, fr, mc)
    } else {
        ## not supported
        warning("lmerNoFit not implemented for this version of lme4, returning regular lmer fit")
        lmer(...)
    }
}

