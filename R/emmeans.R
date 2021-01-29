#' @export
recover_data.rlmerMod <- function(object, ...) {
    fcall <- object@call
    emmeans::recover_data(fcall, delete.response(terms(object)),
                          attr(object@frame, "na.action"),
                          frame = object@frame, ...)
}

#' @export
emm_basis.rlmerMod <- function(object, trms, xlev, grid, vcov., ...) {
    if (missing(vcov.)) {
        V <- as.matrix(vcov(object, correlation = FALSE))
    } else {
        V <- as.matrix(emmeans::.my.vcov(object, vcov.))
    }

    dfargs <- list()
    dffun <- function(k, dfargs) Inf
    contrasts <- attr(object@pp$X, "contrasts")
    m <- model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X <- model.matrix(trms, m, contrasts.arg = contrasts)
    bhat <- fixef(object)
    nbasis <- estimability::all.estble
    mm <- emmeans::.cmpMM(object@pp$X, object@pp$Xwts^2,
                          attr(object@pp$X, "assign"))

    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun,
         dfargs=dfargs, model.matrix = mm)
}

.onLoad <- function(libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE)) {
        emmeans::.emm_register("rlmerMod", pkgname)
    }
}
