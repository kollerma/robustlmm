#' @export
recover_data.rlmerMod <- function(object, ...) {
    fcall <- object@call
    emmeans::recover_data(fcall, delete.response(terms(object)),
                          attr(object@frame, "na.action"),
                          frame = object@frame, ...)
}

#' @export
emm_basis.rlmerMod <- function(object, trms, xlev, grid, vcov., ...) {
    useSatterthwaite <- missing(vcov.)
    if (missing(vcov.)) {
        V <- as.matrix(vcov(object, correlation = FALSE))
    } else {
        V <- as.matrix(emmeans::.my.vcov(object, vcov.))
    }

    ## WS16 step 3: report a finite Satterthwaite df instead of the
    ## z-based Inf. Only for the default vcov (the IF-based
    ## .satterthwaite_df uses sigma^2 unsc(theta); a user-supplied
    ## vcov. would be inconsistent with it, so we keep Inf there).
    ##
    ## emmeans re-environments dffun to baseenv() (to avoid the closure
    ## capturing the model), so dffun cannot see robustlmm's namespace.
    ## Everything it needs therefore travels through dfargs: the fit,
    ## the two worker functions (which still carry their own namespace
    ## environment, so their bodies resolve internals correctly), and a
    ## cache environment. The expensive implicitIF_full is computed
    ## lazily on the first df request and memoised in that environment
    ## (emmeans calls dffun once per estimate). Any failure -- an
    ## unsupported design (> 1 grouping factor or Mallows eta, which
    ## .scoreByCluster rejects), a variance-component boundary, or a
    ## degenerate contrast -- falls back to Inf, i.e. the previous
    ## z-based behaviour.
    if (useSatterthwaite) {
        dfargs <- list(object = object,
                       implicitIF_full = implicitIF_full,
                       satterthwaite_df = .satterthwaite_df,
                       cache = new.env(parent = emptyenv()))
        dffun <- function(k, dfargs) {
            if (!isTRUE(any(k != 0, na.rm = TRUE))) return(Inf)
            cache <- dfargs$cache
            if (is.null(cache$IF) && !isTRUE(cache$failed)) {
                cache$IF <- tryCatch(
                    dfargs$implicitIF_full(dfargs$object),
                    error = function(e) NULL)
                if (is.null(cache$IF)) cache$failed <- TRUE
            }
            if (is.null(cache$IF)) return(Inf)
            res <- tryCatch(suppressWarnings(
                dfargs$satterthwaite_df(dfargs$object,
                                        Lmat = matrix(k, nrow = 1L),
                                        IF = cache$IF)),
                error = function(e) NA_real_)
            ## WS14: a reducible boundary (variance components at 0 but the
            ## remaining covariance regular) still yields a valid
            ## conditional df; only a genuinely singular fit -> Inf.
            if (length(res) != 1L || !is.finite(res) ||
                (isTRUE(attr(res, "boundary")) &&
                 !isTRUE(attr(res, "reducible"))))
                return(Inf)
            max(as.numeric(res), 1)
        }
    } else {
        dfargs <- list()
        dffun <- function(k, dfargs) Inf
    }
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

#' @importFrom Matrix Matrix.Version
.Matrix.Version.abi.on.build <- Matrix.Version()$abi
.onLoad <- function(libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE)) {
        emmeans::.emm_register("rlmerMod", pkgname)
    }
    if (.Matrix.Version.abi.on.build != Matrix.Version()$abi) {
        warning("Package robustlmm was built for another version of ",
                "the Matrix package. Please re-install robustlmm.",
                call. = FALSE, immediate. = TRUE)
    }
}
