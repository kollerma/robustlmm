## confint.rlmerMod: Wald intervals for the fixed effects (closed form
## from vcov), with method = "boot" / "BCa" delegated to
## confintROB::confintROB (Mason, Cantoni & Ghisletta 2021, 2024) so
## users get the published, peer-reviewed bootstrap (parametric or
## wild, with optional bias-correction-and-acceleration) rather than a
## hand-rolled percentile.

##' Confidence intervals for the fixed-effect coefficients of an rlmer fit.
##'
##' Two routes.
##' \describe{
##'   \item{\code{method = "Wald"}}{the default and the only path
##'     implemented in this package. Per-coefficient closed-form CI
##'     \eqn{\hat\beta_k \pm z_{1-\alpha/2} \cdot SE_k} with \eqn{V}
##'     taken from \code{vcov(object, type = vcov_type)} -- so the
##'     \code{"sandwich"} option carries through to the interval.}
##'   \item{\code{method = "boot"} or \code{"BCa"}}{a thin dispatch to
##'     \code{\link[confintROB]{confintROB}} (Mason, Cantoni &
##'     Ghisletta 2021, 2024) with \code{boot.type} forwarded. The
##'     wrapper subsets the returned matrix to the fixed-effect rows so
##'     the shape matches the \code{"Wald"} path; variance-component
##'     CIs from \code{confintROB} are dropped here. \code{vcov_type}
##'     is not honoured on these paths (\code{confintROB} uses its own
##'     internal covariance).}
##' }
##'
##' Guidance (Koller 2014; Mason et al. 2024). The chi-sq-p Wald limit
##' is adequate for \eqn{J \gtrsim 20} groups; the bootstrap earns its
##' (substantial) cost mainly at smaller \eqn{J}. \code{boot.type =
##'   "wild"} (the default, following \pkg{confintROB}) is robust to
##' misspecification of the response covariance, while
##' \code{"parametric"} is exact under the fitted central LMM.
##' \code{method = "BCa"} adds the bias-correction-and-acceleration
##' adjustment to the bootstrap percentile (preferred when the
##' bootstrap distribution is skewed).
##'
##' \strong{Small-J caveat for the sandwich.} \code{vcov_type =
##'   "sandwich"} under-covers at \eqn{J < 20} (in a simulation study:
##' coverage ~0.89 at \eqn{J = 8} vs. nominal 0.95); the sandwich path
##' emits a warning. For Wald CIs at small \eqn{J} prefer \code{vcov_type
##'   = "default"} or use \code{method = "boot"} / \code{"BCa"}.
##'
##' @param object An \code{rlmerMod} object.
##' @param parm Either \code{NULL} (all fixed-effect coefficients), an
##'   integer vector of coefficient indices, or a character vector of
##'   coefficient names.
##' @param level Coverage level; default 0.95.
##' @param method One of \code{"Wald"} (default; closed form),
##'   \code{"boot"} (bootstrap percentile via \pkg{confintROB}), or
##'   \code{"BCa"} (bias-corrected bootstrap via \pkg{confintROB}).
##' @param vcov_type Covariance to use for \eqn{V} when \code{method =
##'   "Wald"}: \code{"default"} (the linearised lme4 vcov; pre-existing
##'   behaviour) or \code{"sandwich"} (the robust cluster-sandwich
##'   \code{\link{vcov_sandwich}}; exact for a single nested grouping
##'   factor, approximate for crossed designs). Ignored when
##'   \code{method = "boot"} or \code{"BCa"}.
##' @param df Critical-value degrees of freedom for \code{method =
##'   "Wald"}. \code{"none"} (default) uses the normal quantile
##'   \eqn{z_{1-\alpha/2}} as before; \code{"satterthwaite"} uses a
##'   per-coefficient Satterthwaite t-quantile (the robust IF-based
##'   df), matching \code{summary(object, df = "satterthwaite")}. It
##'   requires \code{vcov_type = "default"} and a single grouping factor;
##'   otherwise it warns and falls back to the normal quantile.
##' @param boot.type Bootstrap kind passed to \code{confintROB} when
##'   \code{method = "boot"} or \code{"BCa"}: one of
##'   \code{"wild"} (default, the \code{confintROB} recommendation) or
##'   \code{"parametric"}.
##' @param nsim Bootstrap replicates; default 1000.
##' @param seed Optional RNG seed for reproducibility of the bootstrap.
##' @param ... Additional arguments forwarded to \code{confintROB}
##'   (e.g. \code{clusterID} for the wild bootstrap).
##' @return A 2-column matrix with one row per selected fixed-effect
##'   coefficient, columns \code{"<alpha/2> \%"} /
##'   \code{"<1-alpha/2> \%"}. Attributes \code{"method"},
##'   \code{"vcov_type"} (and, for the bootstrap paths,
##'   \code{"boot.type"}) record the options used.
##' @seealso \code{\link[=vcov.rlmerMod]{vcov}},
##'   \code{\link[confintROB]{confintROB}}
##' @references Mason F, Cantoni E, Ghisletta P (2021). \emph{Parametric
##'   and bootstrap-based inference for linear mixed-effects models in
##'   the presence of outliers}. Methodology 17(4): 271--293.
##' @references Mason F, Cantoni E, Ghisletta P (2024). \emph{Bootstrap
##'   confidence intervals for fixed effects in mixed-effects models
##'   with outliers}. Psychological Methods.
##' @importFrom stats confint qnorm qt
##' @method confint rlmerMod
##' @export
confint.rlmerMod <- function(object, parm = NULL, level = 0.95,
                             method = c("Wald", "boot", "BCa"),
                             vcov_type = c("default", "sandwich"),
                             df = c("none", "satterthwaite"),
                             boot.type = c("wild", "parametric"),
                             nsim = 1000L, seed = NULL, ...) {
    stopifnot(is(object, "rlmerMod"))
    method    <- match.arg(method)
    vcov_type <- match.arg(vcov_type)
    df        <- match.arg(df)
    boot.type <- match.arg(boot.type)

    beta <- .fixef(object)
    nms  <- colnames(object@pp$X)
    names(beta) <- nms

    if (method == "Wald") {
        parm_idx <- .resolve_parm(parm, beta, nms)
        a   <- (1 - level) / 2
        V   <- as.matrix(vcov(object, type = vcov_type))
        se  <- sqrt(diag(V))
        ## WS16 step 5: critical value is the normal quantile by default,
        ## or a per-coefficient Satterthwaite t-quantile when df =
        ## "satterthwaite" (default vcov only). Any obstacle -- a
        ## non-default vcov, an unsupported design (>1 grouping factor or
        ## Mallows eta), or a variance-component boundary -- warns and
        ## falls back to the z-quantile, i.e. the previous behaviour.
        crit <- rep(qnorm(1 - a), length(beta))
        if (df == "satterthwaite") {
            if (vcov_type != "default") {
                warning("df = \"satterthwaite\" requires vcov_type = ",
                        "\"default\"; using z-quantiles.", call. = FALSE)
            } else {
                sw <- tryCatch(suppressWarnings(.satterthwaite_df(object)),
                               error = function(e) e)
                if (inherits(sw, "error")) {
                    warning("Satterthwaite df unavailable (",
                            conditionMessage(sw),
                            "); using z-quantiles.", call. = FALSE)
                } else if (isTRUE(attr(sw, "boundary")) &&
                           !isTRUE(attr(sw, "reducible"))) {
                    warning("Variance-component boundary (singular fit): ",
                            "Satterthwaite df not identifiable; using ",
                            "z-quantiles.", call. = FALSE)
                } else {
                    ## WS14: a reducible boundary yields a valid df
                    ## conditional on the zeroed component(s); use it.
                    crit <- qt(1 - a, df = pmax(as.numeric(sw), 1))
                }
            }
        }
        margin <- crit * se
        out <- cbind(beta - margin, beta + margin)[parm_idx, , drop = FALSE]
        colnames(out) <- .pct_cols(a)
        rownames(out) <- nms[parm_idx]
        attr(out, "method")    <- "Wald"
        attr(out, "vcov_type") <- vcov_type
        return(out)
    }

    ## method = "boot" or "BCa": dispatch to confintROB.
    if (!requireNamespace("confintROB", quietly = TRUE))
        stop("method = \"", method, "\" requires the 'confintROB' package; ",
             "install confintROB or use method = \"Wald\".",
             call. = FALSE)
    if (vcov_type != "default")
        warning("vcov_type = \"", vcov_type, "\" is ignored when ",
                "method = \"", method, "\"; confintROB uses its own ",
                "internal covariance.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)

    ## NOTE: do NOT forward parm to confintROB. confintROB treats an
    ## explicit parm = NULL as "no parameters" (returns a 0-row matrix)
    ## while omitting parm means "all parameters". We always ask for
    ## everything and subset rows in this wrapper so confint(fit) and
    ## confint(fit, parm = "Days") both behave consistently.
    ci_full <- confintROB::confintROB(object,
                                       level = level,
                                       method = method,
                                       nsim  = nsim,
                                       boot.type = boot.type,
                                       ...)
    ## confintROB returns CIs for fixed effects + variance components;
    ## subset to fixed effects so the row set matches the Wald path.
    keep <- rownames(ci_full) %in% nms
    out  <- ci_full[keep, , drop = FALSE]
    ## Reorder to match the original fixef order; honour parm.
    parm_idx <- .resolve_parm(parm, beta, nms)
    requested <- nms[parm_idx]
    out <- out[match(requested, rownames(out)), , drop = FALSE]
    attr(out, "method")    <- method
    attr(out, "vcov_type") <- vcov_type
    attr(out, "boot.type") <- boot.type
    out
}

.resolve_parm <- function(parm, beta, nms) {
    if (is.null(parm)) return(seq_along(beta))
    if (is.character(parm)) {
        idx <- match(parm, nms)
        if (anyNA(idx))
            stop("Unknown coefficient name(s): ",
                 paste(parm[is.na(idx)], collapse = ", "),
                 call. = FALSE)
        return(idx)
    }
    as.integer(parm)
}

.pct_cols <- function(a) {
    pct <- function(p) paste0(format(100 * p, trim = TRUE,
                                     scientific = FALSE), " %")
    c(pct(a), pct(1 - a))
}
