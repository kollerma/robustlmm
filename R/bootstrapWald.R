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
##'   \item{\code{method = "Wald"}}{the default: \eqn{\hat\beta_k \pm q\,
##'     SE_k} with \eqn{SE_k} from \code{vcov(object, type = vcov_type)}
##'     and \eqn{q} a normal or \eqn{t} quantile (see \code{df}).}
##'   \item{\code{method = "boot"} or \code{"BCa"}}{delegated to
##'     \code{\link[confintROB]{confintROB}} (Mason, Cantoni and Ghisletta
##'     2021, 2024), returning its fixed-effect rows; \code{vcov_type =
##'     "sandwich"}, \code{df} and \code{cluster} cannot be combined with
##'     it.}
##' }
##'
##' The Wald limit is adequate for about 20 or more groups; the bootstrap
##' earns its cost mainly with fewer (Mason et al. 2024). \code{boot.type =
##' "wild"} is robust to misspecification of the response covariance,
##' \code{"parametric"} relies on the fitted central model, and
##' \code{method = "BCa"} corrects the percentile interval for bias and
##' skewness. The samples of both bootstraps are also available from
##' \code{\link[=simulate.rlmerMod]{simulate}}.
##'
##' \strong{Sandwich covariance.} With few clusters, normal quantiles with
##' the G1 sandwich under-cover, and for between-cluster coefficients the G1
##' sandwich itself is biased down. For \code{vcov_type = "sandwich"} the
##' default (\code{df = "auto"}, i.e. \code{"cluster"}) is therefore the
##' bias-reduced CR2 sandwich of \code{\link{vcov_sandwich}} with a
##' Bell--McCaffrey \eqn{t} quantile per coefficient (the degrees of freedom
##' are in \code{attr(., "df")}); see the section \dQuote{Small-sample
##' references for the cluster sandwich} of \code{vignette("simulationStudies")}.
##' When the random effects are not nested within the clusters, the G1
##' sandwich with \eqn{t_{J-1}} quantiles is used. A coefficient that a
##' single cluster informs on its own (e.g. a factor level present in one
##' cluster) gets an \code{NA} interval, with a warning.
##'
##' @param object An \code{rlmerMod} object.
##' @param parm Either \code{NULL} (all fixed-effect coefficients), an
##'   integer vector of coefficient indices, or a character vector of
##'   coefficient names.
##' @param level Coverage level; default 0.95.
##' @param method One of \code{"Wald"} (default; closed form),
##'   \code{"boot"} (bootstrap percentile via \pkg{confintROB}), or
##'   \code{"BCa"} (bias-corrected bootstrap via \pkg{confintROB}).
##' @param vcov_type Covariance for \code{method = "Wald"}:
##'   \code{"default"} (the model-based covariance) or \code{"sandwich"}
##'   (the cluster sandwich, \code{\link{vcov_sandwich}}).
##' @param df Quantile for \code{method = "Wald"}. \code{"auto"} (default):
##'   \code{"none"} for \code{vcov_type = "default"} and \code{"cluster"}
##'   for \code{vcov_type = "sandwich"}. \code{"none"}: the normal quantile
##'   (with the G1 sandwich for \code{vcov_type = "sandwich"}).
##'   \code{"cluster"}: the CR2 sandwich with Bell--McCaffrey \eqn{t}
##'   quantiles (sandwich only). \code{"satterthwaite"}: robust Satterthwaite
##'   \eqn{t} quantiles as in \code{summary(object, df = "satterthwaite")}
##'   (default covariance and a single grouping factor only; otherwise the
##'   default quantile, with a warning).
##' @param boot.type Bootstrap kind passed to \code{confintROB} when
##'   \code{method = "boot"} or \code{"BCa"}: one of
##'   \code{"wild"} (default, the \code{confintROB} recommendation) or
##'   \code{"parametric"}.
##' @param nsim Bootstrap replicates; default 1000.
##' @param cluster For \code{vcov_type = "sandwich"}, the clusters, passed
##'   to \code{\link{vcov_sandwich}}; needed for models with several
##'   grouping factors.
##' @param seed Optional RNG seed for reproducibility of the bootstrap.
##' @param ... Additional arguments forwarded to \code{confintROB}
##'   (e.g. \code{clusterID} for the wild bootstrap).
##' @return A 2-column matrix with one row per selected fixed-effect
##'   coefficient, columns \code{"<alpha/2> \%"} /
##'   \code{"<1-alpha/2> \%"}. Attributes \code{"method"},
##'   \code{"vcov_type"} (and, for the bootstrap paths,
##'   \code{"boot.type"}) record the options used; the small-J sandwich
##'   reference adds \code{"df"} (per selected coefficient) and
##'   \code{"correction"} (\code{"CR2"}, or \code{"G1"} when CR2 is not
##'   defined).
##' @seealso \code{\link[=vcov.rlmerMod]{vcov}},
##'   \code{\link[confintROB]{confintROB}}
##' @references Bell RM, McCaffrey DF (2002). Bias reduction in standard
##'   errors for linear regression with multi-stage samples. \emph{Survey
##'   Methodology} 28(2), 169--181.
##' @references Mason F, Cantoni E, Ghisletta P (2021). Parametric and
##'   semi-parametric bootstrap-based confidence intervals for robust linear
##'   mixed models. \emph{Methodology} 17(4), 271--295.
##' @references Mason F, Cantoni E, Ghisletta P (2024). Linear mixed models
##'   and latent growth curve models for group comparison studies
##'   contaminated by outliers. \emph{Psychological Methods} 30(1),
##'   155--173.
##' @importFrom stats confint qnorm qt
##' @method confint rlmerMod
##' @section Lifecycle:
##'   Stable.
##' @export
confint.rlmerMod <- function(object, parm = NULL, level = 0.95,
                             method = c("Wald", "boot", "BCa"),
                             vcov_type = c("default", "sandwich"),
                             df = c("auto", "none", "satterthwaite",
                                    "cluster"),
                             boot.type = c("wild", "parametric"),
                             nsim = 1000L, seed = NULL, cluster = NULL, ...) {
    stopifnot(is(object, "rlmerMod"))
    dfGiven   <- !missing(df)
    method    <- match.arg(method)
    vcov_type <- match.arg(vcov_type)
    df        <- match.arg(df)
    boot.type <- match.arg(boot.type)
    .check_level(level)
    if (method == "Wald") {
        if (!is.null(cluster) && vcov_type != "sandwich")
            stop("'cluster' is used only with vcov_type = \"sandwich\"",
                 call. = FALSE)
    } else {
        ## arguments of the Wald interval have no effect on the bootstrap
        if (vcov_type != "default")
            stop("vcov_type = \"", vcov_type, "\" cannot be combined with ",
                 "method = \"", method, "\": the bootstrap of confintROB ",
                 "does not use a covariance of the fixed effects",
                 call. = FALSE)
        if (dfGiven)
            stop("'df' is used only with method = \"Wald\"", call. = FALSE)
        if (!is.null(cluster))
            stop("'cluster' is used only with method = \"Wald\"; the wild ",
                 "bootstrap of confintROB takes 'clusterID'", call. = FALSE)
        .check_number(nsim, "nsim", lower = 1, integer = TRUE)
    }

    beta <- .fixef(object)
    nms  <- colnames(object@pp$X)
    names(beta) <- nms

    if (method == "Wald") {
        parm_idx <- .resolve_parm(parm, beta, nms)
        a   <- (1 - level) / 2
        ## A df choice unavailable for the chosen vcov falls back to that
        ## vcov's default (df = "auto").
        if (df == "satterthwaite" && vcov_type != "default") {
            warning("df = \"satterthwaite\" requires vcov_type = ",
                    "\"default\"; using the sandwich's default (df = ",
                    "\"auto\").", call. = FALSE)
            df <- "auto"
        }
        if (df == "cluster" && vcov_type != "sandwich") {
            warning("df = \"cluster\" applies to vcov_type = \"sandwich\"; ",
                    "using z-quantiles.", call. = FALSE)
            df <- "none"
        }
        if (df == "auto" && vcov_type == "sandwich")
            df <- if (.ddf_auto_sandwich == "cluster") "cluster" else "none"
        ## Sandwich with df = "cluster": the small-J reference of
        ## theory-methods F3 -- the CR2 (bias-reduced) sandwich with a
        ## Bell-McCaffrey t df per coefficient (G1 with t(J - 1) when CR2 is
        ## not defined, i.e. random effects not nested within clusters).
        Vraw <- if (vcov_type == "sandwich" && df == "cluster")
                    .cluster_vcov(object, cluster)
                else vcov(object, type = vcov_type, cluster = cluster)
        V   <- as.matrix(Vraw)
        se  <- sqrt(diag(V))
        ## WS16 step 5: critical value is the normal quantile by default,
        ## or a per-coefficient Satterthwaite t-quantile when df =
        ## "satterthwaite" (default vcov only). Any obstacle -- a
        ## non-default vcov, an unsupported design (>1 grouping factor or
        ## Mallows eta), or a variance-component boundary -- warns and
        ## falls back to the z-quantile, i.e. the previous behaviour.
        crit <- rep(qnorm(1 - a), length(beta))
        dfv  <- NULL
        if (vcov_type == "sandwich" && df == "cluster") {
            dfv <- if (identical(attr(Vraw, "correction"), "CR2"))
                       vapply(seq_along(beta), function(i) {
                           L <- matrix(0, 1L, length(beta)); L[1L, i] <- 1
                           .cr2_df(Vraw, L)
                       }, numeric(1))
                   else rep(attr(Vraw, "n.clusters") - 1, length(beta))
            low <- which(is.finite(dfv) & dfv < 1)
            if (length(low))
                warning("Fewer than one cluster degree of freedom for ",
                        paste(nms[low], collapse = ", "),
                        ": the clusters carry almost no information on ",
                        if (length(low) > 1L) "these coefficients" else
                        "this coefficient", ".", call. = FALSE)
            crit <- qt(1 - a, df = dfv)
        }
        if (df == "satterthwaite") {
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
        margin <- crit * se
        out <- cbind(beta - margin, beta + margin)[parm_idx, , drop = FALSE]
        colnames(out) <- .pct_cols(a)
        rownames(out) <- nms[parm_idx]
        attr(out, "method")    <- "Wald"
        attr(out, "vcov_type") <- vcov_type
        if (!is.null(dfv)) {
            attr(out, "df") <- dfv[parm_idx]
            attr(out, "correction") <- attr(Vraw, "correction")
        }
        return(out)
    }

    ## method = "boot" or "BCa": dispatch to confintROB.
    if (!requireNamespace("confintROB", quietly = TRUE))
        stop("method = \"", method, "\" requires the 'confintROB' package; ",
             "install confintROB or use method = \"Wald\".",
             call. = FALSE)
    ## confintROB refits with lmer(), roughly 1/50 of a robust refit
    .bootMessage(object, nsim, 0.02,
                 sprintf("confint(method = \"%s\")", method),
                 "Reduce 'nsim', or use method = \"Wald\".")
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
