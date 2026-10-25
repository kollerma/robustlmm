## contrast_test(): Wald test of user-specified linear contrasts of the
## fixed effects, with the same references as summary(), anova() and the
## emmeans support (robust Satterthwaite df for the model-based covariance,
## CR2 / approximate-Hotelling df for the cluster sandwich).

## Contrast matrix with one column per fixed effect. A vector is one
## contrast; column names are matched to the coefficient names and
## coefficients not named get a zero.
.contrast_matrix <- function(L, nms) {
    if (!is.numeric(L) || anyNA(L))
        stop("'L' must be a numeric vector or matrix without missing values",
             call. = FALSE)
    if (is.null(dim(L)))
        L <- matrix(L, nrow = 1L, dimnames = list(NULL, names(L)))
    L <- as.matrix(L)
    p <- length(nms)
    if (!is.null(colnames(L))) {
        if (anyDuplicated(colnames(L)))
            stop("'L' has duplicated column names", call. = FALSE)
        bad <- setdiff(colnames(L), nms)
        if (length(bad))
            stop("'L' has columns that are not fixed-effect coefficients: ",
                 paste(bad, collapse = ", "), "; the coefficients are ",
                 paste(nms, collapse = ", "), call. = FALSE)
        full <- matrix(0, nrow(L), p, dimnames = list(rownames(L), nms))
        full[, colnames(L)] <- L
        L <- full
    } else if (ncol(L) != p) {
        stop("'L' has ", ncol(L), " column(s) but the model has ", p,
             " fixed-effect coefficient(s); give one column per coefficient ",
             "or name the columns", call. = FALSE)
    } else {
        colnames(L) <- nms
    }
    if (nrow(L) == 0L) stop("'L' has no rows", call. = FALSE)
    if (any(rowSums(abs(L)) == 0))
        stop("'L' has a row of zeros", call. = FALSE)
    L
}

## Satterthwaite df, or NULL with a warning when unavailable
.satt_or_warn <- function(fun) {
    res <- tryCatch(fun(), error = function(e) e)
    if (inherits(res, "error")) {
        warning("Satterthwaite df unavailable (", conditionMessage(res),
                "); using the normal reference.", call. = FALSE)
        return(NULL)
    }
    res
}

## t test of one contrast (L: 1 x p) with estimate `est` and variance `v`
## under the reference `ddf` ("satterthwaite", "cluster" or "none"); Vraw
## is the covariance the variance was taken from.
.contrast_one <- function(object, L, est, v, ddf, Vraw) {
    se <- sqrt(v)
    df <- Inf
    if (ddf == "satterthwaite") {
        sw <- .satt_or_warn(function() .satterthwaite_df(object, Lmat = L))
        if (!is.null(sw)) {
            if (isTRUE(attr(sw, "boundary")) &&
                !isTRUE(attr(sw, "reducible")))
                warning("Variance-component boundary (singular fit): ",
                        "Satterthwaite df not identifiable; using the ",
                        "normal reference.", call. = FALSE)
            else df <- as.numeric(sw)
        }
    } else if (ddf == "cluster") {
        df <- if (identical(attr(Vraw, "correction"), "CR2"))
                  .cr2_df(Vraw, L)
              else attr(Vraw, "n.clusters") - 1
    }
    tval <- est / se
    data.frame(Estimate = est, `Std. Error` = se, df = df,
               `t value` = tval,
               `Pr(>|t|)` = 2 * stats::pt(abs(tval), df, lower.tail = FALSE),
               check.names = FALSE)
}

##' Wald test of linear contrasts of the fixed effects
##'
##' Tests \eqn{L \beta = rhs} for a fitted \code{\link{rlmer}} model. A
##' single contrast gives a \eqn{t} test, several contrasts a joint \eqn{F}
##' test. The reference distributions are those of \code{summary},
##' \code{\link[=anova.rlmerMod]{anova}} and the \pkg{emmeans} support.
##'
##' @param object an \code{rlmerMod} object.
##' @param L the contrasts: a numeric vector (one contrast) or a matrix with
##'   one contrast per row. Either one column per fixed-effect coefficient,
##'   in the order of \code{fixef(object)}, or named columns that are
##'   matched to the coefficient names; coefficients that are not named get
##'   a zero.
##' @param rhs the hypothesised value of \eqn{L \beta}: one number or one
##'   per contrast.
##' @param ddf the denominator degrees of freedom. \code{"auto"} (default):
##'   \code{"satterthwaite"} for \code{vcov_type = "default"} and
##'   \code{"cluster"} for \code{vcov_type = "sandwich"}.
##'   \code{"satterthwaite"}: the robust Satterthwaite degrees of freedom of
##'   \code{summary(object, df = "satterthwaite")} (model-based covariance
##'   only). \code{"cluster"}: the bias-reduced CR2 sandwich with
##'   Bell--McCaffrey (one contrast) or approximate-Hotelling (several
##'   contrasts) degrees of freedom, as in \code{anova(object, vcov_type =
##'   "sandwich")} (sandwich only). \code{"none"}: the normal and
##'   chi-square references, reported as infinite degrees of freedom.
##' @param vcov_type \code{"default"} (the model-based covariance) or
##'   \code{"sandwich"} (the cluster sandwich, \code{\link{vcov_sandwich}}).
##' @param cluster for \code{vcov_type = "sandwich"}, the clusters, passed
##'   to \code{\link{vcov_sandwich}}; needed for models with several
##'   grouping factors.
##' @param ... not used.
##' @return A data frame of class \code{"anova"}. For one contrast: the
##'   columns \code{Estimate} (of \eqn{L \beta - rhs}), \code{Std. Error},
##'   \code{df}, \code{t value} and \code{Pr(>|t|)}. For several contrasts:
##'   one row with \code{NumDF}, \code{DenDF}, \code{F value} and
##'   \code{Pr(>F)}. The contrast matrix is in \code{attr(., "L")}.
##' @details If the Satterthwaite degrees of freedom are not available for
##'   the fit (for example at a variance-component boundary that leaves
##'   them unidentified), the normal or chi-square reference is used, with
##'   a warning.
##' @seealso \code{\link[=anova.rlmerMod]{anova}},
##'   \code{\link[=confint.rlmerMod]{confint}}, \code{\link{vcov_sandwich}}
##' @examples
##' fm <- rlmer(Reaction ~ Days + (1 | Subject), sleepstudy)
##' ## expected reaction time on day 5
##' contrast_test(fm, c(1, 5))
##' ## is it 300 ms?
##' contrast_test(fm, c("(Intercept)" = 1, Days = 5), rhs = 300)
##' ## both coefficients zero, jointly
##' contrast_test(fm, diag(2))
##' @section Lifecycle:
##'   Stable.
##' @export
contrast_test <- function(object, L, rhs = 0,
                          ddf = c("auto", "satterthwaite", "cluster", "none"),
                          vcov_type = c("default", "sandwich"),
                          cluster = NULL, ...) {
    if (!is(object, "rlmerMod"))
        stop("'object' must be an rlmerMod object", call. = FALSE)
    ddf <- match.arg(ddf)
    vcov_type <- match.arg(vcov_type)
    beta <- .fixef(object)
    nms <- colnames(object@pp$X)
    L <- .contrast_matrix(L, nms)
    k <- nrow(L)
    if (!is.numeric(rhs) || anyNA(rhs) || !length(rhs) %in% c(1L, k))
        stop("'rhs' must be one number or one per contrast (", k, ")",
             call. = FALSE)
    rhs <- rep_len(rhs, k)
    if (ddf == "satterthwaite" && vcov_type != "default")
        stop("ddf = \"satterthwaite\" requires vcov_type = \"default\"; ",
             "use ddf = \"cluster\" with the sandwich", call. = FALSE)
    if (ddf == "cluster" && vcov_type != "sandwich")
        stop("ddf = \"cluster\" requires vcov_type = \"sandwich\"",
             call. = FALSE)
    if (!is.null(cluster) && vcov_type != "sandwich")
        stop("'cluster' is used only with vcov_type = \"sandwich\"",
             call. = FALSE)
    if (ddf == "auto")
        ddf <- if (vcov_type == "sandwich") "cluster" else "satterthwaite"

    Vraw <- if (ddf == "cluster") .cluster_vcov(object, cluster)
            else vcov(object, type = vcov_type, cluster = cluster)
    V <- as.matrix(Vraw)
    est <- drop(L %*% beta) - rhs
    VL <- tcrossprod(L %*% V, L)
    if (k > 1L && qr(VL)$rank < k)
        stop("the contrasts in 'L' are linearly dependent", call. = FALSE)

    reference <- switch(ddf,
        satterthwaite = "Satterthwaite df, vcov_type = \"default\"",
        cluster = sub("^ +", "", .cluster_heading(Vraw)),
        none = sprintf("normal reference, vcov_type = \"%s\"", vcov_type))

    if (k == 1L) {
        out <- .contrast_one(object, L, est, VL[1L, 1L], ddf, Vraw)
        heading <- "Robust Wald t-test of a fixed-effects contrast (rlmer)"
    } else {
        W <- as.numeric(crossprod(est, solve(VL, est)))
        Fval <- W / k
        df2 <- Inf
        if (ddf == "satterthwaite") {
            fd <- .satt_or_warn(function() .satterthwaite_Fdenom(object, L))
            if (!is.null(fd)) df2 <- fd$df2
        } else if (ddf == "cluster") {
            fr <- .cluster_Fref(W, L, Vraw, label = "contrast_test")
            Fval <- fr$F
            df2 <- fr$df2
        }
        out <- data.frame(NumDF = k, DenDF = df2, `F value` = Fval,
                          `Pr(>F)` = stats::pf(Fval, k, df2,
                                               lower.tail = FALSE),
                          check.names = FALSE)
        heading <- sprintf(
            "Robust Wald F-test of %d fixed-effects contrasts (rlmer)", k)
    }
    rownames(out) <- NULL
    structure(out, heading = c(heading, paste0("  ", reference)),
              L = L, rhs = rhs, class = c("anova", "data.frame"))
}
