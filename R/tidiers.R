## broom-style methods for rlmerMod objects: tidy(), glance(), augment().
## rlmerMod is not a merMod, so broom.mixed does not cover it. The methods
## are registered for the generics of the 'generics' package when that is
## loaded (see NAMESPACE); column names follow broom.mixed.

## tibble when available, data.frame otherwise
.as_tbl <- function(df) {
    rownames(df) <- NULL
    if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df)
    else df
}

## Smallest random-effect robustness weight per level of each grouping
## factor: a list with one named vector (names = levels) per factor.
.wt_b_by_factor <- function(object) {
    wb <- getME(object, "w_b")
    fl <- object@flist
    out <- lapply(names(fl), function(f) {
        w <- stats::setNames(rep(1, nlevels(fl[[f]])), levels(fl[[f]]))
        for (i in which(names(wb) == f)) {
            wi <- do.call(pmin, unname(as.list(wb[[i]])))
            w[rownames(wb[[i]])] <- pmin(w[rownames(wb[[i]])], wi)
        }
        w
    })
    names(out) <- names(fl)
    out
}

## Random-effect parameters in broom.mixed's layout
.tidy_ran_pars <- function(x) {
    vc <- VarCorr(x)
    reCovs <- attr(x, "reCovs")
    rows <- list()
    for (i in seq_along(vc)) {
        sdv <- attr(vc[[i]], "stddev")
        nm <- names(sdv)
        rows[[length(rows) + 1L]] <-
            data.frame(group = names(vc)[i], term = paste0("sd__", nm),
                       estimate = unname(sdv), stringsAsFactors = FALSE)
        ## a diagonal structure has no correlation parameters
        isDiag <- !is.null(reCovs) && length(reCovs) >= i &&
            inherits(reCovs[[i]], "Covariance.diag")
        if (length(sdv) > 1L && !isDiag) {
            cr <- attr(vc[[i]], "correlation")
            idx <- which(lower.tri(cr), arr.ind = TRUE)
            rows[[length(rows) + 1L]] <-
                data.frame(group = names(vc)[i],
                           term = paste0("cor__", nm[idx[, "col"]], ".",
                                         nm[idx[, "row"]]),
                           estimate = cr[idx], stringsAsFactors = FALSE)
        }
    }
    rows[[length(rows) + 1L]] <-
        data.frame(group = "Residual", term = "sd__Observation",
                   estimate = .sigma(x), stringsAsFactors = FALSE)
    do.call(rbind, rows)
}

##' Tidying methods for robust linear mixed models
##'
##' Methods for the \code{tidy}, \code{glance} and \code{augment} generics
##' of the \pkg{generics} package (re-exported by \pkg{broom} and
##' \pkg{broom.mixed}) for models fitted by \code{\link{rlmer}}. The column
##' names follow \pkg{broom.mixed}. The methods are available once one of
##' these packages is loaded. They return a tibble if the \pkg{tibble}
##' package is installed and a data frame otherwise.
##'
##' \code{augment} returns the \strong{robustness weights} of the fit next
##' to the data: \code{.wt_e} is the weight of each observation and
##' \code{.wt_b} the weight of the random effect of its cluster. A weight of
##' 1 means the observation or cluster is treated as in the classical fit;
##' a weight well below 1 marks it as outlying and downweighted. With
##' several grouping factors there is one column per factor,
##' \code{.wt_b.<factor>}. If a cluster has several random effects, the
##' smallest of their weights is reported.
##'
##' @section Likelihood-based statistics:
##'   \code{glance} has no \code{logLik}, \code{AIC}, \code{BIC} or
##'   \code{deviance} columns: the robust estimator solves estimating
##'   equations that are not the derivative of a likelihood, so these
##'   quantities are not defined for an \code{rlmerMod}.
##'
##' @param x an \code{rlmerMod} object.
##' @param effects which rows to return: \code{"fixed"} (fixed-effect
##'   coefficients) and / or \code{"ran_pars"} (standard deviations and
##'   correlations of the random effects, and the residual standard
##'   deviation).
##' @param conf.int add confidence limits for the fixed effects
##'   (\code{conf.low}, \code{conf.high}, from
##'   \code{\link[=confint.rlmerMod]{confint}})?
##' @param conf.level their confidence level.
##' @param vcov_type covariance of the fixed effects: \code{"default"}
##'   (model-based) or \code{"sandwich"} (cluster sandwich with the
##'   bias-reduced CR2 correction, as \code{confint} uses).
##' @param df for \code{vcov_type = "default"}, passed to
##'   \code{\link[=rlmerMod-class]{summary}}: the degrees of freedom and
##'   p-values are those of the coefficient table of \code{summary(x, df =
##'   df)}, that is robust Satterthwaite when available and \code{NA}
##'   otherwise.
##' @param cluster for \code{vcov_type = "sandwich"}, the clusters; see
##'   \code{\link{vcov_sandwich}}.
##' @param threshold robustness weights below this value count as
##'   downweighted in \code{glance}.
##' @param data the data to augment; by default the model frame. The
##'   original data can be given if it has one row per observation used in
##'   the fit, or additionally the rows dropped for missing values.
##' @param newdata new data to predict for. If given, \code{augment} returns
##'   it with \code{.fitted} (and \code{.se.fit}) only.
##' @param se_fit add the standard errors of the fitted values
##'   (\code{.se.fit}, from \code{\link[=rlmerMod-class]{predict}} with
##'   \code{interval = "confidence"})?
##' @param influence add Cook's distances (\code{.cooksd}, from
##'   \code{\link[=cooks.distance.rlmerMod]{cooks.distance}})? This needs
##'   the influence function of the fit, which is computed once and cached.
##' @param ... \code{augment}: passed to \code{predict}. Otherwise not
##'   used.
##' @return \describe{
##'   \item{\code{tidy}}{one row per parameter, with columns \code{effect},
##'     \code{group}, \code{term}, \code{estimate}, \code{std.error},
##'     \code{statistic}, \code{df}, \code{p.value} and, if requested,
##'     \code{conf.low} and \code{conf.high}. Random-effect parameters are
##'     named as in \pkg{broom.mixed} (\code{sd__(Intercept)},
##'     \code{cor__(Intercept).x}, \code{sd__Observation}) and have an
##'     estimate only.}
##'   \item{\code{glance}}{one row: \code{nobs}, \code{sigma},
##'     \code{converged} (\code{NA} for an object that was not fitted),
##'     the number of levels of each grouping factor
##'     (\code{ngrps.<factor>}), and the number of observations
##'     (\code{nobs.downweighted}) and of clusters per grouping factor
##'     (\code{ngrps.downweighted.<factor>}) with a robustness weight below
##'     \code{threshold}.}
##'   \item{\code{augment}}{the data with the columns \code{.fitted},
##'     \code{.resid} (response residuals), \code{.wt_e}, \code{.wt_b} (see
##'     above), \code{.hat} (\code{\link[=hatvalues.rlmerMod]{hatvalues}})
##'     and, if requested, \code{.se.fit} and \code{.cooksd}.}
##' }
##' @examples
##' fm <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy)
##' if (requireNamespace("generics", quietly = TRUE)) {
##'   print(generics::tidy(fm, conf.int = TRUE))
##'   print(generics::glance(fm))
##'   ## the least trusted observations and subjects
##'   aug <- generics::augment(fm)
##'   print(aug[order(aug$.wt_e)[1:5], ])
##' }
##' @section Lifecycle:
##'   Stable.
##' @name rlmerMod-tidiers
##' @aliases tidy.rlmerMod glance.rlmerMod augment.rlmerMod
NULL

##' @rdname rlmerMod-tidiers
tidy.rlmerMod <- function(x, effects = c("fixed", "ran_pars"),
                          conf.int = FALSE, conf.level = 0.95,
                          vcov_type = c("default", "sandwich"),
                          df = c("auto", "none", "satterthwaite"),
                          cluster = NULL, ...) {
    effects <- match.arg(effects, several.ok = TRUE)
    vcov_type <- match.arg(vcov_type)
    df <- match.arg(df)
    .check_flag(conf.int, "conf.int")
    .check_level(conf.level, "conf.level")
    out <- list()
    if ("fixed" %in% effects) {
        nms <- colnames(x@pp$X)
        p <- length(nms)
        if (vcov_type == "default") {
            co <- coef(summary(x, df = df))
            hasDf <- "df" %in% colnames(co)
            fixed <- data.frame(
                effect = "fixed", group = NA_character_, term = nms,
                estimate = unname(co[, "Estimate"]),
                std.error = unname(co[, "Std. Error"]),
                statistic = unname(co[, "t value"]),
                df = if (hasDf) unname(co[, "df"]) else NA_real_,
                p.value = if (hasDf) unname(co[, "Pr(>|t|)"]) else NA_real_,
                stringsAsFactors = FALSE)
            ciDf <- if (hasDf && !anyNA(co[, "df"])) "satterthwaite"
                    else "auto"
        } else {
            Vraw <- .cluster_vcov(x, cluster)
            V <- as.matrix(Vraw)
            beta <- .fixef(x)
            rows <- do.call(rbind, lapply(seq_len(p), function(j) {
                L <- matrix(0, 1L, p); L[1L, j] <- 1
                .contrast_one(x, L, beta[j], V[j, j], "cluster", Vraw)
            }))
            fixed <- data.frame(
                effect = "fixed", group = NA_character_, term = nms,
                estimate = rows$Estimate, std.error = rows$`Std. Error`,
                statistic = rows$`t value`, df = rows$df,
                p.value = rows$`Pr(>|t|)`, stringsAsFactors = FALSE)
            ciDf <- "auto"
        }
        if (conf.int) {
            ## same reference as the p-values
            ci <- suppressWarnings(
                confint(x, level = conf.level, vcov_type = vcov_type,
                        df = ciDf, cluster = cluster))
            fixed$conf.low <- unname(ci[, 1L])
            fixed$conf.high <- unname(ci[, 2L])
        }
        out$fixed <- fixed
    }
    if ("ran_pars" %in% effects) {
        rp <- .tidy_ran_pars(x)
        ran <- data.frame(effect = "ran_pars", group = rp$group,
                          term = rp$term, estimate = rp$estimate,
                          std.error = NA_real_, statistic = NA_real_,
                          df = NA_real_, p.value = NA_real_,
                          stringsAsFactors = FALSE)
        if (conf.int) {
            ran$conf.low <- NA_real_
            ran$conf.high <- NA_real_
        }
        out$ran_pars <- ran
    }
    .as_tbl(do.call(rbind, out[effects]))
}

##' @rdname rlmerMod-tidiers
glance.rlmerMod <- function(x, threshold = 0.5, ...) {
    .check_number(threshold, "threshold", lower = 0, upper = 1)
    fl <- x@flist
    ngrps <- vapply(fl, nlevels, integer(1))
    names(ngrps) <- paste0("ngrps.", names(fl))
    low <- vapply(.wt_b_by_factor(x), function(w) sum(w < threshold),
                  integer(1))
    names(low) <- paste0("ngrps.downweighted.", names(fl))
    code <- x@optinfo$conv$opt
    out <- data.frame(nobs = nobs(x), sigma = .sigma(x),
                      converged = if (is.null(code)) NA else code == 0,
                      as.list(ngrps),
                      nobs.downweighted = sum(wgt.e(x) < threshold),
                      as.list(low), check.names = FALSE)
    .as_tbl(out)
}

##' @rdname rlmerMod-tidiers
augment.rlmerMod <- function(x, data = model.frame(x), newdata = NULL,
                             se_fit = FALSE, influence = FALSE, ...) {
    .check_flag(se_fit, "se_fit")
    .check_flag(influence, "influence")
    if (!is.null(newdata)) {
        if (!is.data.frame(newdata))
            stop("'newdata' must be a data frame", call. = FALSE)
        out <- newdata
        if (se_fit) {
            pr <- predict(x, newdata = newdata, interval = "confidence", ...)
            out$.fitted <- pr$fit
            out$.se.fit <- pr$se
        } else {
            out$.fitted <- unname(predict(x, newdata = newdata, ...))
        }
        return(.as_tbl(out))
    }
    if (!is.data.frame(data))
        stop("'data' must be a data frame", call. = FALSE)
    n <- nobs(x)
    cols <- data.frame(.fitted = unname(x@resp$mu),
                       .resid = x@resp$y - x@resp$mu,
                       .wt_e = wgt.e(x))
    wb <- .wt_b_by_factor(x)
    for (f in names(wb)) {
        nm <- if (length(wb) == 1L) ".wt_b" else paste0(".wt_b.", f)
        cols[[nm]] <- unname(wb[[f]][as.character(x@flist[[f]])])
    }
    cols$.hat <- unname(hatvalues(x))
    na.act <- attr(model.frame(x), "na.action")
    if (se_fit) {
        se <- predict(x, interval = "confidence", ...)$se
        ## predict() pads the rows excluded by na.exclude
        if (length(se) == n + length(na.act) && length(na.act))
            se <- se[-na.act]
        cols$.se.fit <- se
    }
    if (influence)
        cols$.cooksd <- unname(cooks.distance(x))
    ## rows dropped for missing values get NA
    if (nrow(data) != n) {
        if (!is.null(na.act) && nrow(data) == n + length(na.act)) {
            full <- cols[rep(NA_integer_, nrow(data)), , drop = FALSE]
            full[-na.act, ] <- cols
            cols <- full
        } else {
            stop("'data' has ", nrow(data), " rows but the model was ",
                 "fitted to ", n, " observations", call. = FALSE)
        }
    }
    .as_tbl(cbind(data, cols))
}
