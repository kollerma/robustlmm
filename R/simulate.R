## simulate.rlmerMod: parametric and wild bootstrap samples of the response.

## One response from the central Gaussian model
##   y = Xb + Z (sig * Ub z_b) + sig * eps / sqrt(w),
## consuming rnorm(q) followed by rnorm(n). Xb is the fixed part
## including any offset; b (optional) holds the random effects fixed.
## Shared with the bootstrap of anova() (.boot_sim_data), whose seeded
## replicates therefore do not change.
.sim_parametric <- function(Xb, Z, Ub, sig, n, q, w = NULL, b = NULL) {
    if (is.null(b)) b <- sig * as.numeric(Ub %*% rnorm(q))
    eps <- rnorm(n, sd = if (is.null(w)) sig else sig / sqrt(w))
    Xb + as.numeric(Z %*% b) + eps
}

## Moore-Penrose inverse by the singular value decomposition (the
## algorithm of MASS::ginv, so the wild bootstrap reproduces confintROB's
## numbers)
.ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
    s <- svd(X)
    pos <- s$d > max(tol * s$d[1L], 0)
    if (all(pos)) s$v %*% (1 / s$d * t(s$u))
    else if (!any(pos)) array(0, dim(X)[2L:1L])
    else s$v[, pos, drop = FALSE] %*%
             ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
}

## The clusters of the wild bootstrap: the grouping factor for one factor,
## the coarsest factor of a nested design, and `cluster` (name of a
## grouping factor or a vector) otherwise. Crossed designs have no default.
.wild_cluster <- function(object, cluster, n) {
    fl <- object@flist
    if (is.null(cluster)) {
        top <- .top_cluster(object)
        if (is.null(top))
            stop("the grouping factors (", paste(names(fl), collapse = ", "),
                 ") are crossed, so the model has no independent clusters; ",
                 "give 'cluster' to choose the factor (or a vector) whose ",
                 "levels get a common wild weight", call. = FALSE)
        return(top$g)
    }
    if (is.character(cluster) && length(cluster) == 1L) {
        if (!cluster %in% names(fl))
            stop("cluster '", cluster, "' is not a grouping factor; ",
                 "available: ", paste(names(fl), collapse = ", "),
                 call. = FALSE)
        return(fl[[cluster]])
    }
    if (length(cluster) != n || anyNA(cluster))
        stop("'cluster' must be the name of a grouping factor or a vector ",
             "with one value per observation (", n, "), without missing ",
             "values", call. = FALSE)
    factor(cluster)
}

## Generator of wild bootstrap responses (Modugno and Giannerini 2015, as
## implemented in confintROB): marginal residuals, rescaled per cluster by
## sqrt(diag((I - P_j)^+)) with P_j the cluster's block of the hat matrix
## of the fixed-effects design, times one two-point weight per cluster.
## Prior weights enter through the whitened design and residuals; the
## offset is not part of the residual and is added back unchanged.
.wild_generator <- function(object, cluster, wild.weights) {
    y <- object@resp$y
    n <- length(y)
    off <- object@resp$offset
    w <- object@resp$weights
    X <- as.matrix(.X(object))
    bet <- unname(.fixef(object))
    id <- .wild_cluster(object, cluster, n)
    ids <- unique(id)
    TT <- length(ids)
    Xbet <- drop(X %*% bet)
    rt <- unname(y - off - Xbet)
    weighted <- any(w != 1)
    sw <- sqrt(w)
    Xw <- if (weighted) X * sw else X
    ## (X'X)^-1 on the estimable columns
    qrX <- qr(Xw)
    if (qrX$rank < ncol(Xw))
        Xw <- Xw[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE]
    tXX <- solve(crossprod(Xw))
    rw <- if (weighted) sw * rt else rt
    scaled <- numeric(n)
    idxMap <- integer(n)
    for (i in seq_len(TT)) {
        idx <- id == ids[i]
        Xt <- Xw[idx, , drop = FALSE]
        Pt <- tcrossprod(Xt %*% tXX, Xt)
        scaled[idx] <- sqrt(diag(.ginv(diag(sum(idx)) - Pt))) * rw[idx]
        idxMap[idx] <- i
    }
    if (weighted) scaled <- scaled / sw
    mean0 <- if (any(off != 0)) Xbet + off else Xbet
    if (wild.weights == "mammen") {
        v <- c(-(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2)
        p1 <- (sqrt(5) + 1) / (2 * sqrt(5))
        prob <- c(p1, 1 - p1)
        function() {
            wt <- base::sample(v, TT, replace = TRUE, prob = prob)
            mean0 + scaled * wt[idxMap]
        }
    } else {
        function() {
            wt <- base::sample(c(-1, 1), TT, replace = TRUE)
            mean0 + scaled * wt[idxMap]
        }
    }
}

##' Simulate responses from a fitted robust linear mixed model
##'
##' Generates new responses for the design of a model fitted by
##' \code{\link{rlmer}}: parametric samples from the fitted model, or wild
##' bootstrap samples built from its residuals. These are the samples of
##' the parametric and wild bootstraps of \pkg{confintROB}; use them to run
##' a bootstrap of your own.
##'
##' \describe{
##'   \item{\code{type = "parametric"}}{draws from the central Gaussian model
##'     at the robust estimates: random effects with the fitted covariance
##'     and errors with variance \eqn{\hat\sigma^2 / w}{sigma^2 / w} for
##'     prior weights \eqn{w}, added to the fitted fixed part (including any
##'     offset). The samples follow the model the robust fit describes for
##'     the bulk of the data; they do \emph{not} reproduce contamination, so
##'     they contain no outliers. With \code{use.u = TRUE} the fitted
##'     random effects are kept and only the errors are redrawn. This leaves
##'     out the between-cluster variation and is not suitable for confidence
##'     intervals of the model parameters.}
##'   \item{\code{type = "wild"}}{the cluster wild bootstrap of Modugno and
##'     Giannerini (2015). The marginal residuals \eqn{y - X
##'     \hat\beta}{y - X beta} are rescaled within each cluster for the
##'     leverage of the fixed-effects design, multiplied by one random
##'     two-point weight per cluster, and added back to the fixed part. It
##'     makes no assumption on the distribution or the covariance of the
##'     response within a cluster.}
##' }
##'
##' @section Clusters of the wild bootstrap:
##'   All observations of a cluster get the same random weight, so the
##'   clusters have to be independent of each other. With one grouping
##'   factor these are its levels. For nested grouping factors the default
##'   is the outermost (coarsest) factor. Crossed grouping factors have no
##'   independent clusters; \code{cluster} must then be given, and the
##'   dependence across its levels is ignored.
##'
##' @param object an \code{rlmerMod} object.
##' @param nsim number of responses to simulate.
##' @param seed \code{NULL} or an integer passed to \code{\link{set.seed}}
##'   before simulating; the previous state of the random number generator
##'   is then restored on exit.
##' @param type \code{"parametric"} (default) or \code{"wild"}, see Details.
##' @param use.u for \code{type = "parametric"}: keep the fitted random
##'   effects (\code{TRUE}) or draw new ones (\code{FALSE}, default)?
##' @param cluster for \code{type = "wild"}: the clusters, as the name of a
##'   grouping factor of the model or a vector with one value per
##'   observation. See the section on clusters.
##' @param wild.weights for \code{type = "wild"}: the two-point distribution
##'   of the cluster weights, \code{"mammen"} (default; Mammen 1993) or
##'   \code{"rademacher"} (\eqn{\pm 1} with equal probability).
##' @param ... not used.
##' @return A data frame with \code{nsim} columns \code{sim_1}, \ldots, each
##'   a simulated response for the observations used in the fit. The seed
##'   is stored in \code{attr(, "seed")}, as by
##'   \code{\link[lme4:simulate.merMod]{simulate}} for \code{lmer} fits.
##' @references Mammen E (1993). Bootstrap and wild bootstrap for high
##'   dimensional linear models. \emph{The Annals of Statistics} 21(1),
##'   255--285.
##' @references Mason F, Cantoni E, Ghisletta P (2021). Parametric and
##'   semi-parametric bootstrap-based confidence intervals for robust linear
##'   mixed models. \emph{Methodology} 17(4), 271--295.
##' @references Modugno L, Giannerini S (2015). The wild bootstrap for
##'   multilevel models. \emph{Communications in Statistics -- Theory and
##'   Methods} 44(22), 4812--4825.
##' @seealso \code{\link[=confint.rlmerMod]{confint}} with \code{method =
##'   "boot"}, which runs the bootstrap of \pkg{confintROB}.
##' @examples
##' fm <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy)
##' ysim <- simulate(fm, nsim = 3, seed = 1)
##' head(ysim)
##'
##' \dontrun{
##'   ## a wild bootstrap of the slope by hand: refit each sample.
##'   ## confintROB refits with lmer(), which is much faster:
##'   ystar <- simulate(fm, nsim = 200, seed = 2, type = "wild")
##'   slopes <- sapply(ystar, function(y)
##'       fixef(lmer(y ~ Days + (Days | Subject), sleepstudy))["Days"])
##'   quantile(slopes, c(0.025, 0.975))
##'   ## the same with robust refits
##'   slopes <- sapply(ystar[1:20], function(y) {
##'       d <- sleepstudy; d$Reaction <- y
##'       fixef(update(fm, data = d))["Days"]
##'   })
##' }
##' @importFrom stats simulate rnorm runif
##' @method simulate rlmerMod
##' @section Lifecycle:
##'   Stable.
##' @export
simulate.rlmerMod <- function(object, nsim = 1, seed = NULL,
                              type = c("parametric", "wild"),
                              use.u = FALSE, cluster = NULL,
                              wild.weights = c("mammen", "rademacher"),
                              ...) {
    type <- match.arg(type)
    wild.weights <- match.arg(wild.weights)
    .check_number(nsim, "nsim", lower = 1, integer = TRUE)
    .check_flag(use.u, "use.u")
    if (type == "parametric" && !is.null(cluster))
        stop("'cluster' is used only with type = \"wild\"", call. = FALSE)
    if (type == "wild" && !missing(use.u))
        stop("'use.u' is used only with type = \"parametric\"",
             call. = FALSE)
    ## seed handling as in stats:::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    n <- nobs(object)
    if (type == "parametric") {
        pp <- object@pp
        Xb <- as.numeric(.X(object) %*% .fixef(object)) + object@resp$offset
        Z <- t(.Zt(object))
        w <- object@resp$weights
        if (all(w == 1)) w <- NULL
        b <- if (use.u) .b(object) else NULL
        gen <- function()
            .sim_parametric(Xb, Z, pp$U_b, .sigma(object), n, pp$q, w, b)
    } else {
        gen <- .wild_generator(object, cluster, wild.weights)
    }
    val <- as.data.frame(replicate(nsim, gen(), simplify = FALSE),
                         col.names = paste0("sim_", seq_len(nsim)))
    rownames(val) <- rownames(model.frame(object))
    attr(val, "seed") <- RNGstate
    val
}

##' Refit a robust linear mixed model to a new response
##'
##' Fits the model of \code{object} again with the response replaced by
##' \code{newresp}, keeping the psi-functions, the method, the prior weights
##' and the offset stored in \code{object}; the original data is not needed.
##' The fit is started from the estimates of \code{object}. This is the \code{\link[lme4]{refit}}
##' generic of \pkg{lme4}; together with
##' \code{\link[=simulate.rlmerMod]{simulate}} it gives a bootstrap with
##' robust refits.
##'
##' @param object an \code{rlmerMod} object.
##' @param newresp the new response: a numeric vector with one value per
##'   observation used in the fit, or a one-column data frame such as a
##'   column of \code{simulate(object)}.
##' @param ... not used.
##' @return An \code{rlmerMod} object. Its call is that of \code{object}.
##' @note The \pkg{generics} package has a generic of the same name. If it
##'   is attached after \pkg{lme4}, call \code{lme4::refit}.
##' @examples
##' fm <- rlmer(Reaction ~ Days + (1 | Subject), sleepstudy)
##' ystar <- simulate(fm, seed = 1)
##' fixef(refit(fm, ystar$sim_1))
##' @importFrom lme4 refit
##' @method refit rlmerMod
##' @section Lifecycle:
##'   Stable.
##' @export
refit.rlmerMod <- function(object, newresp, ...) {
    if (is.list(newresp) && length(newresp) == 1L) newresp <- newresp[[1L]]
    n <- nobs(object)
    if (!is.numeric(newresp) || length(newresp) != n || anyNA(newresp))
        stop("'newresp' must be a numeric vector with one value per ",
             "observation (", n, "), without missing values", call. = FALSE)
    ## a copy of the fit with the new response; the model, the
    ## psi-functions, the prior weights and the offset are those stored in
    ## the object, so nothing of the original call is evaluated again
    fit <- object
    fit@pp <- object@pp$copy()
    fit@resp <- object@resp$copy()
    fit@resp$y <- as.numeric(newresp)
    fit@frame[[1L]] <- as.numeric(newresp)
    ## re-solve the effects and sigma for the new response at the old
    ## theta: the iterations test the estimating equations at the current
    ## effects and would otherwise stop at once
    setTheta(fit, theta(fit), fit.effects = TRUE, update.sigma = TRUE)
    oi <- object@optinfo
    rel.tol <- if (is.null(oi$rel.tol)) 1e-8 else oi$rel.tol
    max.iter <- if (is.null(oi$max.iter)) 40 * (len(fit, "theta") + 1)^2
                else oi$max.iter
    .rlmer(fit, rel.tol, max.iter, verbose = 0, doFit = TRUE,
           size_obr = isTRUE(oi$size_obr))
}
