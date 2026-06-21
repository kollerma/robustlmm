##' RANSAC-style random subsample initial estimator for linear
##' mixed-effects models.
##'
##' For \code{K} random subsamples of the data, fit a classical
##' \code{\link[lme4]{lmer}} on each subsample, score by a robust
##' scale of residuals computed on the \emph{full} data, and return
##' the lmer fit minimising that score.
##'
##' The motivation is that for redescending psi-functions
##' (e.g. \code{\link{bisquarePsi}}) the \code{\link{rlmer}}
##' optimiser benefits from a starting value close to the true
##' parameters. A bad initial estimate can produce phony local
##' minima (e.g. random-effects correlation pinned at +/- 1; see
##' Koller and Stahel 2022, Section 4.4). RANSAC is a classical
##' way of generating a high-breakdown-point initial estimate by
##' subsampling.
##'
##' @title RANSAC initial estimator for LMM
##' @param formula model formula in \code{\link[lme4]{lmer}} syntax.
##' @param data full data frame.
##' @param K number of random subsamples (default 200).
##' @param sub_frac fraction of the data per subsample (default 0.5).
##' @param scale_fn function from a numeric residual vector to a
##'   scalar scale. Default \code{\link[robustbase]{Qn}}.
##' @param seed optional RNG seed for reproducibility.
##' @param verbose logical; print progress every 50 subsamples.
##' @return list with \code{fit} (best \code{lmerMod}), \code{scale}
##'   (its score), \code{subset} (its row indices), \code{scales}
##'   (vector of all K scores), \code{K}, \code{n_sub}.
##' @examples
##'   set.seed(1)
##'   res <- ransac_lme4(Reaction ~ Days + (Days | Subject),
##'                       data = sleepstudy, K = 30)
##'   res$scale
##' @importFrom robustbase Qn
##' @export
ransac_lme4 <- function(formula, data, K = 200L, sub_frac = 0.5,
                        scale_fn = robustbase::Qn,
                        seed = NULL, verbose = FALSE) {
    if (!is.null(seed)) set.seed(seed)
    if (!is.numeric(K) || length(K) != 1L || K < 1L)
        stop("'K' must be a positive integer scalar")
    if (!is.numeric(sub_frac) || length(sub_frac) != 1L ||
        sub_frac <= 0 || sub_frac > 1)
        stop("'sub_frac' must be in (0, 1]")
    n <- nrow(data)
    n_sub <- ceiling(sub_frac * n)
    yname <- as.character(formula[[2]])

    best_scale  <- Inf
    best_fit    <- NULL
    best_subset <- NULL
    scales <- rep(NA_real_, K)

    for (k in seq_len(K)) {
        idx <- sample.int(n, n_sub)
        sub_data <- data[idx, , drop = FALSE]
        fit_sub <- tryCatch(
            suppressMessages(suppressWarnings(
                lme4::lmer(formula, data = sub_data, REML = TRUE))),
            error = function(e) NULL)
        if (is.null(fit_sub)) next
        pred_full <- tryCatch(stats::predict(fit_sub, newdata = data,
                                              allow.new.levels = TRUE),
                              error = function(e) NULL)
        if (is.null(pred_full)) next
        resid_full <- data[[yname]] - pred_full
        s <- scale_fn(resid_full)
        scales[k] <- s
        if (!is.na(s) && s < best_scale) {
            best_scale  <- s
            best_fit    <- fit_sub
            best_subset <- idx
        }
        if (verbose && k %% 50 == 0)
            message(sprintf("ransac_lme4: k = %d, best scale so far = %.4f",
                            k, best_scale))
    }
    list(fit = best_fit, scale = best_scale, subset = best_subset,
         scales = scales, K = K, n_sub = n_sub)
}

##' Fit \code{\link{rlmer}} with a RANSAC-derived initial estimator.
##'
##' Convenience wrapper that calls \code{\link{ransac_lme4}} to obtain
##' a starting value and then passes it to \code{\link{rlmer}}'s
##' \code{init} argument.
##'
##' @title rlmer with RANSAC initial estimator
##' @param formula,data passed to \code{\link{rlmer}}.
##' @param K,sub_frac,seed passed to \code{\link{ransac_lme4}}.
##' @param ... other arguments to \code{\link{rlmer}} (e.g.
##'   \code{rho.e}, \code{rho.b}, \code{method}).
##' @return rlmerMod object.
##' @examples
##'   \donttest{
##'   fit <- rlmer_ransac(Reaction ~ Days + (Days | Subject),
##'                        data = sleepstudy, K = 30)
##'   }
##' @export
rlmer_ransac <- function(formula, data, K = 200L, sub_frac = 0.5,
                         seed = NULL, ...) {
    init <- ransac_lme4(formula, data, K = K, sub_frac = sub_frac,
                        seed = seed)
    if (is.null(init$fit))
        stop("ransac_lme4: no successful lmer fit across ", K, " subsamples")
    rlmer(formula, data = data, init = init$fit, ...)
}
