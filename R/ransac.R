##' RANSAC-style random subsample initial estimator for linear
##' mixed-effects models.
##'
##' For \code{K} random subsamples of the data, fit a classical
##' \code{\link[lme4]{lmer}} on each subsample, score by a robust
##' scale of residuals computed on the \emph{full} data, and return
##' the lmer fit minimising that score.
##'
##' The motivation is that for redescending psi-functions
##' (e.g. \code{\link{lqqPsi}}, the recommended redescender, or the
##' faster-redescending \code{\link{bisquarePsi}}) the \code{\link{rlmer}}
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
##' @param K maximum number of random subsamples (default 200). With
##'   \code{adaptive = TRUE} the search stops early once the best score
##'   plateaus; \code{K} is then an upper budget rather than a fixed
##'   count.
##' @param sub_frac fraction of the data per subsample (default 0.5).
##' @param scale_fn function from a numeric residual vector to a
##'   scalar scale. Default \code{\link[robustbase]{Qn}}.
##' @param adaptive logical (default \code{TRUE}); stop drawing
##'   subsamples once the best score has not improved by more than
##'   \code{tol} (relative) for \code{patience} consecutive draws, after
##'   at least \code{K_min} draws. \code{FALSE} runs exactly \code{K}
##'   draws (the previous behaviour).
##' @param patience number of consecutive non-improving draws that
##'   triggers the adaptive stop (default 50).
##' @param K_min minimum draws before the adaptive stop can fire
##'   (default 50).
##' @param tol relative-improvement threshold for the adaptive stop
##'   (default 1e-3).
##' @param seed optional RNG seed for reproducibility.
##' @param verbose logical; print progress every 50 subsamples.
##' @return list with \code{fit} (best \code{lmerMod}), \code{scale}
##'   (its score), \code{subset} (its row indices), \code{scales}
##'   (length-\code{K} vector of scores; \code{NA} for draws not run
##'   under an adaptive early stop), \code{K} (the requested cap),
##'   \code{K_used} (draws actually run), \code{n_sub}, \code{n_singular}
##'   (total number of collinear candidate observations skipped by the
##'   nonsingular-subsampling draw across all subsamples; see below).
##' @examples
##'   set.seed(1)
##'   res <- ransac_lme4(Reaction ~ Days + (Days | Subject),
##'                       data = sleepstudy, K = 30)
##'   res$scale
##' @param stratify logical (default \code{TRUE}); draw each subsample
##'   stratified by the random-effects grouping factor (the one with the
##'   most levels, parsed from the formula), taking
##'   \code{ceiling(sub_frac * n_g)} (at least one) rows within each level
##'   so every grouping level is represented. This keeps the per-subsample
##'   \code{lmer} fittable on designs with many small clusters, where plain
##'   random subsampling can leave a level empty. \code{FALSE} (or no
##'   grouping factor) falls back to simple random subsampling.
##' @param n_keep number of distinct best-scoring candidate starts to
##'   return in \code{$candidates} (default 1). The multi-start consensus
##'   of \code{\link{rlmer_ransac}} (\code{n_starts > 1}) uses these to
##'   sample several basins of a redescending psi.
##' @section Nonsingular subsampling:
##'   With categorical predictors a rank-deficient subsample often does
##'   \emph{not} make \code{\link[lme4]{lmer}} error --- it silently drops
##'   the aliased (all-zero) fixed-effect columns of a dropped factor
##'   level --- so a degenerate candidate would otherwise enter the scale
##'   competition with the wrong number of parameters. Rather than draw and
##'   then repair, each subsample is drawn \emph{nonsingular by
##'   construction} using the nonsingular-subsampling algorithm of Koller
##'   and Stahel (2017, Algorithm 1; the method behind
##'   \code{robustbase::\link[robustbase]{lmrob.control}(setting =
##'   "KS2014")}). The draw units (clusters when a grouping factor is
##'   available, else single observations) are randomly permuted; a
##'   Gaxpy-variant LU decomposition with partial pivoting and column
##'   skipping (implemented in C++) then walks the permuted observations and
##'   greedily selects the first \eqn{p} that are linearly independent,
##'   \emph{skipping} any observation collinear with those already chosen.
##'   This yields a full-rank fixed-effects core whenever the full design
##'   has full column rank. Whole clusters carrying that core are the
##'   mandatory seed; further clusters are then added in the permuted order
##'   until the retained data are identifiable --- (i) a full-rank
##'   fixed-effects design [guaranteed by the core], (ii) at least two
##'   observed levels of every random-effects grouping factor, and (iii)
##'   non-constant numeric random-slope variables --- and finally up to the
##'   target subsample size. Because adding clusters never lowers the rank
##'   or removes a level this is monotone and terminates. \code{n_singular}
##'   counts the collinear candidate observations skipped by the LU across
##'   all draws (a measure of the collinearity encountered); it is
##'   \code{0} on a purely continuous, full-rank design, where the algorithm
##'   selects exactly the first \eqn{p} permuted observations, so the draw
##'   is a uniform random cluster subsample --- statistically identical to
##'   plain random subsampling. Note that a nonsingular \emph{start} does
##'   not guarantee the fit stays nonsingular through the \code{\link{rlmer}}
##'   refinement: a redescending psi can zero-weight observations back into
##'   a rank-deficient design (Koller and Stahel 2017, Remark 2), which
##'   \code{\link{rlmer_ransac}} checks for and warns about.
##' @references
##' Koller, M. and Stahel, W. A. (2017) Nonsingular subsampling for
##' regression S estimators with categorical predictors.
##' \emph{Computational Statistics} \bold{32}(2), 631--646.
##' @importFrom robustbase Qn
##' @export
## Nonsingular subsampling (Koller and Stahel 2017). Static design
## information used to detect rank-deficient subsamples: the fixed-effects
## terms (to build the model matrix), the factor predictors among them,
## the random-effects grouping expressions, the numeric random-slope
## variables, and the expected number of fixed-effect coefficients on the
## full (assumed full-rank) data.
.ransac_design_info <- function(formula, data) {
    info <- list(fe_terms = NULL, fe_factors = character(0),
                 groupers = list(), slope_numeric = character(0),
                 p_expected = NA_integer_)
    fe_form <- tryCatch(reformulas::nobars(formula), error = function(e) NULL)
    if (is.null(fe_form)) return(info)
    fe_terms <- tryCatch(stats::delete.response(stats::terms(fe_form, data = data)),
                         error = function(e) NULL)
    if (is.null(fe_terms)) return(info)
    info$fe_terms <- fe_terms
    fe_vars <- all.vars(fe_terms)
    info$fe_factors <- fe_vars[vapply(fe_vars, function(v)
        v %in% names(data) && (is.factor(data[[v]]) || is.character(data[[v]])),
        logical(1))]
    info$p_expected <- tryCatch(ncol(stats::model.matrix(fe_terms, data)),
                                error = function(e) NA_integer_)
    bars <- tryCatch(reformulas::findbars(formula), error = function(e) NULL)
    if (!is.null(bars) && length(bars)) {
        info$groupers <- lapply(bars, function(b) b[[3L]])
        info$slope_numeric <- unique(unlist(lapply(bars, function(b) {
            vs <- all.vars(b[[2L]])
            vs[vapply(vs, function(v)
                v %in% names(data) && is.numeric(data[[v]]), logical(1))]
        }), use.names = FALSE))
    }
    info
}

## TRUE/FALSE whether the subsample `idx` is nonsingular, plus the factor
## predictors that lost a level (used to steer the repair). Consumes no
## RNG, so a full-rank draw leaves the subsample and RNG stream untouched.
.ransac_check_subsample <- function(idx, data, info) {
    sub <- data[idx, , drop = FALSE]
    ok <- TRUE
    missing_fac <- character(0)
    if (!is.null(info$fe_terms)) {
        X <- tryCatch(stats::model.matrix(info$fe_terms, sub),
                      error = function(e) NULL)
        if (is.null(X) || ncol(X) == 0L) {
            ok <- FALSE
        } else {
            r <- tryCatch(qr(X)$rank, error = function(e) NA_integer_)
            if (is.na(r) || r < ncol(X)) ok <- FALSE
        }
    }
    for (v in info$fe_factors) {
        if (is.factor(data[[v]]) &&
            length(unique(sub[[v]])) < nlevels(data[[v]])) {
            missing_fac <- c(missing_fac, v)
            ok <- FALSE
        }
    }
    for (g in info$groupers) {
        gv <- tryCatch(eval(g, envir = sub), error = function(e) NULL)
        if (is.null(gv) || length(unique(gv)) < 2L) ok <- FALSE
    }
    for (v in info$slope_numeric) {
        if (v %in% names(sub) && is.numeric(sub[[v]]) &&
            length(unique(sub[[v]])) < 2L) ok <- FALSE
    }
    list(ok = ok, missing_fac = missing_fac)
}

## Draw a nonsingular subsample by construction (Koller and Stahel 2017,
## Algorithm 1) rather than draw-then-check-then-repair. Consumes the
## ambient RNG set by the caller's `set.seed`; does NOT set a seed itself.
##
## Draw units are clusters when a grouping factor is available (`strata`,
## the finest grouping factor from .ransac_strata), else single
## observations. The clusters (or observations) are randomly permuted; the
## observation walk-order is the observations grouped in permuted-cluster
## order. On that walk order the Gaxpy-LU column-skipping routine
## `nonsingularSubsampleLU` selects a full-rank fixed-effects CORE (the
## first p linearly independent observations); the clusters carrying that
## core are the mandatory seed set. Whole clusters are then added in
## permuted order until the model is identifiable on the retained data --
## (a) X full column rank [guaranteed by the core], (b) every grouping
## factor has >= 2 retained levels, (c) every random-slope covariate varies
## -- and finally up to the same target size the old draw used. Because
## adding clusters is monotone (never lowers rank or removes a level) this
## terminates whenever the full design is identifiable. When there is no
## collinearity the core is the first p observations, so the result is the
## first ceil(sub_frac * .) clusters of the permutation, i.e. a uniform
## random cluster subsample -- statistically identical to plain random
## subsampling (the paper's guarantee).
.ransac_nonsingular_subsample <- function(formula, data, sub_frac, strata,
                                          design_info) {
    n <- nrow(data)

    ## -- draw units and the observation walk order (consumes RNG) --------
    if (!is.null(strata)) {
        cl <- as.factor(strata)
        units <- split(seq_len(n), cl)              # obs indices per cluster
        perm  <- sample(names(units))               # permuted cluster ids
        tab    <- as.integer(table(cl))
        target <- sum(pmax(1L, ceiling(sub_frac * tab)))
    } else {
        units <- stats::setNames(as.list(seq_len(n)), as.character(seq_len(n)))
        perm  <- as.character(sample.int(n))        # permuted observations
        target <- ceiling(sub_frac * n)
    }
    walk <- unlist(units[perm], use.names = FALSE)  # observation walk order

    ## -- nonsingular fixed-effects core (Algorithm 1) --------------------
    ## Skip the LU when the design cannot be collinear across subsamples
    ## (intercept-only or purely numeric, full-rank almost surely): every
    ## subsample is then trivially full-rank and the core is empty.
    n_skipped <- 0L
    seed_units <- character(0)
    if (!is.null(design_info$fe_terms) && length(design_info$fe_factors)) {
        Xt <- tryCatch(t(stats::model.matrix(design_info$fe_terms, data)),
                       error = function(e) NULL)
        if (!is.null(Xt) && nrow(Xt) >= 1L) {
            lu <- nonsingularSubsampleLU(Xt, as.integer(walk - 1L), 1e-7)
            n_skipped <- as.integer(lu$n_skipped)
            if (isTRUE(lu$singular)) {
                ## full design rank-deficient even using all data: the
                ## model is not identifiable (should not happen for a
                ## validly-specified model). Report failure.
                return(list(idx = sort(walk), n_singular = n_skipped,
                            ok = FALSE))
            }
            core_obs <- lu$selected + 1L            # 1-based core observations
            ## map core observations to their draw units (clusters)
            if (!is.null(strata)) {
                seed_units <- unique(as.character(cl[core_obs]))
            } else {
                seed_units <- as.character(core_obs)
            }
        }
    }

    ## -- assemble: mandatory seed clusters first, then permuted order ----
    is_seed  <- perm %in% seed_units
    seed_ord <- perm[is_seed]                       # in permuted order
    rest_ord <- perm[!is_seed]

    idx <- integer(0)
    ## phase 1: all seed clusters (guarantees the full-rank core)
    for (u in seed_ord) idx <- c(idx, units[[u]])
    ## phase 2: add whole clusters until identifiable AND at target size
    for (u in rest_ord) {
        if (length(idx) >= target &&
            .ransac_check_subsample(idx, data, design_info)$ok)
            break
        idx <- c(idx, units[[u]])
    }
    ok <- .ransac_check_subsample(idx, data, design_info)$ok
    list(idx = sort(idx), n_singular = n_skipped, ok = ok)
}

ransac_lme4 <- function(formula, data, K = 200L, sub_frac = 0.5,
                        scale_fn = robustbase::Qn,
                        adaptive = TRUE, patience = 50L, K_min = 50L,
                        tol = 1e-3, stratify = TRUE, n_keep = 1L,
                        seed = NULL, verbose = FALSE) {
    if (!is.null(seed)) set.seed(seed)
    if (!is.numeric(K) || length(K) != 1L || K < 1L)
        stop("'K' must be a positive integer scalar")
    if (!is.numeric(sub_frac) || length(sub_frac) != 1L ||
        sub_frac <= 0 || sub_frac > 1)
        stop("'sub_frac' must be in (0, 1]")
    if (!is.numeric(n_keep) || length(n_keep) != 1L || n_keep < 1L)
        stop("'n_keep' must be a positive integer scalar")
    n_keep <- as.integer(n_keep)
    n <- nrow(data)
    yname <- as.character(formula[[2]])

    ## WS12: stratify subsamples by the (finest) grouping factor so every
    ## random-effect level is represented in each draw.
    strata <- if (stratify) .ransac_strata(formula, data) else NULL
    use_stratify <- !is.null(strata)
    if (use_stratify) {
        tab   <- as.integer(table(strata))
        n_sub <- sum(pmax(1L, ceiling(sub_frac * tab)))
    } else {
        n_sub <- ceiling(sub_frac * n)
    }

    best_scale   <- Inf
    best_fit     <- NULL
    best_subset  <- NULL
    scales       <- rep(NA_real_, K)
    last_improve <- 0L     # draw index of the last > tol improvement
    K_used       <- K
    topk         <- list()  # WS12-D: top-n_keep candidates (best first)
    n_singular   <- 0L      # KS2017: draws skipped/discarded as deficient

    ## Nonsingular-subsampling design info (Koller and Stahel 2017): the
    ## fixed-effects terms / factors / grouping structure used to draw a
    ## full-rank subsample by construction.
    design_info <- .ransac_design_info(formula, data)

    for (k in seq_len(K)) {
        ## KS2017 Algorithm 1: draw a subsample that is nonsingular by
        ## construction, so no rank-deficient (degenerate) candidate enters
        ## the scale competition. `n_singular` accumulates the collinear
        ## candidate observations the LU skipped while building the draw.
        draw <- .ransac_nonsingular_subsample(formula, data, sub_frac,
                                              strata, design_info)
        idx <- draw$idx
        n_singular <- n_singular + draw$n_singular
        if (!draw$ok) next
        sub_data <- data[idx, , drop = FALSE]
        fit_sub <- tryCatch(
            suppressMessages(suppressWarnings(
                lme4::lmer(formula, data = sub_data, REML = TRUE))),
            error = function(e) NULL)
        ## KS2017: reject (and count) a fit whose lmer nonetheless dropped
        ## a fixed-effect coefficient (aliased column silently removed).
        if (!is.null(fit_sub) && !is.na(design_info$p_expected) &&
            length(lme4::fixef(fit_sub)) != design_info$p_expected) {
            n_singular <- n_singular + 1L
            fit_sub <- NULL
        }
        if (!is.null(fit_sub)) {
            pred_full <- tryCatch(stats::predict(fit_sub, newdata = data,
                                                 allow.new.levels = TRUE),
                                  error = function(e) NULL)
            if (!is.null(pred_full)) {
                resid_full <- data[[yname]] - pred_full
                s <- scale_fn(resid_full)
                scales[k] <- s
                if (!is.na(s) && s < best_scale * (1 - tol)) {
                    last_improve <- k          # a real (> tol) improvement
                }
                if (!is.na(s) && s < best_scale) {
                    best_scale  <- s
                    best_fit    <- fit_sub
                    best_subset <- idx
                }
                ## WS12-D: keep the top n_keep candidates for multi-start
                ## consensus (bounded buffer; only when more than one is
                ## requested).
                if (n_keep > 1L && !is.na(s)) {
                    topk[[length(topk) + 1L]] <-
                        list(scale = s, fit = fit_sub, subset = idx)
                    if (length(topk) > n_keep) {
                        ord  <- order(vapply(topk, `[[`, numeric(1), "scale"))
                        topk <- topk[ord[seq_len(n_keep)]]
                    }
                }
            }
        }
        if (verbose && k %% 50 == 0)
            message(sprintf(paste0("ransac_lme4: k = %d, best scale so far ",
                                   "= %.4f, n_singular = %d"),
                            k, best_scale, n_singular))
        ## adaptive early stop: best score has plateaued
        if (adaptive && k >= K_min && is.finite(best_scale) &&
            (k - last_improve) >= patience) {
            K_used <- k
            break
        }
    }
    basin_radius <- if (!is.null(best_fit))
        tryCatch(as.numeric(ransac_basin_radius(best_fit)),
                 error = function(e) NA_real_) else NA_real_
    ## WS12-D: assemble the candidate list (best first), dropping
    ## near-duplicate starts (same fixed effects) so the consensus samples
    ## distinct starts. n_keep == 1 -> just the best fit.
    candidates <- if (n_keep > 1L && length(topk)) {
        ord  <- order(vapply(topk, `[[`, numeric(1), "scale"))
        .dedup_by_fixef(lapply(topk[ord], `[[`, "fit"))
    } else if (!is.null(best_fit)) list(best_fit) else list()
    if (verbose)
        message(sprintf("ransac_lme4: done, K_used = %d, n_singular = %d",
                        K_used, n_singular))
    list(fit = best_fit, scale = best_scale, subset = best_subset,
         scales = scales, K = K, K_used = K_used, n_sub = n_sub,
         n_singular = n_singular, stratified = use_stratify,
         basin_radius = basin_radius, candidates = candidates)
}

## WS12-D: drop lmer fits whose fixed effects duplicate an already-kept
## one (relative L2 distance < eps), keeping the first (best-scoring).
.dedup_by_fixef <- function(fits) {
    keep <- list(); betas <- list()
    for (f in fits) {
        b <- tryCatch(as.numeric(lme4::fixef(f)), error = function(e) NULL)
        if (is.null(b)) next
        dup <- any(vapply(betas, function(bb)
            sqrt(sum((b - bb)^2)) <= 1e-6 * (1 + sqrt(sum(bb^2))),
            logical(1)))
        if (!dup) { keep[[length(keep) + 1L]] <- f
                    betas[[length(betas) + 1L]] <- b }
    }
    keep
}

## WS12: the grouping factor to stratify RANSAC subsamples by. Parses the
## random-effects bars of the formula, evaluates each grouping expression
## in the data (so nesting / interactions like school:class work), and
## returns the factor with the most levels (the finest -- hardest to cover
## by chance), or NULL when there is no usable grouping factor.
.ransac_strata <- function(formula, data) {
    bars <- tryCatch(reformulas::findbars(formula), error = function(e) NULL)
    if (is.null(bars) || !length(bars)) return(NULL)
    facs <- lapply(bars, function(b)
        tryCatch(as.factor(eval(b[[3L]], envir = data)),
                 error = function(e) NULL))
    keep <- vapply(facs, function(f)
        !is.null(f) && length(f) == nrow(data) && nlevels(f) >= 2L,
        logical(1))
    facs <- facs[keep]
    if (!length(facs)) return(NULL)
    facs[[which.max(vapply(facs, nlevels, integer(1)))]]
}

##' Support-preservation (basin) radius for a fixed-effects start.
##'
##' Computes the radius \eqn{r^\star(c) = c\,\sigma / (2 \max_j
##' \|x_j\|)} of the ball around the initial fixed-effects estimate
##' within which every observation keeps its redescender-support
##' membership, so the population Hessian stays positive definite (Koller
##' and Stahel; the RANSAC-RSE basin theorem). Here \eqn{c} is the
##' \emph{rejection point} of the redescending \eqn{\psi} --- the smallest
##' \eqn{x > 0} with \eqn{\psi(x) = 0}. For the bisquare this is the tuning
##' cutoff (default \code{4.685}); the geometry generalises to any
##' finite-rejection-point redescender (e.g. \code{\link{lqqPsi}}) by
##' finding its rejection point numerically. A redescending \eqn{\psi} is
##' safe to engage from a start only if the eventual estimate stays within
##' this radius of the (high-breakdown) start --- otherwise the iteration
##' may leave the basin and converge to a phony solution.
##'
##' @title Redescender basin (support-preservation) radius
##' @param object a fitted \code{lmerMod} or \code{rlmerMod} supplying
##'   the error scale \eqn{\sigma} and the fixed-effects design matrix.
##' @param cc the rejection point \eqn{c} of the redescender. If
##'   \code{NULL} and \code{rho} is supplied, \eqn{c} is the rejection
##'   point of \code{rho} found numerically; if both are \code{NULL}, the
##'   bisquare default \code{4.685}.
##' @param rho an optional redescending \code{psi_func_rcpp} (e.g.
##'   \code{\link{bisquarePsi}}, \code{\link{lqqPsi}}) whose rejection
##'   point supplies \eqn{c} when \code{cc} is \code{NULL}.
##' @return the basin radius \eqn{r^\star(c)} in the units of
##'   \eqn{\beta} (a scalar), with attribute \code{"max_xnorm"}.
##' @export
ransac_basin_radius <- function(object, cc = NULL, rho = NULL) {
    if (is.null(cc)) {
        rp <- .psi_rejection_point(rho)
        cc <- if (!is.null(rp) && is.finite(rp)) rp else 4.685
    }
    sigma_e <- tryCatch(as.numeric(sigma(object)), error = function(e) NA_real_)
    X <- tryCatch(as.matrix(getME(object, "X")), error = function(e) NULL)
    if (is.null(X) || !is.finite(sigma_e))
        stop("ransac_basin_radius: cannot extract sigma / design matrix ",
             "from 'object'.")
    max_xnorm <- max(sqrt(rowSums(X^2)))
    r <- cc * sigma_e / (2 * max_xnorm)
    attr(r, "max_xnorm") <- max_xnorm
    r
}

## The rejection point of a redescending psi: the smallest x > 0 with
## psi(x) == 0 (beyond the psi peak). Found numerically on a grid, then
## refined by bisection. Returns NULL when rho is not a redescender or has
## no finite rejection point (monotone / unbounded psi) within `upper`.
.psi_rejection_point <- function(rho, upper = 100, tol = 1e-8) {
    if (is.null(rho) || !methods::is(rho, "psi_func_rcpp")) return(NULL)
    f <- function(x) abs(tryCatch(rho@psi(x), error = function(e) NA_real_))
    xs <- seq(0, upper, length.out = 20001L)
    v  <- f(xs)
    pos <- which(v > tol)
    if (!length(pos)) return(NULL)
    last_pos <- max(pos)
    if (last_pos >= length(xs)) return(NULL)   # no finite rejection point
    lo <- xs[last_pos]; hi <- xs[last_pos + 1L]
    for (i in seq_len(60L)) {
        mid <- (lo + hi) / 2
        if (f(mid) > tol) lo <- mid else hi <- mid
    }
    hi
}

## TRUE if a psi_func is redescending (psi -> 0 in the tails), as opposed
## to monotone-bounded (Huber: psi -> k) or unbounded (classical).
.is_redescending <- function(rho) {
    if (is.null(rho) || !methods::is(rho, "psi_func_rcpp")) return(FALSE)
    big <- tryCatch(abs(rho@psi(1e3)), error = function(e) NA_real_)
    isTRUE(big < 1e-8)
}

## the redescender tuning constant c of a psi_func (tDefs "c", else "k",
## else the first), or NULL.
.rho_tuning_c <- function(rho) {
    if (is.null(rho) || !methods::is(rho, "psi_func_rcpp")) return(NULL)
    td <- rho@tDefs
    if (!length(td)) return(NULL)
    if (!is.null(names(td)) && "c" %in% names(td)) return(unname(td["c"]))
    if (!is.null(names(td)) && "k" %in% names(td)) return(unname(td["k"]))
    unname(td[1L])
}

## WS12 basin gate: when a redescending rho.e is engaged from a RANSAC
## start, warn if the converged FIXED EFFECTS left the bisquare support-
## preservation radius r*(c) of the start. NB this is the beta-side
## (Hessian support-preservation) condition only; it is a WEAK detector of
## the random-effects phony correlation (P(warn | phony) ~ 0.17 in
## inst/simulationStudy/ransacBasin.R) -- the rho -> +-1 attractor lives
## in a different (RE-covariance) basin and is checked separately by
## .ransac_check_phony.
.ransac_check_basin <- function(fit, start_fit, rho_e) {
    if (!.is_redescending(rho_e)) return(invisible(NULL))
    ## Use the rejection point of rho_e (the general geometry); for the
    ## bisquare this is its cutoff c, for lqq etc. it is found numerically.
    r  <- tryCatch(ransac_basin_radius(start_fit, rho = rho_e),
                   error = function(e) NULL)
    if (is.null(r)) return(invisible(NULL))
    b0 <- tryCatch(lme4::fixef(start_fit), error = function(e) NULL)
    b1 <- tryCatch(.fixef(fit), error = function(e) NULL)
    if (is.null(b0) || is.null(b1) || length(b0) != length(b1))
        return(invisible(NULL))
    d <- sqrt(sum((b1 - b0)^2))
    if (d > r)
        warning(sprintf(
            paste0("RANSAC + redescending rho: the fitted fixed effects ",
                   "moved %.3g from the high-breakdown start, beyond the ",
                   "support-preservation radius r* = %.3g (from the ",
                   "rejection point of rho.e), so the redescender support ",
                   "membership of every observation is no ",
                   "longer guaranteed preserved (the population Hessian may ",
                   "be indefinite on the path from the start). This is the ",
                   "fixed-effect (beta) basin condition; it does NOT by ",
                   "itself certify the random-effects covariance (see the ",
                   "separate phony-correlation check). Consider a monotone ",
                   "rho (e.g. smoothPsi), more RANSAC draws (raise K), or ",
                   "verify against confint(., method = \"boot\")."),
            d, as.numeric(r)), call. = FALSE)
    invisible(list(distance = d, radius = as.numeric(r), in_basin = d <= r))
}

## WS12 follow-up (b): direct detector of the random-effects phony
## correlation -- the rho -> +-1 boundary attractor (Koller and Stahel
## 2022, Section 4.4) that the beta support-preservation radius does NOT
## detect. Returns the maximum absolute random-effect correlation of the
## fit; NA when there is no correlation parameter (diagonal V_b).
.re_max_abscor <- function(fit) {
    vc <- tryCatch(VarCorr(fit), error = function(e) NULL)
    if (is.null(vc)) return(NA_real_)
    cors <- unlist(lapply(vc, function(b) {
        cm <- attr(b, "correlation")
        if (is.null(cm) || nrow(cm) < 2L) return(numeric(0))
        abs(cm[lower.tri(cm)])
    }), use.names = FALSE)
    if (!length(cors)) return(NA_real_)
    max(cors)
}

## Warn when the fitted random-effects covariance is (near-)singular --
## |rho-hat| close to 1 -- i.e. the phony-correlation failure, regardless
## of the rho-function. This is the RE-side check the beta basin radius
## cannot supply; fires for any rlmer fit, not only redescenders.
.ransac_check_phony <- function(fit, threshold = 0.99) {
    rho <- .re_max_abscor(fit)
    if (!is.finite(rho)) return(invisible(NULL))
    if (rho > threshold)
        warning(sprintf(
            paste0("Phony random-effects correlation: the fitted RE ",
                   "covariance is near-singular (max |rho-hat| = %.4f > ",
                   "%.2f), the rho -> +-1 boundary attractor of Koller and ",
                   "Stahel (2022, Sec. 4.4). The variance components and ",
                   "any inference depending on them are unreliable; refit ",
                   "from a different start (raise the RANSAC K), use a ",
                   "monotone rho, or drop the random-effects correlation ",
                   "term."), rho, threshold), call. = FALSE)
    invisible(list(max_abscor = rho, phony = rho > threshold,
                   threshold = threshold))
}

## Refinement-singularity guard (Koller and Stahel 2017, Remark 2, flavour
## (i)). A nonsingular subsample only guarantees a full-rank design at the
## START. During the rlmer refinement a redescending psi can assign
## (near-)zero robustness weight to observations and so collapse the
## effective (positive-weight) fixed-effects design to less than full rank
## -- the fit is then no longer identified on the observations it actually
## uses. This inspects the e-side robustness weights of a converged fit and
## reports the rank of the design restricted to the observations with
## weight > eps.
.ransac_refinement_rank <- function(fit, eps = 1e-4) {
    w <- tryCatch(as.numeric(wgt.e(fit)), error = function(e) NULL)
    X <- tryCatch(as.matrix(getME(fit, "X")), error = function(e) NULL)
    if (is.null(w) || is.null(X) || length(w) != nrow(X))
        return(NULL)
    p   <- ncol(X)
    pos <- is.finite(w) & w > eps
    Xp  <- X[pos, , drop = FALSE]
    r   <- if (nrow(Xp) == 0L) 0L
           else tryCatch(qr(Xp)$rank, error = function(e) NA_integer_)
    list(rank = r, p = p, n_pos = sum(pos), n = length(w),
         singular = is.na(r) || r < p)
}

## Warn (and return the diagnosis) when the refinement collapsed the
## positive-weight design to rank-deficiency. Called on the converged fit.
.ransac_check_refinement <- function(fit, eps = 1e-4, tries = NULL) {
    d <- .ransac_refinement_rank(fit, eps)
    if (is.null(d)) return(invisible(NULL))
    if (isTRUE(d$singular)) {
        tried <- if (is.null(tries) || tries <= 1L) "" else
                 sprintf(" %d re-seeded RANSAC starts all collapsed the same way;",
                         tries)
        warning(sprintf(
            paste0("RANSAC + redescending rho: the refinement zero-weighted ",
                   "observations into a rank-deficient design -- only %d of ",
                   "%d observations keep a robustness weight > %g and their ",
                   "fixed-effects design has rank %s < %d (Koller and Stahel ",
                   "2017, Remark 2).%s the fit is not identified on the ",
                   "observations it uses; its fixed effects and any inference ",
                   "are unreliable. Raise max_tries / K / n_starts, or use a ",
                   "psi-function that always gives positive weights (e.g. ",
                   "rho.e = smoothPsi) so no observation is fully rejected."),
            d$n_pos, d$n, eps, as.character(d$rank), d$p, tried), call. = FALSE)
    }
    invisible(d)
}

##' Fit \code{\link{rlmer}} with a RANSAC-derived initial estimator.
##'
##' Convenience wrapper that calls \code{\link{ransac_lme4}} to obtain
##' a starting value and then passes it to \code{\link{rlmer}}'s
##' \code{init} argument.
##'
##' With \code{n_starts > 1} it runs a multi-start consensus:
##' \code{rlmer} is fitted from each of the \code{n_starts} best distinct
##' RANSAC candidate starts, and the returned fit is the
##' lowest-residual-scale one whose random-effects covariance is interior
##' (\eqn{|\hat\rho| \le} \code{phony_threshold}). This samples several
##' basins of a redescending \eqn{\psi} and so recovers the interior
##' solution when the single best start happens to fall into the phony
##' \eqn{|\hat\rho| \to 1} attractor. If every start lands phony, the
##' best-scoring fit is returned with a warning. The per-start summary is
##' attached as \code{attr(fit, "consensus")}.
##'
##' The nonsingular subsample guarantees a full-rank fixed-effects design
##' at the \emph{start} only; a redescending \eqn{\psi} can zero-weight
##' observations during the refinement and collapse the positive-weight
##' design to rank-deficiency (Koller and Stahel 2017, Remark 2). After the
##' fit converges its e-side robustness weights are inspected and a warning
##' is issued if the design restricted to the positively-weighted
##' observations is rank-deficient, suggesting a positive-weight
##' \eqn{\psi} (e.g. \code{rho.e = smoothPsi}) or a different start.
##'
##' @title rlmer with RANSAC initial estimator
##' @param formula,data passed to \code{\link{rlmer}}.
##' @param K,sub_frac,seed passed to \code{\link{ransac_lme4}}.
##' @param n_starts number of distinct RANSAC starts for the multi-start
##'   consensus (default 1 = single best start, the previous behaviour).
##' @param phony_threshold a fit is treated as phony (non-interior) when
##'   its maximum \eqn{|\hat\rho|} exceeds this (default 0.99).
##' @param max_tries maximum number of RANSAC re-seeds when a redescending
##'   \code{rho.e} zero-weights the refinement into a rank-deficient
##'   positive-weight design (Koller and Stahel 2017, Remark 2). The fit is
##'   re-drawn and re-fitted from a fresh nonsingular start until its
##'   positive-weight fixed-effects design is full rank, up to
##'   \code{max_tries} (default 5); if all attempts still collapse the last
##'   fit is returned with a warning. With a fixed \code{seed}, attempt
##'   \eqn{t} uses \code{seed + t - 1}, so the first attempt is reproducible
##'   and the retries are deterministic.
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
                         n_starts = 1L, phony_threshold = 0.99,
                         seed = NULL, max_tries = 5L, ...) {
    if (!is.numeric(n_starts) || length(n_starts) != 1L || n_starts < 1L)
        stop("'n_starts' must be a positive integer scalar")
    if (!is.numeric(max_tries) || length(max_tries) != 1L || max_tries < 1L)
        stop("'max_tries' must be a positive integer scalar")
    max_tries <- as.integer(max_tries)
    rho_e <- list(...)[["rho.e"]]

    ## One RANSAC -> rlmer attempt from a fresh draw. Returns
    ## list(fit, start, all_phony, n_fits, min_rho) or NULL, and attaches
    ## the consensus table for the multi-start path. Emits no informational
    ## warnings so that re-seeded retries below do not spam them.
    one_try <- function(try_seed) {
        init <- ransac_lme4(formula, data, K = K, sub_frac = sub_frac,
                            n_keep = as.integer(n_starts), seed = try_seed)
        if (is.null(init$fit)) return(NULL)
        if (n_starts == 1L) {
            fit <- tryCatch(rlmer(formula, data = data, init = init$fit, ...),
                            error = function(e) NULL)
            if (is.null(fit)) return(NULL)
            return(list(fit = fit, start = init$fit, all_phony = FALSE))
        }
        starts <- init$candidates
        fits <- lapply(starts, function(st)
            tryCatch(suppressWarnings(rlmer(formula, data = data, init = st, ...)),
                     error = function(e) NULL))
        ok <- !vapply(fits, is.null, logical(1))
        if (!any(ok)) return(NULL)
        fits <- fits[ok]; starts <- starts[ok]
        rho   <- vapply(fits, .re_max_abscor, numeric(1))
        scale <- vapply(fits, function(f)
            tryCatch(robustbase::Qn(residuals(f)), error = function(e) NA_real_),
            numeric(1))
        interior <- is.finite(rho) & rho <= phony_threshold
        ## prefer an interior solution; among the eligible, the lowest robust
        ## residual scale. If none interior, fall back to best scale overall.
        pool <- if (any(interior)) which(interior) else seq_along(fits)
        chosen <- pool[which.min(scale[pool])]
        fit <- fits[[chosen]]
        consensus <- data.frame(start = seq_along(fits), max_abscor = rho,
                                resid_scale = scale, interior = interior,
                                chosen = FALSE)
        consensus$chosen[chosen] <- TRUE
        attr(fit, "consensus") <- consensus
        list(fit = fit, start = starts[[chosen]], all_phony = !any(interior),
             n_fits = length(fits), min_rho = min(rho, na.rm = TRUE))
    }

    ## KS2017 Remark 2 (i): the redescender may zero-weight observations
    ## into a rank-deficient design during refinement. Re-seed and refit
    ## from a fresh nonsingular RANSAC draw until the positive-weight design
    ## is full rank, up to max_tries. A given seed makes try t reproducible
    ## (seed + t - 1); with seed = NULL each try draws fresh from the RNG.
    last <- NULL; used <- 0L
    for (t in seq_len(max_tries)) {
        used <- t
        try_seed <- if (is.null(seed)) NULL else as.integer(seed) + (t - 1L)
        cand <- tryCatch(one_try(try_seed), error = function(e) NULL)
        if (is.null(cand)) next
        last <- cand
        d <- .ransac_refinement_rank(cand$fit)
        if (is.null(d) || !isTRUE(d$singular)) break  # identified: accept
    }
    if (is.null(last))
        stop("ransac_lme4: no successful rlmer fit across ", K, " subsamples",
             if (max_tries > 1L) paste0(" in ", used, " tries") else "", ".")

    fit <- last$fit
    ## informational post-fit checks on the returned fit (basin radius for a
    ## redescending rho.e; phony near-singular RE correlation).
    .ransac_check_basin(fit, last$start, rho_e)
    if (n_starts == 1L) {
        .ransac_check_phony(fit)
    } else if (isTRUE(last$all_phony)) {
        warning(sprintf(
            paste0("rlmer_ransac consensus: all %d RANSAC starts converged ",
                   "to a phony random-effects correlation (min |rho-hat| = ",
                   "%.4f > %.2f). The data may not support the random-effects ",
                   "correlation; consider dropping it or a monotone rho."),
            last$n_fits, last$min_rho, phony_threshold), call. = FALSE)
    }
    ## refinement guard: warn only if still rank-deficient after all retries.
    .ransac_check_refinement(fit, tries = used)
    fit
}
