## Breakdown ceiling of the fixed-effects design.
##
## For a regression-equivariant estimator the finite-sample breakdown
## point is at most floor((n - k + 1)/2)/n, where k is the largest number
## of carrier rows x_i lying on one hyperplane through the origin
## (Mili and Coakley 1996). In general position k = p - 1; categorical
## predictors make k large: with an intercept and a balanced dummy, the
## rows of either level lie on one hyperplane, k = n/2, and no estimator
## can exceed breakdown 1/4.
##
## Yohai's MM theorem inherits the breakdown of the initial estimate only
## up to min(b, 1 - b - k/n) (b the constant of the S-scale; the cap is
## machine-checked in the robustlmm theory notes, Topic3.YohaiMM,
## yohai_bound), which is vacuous at b = 1/2 once k/n >= 1/2; it is
## maximised at b = (n - k)/(2n), where it matches the ceiling above.
##
## k is found by enumerating the hyperplanes spanned by p - 1 linearly
## independent distinct carrier rows: when X has full rank the maximum is
## attained on such a hyperplane (extend the spanning set of a maximising
## hyperplane by further carriers). Factor-heavy designs, the only ones
## where k is large, have few distinct rows, so this is exact there. Above
## `max_subsets` subsets a random sample of them is used and k is only a
## lower bound (the ceilings are then upper bounds); the result says so.

.design_k <- function(X, max_subsets = 2e5, tol = 1e-8) {
    X <- as.matrix(X)
    n <- nrow(X); p <- ncol(X)
    key <- apply(signif(X, 12), 1L, paste, collapse = "\r")
    ukey <- unique(key)
    U <- X[match(ukey, key), , drop = FALSE]
    cnt <- as.numeric(table(factor(key, levels = ukey)))
    nrm <- sqrt(rowSums(U^2))
    zero <- nrm <= tol
    if (p == 1L) return(list(k = sum(cnt[zero]), exact = TRUE))
    D <- nrow(U)
    nsub <- choose(D, p - 1L)
    exact <- nsub <= max_subsets
    count_on <- function(rows) {
        B <- U[rows, , drop = FALSE]
        if (qr(B)$rank < p - 1L) return(-1)
        v <- qr.Q(qr(t(B)), complete = TRUE)[, p]     # normal of span(B)
        on <- abs(U %*% v) <= tol * pmax(nrm, 1)
        sum(cnt[on])
    }
    best <- sum(cnt[zero])
    if (exact) {
        cmb <- utils::combn(D, p - 1L)
        for (s in seq_len(ncol(cmb))) best <- max(best, count_on(cmb[, s]))
    } else {
        ## random subsets from a fixed seed, leaving the caller's RNG stream
        ## untouched (the result is then also reproducible)
        had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
        if (had) old <- get(".Random.seed", envir = globalenv())
        on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
                else rm(".Random.seed", envir = globalenv()))
        set.seed(1L)
        for (s in seq_len(max_subsets))
            best <- max(best, count_on(sample.int(D, p - 1L)))
    }
    list(k = best, exact = exact)
}

##' Breakdown ceiling implied by the fixed-effects design.
##'
##' Computes \eqn{k}, the largest number of rows of the fixed-effects
##' model matrix that lie on one hyperplane through the origin, and from it
##' \itemize{
##'   \item \code{max_bp}: the largest finite-sample breakdown point any
##'     regression-equivariant estimator of the fixed effects can have,
##'     \eqn{\lfloor (n - k + 1)/2 \rfloor / n}{floor((n - k + 1)/2)/n}
##'     (Mili and Coakley 1996);
##'   \item \code{mm_bp}: the breakdown point that an MM-type refinement is
##'     guaranteed to inherit from its start, \eqn{\min(b, 1 - b - k/n)}
##'     at \eqn{b = 1/2} (following Yohai 1987), and \code{b_opt}
##'     \eqn{= (n - k)/(2n)}, the scale constant at which that guarantee
##'     reaches \code{max_bp}.
##' }
##' In general position \eqn{k = p - 1} and the ceiling is about 1/2.
##' Categorical predictors raise \eqn{k}: with an intercept and a balanced
##' two-level factor, \eqn{k = n/2}, the ceiling is about 1/4 and the MM
##' guarantee at \eqn{b = 1/2} is empty. No start, including
##' \code{\link{rlmer_ransac}}, can do better than this ceiling for the
##' fixed effects. \code{b_opt} is not used by \code{\link{ransac_lme4}};
##' see the section \dQuote{Ranking RANSAC candidates} of
##' \code{vignette("simulationStudies")}.
##'
##' @title Breakdown ceiling of the fixed-effects design
##' @param object an \code{rlmerMod} or \code{merMod} fit, or a numeric
##'   model matrix.
##' @param max_subsets budget for the exact enumeration of hyperplanes;
##'   above it a random search is used and \code{k} is a lower bound
##'   (\code{exact = FALSE}), so the reported ceilings are upper bounds.
##' @return a list with \code{n}, \code{p}, \code{k}, \code{k_share}
##'   (\eqn{k/n}), \code{max_bp}, \code{mm_bp}, \code{b_opt} and
##'   \code{exact}.
##' @references Mili L, Coakley CW (1996). Robust estimation in structured
##'   linear regression. \emph{The Annals of Statistics} 24(6), 2593--2607.
##' @references Yohai VJ (1987). High breakdown-point and high efficiency
##'   robust estimates for regression. \emph{The Annals of Statistics}
##'   15(2), 642--656.
##' @examples
##' d <- data.frame(y = rnorm(40), g = gl(2, 20), x = rnorm(40))
##' design_breakdown(model.matrix(~ g, d))       # dummy: ceiling ~ 1/4
##' design_breakdown(model.matrix(~ x, d))       # general position
##' @section Lifecycle:
##'   \strong{Experimental.} The interface and the results may change in a
##'   later version.
##' @export
design_breakdown <- function(object, max_subsets = 2e5) {
    X <- if (is.matrix(object)) object
         else if (is(object, "rlmerMod")) object@pp$X
         else if (is(object, "merMod")) lme4::getME(object, "X")
         else stop("'object' must be an rlmerMod or merMod fit or a matrix")
    n <- nrow(X); p <- ncol(X)
    kk <- .design_k(X, max_subsets = max_subsets)
    k <- kk$k
    list(n = n, p = p, k = k, k_share = k / n,
         max_bp = floor((n - k + 1) / 2) / n,
         mm_bp = max(0, min(0.5, 1 - 0.5 - k / n)),
         b_opt = (n - k) / (2 * n),
         exact = kk$exact)
}
