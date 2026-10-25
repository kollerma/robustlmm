## Fragility of the variance-component parameters under cluster
## contamination: the first-order contamination fraction that drives a
## bounded parameter to its bound (paper 2 thm:threshold, whose von Mises
## expansion is machine-checked in IF-thread1 VonMisesThreshold).
##
## Along F_eps = (1 - eps) F + eps G the parameter moves as
## theta_k(F_eps) = theta_k(F) + eps IF_k(G) + O(eps^2), so a contaminating
## cluster G with IF_k(G) < 0 drives theta_k to zero at
## eps* = theta_k / (-IF_k(G)) to first order. At the empirical
## distribution the case-weight path is Hampel's influence function
## (paper 7), and for observed cluster j
##     IF(j) = -J * Jpar^{-1} s_j,
## with s_j the cluster's score contribution (.scoreByCluster) and Jpar the
## parameter Jacobian of the full score (implicitIF_full). The worst
## observed cluster gives eps*_k = theta_k / max_j(-IF_k(j)). The
## theoretical threshold uses the supremum over all conceivable clusters,
## which can only be larger than the observed maximum, so the value
## reported here is an upper bound on it: a small value flags a fragile
## parameter, a large one does not certify robustness.

.fragility_table <- function(fit, IF = NULL) {
    tc <- .top_cluster(fit)
    if (is.null(tc))
        stop("fragility() needs a single or nested grouping structure; ",
             "this fit has crossed grouping factors.", call. = FALSE)
    if (is.null(IF)) IF <- implicitIF_full(fit)
    S  <- .scoreByCluster(fit)
    J  <- ncol(S)
    Xi <- .applyJparInv(IF, S)      # = -Jpar^{-1} s_j (per-cluster delta)
    p  <- fit@pp$p; q <- fit@pp$q
    th <- getME(fit, "theta")
    lower <- getME(fit, "lower")
    L  <- length(th)
    IFth <- J * Xi[p + q + 1L + seq_len(L), , drop = FALSE]   # Hampel IF
    bounded <- which(is.finite(lower) & lower == 0)
    if (!length(bounded)) return(NULL)
    ## which diagonal position each bounded theta occupies in its block
    ## (1 = first component of the block: its relative SD; > 1: the
    ## conditional SD, whose bound is a singular covariance)
    pos <- integer(L)
    for (B in fit@blocks) {
        B <- as.matrix(B)
        d <- diag(B)
        pos[d[d > 0]] <- seq_along(d)[d > 0]
    }
    levs <- levels(tc$g)
    rows <- lapply(bounded, function(k) {
        down <- -IFth[k, ]
        jmax <- which.max(down)
        gmax <- if (length(jmax)) down[jmax] else NA_real_
        eps <- if (th[k] <= 1e-8) 0
               else if (is.finite(gmax) && gmax > 0) th[k] / gmax
               else Inf
        data.frame(parameter = names(th)[k],
                   bound = if (pos[k] <= 1L) "zero variance"
                           else "singular covariance",
                   theta = unname(th[k]),
                   eps = eps,
                   ## m added copies of the driver are a fraction m / (J + m)
                   clusters = if (eps < 1) eps * J / (1 - eps) else Inf,
                   driver = if (is.finite(gmax) && gmax > 0) levs[jmax]
                            else NA_character_,
                   stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    attr(out, "J") <- J
    attr(out, "grouping") <- tc$name
    out
}

##' First-order contamination fraction that would drive each variance
##' parameter of a robust LMM to its bound.
##'
##' For every variance parameter \eqn{\theta_k} with a lower bound of zero
##' (the diagonal entries of the relative Cholesky factor), \code{fragility}
##' reports
##' \deqn{\epsilon^*_k = \theta_k / \max_j(-\mathrm{IF}_k(j)),}{eps*_k =
##'   theta_k / max_j(-IF_k(j)),}
##' where \eqn{\mathrm{IF}_k(j)} is the influence of cluster \eqn{j} on
##' \eqn{\theta_k} at the fitted model. Replacing a fraction \eqn{\epsilon}
##' of the clusters by copies of the most influential one moves
##' \eqn{\theta_k} by \eqn{\epsilon\,\mathrm{IF}_k + O(\epsilon^2)}, so it
##' reaches its bound at \eqn{\epsilon^*_k} to first order: a zero variance
##' for the first component of a random-effects block, a singular
##' covariance matrix (for a \eqn{2 \times 2} block a correlation of
##' \eqn{\pm 1}) for the later ones.
##'
##' Only the observed clusters are searched, so \eqn{\epsilon^*} is an upper
##' bound on the worst case: a small value flags a fragile parameter, a
##' large one does not certify robustness, and values above 1/2 are only
##' indicative. \code{clusters}, \eqn{\epsilon^* J / (1 - \epsilon^*)}, is
##' the number of copies of the driving cluster that, added to the \eqn{J}
##' observed clusters, make up the fraction \eqn{\epsilon^*}. Crossed
##' designs are not supported.
##'
##' @title Fragility of variance parameters under contamination
##' @param object an \code{rlmerMod} object.
##' @param ... ignored.
##' @return a data frame with one row per bounded variance parameter:
##'   \code{parameter}, \code{bound}, \code{theta} (estimate, relative to
##'   the residual scale), \code{eps} (\eqn{\epsilon^*}), \code{clusters}
##'   (\eqn{\epsilon^* J / (1 - \epsilon^*)}) and \code{driver} (the most
##'   influential cluster).
##' @examples
##' fit <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
##'              method = "DASvar")
##' fragility(fit)
##' @section Lifecycle:
##'   \strong{Experimental.} The interface and the results may change in a
##'   later version.
##' @export
fragility <- function(object, ...) UseMethod("fragility")

##' @rdname fragility
##' @export
fragility.rlmerMod <- function(object, ...) .fragility_table(object)
