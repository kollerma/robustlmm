## Influence-function substrate for rlmerMod: case-weight IF and the
## robust cluster-sandwich covariance of the fixed effects.
##
## Port of IF-thread1/implicitIF.R (the partial version: sigma-hat and
## theta-hat held fixed). The full (sigma, theta) IF used by
## cooks.distance.rlmerMod lives in influence_full.R.

##' Implicit (numerical) influence function building blocks for an rlmer fit.
##'
##' Computes the local-shift sensitivity \eqn{\partial \hat{\beta}/\partial y}
##' (and the analogous quantity for the spherical random effects \eqn{\hat{u}})
##' by implicit differentiation of the joint \eqn{(\beta, u)}-scoring
##' equations at the fit. The returned list is a building block for
##' \code{\link{caseweightIF}} and \code{\link{vcov_sandwich}}; users who
##' only want robust standard errors should call \code{vcov(object, type = "sandwich")}.
##'
##' Prior weights and offsets are supported: the score is evaluated in
##' the same prior-weight-whitened coordinates the estimator iterates
##' in (\eqn{U_e = \mathrm{diag}(\sqrt{w})}, residuals offset-corrected).
##'
##' This is the partial (\eqn{\sigma, \theta} held fixed) version. It is
##' \emph{not} the Hampel influence function; the Hampel object is
##' \code{\link{caseweightIF}}, which reuses the same Jacobian but with the
##' score \emph{value} \eqn{\psi} on the right-hand side rather than
##' \eqn{\partial_y \psi}. See Koller (2026, Paper 1) for a discussion.
##'
##' @param fit An \code{rlmerMod} object.
##' @param use.expected If \code{TRUE}, use \eqn{E[\psi']} for the diagonal
##'   blocks of the Jacobian; the default (\code{FALSE}) uses the actual
##'   \eqn{\psi'} values at the fitted residuals (the empirical version).
##' @return A list with components \code{IF_beta}, \code{IF_u}, \code{IF_b},
##'   the Jacobian blocks \code{Jbb}, \code{Jbu}, \code{Jub}, \code{S},
##'   and intermediate rotated-design matrices used by the cluster
##'   sandwich. \code{IF_beta} is the local-shift sensitivity matrix
##'   \eqn{p \times n}.
##' @keywords internal
implicitIF <- function(fit, use.expected = FALSE) {
    stopifnot(is(fit, "rlmerMod"))

    pp     <- fit@pp
    rho_e  <- fit@rho.e
    rho_b  <- fit@rho.b
    sig    <- .sigma(fit)
    lambda_e <- rho_e@EDpsi()
    lambda_b_q <- pp$D_b@x

    X     <- pp$X
    Z     <- t(pp$Zt)
    U_b   <- pp$U_b
    U_e   <- pp$U_e
    n     <- pp$n
    p     <- pp$p
    q     <- pp$q

    beta  <- pp$beta
    u     <- pp$b.s
    y     <- fit@resp$y
    off   <- fit@resp$offset

    ## The estimator works in prior-weight-whitened coordinates: with
    ## prior weights w (Var(e_i) = sigma^2 / w_i), U_e = diag(sqrt(w))
    ## and the psi functions are evaluated at U_e r / sigma (equals
    ## resp$wtres / sigma, see std.e() in helpers.R). For unit weights
    ## U_e = I.
    r_std <- as.numeric(U_e %*% (y - off - X %*% beta - Z %*% (U_b %*% u))) / sig
    u_std <- as.numeric(u) / sig

    Dpsi_e <- if (use.expected) rep(lambda_e, n) else rho_e@Dpsi(r_std)
    D_e_act <- Diagonal(x = Dpsi_e)
    ## WS11: the e-side score is eta-weighted (X~' H psi_e, Z~' H psi_e),
    ## so its curvature / y-derivative carry H D_e_act = diag(eta *
    ## Dpsi_e). The b-side penalty (Lambda_ratio D_b_act) is unweighted.
    ## eta == 1 gives D_e_eta == D_e_act (bit-identical).
    eta     <- pp$eta
    D_e_eta <- Diagonal(x = eta * Dpsi_e)

    if (use.expected) {
        D_b_act <- Diagonal(x = lambda_b_q)
    } else {
        D_b_act <- .buildBlockDpsi_b(fit, u_std, rho_b)
    }

    ## Whitened design matrices, identical to the ones the estimator
    ## iterates with (pp$U_eX, pp$U_eZU_b).
    A_X     <- as.matrix(U_e %*% X)
    A_ZUb   <- as.matrix(U_e %*% (Z %*% U_b))

    Jbb <- crossprod(A_X, D_e_eta %*% A_X)
    Jbu <- crossprod(A_X, D_e_eta %*% A_ZUb)
    Lambda_ratio <- Diagonal(x = lambda_e / lambda_b_q)
    Juu <- crossprod(A_ZUb, D_e_eta %*% A_ZUb) + Lambda_ratio %*% D_b_act
    Juu <- (Juu + t(Juu)) / 2
    Jub <- t(Jbu)

    ## d r_std / d y = U_e / sigma, so the y-rows of the Jacobian carry
    ## another factor U_e (the sigma factors cancel against the ones
    ## implicit in the parameter rows). The e-side weight H rides along.
    Gby <- crossprod(A_X, D_e_eta %*% U_e)
    Guy <- crossprod(A_ZUb, D_e_eta %*% U_e)

    Jbb_inv_Gby <- .solveSafe(Jbb, Gby)
    Jbb_inv_Jbu <- .solveSafe(Jbb, Jbu)
    S           <- Juu - Jub %*% Jbb_inv_Jbu
    IF_u        <- as.matrix(.solveSafe(S, Guy - Jub %*% Jbb_inv_Gby))
    IF_beta     <- as.matrix(Jbb_inv_Gby - Jbb_inv_Jbu %*% IF_u)
    IF_b        <- as.matrix(U_b %*% IF_u)

    list(IF_beta = IF_beta, IF_u = IF_u, IF_b = IF_b,
         D_e = Dpsi_e, D_b_act = D_b_act,
         lambda_e = lambda_e, lambda_b_q = lambda_b_q, sigma = sig,
         Jbb = Jbb, Jbu = Jbu, Jub = Jub, S = S,
         Jbb_inv_Jbu = Jbb_inv_Jbu,
         A_X = A_X, A_ZUb = A_ZUb, U_b = U_b, U_e = U_e,
         Lambda_ratio = Lambda_ratio,
         r_std = r_std, u_std = u_std)
}

## solve() with a Moore-Penrose fallback at the variance-component
## boundary, where Jbb / Schur can become rank-deficient.
.solveSafe <- function(A, b) {
    tryCatch(solve(A, b), error = function(e) {
        warning("Jacobian rank-deficient (variance-component boundary?); ",
                "using pseudo-inverse fallback.", call. = FALSE)
        sv <- svd(as.matrix(A))
        tol <- max(dim(A)) * .Machine$double.eps * sv$d[1L]
        Ainv <- tcrossprod(sv$v %*% diag(ifelse(sv$d > tol, 1 / sv$d, 0)), sv$u)
        Ainv %*% as.matrix(b)
    })
}

##' Case-weight (Hampel) influence function of the fixed effects of a
##' fitted \code{rlmerMod} object.
##'
##' Returns \eqn{IF(y_i) = -\hat{J}^{-1} \psi(\hat{par}; y_i)} (the Hampel
##' object built from the score \emph{value} \eqn{\psi}), reusing the
##' Jacobian \eqn{\hat{J}} formed by \code{\link{implicitIF}}. This is the
##' right object for robustness quantities (gross-error sensitivity,
##' breakdown), distinct from the local-shift sensitivity returned by
##' \code{\link{implicitIF}} for a redescending \eqn{\psi}.
##'
##' \eqn{\hat{\sigma}, \hat{\theta}} are held fixed (partial IF).
##'
##' @param fit An \code{rlmerMod} object.
##' @param idx Optional integer indices selecting observations to compute
##'   (default: all observations).
##' @param use.expected Passed to \code{\link{implicitIF}}; default
##'   \code{FALSE} (empirical Jacobian).
##' @return A list with \code{IF_beta} (\eqn{p \times |idx|}), \code{IF_u}
##'   (\eqn{q \times |idx|}), \code{IF_b}, and the per-observation
##'   score-contribution columns \code{Gb_cw}, \code{Gu_cw}.
##' @seealso \code{\link{vcov_sandwich}}, \code{\link[=vcov.rlmerMod]{vcov}}
##' @export
caseweightIF <- function(fit, idx = NULL, use.expected = FALSE) {
    stopifnot(is(fit, "rlmerMod"))
    parts <- implicitIF(fit, use.expected = use.expected)

    pp    <- fit@pp
    rho_e <- fit@rho.e
    n     <- pp$n
    if (is.null(idx)) idx <- seq_len(n)

    sig <- parts$sigma
    eta <- pp$eta
    if (length(eta) != n) eta <- rep.int(1, n)
    ## Mallows design weights multiply each observation's score
    ## contribution, as in the estimating equations.
    psi_e_vals <- eta * rho_e@psi(parts$r_std)

    W     <- Diagonal(x = psi_e_vals[idx])
    Gb_cw <- as.matrix(crossprod(parts$A_X[idx, , drop = FALSE], W))
    Gu_cw <- as.matrix(crossprod(parts$A_ZUb[idx, , drop = FALSE], W))

    Jbb_inv_Gb <- .solveSafe(parts$Jbb, Gb_cw)
    IF_u    <- as.matrix(.solveSafe(parts$S, Gu_cw - parts$Jub %*% Jbb_inv_Gb))
    IF_beta <- as.matrix(Jbb_inv_Gb - parts$Jbb_inv_Jbu %*% IF_u)
    IF_beta <- sig * IF_beta
    IF_u    <- sig * IF_u

    colnames(IF_beta) <- paste0("y", idx)
    list(IF_beta = IF_beta, IF_u = IF_u,
         IF_b = as.matrix(parts$U_b %*% IF_u),
         Gb_cw = Gb_cw, Gu_cw = Gu_cw, sigma = sig)
}

##' Resolve a cluster specification for a fitted \code{rlmerMod} object.
##'
##' Used by \code{\link{vcov_sandwich}}. When every grouping factor is
##' nested within the clusters (a single grouping factor, or the top level
##' of a nested design) the cluster-robust sandwich is exact; for crossed
##' factors no single clustering nests every random-effect block, so a
##' warning is issued and the sandwich is approximate.
##'
##' @param fit \code{rlmerMod} object.
##' @param cluster \code{NULL} (auto-detect: use the sole grouping factor,
##'   or error if the design is crossed), a character string naming a
##'   grouping factor of the model, or a length-\code{n} vector of cluster
##'   memberships.
##' @param n Number of observations.
##' @return A length-\code{n} factor of cluster memberships.
##' @keywords internal
resolveCluster <- function(fit, cluster, n) {
    fl <- getME(fit, "flist")
    if (is.null(cluster)) {
        if (length(fl) == 1L) return(fl[[1L]])
        stop("Design has multiple grouping factors (",
             paste(names(fl), collapse = ", "),
             "); specify `cluster` explicitly. For crossed factors the ",
             "cluster-robust sandwich is only approximate.")
    }
    if (is.character(cluster) && length(cluster) == 1L) {
        if (!cluster %in% names(fl))
            stop("cluster '", cluster, "' is not a grouping factor; ",
                 "available: ", paste(names(fl), collapse = ", "))
        g <- fl[[cluster]]
    } else {
        if (length(cluster) != n)
            stop("cluster vector must have length n = ", n)
        g <- factor(cluster)
    }
    ## a grouping factor nested within the clusters (the lower levels of a
    ## nested design) keeps the sandwich exact; one that is not (crossed
    ## factors) leaves per-cluster u-scores that do not vanish
    crossed <- names(fl)[!vapply(fl, .nested_within, logical(1), g = g)]
    if (length(crossed))
        warning("Grouping factor(s) ", paste(crossed, collapse = ", "),
                " not nested within the clusters: their per-cluster ",
                "u-scores do not vanish, so the sandwich is approximate.",
                call. = FALSE)
    g
}

## Is every level of factor f contained in a single level of g?
.nested_within <- function(f, g) {
    f <- as.integer(factor(f)); g <- as.integer(factor(g))
    ok <- !is.na(f) & !is.na(g)
    first <- g[ok][match(f[ok], f[ok])]
    all(first == g[ok])
}

##' Robust cluster-sandwich covariance of the fixed effects of a fitted
##' \code{rlmerMod} object.
##'
##' Computes the robust score sandwich \eqn{\hat{V}_{IF} = \hat{A}^{-1}
##' \hat{B} \hat{A}^{-T}}, where \eqn{\hat{A}} is the Schur-complement
##' (marginal) Jacobian of the profiled \eqn{\beta}-score and \eqn{\hat{B}
##' = \sum_j s_j s_j^T} sums the per-cluster \eqn{\beta}-score
##' contributions \eqn{s_j = \sum_{i \in j} \eta_i x_i \psi_e(\hat{r}_i)}
##' (\eqn{\eta_i} the Mallows design weights, 1 unless
##' \code{design.weights} is used). Equal to the user-facing
##' \code{vcov(object, type = "sandwich")}.
##'
##' Exact when every grouping factor is nested within the clusters (a
##' single grouping factor, or \code{cluster} naming the top level of a
##' nested design); approximate for crossed factors, with a warning.
##'
##' \strong{Small-sample corrections.} \code{"G1"} (default) scales by
##' \eqn{J/(J-1)}. This removes the small-sample bias for the intercept but
##' not for between-cluster coefficients, for which the sandwich stays
##' biased down when there are few clusters. \code{"CR2"} is the
##' bias-reduced sandwich of Bell and McCaffrey (2002) adapted to the robust
##' estimating equations: with \eqn{\hat{A} = \sum_j \hat{A}_j} split into
##' cluster contributions, each \eqn{s_j} is replaced by
##' \eqn{(I - \hat{A}_j\hat{A}^{-1})^{-1/2} s_j}. For a Gaussian fit it is
##' unbiased under the fitted model. It takes the symmetric square root
##' after whitening by the clusters' marginal covariances, whereas the CR2
##' of Pustejovsky and Tipton (2018) uses their Cholesky factors; the two
##' coincide when these covariances are multiples of the identity and
##' otherwise differ slightly. \code{"CR3"} uses
##' \eqn{(I - \hat{A}_j\hat{A}^{-1})^{-1}}: \eqn{J/(J-1)} times the one-step
##' leave-one-cluster-out jackknife. CR2 and CR3 need every random effect
##' nested within one
##' cluster. \code{\link[=anova.rlmerMod]{anova}} and
##' \code{\link[=confint.rlmerMod]{confint}} use CR2 with Bell--McCaffrey /
##' approximate Hotelling degrees of freedom (Pustejovsky and Tipton 2018)
##' for \code{vcov_type = "sandwich"}.
##'
##' \eqn{\hat{\sigma}} and \eqn{\hat{\theta}} are held fixed (partial
##' sandwich).
##'
##' @param fit \code{rlmerMod} object.
##' @param cluster Cluster specification; see \code{\link{resolveCluster}}.
##' @param correction One of \code{"G1"} (default, applies \eqn{J/(J-1)}),
##'   \code{"none"}, \code{"CR2"} or \code{"CR3"}; see Details.
##' @return A \eqn{p \times p} covariance matrix for \eqn{\hat{\beta}},
##'   with dimnames from the fixed-effect coefficient names and attribute
##'   \code{"n.clusters"}.
##' @references Bell RM, McCaffrey DF (2002). Bias reduction in standard
##'   errors for linear regression with multi-stage samples. \emph{Survey
##'   Methodology} 28(2), 169--181.
##' @references Pustejovsky JE, Tipton E (2018). Small-sample methods for
##'   cluster-robust variance estimation and hypothesis testing in fixed
##'   effects models. \emph{Journal of Business & Economic Statistics}
##'   36(4), 672--683.
##' @seealso \code{\link[=vcov.rlmerMod]{vcov}}, \code{\link{caseweightIF}}
##' @section Lifecycle:
##'   Stable.
##' @export
vcov_sandwich <- function(fit, cluster = NULL,
                          correction = c("G1", "none", "CR2", "CR3")) {
    stopifnot(is(fit, "rlmerMod"))
    correction <- match.arg(correction)
    V <- .vcov_sandwich(fit, cluster = cluster, correction = correction)
    J <- attr(V, "n.clusters")
    if (J < 20L && correction %in% c("G1", "none"))
        warning("Cluster sandwich with J = ", J, " < 20 clusters: with ",
                "correction = \"", correction, "\" it is biased down for ",
                "between-cluster coefficients, and chi-square / normal ",
                "references are anti-conservative. anova() and confint() ",
                "use correction = \"CR2\" with Bell-McCaffrey degrees of ",
                "freedom for vcov_type = \"sandwich\".", call. = FALSE)
    V
}

## Workhorse of vcov_sandwich() without the small-J warning (anova() and
## confint() call it with the correction they then calibrate).
.vcov_sandwich <- function(fit, cluster = NULL, correction = "G1") {
    parts <- implicitIF(fit)
    rho_e <- fit@rho.e
    n     <- fit@pp$n
    p     <- fit@pp$p
    sig   <- parts$sigma
    eta   <- fit@pp$eta
    if (length(eta) != n) eta <- rep.int(1, n)

    ## The e-side estimating equation is sum_i eta_i x_i psi_e(r_i): the
    ## Mallows design weights multiply each score contribution.
    psi_e <- rho_e@psi(parts$r_std)
    Gb    <- t(parts$A_X * (eta * psi_e))

    g       <- resolveCluster(fit, cluster, n)
    groups  <- split(seq_len(n), g)
    J       <- length(groups)

    Smat <- vapply(groups,
                   function(ix) rowSums(Gb[, ix, drop = FALSE]),
                   numeric(p))
    if (is.null(dim(Smat))) Smat <- matrix(Smat, nrow = p)

    cr <- NULL
    if (correction %in% c("CR2", "CR3")) {
        Aj <- .cluster_jacobians(fit, parts, g, eta)
        ## Whitening root B with B B' = A^{-1}: the Cholesky factor, not a
        ## truncated eigen-root. CR2 does not depend on the root chosen,
        ## A^{-1} (I - A_j A^{-1})^{-1/2} = B (I - M_j)^{-1/2} B' with
        ## M_j = B' A_j B, but a relative eigenvalue cut-off on A silently
        ## dropped directions of badly scaled designs (a raw calendar year:
        ## cond(A) ~ 1e10).
        B  <- .whitening_root(Reduce(`+`, Aj))
        Mj <- lapply(Aj, function(a) {
            M <- crossprod(B, a %*% B)
            (M + t(M)) / 2
        })
        ## A direction of M_j with eigenvalue 1 is informed by cluster j
        ## alone (e.g. a factor level present in one cluster): (I - M_j)^{-1/2}
        ## is infinite there and the sandwich cannot estimate its variance.
        ## Such directions are recorded (and contrasts loading on them get
        ## no reference distribution, see .cr2_df) instead of silently
        ## dropped by the pseudo-inverse.
        degen <- .cr_degenerate_directions(Mj)
        if (ncol(degen))
            warning("The cluster sandwich cannot estimate the variance in ",
                    ncol(degen), " direction(s) of the fixed effects that a ",
                    "single cluster informs on its own (for example a factor ",
                    "level present in only one cluster); contrasts involving ",
                    "them get no small-sample reference (NA). Use vcov_type = ",
                    "\"default\" for such terms.", call. = FALSE)
        Rj <- lapply(Mj, function(M)
            .msym_pow(diag(p) - M, if (correction == "CR2") -1/2 else -1))
        IF <- vapply(seq_len(J), function(j)
            sig * as.numeric(B %*% (Rj[[j]] %*% crossprod(B, Smat[, j]))),
            numeric(p))
        if (is.null(dim(IF))) IF <- matrix(IF, nrow = p)
        V <- tcrossprod(IF)
        if (correction == "CR2") {
            ## kept in an environment so that printing V stays short
            cr <- new.env(parent = emptyenv())
            cr$B <- B; cr$Mj <- Mj; cr$Rj <- Rj; cr$degen <- degen
        }
    } else {
        Jbb_inv_S <- .solveSafe(parts$Jbb, Smat)
        xu        <- .solveSafe(parts$S, -parts$Jub %*% Jbb_inv_S)
        xb        <- Jbb_inv_S - parts$Jbb_inv_Jbu %*% xu
        IF        <- sig * as.matrix(xb)
        V <- tcrossprod(IF)
        if (correction == "G1" && J > 1L) V <- V * (J / (J - 1))
    }
    V <- (V + t(V)) / 2
    nm <- colnames(fit@pp$X)
    dimnames(V) <- list(nm, nm)
    attr(V, "n.clusters") <- J
    attr(V, "correction") <- correction
    if (!is.null(cr)) attr(V, "cr2") <- cr
    V
}

## Symmetric matrix power via the eigendecomposition; eigenvalues below
## tol (relative) are treated as zero (Moore-Penrose convention, as for
## the CR2 adjustment matrices of Pustejovsky and Tipton 2018). Used only
## for I - M_j, whose eigenvalues lie in [0, 1].
.msym_pow <- function(M, pw, tol = 1e-10) {
    e <- eigen((M + t(M)) / 2, symmetric = TRUE)
    lam <- e$values
    keep <- lam > tol * max(1, abs(lam[1L]))
    v <- numeric(length(lam))
    v[keep] <- lam[keep]^pw
    e$vectors %*% (v * t(e$vectors))
}

## B with B B' = M^{-1} for a symmetric positive definite M: B = R^{-1} for
## the Cholesky factor M = R'R. Scale-equivariant, so no direction of a
## badly scaled design is lost. A singular M (fixed effects without
## information) falls back to a pseudo-inverse root with a warning.
.whitening_root <- function(M) {
    M <- (M + t(M)) / 2
    R <- tryCatch(chol(M), error = function(e) NULL)
    if (!is.null(R)) return(backsolve(R, diag(nrow(M))))
    warning("CR2: the fixed-effects Jacobian is singular; using a ",
            "pseudo-inverse.", call. = FALSE)
    .msym_pow(M, -1/2)
}

## Directions (columns, in the whitened coordinates of M_j) in which some
## M_j has eigenvalue 1 up to tol: informed by that cluster alone.
.cr_degenerate_directions <- function(Mj, tol = 1e-8) {
    p <- nrow(Mj[[1L]])
    vs <- lapply(Mj, function(M) {
        e <- eigen(M, symmetric = TRUE)
        e$vectors[, e$values >= 1 - tol, drop = FALSE]
    })
    V <- do.call(cbind, c(list(matrix(0, p, 0)), vs))
    if (!ncol(V)) return(V)
    q <- qr(V)
    qr.Q(q)[, seq_len(q$rank), drop = FALSE]
}

## Cluster contributions A_j of the marginal (u-profiled) beta-Jacobian,
## sum_j A_j = A = Jbb - Jbu Juu^{-1} Jub. Requires every random effect to
## load on a single cluster, so that Juu is block-diagonal by cluster.
.cluster_jacobians <- function(fit, parts, g, eta) {
    gi <- as.integer(g)
    J  <- nlevels(g)
    De <- eta * parts$D_e
    AX <- parts$A_X
    AZ <- parts$A_ZUb
    q  <- ncol(AZ)
    Pen <- as.matrix(parts$Lambda_ratio %*% parts$D_b_act)
    Pen <- (Pen + t(Pen)) / 2
    ## random effect -> cluster
    Zt_t <- as(fit@pp$Zt, "TsparseMatrix")
    uc <- unique(cbind(Zt_t@i + 1L, gi[Zt_t@j + 1L]))
    if (anyDuplicated(uc[, 1L]))
        stop(structure(class = c("rlmer_cr_not_nested", "error", "condition"),
            list(message = paste0(
                "correction = \"CR2\" / \"CR3\" needs every random effect ",
                "nested within one cluster; here some load on several ",
                "clusters (crossed factors, or a cluster variable finer than ",
                "a grouping factor). Use correction = \"G1\"."),
                call = NULL)))
    u_cl <- integer(q)
    u_cl[uc[, 1L]] <- uc[, 2L]
    ## random effects without data (all-zero covariate in their level)
    ## join the cluster of the block they are coupled to
    for (it in seq_len(3L)) {
        zero <- which(u_cl == 0L)
        if (!length(zero)) break
        for (k in zero) {
            nb <- which(Pen[k, ] != 0 & u_cl > 0L)
            if (length(nb)) u_cl[k] <- u_cl[nb[1L]]
        }
    }
    lapply(seq_len(J), function(j) {
        ix <- which(gi == j)
        uj <- which(u_cl == j)
        Xj <- AX[ix, , drop = FALSE]
        dj <- De[ix]
        Aj <- crossprod(Xj, dj * Xj)
        if (length(uj)) {
            Zj  <- AZ[ix, uj, drop = FALSE]
            Bxu <- crossprod(Xj, dj * Zj)
            Buu <- crossprod(Zj, dj * Zj) + Pen[uj, uj, drop = FALSE]
            Aj  <- Aj - Bxu %*% .solveSafe(Buu, t(Bxu))
        }
        as.matrix((Aj + t(Aj)) / 2)
    })
}

## Degrees of freedom for the CR2 Wald test of L beta = 0 (L: k x p): the
## approximate Hotelling T^2 of Pustejovsky and Tipton (2018), eta =
## k (k + 1) / sum_{s,t} Var(Z_st), Z the CR2 estimate of L V L' whitened
## by its working-model value, under the working model Cov(s_j) ~ A_j with
## independent clusters. For k = 1 this is the Bell and McCaffrey (2002)
## Satterthwaite df. The test is F = (eta - k + 1) / (eta k) * W on
## (k, eta - k + 1) df, W the Wald statistic.
.cr2_df <- function(V, L) {
    cr <- attr(V, "cr2")
    if (is.null(cr)) return(NA_real_)
    L  <- rbind(L)
    k  <- nrow(L)
    B  <- cr$B
    p  <- ncol(B)
    J  <- length(cr$Mj)
    ## whitened contrasts: rows of Lt satisfy Lt A^{-1} Lt' = I
    Lw <- L %*% B                                       # k x p
    Lt <- .whitening_root(tcrossprod(Lw))
    ct <- t(crossprod(Lt, Lw))                          # p x k, B' Lt'
    ## a contrast loading on a direction informed by one cluster alone has
    ## no estimable sandwich variance (see .cr_degenerate_directions)
    if (ncol(cr$degen) &&
        max(abs(crossprod(cr$degen, ct))) > 1e-6)
        return(NA_real_)
    ## Per cluster j and whitened contrast s: a_js = R_j c_s, u_js = M_j a_js,
    ## d_j^{st} = a_js' u_jt. With U_s = [u_1s .. u_Js] (p x J) the working-
    ## model covariances of the cluster terms are G^{st} = diag(d^{st}) -
    ## U_s' U_t, and Var(Z_st) = <G^{ss}, G^{tt}> + tr(G^{st} G^{st}); both
    ## are expanded below so that only p x p products are formed (O(J p^2)).
    U <- array(0, c(p, J, k))
    d <- array(0, c(J, k, k))
    for (j in seq_len(J)) {
        a  <- cr$Rj[[j]] %*% ct
        uu <- cr$Mj[[j]] %*% a
        U[, j, ] <- uu
        d[j, , ] <- crossprod(a, uu)
    }
    Us  <- lapply(seq_len(k), function(s) matrix(U[, , s], p))
    nrm <- vapply(Us, function(u) colSums(u^2), numeric(J))
    nrm <- matrix(nrm, J)
    tot <- 0
    for (s in seq_len(k)) for (t in seq_len(k)) {
        dss <- d[, s, s]; dtt <- d[, t, t]; dst <- d[, s, t]
        UsUt <- tcrossprod(Us[[s]], Us[[t]])            # p x p
        UtUs <- t(UsUt)
        inner_st <- colSums(Us[[s]] * Us[[t]])          # u_js . u_jt
        g1 <- sum(dss * dtt) - sum(dss * nrm[, t]) - sum(dtt * nrm[, s]) +
            sum(UsUt^2)
        g2 <- sum(dst^2) - 2 * sum(dst * inner_st) + sum(UtUs * t(UtUs))
        tot <- tot + g1 + g2
    }
    k * (k + 1) / tot
}

## Model-based delta-method covariance of beta-hat (internal comparator).
## NOT a robust sandwich: sandwiches the local-shift sensitivity between
## the assumed Gaussian Cov(y). Equals vcov_sandwich() only under
## correctly-specified Gaussian y with the true V_b.
.vcov_model_delta <- function(fit, parts = NULL) {
    if (is.null(parts)) parts <- implicitIF(fit)
    pp  <- fit@pp
    Z_full <- crossprod(pp$Zt, pp$U_b)
    n      <- pp$n
    ## Var(y) / sigma^2 = Z U_b U_b' Z' + W^{-1}, W = diag(prior weights)
    cov_y_unit <- as.matrix(tcrossprod(Z_full) +
                            Diagonal(x = 1 / diag(pp$U_e)^2))
    V <- parts$sigma^2 * tcrossprod(parts$IF_beta %*% cov_y_unit, parts$IF_beta)
    V <- (V + t(V)) / 2
    nm <- colnames(pp$X)
    dimnames(V) <- list(nm, nm)
    V
}

## Block-diagonal Jacobian of psi_b wrt u, evaluated at u_std. For block
## size 1 it collapses to a diagonal; for size > 1 each block k of size
## m is w(d_k) I_m + 2 w'(d_k) u_k u_k', d_k = ||u_k||^2.
.buildBlockDpsi_b <- function(fit, u_std, rho_b) {
    q <- length(u_std)
    is_diag_only <- all(vapply(fit@idx, nrow, integer(1)) == 1L)
    if (is_diag_only) {
        Dpsi_b <- numeric(q)
        for (bt in seq_along(rho_b)) {
            idx <- which(fit@ind == bt)
            Dpsi_b[idx] <- rho_b[[bt]]@Dpsi(u_std[idx])
        }
        return(Diagonal(x = Dpsi_b))
    }
    blocks_list <- vector("list", 0)
    for (bt in seq_along(rho_b)) {
        idx_mat <- fit@idx[[bt]]
        m       <- nrow(idx_mat)
        wfun    <- rho_b[[bt]]@wgt
        Dwfun   <- rho_b[[bt]]@Dwgt
        for (k in seq_len(ncol(idx_mat))) {
            j_set <- as.vector(idx_mat[, k])
            u_k   <- u_std[j_set]
            d_k   <- sum(u_k * u_k)
            blockMat <- wfun(d_k) * diag(m) + 2 * Dwfun(d_k) * tcrossprod(u_k)
            blocks_list[[length(blocks_list) + 1L]] <- list(idx = j_set,
                                                            mat = blockMat)
        }
    }
    i_idx <- unlist(lapply(blocks_list,
                           function(b) rep(b$idx, length(b$idx))))
    j_idx <- unlist(lapply(blocks_list,
                           function(b) rep(b$idx, each = length(b$idx))))
    x_val <- unlist(lapply(blocks_list, function(b) as.vector(b$mat)))
    sparseMatrix(i = i_idx, j = j_idx, x = x_val, dims = c(q, q),
                 symmetric = FALSE)
}

## Reconstruct the random-effects design matrix Z for newdata, mapped
## to the fit's column ordering. Newdata levels not seen in the fit get
## a zero row (consistent with predict.merMod's allow.new.levels=TRUE
## behaviour). Returns an n_new x q dense matrix in the fit's column
## space.
.reconstruct_Z_new <- function(fit, newdata) {
    pp_fit   <- fit@pp
    cnms     <- fit@cnms
    fit_flist <- fit@flist
    q_fit    <- pp_fit$q
    n_new    <- nrow(newdata)
    Z_new    <- matrix(0, n_new, q_fit)
    Gp <- fit@Gp
    for (gname in names(cnms)) {
        terms_g <- cnms[[gname]]
        rhs <- paste(terms_g, collapse = " + ")
        if ("(Intercept)" %in% terms_g) {
            rhs <- if (length(terms_g) == 1L) "1"
                   else paste(setdiff(terms_g, "(Intercept)"),
                              collapse = " + ")
        }
        fml <- as.formula(paste("~", rhs))
        mm  <- model.matrix(fml, data = newdata)
        new_g       <- as.factor(newdata[[gname]])
        fit_levels  <- levels(fit_flist[[gname]])
        col_idx     <- match(as.character(new_g), fit_levels)
        n_terms     <- ncol(mm)
        gp_start    <- Gp[which(names(cnms) == gname)] + 1L
        for (i in seq_len(n_new)) {
            li <- col_idx[i]
            if (is.na(li)) next
            cols <- gp_start + (li - 1L) * n_terms + seq_len(n_terms) - 1L
            Z_new[i, cols] <- mm[i, ]
        }
    }
    Z_new
}

## SE of new predictions. Fixed-effect part from vcov_sandwich (robust);
## random-effects part from the partial-IF model-delta variance of b-hat
## (sigma^2 * IF_b %*% Cov(y) %*% IF_b'). No cross term -- a deliberate
## conservative choice; the spec ("RE term via IF_u" in the prediction
## section) does not specify one and the resulting SE is monotone in the
## two PSD pieces.
##
## prediction = TRUE adds sigma-hat^2 for the new observation's residual
## variance.
##
## Returns the length-n_new SE vector.
.predictSE <- function(fit, X_new, Z_new = NULL,
                       cluster = NULL, correction = "G1",
                       prediction = FALSE) {
    V_beta <- vcov_sandwich(fit, cluster = cluster, correction = correction)
    X_new  <- as.matrix(X_new)
    se2 <- rowSums((X_new %*% V_beta) * X_new)
    if (!is.null(Z_new)) {
        parts   <- implicitIF(fit)
        sig2    <- parts$sigma^2
        U_b     <- parts$U_b
        Z_full  <- crossprod(fit@pp$Zt, U_b)
        cov_y_u <- as.matrix(tcrossprod(Z_full) +
                             Diagonal(x = 1 / diag(fit@pp$U_e)^2))
        IF_b    <- as.matrix(U_b %*% parts$IF_u)
        V_b     <- sig2 * tcrossprod(IF_b %*% cov_y_u, IF_b)
        Z_new   <- as.matrix(Z_new)
        se2 <- se2 + rowSums((Z_new %*% V_b) * Z_new)
    }
    if (prediction) se2 <- se2 + .sigma(fit)^2
    sqrt(pmax(se2, 0))
}

## Block-aware psi_b applied at u_std. For block k of size m,
## psi_b_k(u_k_std) = w_k(||u_k_std||^2) * u_k_std.
.applyBlockPsi_b <- function(fit, u_std, rho_b) {
    q  <- length(u_std)
    out <- numeric(q)
    is_diag_only <- all(vapply(fit@idx, nrow, integer(1)) == 1L)
    if (is_diag_only) {
        for (bt in seq_along(rho_b)) {
            idx <- which(fit@ind == bt)
            out[idx] <- rho_b[[bt]]@psi(u_std[idx])
        }
        return(out)
    }
    for (bt in seq_along(rho_b)) {
        idx_mat <- fit@idx[[bt]]
        m       <- nrow(idx_mat)
        wfun    <- rho_b[[bt]]@wgt
        for (k in seq_len(ncol(idx_mat))) {
            j_set <- as.vector(idx_mat[, k])
            u_k   <- u_std[j_set]
            out[j_set] <- wfun(sum(u_k * u_k)) * u_k
        }
    }
    out
}
