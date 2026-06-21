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

    r_std <- as.numeric(solve(U_e, y - X %*% beta - Z %*% (U_b %*% u))) / sig
    u_std <- as.numeric(u) / sig

    Dpsi_e <- if (use.expected) rep(lambda_e, n) else rho_e@Dpsi(r_std)
    D_e_act <- Diagonal(x = Dpsi_e)

    if (use.expected) {
        D_b_act <- Diagonal(x = lambda_b_q)
    } else {
        D_b_act <- .buildBlockDpsi_b(fit, u_std, rho_b)
    }

    Ue_inv  <- Diagonal(x = 1 / diag(U_e))
    A_X     <- as.matrix(Ue_inv %*% X)
    A_ZUb   <- as.matrix(Ue_inv %*% (Z %*% U_b))

    Jbb <- t(A_X) %*% D_e_act %*% A_X
    Jbu <- t(A_X) %*% D_e_act %*% A_ZUb
    Lambda_ratio <- Diagonal(x = lambda_e / lambda_b_q)
    Juu <- t(A_ZUb) %*% D_e_act %*% A_ZUb + Lambda_ratio %*% D_b_act
    Juu <- (Juu + t(Juu)) / 2
    Jub <- t(Jbu)

    Gby <- t(A_X) %*% D_e_act %*% Ue_inv
    Guy <- t(A_ZUb) %*% D_e_act %*% Ue_inv

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
         A_X = A_X, A_ZUb = A_ZUb, U_b = U_b, Ue_inv = Ue_inv,
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
        Ainv <- sv$v %*% diag(ifelse(sv$d > tol, 1 / sv$d, 0)) %*% t(sv$u)
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
    psi_e_vals <- rho_e@psi(parts$r_std)

    W     <- Diagonal(x = psi_e_vals[idx])
    Gb_cw <- as.matrix(t(parts$A_X[idx, , drop = FALSE]) %*% W)
    Gu_cw <- as.matrix(t(parts$A_ZUb[idx, , drop = FALSE]) %*% W)

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
##' Used by \code{\link{vcov_sandwich}}. For a single grouping factor
##' (nested design) the cluster-robust sandwich is exact; for crossed
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
        if (length(fl) > 1L)
            warning("Crossed/multiple grouping factors: the per-cluster ",
                    "u-scores of random effects NOT nested within '",
                    cluster, "' do not vanish, so the sandwich is ",
                    "approximate.", call. = FALSE)
        return(fl[[cluster]])
    }
    if (length(cluster) != n)
        stop("cluster vector must have length n = ", n)
    factor(cluster)
}

##' Robust cluster-sandwich covariance of the fixed effects of a fitted
##' \code{rlmerMod} object.
##'
##' Computes the robust score sandwich \eqn{\hat{V}_{IF} = \hat{A}^{-1}
##' \hat{B} \hat{A}^{-T}}, where \eqn{\hat{A}} is the Schur-complement
##' (marginal) Jacobian of the profiled \eqn{\beta}-score and \eqn{\hat{B}
##' = \sum_j s_j s_j^T} sums the per-cluster \eqn{\beta}-score
##' contributions \eqn{s_j = \sum_{i \in j} x_i \psi_e(\hat{r}_i)}. Equal
##' to the user-facing \code{vcov(object, type = "sandwich")}.
##'
##' Exact for a single (nested) grouping factor; approximate for crossed
##' factors (a warning is issued via \code{\link{resolveCluster}}). With
##' few clusters, set \code{correction = "G1"} (default) for the
##' \eqn{J/(J-1)} small-sample scaling.
##'
##' \eqn{\hat{\sigma}, \hat{\theta}} are held fixed (partial sandwich); the
##' returned variance is the leading-order fixed-effects covariance.
##'
##' @param fit \code{rlmerMod} object.
##' @param cluster Cluster specification; see \code{\link{resolveCluster}}.
##' @param correction One of \code{"G1"} (default, applies \eqn{J/(J-1)})
##'   or \code{"none"}.
##' @return A \eqn{p \times p} covariance matrix for \eqn{\hat{\beta}},
##'   with dimnames from the fixed-effect coefficient names and attribute
##'   \code{"n.clusters"}.
##' @seealso \code{\link[=vcov.rlmerMod]{vcov}}, \code{\link{caseweightIF}}
##' @export
vcov_sandwich <- function(fit, cluster = NULL, correction = c("G1", "none")) {
    stopifnot(is(fit, "rlmerMod"))
    correction <- match.arg(correction)
    parts <- implicitIF(fit)
    rho_e <- fit@rho.e
    n     <- fit@pp$n
    p     <- fit@pp$p
    sig   <- parts$sigma

    psi_e <- rho_e@psi(parts$r_std)
    Gb    <- t(parts$A_X * psi_e)

    g       <- resolveCluster(fit, cluster, n)
    groups  <- split(seq_len(n), g)
    J       <- length(groups)

    Smat <- vapply(groups,
                   function(ix) rowSums(Gb[, ix, drop = FALSE]),
                   numeric(p))
    if (is.null(dim(Smat))) Smat <- matrix(Smat, nrow = p)

    Jbb_inv_S <- .solveSafe(parts$Jbb, Smat)
    xu        <- .solveSafe(parts$S, -parts$Jub %*% Jbb_inv_S)
    xb        <- Jbb_inv_S - parts$Jbb_inv_Jbu %*% xu
    IF        <- sig * as.matrix(xb)

    V <- IF %*% t(IF)
    if (correction == "G1" && J > 1L) V <- V * (J / (J - 1))
    V <- (V + t(V)) / 2
    nm <- colnames(fit@pp$X)
    dimnames(V) <- list(nm, nm)
    attr(V, "n.clusters") <- J
    V
}

## Model-based delta-method covariance of beta-hat (internal comparator).
## NOT a robust sandwich: sandwiches the local-shift sensitivity between
## the assumed Gaussian Cov(y). Equals vcov_sandwich() only under
## correctly-specified Gaussian y with the true V_b.
.vcov_model_delta <- function(fit, parts = NULL) {
    if (is.null(parts)) parts <- implicitIF(fit)
    pp  <- fit@pp
    Z_full <- t(pp$Zt) %*% pp$U_b
    n      <- pp$n
    cov_y_unit <- as.matrix(Z_full %*% t(Z_full) + Diagonal(n))
    V <- parts$sigma^2 * (parts$IF_beta %*% cov_y_unit %*% t(parts$IF_beta))
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
        Z_full  <- t(fit@pp$Zt) %*% U_b
        cov_y_u <- as.matrix(Z_full %*% t(Z_full) + Diagonal(fit@pp$n))
        IF_b    <- as.matrix(U_b %*% parts$IF_u)
        V_b     <- sig2 * (IF_b %*% cov_y_u %*% t(IF_b))
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
