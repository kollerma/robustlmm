## Full implicit influence function (beta, u, sigma, theta) and the
## joint Cook's-distance equivalent built on top of it.
##
## The (beta, u) machinery in R/influence.R is the partial IF used by
## the cluster sandwich and the case-weight IF; this file adds the
## sigma- and theta-DAS scoring equations to obtain the joint IF, by
## differentiating the hand-coded score function numerically with
## numDeriv. The Jacobian wrt y is a finite-difference column-by-column
## pass.
##
## Supported tau methods:
##   - "DASvar": closed-form quadratic-in-theta tau_e, tau_b. Fast.
##   - "DAStau", diagonal V_b: tau_e iteratively refined by calcTau
##     via 2D Gauss-Hermite quadrature; tau_b via per-block calcTau.
##   - "DAStau", block-diagonal V_b: tau_b via calcTau.nondiag proper
##     (4D quadrature per block size 2). Block sizes > 2 are not
##     supported by robustlmm itself.
##
## Each theta perturbation in the numerical Jacobian invalidates the
## .tau_e and .Tbk caches so the perturbed tau is reconverged from a
## clean DASvar starting point.

## Full RSE score: (F_beta, F_u, F_sigma, F_theta).
.scoreVec <- function(par, fit, y = fit@resp$y) {
    pp     <- fit@pp
    rho_e  <- fit@rho.e
    rho_b  <- fit@rho.b
    rho_se <- fit@rho.sigma.e
    rho_sb <- fit@rho.sigma.b
    L      <- length(getME(fit, "theta"))
    p      <- pp$p; q <- pp$q

    beta  <- par[1:p]
    u     <- par[p + 1:q]
    sigma <- exp(par[p + q + 1L])
    theta <- par[p + q + 1L + 1:L]

    pp$setTheta(theta)
    if (fit@method == "DAStau") {
        pp$.tau_e    <- numeric(0)
        pp$.setTau_e <- FALSE
        pp$.Tbk      <- list()
        pp$.setTbk   <- FALSE
    }

    U_e    <- pp$U_e
    U_b    <- pp$U_b
    X      <- pp$X
    Z      <- t(pp$Zt)
    Ue_inv <- Diagonal(x = 1 / diag(U_e))

    r     <- as.numeric(y - X %*% beta - Z %*% (U_b %*% u))
    r_std <- as.numeric(Ue_inv %*% r) / sigma
    u_std <- u / sigma

    psi_e       <- rho_e@psi(r_std)
    lambda_e    <- rho_e@EDpsi()
    lambda_b_q  <- pp$D_b@x
    kappa_sigma <- pp$kappa_e

    F_beta <- as.numeric(t(X) %*% (Ue_inv %*% psi_e))

    psi_b  <- .applyBlockPsi_b(fit, u_std, rho_b)
    F_u    <- as.numeric(t(U_b) %*% (t(Z) %*% (Ue_inv %*% psi_e))) -
              (lambda_e / lambda_b_q) * psi_b

    tau_e <- pp$tau_e()
    if (any(is.na(tau_e))) stop("tau_e contains NA at perturbed theta.")

    arg_sigma <- r_std / tau_e
    w_sigma   <- rho_se@wgt(arg_sigma)
    F_sigma   <- sum(tau_e^2 * w_sigma * (arg_sigma^2 - kappa_sigma))

    block_sizes <- vapply(fit@idx, nrow, integer(1))
    if (all(block_sizes == 1L)) {
        F_theta <- .scoreTheta_diag(fit, u_std, sigma, pp, rho_b, rho_sb)
    } else {
        F_theta <- .scoreTheta_block(fit, u_std, sigma, pp, rho_sb)
    }
    c(F_beta, F_u, F_sigma, F_theta)
}

## tau_b per RE for the diagonal-V_b case, matching updateThetaTau.
.getTauB_diag <- function(fit, pp, rho_b, rho_sb) {
    q <- pp$q
    tau_b <- numeric(q)
    if (fit@method == "DASvar") {
        tau_b[] <- sqrt(diag(pp$Tb()))
        return(tau_b)
    }
    skbs <- .s(fit, theta = TRUE, pp = pp)
    kappas_per_type <- unlist(pp$kappa_b)
    L_diag <- diag(pp$L)
    for (type in seq_along(rho_b)) {
        bidx_v <- as.vector(fit@idx[[type]])
        tau_per_re <- calcTau(L_diag[bidx_v],
                              skbs[bidx_v],
                              rho_b[[type]],
                              rho_sb[[type]],
                              pp,
                              kappas_per_type[type])
        tau_b[bidx_v] <- tau_per_re
    }
    tau_b
}

## Diagonal-V_b theta-scoring: same equations rlmer's updateThetaTau
## iterates against. Returns L = length(theta) values, one per block-type.
.scoreTheta_diag <- function(fit, u_std, sigma, pp, rho_b, rho_sb) {
    tau_b  <- .getTauB_diag(fit, pp, rho_b, rho_sb)
    arg_th <- u_std / tau_b
    L      <- length(rho_sb)
    F_theta <- numeric(L)
    kappa_per_block <- pp$kappa_b
    for (l in seq_len(L)) {
        is_in  <- (fit@ind == l)
        w_th   <- rho_sb[[l]]@wgt(arg_th[is_in])
        F_theta[l] <- sum(w_th * u_std[is_in]^2) -
                      kappa_per_block[[l]] * sum(w_th * tau_b[is_in]^2)
    }
    F_theta
}

## Block-diagonal-V_b theta-scoring: matches rlmer.fit.DAS.nondiag.
.scoreTheta_block <- function(fit, u_std, sigma, pp, rho_sb) {
    Tb_full <- pp$Tb()
    kappa_per_block <- pp$kappa_b
    F_theta <- numeric(0)
    for (type in seq_along(fit@idx)) {
        idx_mat <- fit@idx[[type]]
        s <- nrow(idx_mat); K <- ncol(idx_mat)
        lhs <- matrix(0, s, s); rhs <- matrix(0, s, s)
        psi_fn  <- rho_sb[[type]]@psi
        wgt_fn  <- rho_sb[[type]]@wgt
        kappa_t <- kappa_per_block[[type]]
        for (k in seq_len(K)) {
            bk_idx  <- as.vector(idx_mat[, k])
            b_block <- u_std[bk_idx]
            Tk      <- as.matrix(Tb_full[bk_idx, bk_idx, drop = FALSE])
            Lk      <- t(chol(Tk))
            T_b     <- forwardsolve(Lk, b_block)
            d_k     <- sum(T_b * T_b)
            if (s > 1L) {
                w_eta   <- wgt_fn(d_k)
                w_delta <- (psi_fn(d_k) -
                            psi_fn(d_k - s * kappa_t)) / s
            } else {
                w_eta   <- wgt_fn(d_k)
                w_delta <- w_eta * kappa_t
            }
            lhs <- lhs + w_eta   * tcrossprod(b_block)
            rhs <- rhs + w_delta * Tk
        }
        lhs  <- lhs / K
        rhs  <- rhs / K
        diff <- lhs - rhs
        pat  <- fit@blocks[[type]] != 0
        Lind <- fit@blocks[[type]][pat]
        entries <- diff[pat]
        ord <- order(Lind)
        F_theta <- c(F_theta, entries[ord])
    }
    F_theta
}

##' Full implicit influence function (beta, u, sigma, theta) for a
##' fitted \code{rlmerMod} object.
##'
##' Extends \code{\link{implicitIF}} by adding the \eqn{\sigma}-DAS and
##' \eqn{\theta}-DAS scoring equations to the implicit-derivative linear
##' system. The Jacobian of the full score wrt \eqn{(\beta, u, \log
##' \sigma, \theta)} is computed numerically via \code{\link[numDeriv:jacobian]{numDeriv::jacobian}}
##' with Richardson extrapolation; the Jacobian wrt the response is a
##' central finite difference column by column. Only \code{DASvar} and
##' \code{DAStau} methods are supported (DAStau additionally requires
##' block sizes \eqn{\le 2}).
##'
##' This function is the engine behind \code{\link[=cooks.distance.rlmerMod]{cooks.distance}}
##' / \code{\link[=influence.rlmerMod]{influence}}; users who only want
##' Cook's distance should call those S3 methods.
##'
##' The numerical Jacobian is slow (a few seconds on Penicillin); each
##' call recomputes it from scratch.
##'
##' @param fit An \code{rlmerMod} object.
##' @param eps Numerical step for \code{numDeriv::jacobian} (default
##'   \code{1e-6}).
##' @return A list with components \code{IF_beta} (\eqn{p \times n}),
##'   \code{IF_u} (\eqn{q \times n}), \code{IF_sigma} (\eqn{1 \times n}),
##'   \code{IF_theta} (\eqn{L \times n}), and the model-based delta-method
##'   covariance \code{vcov_model_delta}.
##' @seealso \code{\link{implicitIF}},
##'   \code{\link[=cooks.distance.rlmerMod]{cooks.distance}}
##' @export
implicitIF_full <- function(fit, eps = 1e-6) {
    stopifnot(is(fit, "rlmerMod"))
    if (!fit@method %in% c("DASvar", "DAStau"))
        stop("Only DASvar and DAStau methods supported (got '",
             fit@method, "').")
    block_sizes <- vapply(fit@idx, nrow, integer(1))
    if (fit@method == "DAStau" && any(block_sizes > 2L))
        stop("DAStau is not defined for block sizes > 2 in robustlmm ",
             "itself; refit with DASvar.")

    pp <- fit@pp
    p  <- pp$p; q <- pp$q; n <- pp$n
    L  <- length(getME(fit, "theta"))

    beta0  <- .fixef(fit)
    u0     <- pp$b.s
    sigma0 <- .sigma(fit)
    theta0 <- getME(fit, "theta")
    y0     <- fit@resp$y

    par0   <- c(beta0, u0, log(sigma0), theta0)
    names(par0) <- c(paste0("beta",  seq_len(p)),
                     paste0("u",     seq_len(q)),
                     "log_sigma",
                     paste0("theta", seq_len(L)))

    saved_theta <- theta0
    on.exit(pp$setTheta(saved_theta))

    Jpar <- numDeriv::jacobian(
        function(par) .scoreVec(par, fit, y0),
        par0,
        method = "Richardson",
        method.args = list(eps = eps, d = 1e-4, r = 4))

    Jy <- matrix(NA_real_, nrow = p + q + 1L + L, ncol = n)
    h  <- 1e-4
    for (i in seq_len(n)) {
        y_plus  <- y0; y_plus[i]  <- y_plus[i]  + h
        y_minus <- y0; y_minus[i] <- y_minus[i] - h
        Jy[, i] <- (.scoreVec(par0, fit, y_plus) -
                    .scoreVec(par0, fit, y_minus)) / (2 * h)
    }

    d_full     <- p + q + 1L + L
    solve_info <- list(rank = d_full, full = d_full, singular = FALSE,
                       min_kept_sv = NA_real_, max_dropped_sv = NA_real_)
    IF_par <- tryCatch(
        -solve(Jpar, Jy),
        error = function(e) {
            sv  <- svd(Jpar)
            tol <- max(dim(Jpar)) * .Machine$double.eps * max(sv$d)
            pos <- sv$d > tol
            rk  <- sum(pos)
            solve_info$rank           <<- rk
            solve_info$singular       <<- TRUE
            solve_info$min_kept_sv    <<- if (rk > 0L) min(sv$d[pos]) else NA_real_
            solve_info$max_dropped_sv <<- if (any(!pos)) max(sv$d[!pos]) else 0
            warning(sprintf(
                paste0("Jpar rank-deficient (rank %d of %d; smallest kept ",
                       "sv %.2e, largest dropped sv %.2e) -- likely a ",
                       "variance-component boundary. Using pseudo-inverse; ",
                       "IF is the minimum-norm solution in the identifiable ",
                       "subspace, unidentifiable direction(s) zeroed."),
                rk, d_full, solve_info$min_kept_sv, solve_info$max_dropped_sv),
                call. = FALSE)
            dinv <- numeric(length(sv$d)); dinv[pos] <- 1 / sv$d[pos]
            -((sv$v %*% (dinv * t(sv$u))) %*% Jy)
        })
    rownames(IF_par) <- names(par0)

    IF_beta   <- IF_par[1:p, , drop = FALSE]
    IF_u      <- IF_par[p + 1:q, , drop = FALSE]
    IF_logsig <- IF_par[p + q + 1L, , drop = FALSE]
    IF_sigma  <- sigma0 * IF_logsig
    IF_theta  <- IF_par[p + q + 1L + 1:L, , drop = FALSE]

    U_b    <- pp$U_b
    Z      <- t(pp$Zt)
    Z_full <- Z %*% U_b
    cov_y  <- as.matrix(sigma0^2 *
                        (Z_full %*% t(Z_full) + Diagonal(n)))

    list(IF_beta  = IF_beta,
         IF_u     = IF_u,
         IF_sigma = IF_sigma,
         IF_theta = IF_theta,
         vcov_model_delta = IF_beta %*% cov_y %*% t(IF_beta),
         Jpar       = Jpar,
         par0       = par0,
         solve_info = solve_info)
}

##' Per-observation Cook's-distance equivalent for an rlmerMod fit.
##'
##' Joint Mahalanobis influence of each observation on the fitted
##' \eqn{(\hat{\beta}, \hat{\sigma}, \hat{\theta})}. Stacks the per-
##' observation influence vectors into a \eqn{(p + 1 + L) \times n}
##' matrix, computes the empirical variance \eqn{V = (1/n) IF\,IF^T},
##' and returns \eqn{\sqrt{x_i^T V^{-1} x_i}} per observation. If
##' \eqn{V} is singular (variance-component boundary), the Moore-Penrose
##' pseudo-inverse is used.
##'
##' The full IF computation is the expensive part; pre-compute it once
##' via \code{IF = implicitIF_full(fit)} and pass it in if you need
##' \code{cooks.distance} and \code{influence} together.
##'
##' @param model An \code{rlmerMod} object.
##' @param IF Optional pre-computed \code{\link{implicitIF_full}(model)};
##'   computed on demand if \code{NULL}.
##' @param ... Currently unused.
##' @return Numeric vector of length \code{n}.
##' @seealso \code{\link{implicitIF_full}},
##'   \code{\link[=influence.rlmerMod]{influence}}
##' @method cooks.distance rlmerMod
##' @export
cooks.distance.rlmerMod <- function(model, IF = NULL, ...) {
    stopifnot(is(model, "rlmerMod"))
    if (is.null(IF)) IF <- implicitIF_full(model)
    IF_all <- rbind(IF$IF_beta, IF$IF_sigma, IF$IF_theta)
    n      <- ncol(IF_all)
    V      <- (IF_all %*% t(IF_all)) / n
    V      <- (V + t(V)) / 2
    V_inv  <- tryCatch(solve(V), error = function(e) {
        sv  <- svd(V); tol <- max(dim(V)) * .Machine$double.eps * sv$d[1L]
        sv$v %*% diag(ifelse(sv$d > tol, 1 / sv$d, 0)) %*% t(sv$u)
    })
    d <- sqrt(colSums(IF_all * (V_inv %*% IF_all)))
    names(d) <- paste0("obs", seq_len(n))
    d
}

##' Per-observation influence on \eqn{(\hat{\beta}, \hat{\sigma},
##' \hat{\theta})} for an rlmerMod fit.
##'
##' Thin wrapper that returns the stacked influence matrix used by
##' \code{\link[=cooks.distance.rlmerMod]{cooks.distance}}. Cheaper to
##' call this once and pass the result back to \code{cooks.distance} via
##' \code{IF = .} than to recompute the numerical Jacobian twice.
##'
##' @param model An \code{rlmerMod} object.
##' @param do.coef Ignored (kept for compatibility with the
##'   \code{\link[stats]{influence}} generic).
##' @param ... Currently unused.
##' @return A list with \code{IF_beta} (\eqn{p \times n}),
##'   \code{IF_sigma} (\eqn{1 \times n}), \code{IF_theta} (\eqn{L \times
##'   n}), and the full \code{implicitIF_full} object as \code{full}.
##' @seealso \code{\link{implicitIF_full}},
##'   \code{\link[=cooks.distance.rlmerMod]{cooks.distance}}
##' @method influence rlmerMod
##' @export
influence.rlmerMod <- function(model, do.coef = TRUE, ...) {
    stopifnot(is(model, "rlmerMod"))
    IF <- implicitIF_full(model)
    list(IF_beta  = IF$IF_beta,
         IF_sigma = IF$IF_sigma,
         IF_theta = IF$IF_theta,
         full     = IF)
}
