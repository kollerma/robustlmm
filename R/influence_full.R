## Full implicit influence function (beta, u, sigma, theta) and the
## joint Cook's-distance equivalent built on top of it.
##
## The (beta, u) machinery in R/influence.R is the partial IF used by
## the cluster sandwich and the case-weight IF; this file adds the
## sigma- and theta-DAS scoring equations to obtain the joint IF.
## The parameter Jacobian (d score / d par) uses numDeriv; the
## response Jacobian (d score / d y) is closed-form. WS22 speedups
## (~25x at n=240, exact -- see PLAN-WS22): the analytic d score / d y
## replaces a 2n-column finite-difference pass, and on the non-theta
## columns of the parameter Jacobian the converged tau_e and tau_b are
## reused (both are theta-only), so the DAStau (tau | s) alternation
## runs only on the L theta columns.
##
## Supported tau methods:
##   - "DASvar": closed-form quadratic-in-theta tau_e, tau_b. Fast.
##   - "DAStau", diagonal V_b: tau_e iteratively refined by calcTau
##     via 2D Gauss-Hermite quadrature; tau_b via per-block calcTau.
##   - "DAStau", block-diagonal V_b: tau_b via calcTau.nondiag proper
##     (4D quadrature per block of dimension 2; dimensions > 2 fall
##     back to the DASvar T_b inside calcTau.nondiag, mirroring the
##     estimator).
##
## Each theta perturbation in the numerical Jacobian invalidates the
## .tau_e and .Tbk caches so the perturbed tau is reconverged from a
## clean DASvar starting point.

## Full RSE score: (F_beta, F_u, F_sigma, F_theta).
.scoreVec <- function(par, fit, y = fit@resp$y, tau.fixed = NULL) {
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

    ## WS22: setTheta recomputes the (theta-only) geometry M / diagA /
    ## U_eZU_b and clears the M/unsc caches -- redundant when theta is
    ## unchanged (the numDeriv beta/u/sigma columns and any .scoreVec
    ## call at the fitted theta). Skip it then; the geometry already
    ## matches.
    if (any(as.numeric(theta) != as.numeric(pp$theta)))
        pp$setTheta(theta)
    if (fit@method == "DAStau" && !is.null(tau.fixed)) {
        ## WS22 fast path: theta is unchanged from the fit (the caller
        ## perturbed only beta/u/sigma), so tau_e is unchanged -- it
        ## depends on theta, the geometry and the rho-functions, not on
        ## beta/u/sigma. Reuse the converged tau_e instead of re-running
        ## the (tau | s) alternation. This is EXACT (not an
        ## approximation) and removes the alternation from the p+q+1
        ## non-theta columns of the parameter Jacobian. tau.fixed is a
        ## list(e = tau_e, b = tau_b): both are theta-only and reused.
        pp$.tau_e    <- tau.fixed$e
        pp$.setTau_e <- TRUE
    } else if (fit@method == "DAStau") {
        ## Each tau_e() call performs ONE step of the alternating
        ## (tau | s(tau)) scheme: calcTau converges tau given s, but s
        ## itself is computed from the cached .tau_e. The estimator
        ## reaches the joint fixed point only through repeated calls
        ## across its fit iterations, so a single call here evaluates
        ## the score at the wrong tau -- up to 21% off for weighted
        ## fits if restarted cold, because the DASvar initialization
        ## is far from the fixed point there (clean unweighted fits
        ## start ~at the fixed point, which is how this survived the
        ## original validation). Iterate to the joint fixed point the
        ## estimator solves; the cached .tau_e is kept as warm start
        ## (correctness comes from the convergence test, not the
        ## start), and an empty cache falls back to the DASvar
        ## initialization inside tau_e() itself.
        pp$.setTau_e <- FALSE
        pp$.Tbk      <- list()
        pp$.setTbk   <- FALSE
        tau_prev <- pp$tau_e()
        for (k in seq_len(50L)) {
            pp$.setTau_e <- FALSE   # keep .tau_e as the warm start
            tau_cur <- pp$tau_e()
            if (max(abs(tau_cur - tau_prev)) <
                1e-10 * max(1, max(abs(tau_prev)))) break
            tau_prev <- tau_cur
        }
        if (k == 50L)
            warning("tau_e alternating iteration did not converge ",
                    "within 50 steps in .scoreVec")
    }

    U_e    <- pp$U_e
    U_b    <- pp$U_b
    X      <- pp$X
    Z      <- t(pp$Zt)
    off    <- fit@resp$offset

    ## Prior-weight-whitened residuals: U_e = diag(sqrt(w)), matching
    ## resp$wtres / sigma (see std.e() in helpers.R). Offset-corrected.
    r     <- as.numeric(y - off - X %*% beta - Z %*% (U_b %*% u))
    r_std <- as.numeric(U_e %*% r) / sigma
    u_std <- u / sigma

    psi_e       <- rho_e@psi(r_std)
    lambda_e    <- rho_e@EDpsi()
    lambda_b_q  <- pp$D_b@x
    kappa_sigma <- pp$kappa_e
    ## WS11: the e-side score is eta-weighted (X~' H psi_e, Z~' H psi_e,
    ## and the H-weighted scale equation), matching the estimator
    ## (fitEffects + updateSigma). eta == 1 reproduces the unweighted
    ## score exactly. The b-side penalty (psi_b) and the theta score
    ## are not eta-weighted (eta enters theta only via tau_b).
    eta <- pp$eta

    F_beta <- as.numeric(crossprod(X, U_e %*% (eta * psi_e)))

    psi_b  <- .applyBlockPsi_b(fit, u_std, rho_b)
    F_u    <- as.numeric(crossprod(U_b, crossprod(Z, U_e %*% (eta * psi_e)))) -
              (lambda_e / lambda_b_q) * psi_b

    tau_e <- pp$tau_e()
    if (any(is.na(tau_e)))
        stop("tau_e contains NA (", sum(is.na(tau_e)), " of ", length(tau_e),
             ") at perturbed theta = (",
             paste(signif(theta, 6), collapse = ", "), ").")

    arg_sigma <- r_std / tau_e
    w_sigma   <- rho_se@wgt(arg_sigma)
    F_sigma   <- sum(eta * tau_e^2 * w_sigma * (arg_sigma^2 - kappa_sigma))

    block_sizes <- vapply(fit@idx, nrow, integer(1))
    if (all(block_sizes == 1L)) {
        ## WS22: reuse the converged tau_b on the fast path (theta-only).
        F_theta <- .scoreTheta_diag(fit, u_std, sigma, pp, rho_b, rho_sb,
                                    tau_b = tau.fixed$b)
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
.scoreTheta_diag <- function(fit, u_std, sigma, pp, rho_b, rho_sb,
                             tau_b = NULL) {
    ## WS22: tau_b is theta-only; the caller may supply the converged
    ## value (fast path) instead of recomputing it via calcTau.
    if (is.null(tau_b)) tau_b <- .getTauB_diag(fit, pp, rho_b, rho_sb)
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

## Per-block T_b for block-diagonal V_b, dispatching on the tau method
## exactly like rlmer.fit.DAS.nondiag does (rlmer.R, switch on method):
## DASvar uses the closed-form pp$Tb(); DAStau iterates calcTau.nondiag
## (4D Gauss-Hermite quadrature for blocks of dimension 2; dimensions
## > 2 fall back to DASvar inside calcTau.nondiag, with a warning).
## The .Tbk cache has been cleared by .scoreVec, so calcTau.nondiag
## reconverges from a clean DASvar starting point at the current theta.
.getTbFull_block <- function(fit, pp) {
    if (fit@method == "DASvar") return(pp$Tb())
    ghZ   <- as.matrix(expand.grid(pp$ghz, pp$ghz, pp$ghz, pp$ghz))
    ghZ12 <- ghZ[, 1:2]
    ghZ34 <- ghZ[, 3:4]
    ghw   <- apply(as.matrix(expand.grid(pp$ghw, pp$ghw, pp$ghw, pp$ghw)),
                   1, prod)
    calcTau.nondiag(fit, ghZ12, ghZ34, ghw, .S(fit), .kappa_b(fit),
                    max.iter = 200L, rel.tol = 1e-6)
}

## Block-diagonal-V_b theta-scoring: matches rlmer.fit.DAS.nondiag.
.scoreTheta_block <- function(fit, u_std, sigma, pp, rho_sb) {
    Tb_full <- .getTbFull_block(fit, pp)
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
##' \sigma, \theta)} is computed via \code{\link[numDeriv:jacobian]{numDeriv::jacobian}}
##' with Richardson extrapolation, reusing the converged DAS scales on
##' the non-\eqn{\theta} columns (they depend on \eqn{\theta} only); the
##' Jacobian wrt the response is closed-form. Only \code{DASvar} and
##' \code{DAStau} methods are supported (DAStau additionally requires
##' block sizes \eqn{\le 2}).
##'
##' This function is the engine behind \code{\link[=cooks.distance.rlmerMod]{cooks.distance}}
##' / \code{\link[=influence.rlmerMod]{influence}}; users who only want
##' Cook's distance should call those S3 methods.
##'
##' The result is cached on the fit (keyed by \eqn{\theta}) when
##' \code{use.cache = TRUE}, so repeated calls -- and the consumers
##' \code{\link[=cooks.distance.rlmerMod]{cooks.distance}}, the sandwich
##' \code{vcov}, and the Satterthwaite \code{df} of \code{summary} /
##' \code{anova} / \code{emmeans} / \code{confint} -- share a single
##' computation.
##'
##' @param fit An \code{rlmerMod} object.
##' @param eps Numerical step for \code{numDeriv::jacobian} (default
##'   \code{1e-6}).
##' @param use.cache If \code{TRUE} (default), return a result cached on
##'   the fit when one is stored for the current \eqn{\theta}, and cache
##'   the result otherwise.
##' @return A list with components \code{IF_beta} (\eqn{p \times n}),
##'   \code{IF_u} (\eqn{q \times n}), \code{IF_sigma} (\eqn{1 \times n}),
##'   \code{IF_theta} (\eqn{L \times n}), and the model-based delta-method
##'   covariance \code{vcov_model_delta}.
##' @seealso \code{\link{implicitIF}},
##'   \code{\link[=cooks.distance.rlmerMod]{cooks.distance}}
##' @export
implicitIF_full <- function(fit, eps = 1e-6, use.cache = TRUE) {
    stopifnot(is(fit, "rlmerMod"))
    if (!fit@method %in% c("DASvar", "DAStau"))
        stop("Only DASvar and DAStau methods supported (got '",
             fit@method, "').")

    ## WS16: the full IF is the documented bottleneck and is shared by
    ## cooks.distance / vcov_sandwich / the Satterthwaite df. Cache it on
    ## the (mutable) predictor module, keyed by theta -- the IF is a
    ## property of the converged fit, and the only in-place mutation of pp
    ## is the transient setTheta() in the df finite-difference loop, which
    ## restores theta on exit. A theta-keyed cache therefore stays valid
    ## across calls (summary / anova / emmeans / confint share one
    ## computation) without being invalidated by that FD.
    ppc <- fit@pp
    theta_now <- as.numeric(getME(fit, "theta"))
    if (use.cache && length(ppc$cache.IFfull) &&
        length(ppc$cache.IFfull.theta) == length(theta_now) &&
        isTRUE(all.equal(ppc$cache.IFfull.theta, theta_now,
                         tolerance = 1e-10)))
        return(ppc$cache.IFfull)
    block_sizes <- vapply(fit@idx, nrow, integer(1))
    if (fit@method == "DAStau" && any(block_sizes > 2L))
        warning("DAStau uses the DASvar T_b for blocks of dimension > 2 ",
                "(both in the fit and in this influence function).",
                call. = FALSE)

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
    ## snapshot the converged tau caches as the estimator left them, so
    ## the fit object is restored EXACTLY on exit. (Merely emptying the
    ## caches is not enough: tau_e() performs a single alternating
    ## (tau | s(tau)) step per call, so a later cold call would land off
    ## the joint fixed point -- materially so for weighted fits.)
    saved_tau_e    <- pp$.tau_e
    saved_setTau_e <- pp$.setTau_e
    saved_Tbk      <- pp$.Tbk
    saved_setTbk   <- pp$.setTbk
    on.exit({
        pp$setTheta(saved_theta)
        if (fit@method == "DAStau") {
            pp$.tau_e    <- saved_tau_e
            pp$.setTau_e <- saved_setTau_e
            pp$.Tbk      <- saved_Tbk
            pp$.setTbk   <- saved_setTbk
        }
    })

    ## WS22: numDeriv perturbs one coordinate per column. For the
    ## p+q+1 non-theta columns theta is unchanged, so reuse the
    ## converged tau_e (exact); only the L theta columns re-run the
    ## (tau | s) alternation. For DASvar tau_e is closed-form, so no
    ## fast path is needed.
    theta_pos <- p + q + 1L + seq_len(L)
    ## converged tau_e and tau_b at the fit (both theta-only); reused on
    ## every non-theta column. tau_b only for the diagonal-V_b case.
    tau_ref <- NULL
    if (fit@method == "DAStau" && all(block_sizes == 1L)) {
        tau_b0  <- .getTauB_diag(fit, pp, fit@rho.b, fit@rho.sigma.b)
        tau_ref <- list(e = saved_tau_e, b = tau_b0)
    }
    scoreFun  <- function(par) {
        tf <- if (!is.null(tau_ref) &&
                  max(abs(par[theta_pos] - theta0)) < 1e-12) tau_ref else NULL
        .scoreVec(par, fit, y0, tau.fixed = tf)
    }
    Jpar <- numDeriv::jacobian(
        scoreFun, par0,
        method = "Richardson",
        method.args = list(eps = eps, d = 1e-4, r = 4))

    ## Jy = d score / d y, ANALYTIC (WS22). The score depends on y only
    ## through the whitened residual r_std = U_e (y - offset - X beta -
    ## Z U_b u) / sigma. At fixed par: tau is y-independent and F_theta
    ## does not depend on y at all (zero rows). This replaces the former
    ## 2n .scoreVec finite-difference calls -- each of which re-ran the
    ## DAStau (tau | s) alternation redundantly -- with closed-form
    ## matrix products plus one vectorized scalar FD for the scale row.
    ## numDeriv (Jpar above) left pp at a perturbed par; restore it to
    ## par0 (also reconverges tau) before reading the fit-state pieces.
    .scoreVec(par0, fit, y0)
    Dr      <- diag(pp$U_e) / sigma0          # d r_std_i / d y_i
    r_std0  <- as.numeric(pp$U_e %*% (y0 - fit@resp$offset -
                          pp$X %*% beta0 -
                          crossprod(pp$Zt, pp$U_b %*% u0))) / sigma0
    eta_v   <- pp$eta
    w_y     <- eta_v * fit@rho.e@Dpsi(r_std0) * Dr   # e-side y-weight
    Jy <- matrix(0, nrow = p + q + 1L + L, ncol = n)
    Jy[1:p, ]       <- t(as.matrix(pp$U_eX)    * w_y)   # d F_beta / d y
    Jy[p + 1:q, ]   <- t(as.matrix(pp$U_eZU_b) * w_y)   # d F_u    / d y
    ## scale row: d/dy of  eta tau^2 w_sigma(r/tau) ((r/tau)^2 - kappa);
    ## per-observation, so a vectorized central difference in r_std is exact.
    tau_e0 <- pp$tau_e(); kap <- pp$kappa_e; rse <- fit@rho.sigma.e
    s_of <- function(r) { a <- r / tau_e0
        eta_v * tau_e0^2 * rse@wgt(a) * (a^2 - kap) }
    dlt <- 1e-6
    Jy[p + q + 1L, ] <- (s_of(r_std0 + dlt) - s_of(r_std0 - dlt)) /
                        (2 * dlt) * Dr
    ## F_theta rows stay 0 (theta score is y-independent at fixed par).

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
    ## Var(y) = sigma^2 (Z U_b U_b' Z' + W^{-1}), W = diag(prior weights)
    cov_y  <- as.matrix(sigma0^2 *
                        (tcrossprod(Z_full) +
                         Diagonal(x = 1 / diag(pp$U_e)^2)))

    res <- list(IF_beta  = IF_beta,
                IF_u     = IF_u,
                IF_sigma = IF_sigma,
                IF_theta = IF_theta,
                vcov_model_delta = tcrossprod(IF_beta %*% cov_y, IF_beta),
                Jpar       = Jpar,
                par0       = par0,
                solve_info = solve_info)
    if (use.cache) {
        ppc$cache.IFfull       <- res
        ppc$cache.IFfull.theta <- theta_now
    }
    res
}

## WS16: the cached implicitIF_full(fit), or NULL if none is stored for
## the current theta. Lets summary(df = "auto") reuse an IF that another
## call (cooks.distance / vcov_sandwich / confint / a prior summary)
## already paid for, regardless of fit size.
.IFfullCached <- function(fit) {
    pp <- fit@pp
    theta_now <- as.numeric(getME(fit, "theta"))
    if (length(pp$cache.IFfull) &&
        length(pp$cache.IFfull.theta) == length(theta_now) &&
        isTRUE(all.equal(pp$cache.IFfull.theta, theta_now, tolerance = 1e-10)))
        pp$cache.IFfull
    else NULL
}

## WS15: Mahalanobis norm of a (d x k) influence matrix W under the
## empirical metric V = (W W')/k (svd pseudo-inverse at a boundary).
## Shared by the per-observation and per-cluster Cook's distances.
.cooksMahalanobis <- function(W) {
    k <- ncol(W)
    V <- tcrossprod(W) / k
    V <- (V + t(V)) / 2
    V_inv <- tryCatch(solve(V), error = function(e) {
        sv  <- svd(V); tol <- max(dim(V)) * .Machine$double.eps * sv$d[1L]
        tcrossprod(sv$v %*% diag(ifelse(sv$d > tol, 1 / sv$d, 0)), sv$u)
    })
    sqrt(colSums(W * (V_inv %*% W)))
}

##' Cook's-distance equivalent for an rlmerMod fit (per observation or
##' per cluster).
##'
##' Joint Mahalanobis influence on the fitted \eqn{(\hat{\beta},
##' \hat{\sigma}, \hat{\theta})}. With \code{groups = NULL} (default) the
##' unit is the observation: the per-observation influence vectors are
##' stacked into a \eqn{(p + 1 + L) \times n} matrix \eqn{W}, and the
##' result is \eqn{\sqrt{w_i^T V^{-1} w_i}} with \eqn{V = (1/n) W W^T}.
##' With \code{groups} set the unit is the cluster: the
##' per-cluster influence functions \eqn{\Xi = -J_{par}^{-1} S} (\eqn{S}
##' the per-cluster score contributions from \code{.scoreByCluster}),
##' restricted to the \eqn{(\beta, \sigma, \theta)} rows, are scored the
##' same way with \eqn{V = (1/J) \Xi \Xi^T}. This flags whole-cluster
##' outliers that no single observation reveals -- e.g. before relying on
##' the bootstrap variance-component test of \code{\link{anova}}, which
##' is anti-conservative under group contamination. If \eqn{V} is
##' singular (a variance-component boundary) the Moore-Penrose
##' pseudo-inverse is used.
##'
##' The full IF computation is the expensive part; pre-compute it once
##' via \code{IF = implicitIF_full(fit)} and pass it in if you need
##' \code{cooks.distance} and \code{influence} together.
##'
##' @param model An \code{rlmerMod} object.
##' @param groups Cluster grouping for cluster-level Cook's distance:
##'   \code{NULL} (default) for per-observation; \code{TRUE} for the
##'   top-level grouping factor; or a factor / grouping-factor name to
##'   aggregate over. Cluster-level requires a single or nested grouping
##'   structure (crossed designs error, as for the variance-parameter
##'   covariance).
##' @param IF Optional pre-computed \code{\link{implicitIF_full}(model)};
##'   computed on demand if \code{NULL}.
##' @param ... Currently unused.
##' @return Named numeric vector: one entry per observation
##'   (\code{groups = NULL}) or per cluster level.
##' @seealso \code{\link{implicitIF_full}},
##'   \code{\link[=influence.rlmerMod]{influence}},
##'   \code{\link[=hatvalues.rlmerMod]{hatvalues}}
##' @method cooks.distance rlmerMod
##' @export
cooks.distance.rlmerMod <- function(model, groups = NULL, IF = NULL, ...) {
    stopifnot(is(model, "rlmerMod"))
    if (is.null(IF)) IF <- implicitIF_full(model)
    if (is.null(groups)) {
        IF_all <- rbind(IF$IF_beta, IF$IF_sigma, IF$IF_theta)
        d <- .cooksMahalanobis(IF_all)
        names(d) <- paste0("obs", seq_len(ncol(IF_all)))
        return(d)
    }
    ## cluster-level: per-cluster influence Xi = -Jpar^{-1} S,
    ## restricted to (beta, sigma, theta) -- mirrors the per-observation
    ## version, which drops the u rows.
    S   <- .scoreByCluster(model)           # errors on crossed designs
    Xi  <- .applyJparInv(IF, S)             # d_full x J (top-cluster cols)
    pp  <- model@pp; p <- pp$p; q <- pp$q
    L   <- length(getME(model, "theta"))
    sel <- c(seq_len(p), p + q + 1L, p + q + 1L + seq_len(L))
    W   <- Xi[sel, , drop = FALSE]
    colnames(W) <- colnames(S)              # top-cluster level labels
    ## TRUE or the top-cluster factor -> columns as-is; a coarser factor
    ## -> sum the independent top-cluster columns into it (still valid); a
    ## finer factor would break independence -> error.
    if (!isTRUE(groups)) {
        top  <- .top_cluster(model)$g
        f    <- .resolveGrouping(model, groups)
        ## each top cluster must fall in exactly one requested group
        ## (i.e. `groups` is the top cluster or coarser)
        if (any(rowSums(table(top, f) > 0L) != 1L))
            stop("`groups` must be the top-level cluster or coarser ",
                 "(constant within each '", .top_cluster(model)$name,
                 "' cluster); a finer grouping breaks the independence ",
                 "the cluster influence relies on.")
        fmap <- f[match(colnames(W), as.character(top))]   # group per column
        if (nlevels(droplevels(as.factor(fmap))) < ncol(W))
            W <- t(rowsum(t(W), fmap))
    }
    d <- .cooksMahalanobis(W)
    names(d) <- colnames(W)
    d
}

##' Robust leverage (hat values) for an rlmerMod fit.
##'
##' The self-leverage \eqn{A_{ii}} of each observation in the robust,
##' random-effect-whitened convolution that the estimator uses internally
##' -- the robust analogue of the linear mixed-model hat value. It
##' reduces to the classical \code{\link[lme4]{lmer}} leverage at the
##' non-robust limit (\code{rho = cPsi}); at the robust default the
##' effective degrees of freedom \eqn{\sum_i A_{ii}} differ as
##' downweighting changes each observation's pull. With \code{groups} set,
##' the per-observation leverages are summed within each cluster, giving
##' the cluster's leverage (its effective-df contribution).
##'
##' @param model An \code{rlmerMod} object.
##' @param groups Leverage aggregation: \code{NULL} (default) for
##'   per-observation; \code{TRUE} for the top-level grouping factor; or a
##'   factor / grouping-factor name to aggregate over.
##' @param ... Currently unused.
##' @return Named numeric vector: one entry per observation
##'   (\code{groups = NULL}) or per cluster level.
##' @seealso \code{\link[=cooks.distance.rlmerMod]{cooks.distance}},
##'   \code{\link{rlmer}} (the \code{design.weights} argument bounds
##'   high-leverage design points)
##' @method hatvalues rlmerMod
##' @export
hatvalues.rlmerMod <- function(model, groups = NULL, ...) {
    stopifnot(is(model, "rlmerMod"))
    h <- as.numeric(model@pp$diagA)
    n <- model@pp$n
    if (length(h) != n) h <- rep.int(NA_real_, n)  # no-RE / unset guard
    if (is.null(groups)) {
        names(h) <- paste0("obs", seq_len(n))
        return(h)
    }
    f <- .resolveGrouping(model, groups)
    out <- tapply(h, f, sum)
    out[is.na(out)] <- 0
    c(out)
}

## WS15: resolve the `groups` argument of the cluster-level diagnostics to
## a factor of length n. TRUE -> the top-level cluster factor; a string ->
## the named grouping factor; a factor/vector -> used as given.
.resolveGrouping <- function(model, groups) {
    if (isTRUE(groups)) {
        tc <- .top_cluster(model)
        if (is.null(tc))
            stop("no top-level cluster (crossed grouping factors); pass an ",
                 "explicit grouping factor to `groups`.")
        return(tc$g)
    }
    if (is.character(groups) && length(groups) == 1L) {
        if (!groups %in% names(model@flist))
            stop("grouping factor '", groups, "' not found; available: ",
                 paste(names(model@flist), collapse = ", "), ".")
        return(model@flist[[groups]])
    }
    g <- as.factor(groups)
    if (length(g) != model@pp$n)
        stop("`groups` must have length n = ", model@pp$n, ".")
    g
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
