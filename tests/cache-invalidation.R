## Regression guard for the .tau_e / .Tbk cache-invalidation in
## .scoreVec() (R/influence_full.R). The full IF uses
## numDeriv::jacobian() to differentiate .scoreVec() with respect to
## (beta, u, log sigma, theta) -- the result is a valid derivative only
## if .scoreVec(par, fit, y) depends solely on (par, y). For DAStau the
## tau-DAS auxiliary scales are iteratively solved by calcTau seeded
## from a cached .tau_e (and .Tbk for block-diagonal V_b), so a stale
## cache from a previous theta would silently bias the Jacobian.
## .scoreVec() guards against this by clearing both caches on every
## DAStau invocation; this test fails (stop()) if the guard is removed.
##
## Two checks per fit (diagonal V_b on Dyestuff, block-diagonal V_b on
## sleepstudy):
##   (1) round-trip determinism: .scoreVec(par0) before and after a
##       theta excursion;
##   (2) cache-poison resistance: after corrupting .tau_e / .Tbk with
##       a 1e3 offset, .scoreVec(par0) still recomputes the clean value.
##
## Ported from IF-thread1/thread1_validation/validate_cache_invalidation.R.

require(robustlmm)
suppressMessages(require(lme4))
suppressMessages(require(Matrix))

.scoreVec <- robustlmm:::.scoreVec

TOL_DET    <- 1e-8       # path-independence: exact up to FP noise
TOL_POISON <- 1e-6       # calcTau reconverges to the same fixed point

build_par0 <- function(fit) {
    pp <- fit@pp
    p <- pp$p; q <- pp$q; L <- length(getME(fit, "theta"))
    par0 <- c(robustlmm:::.fixef(fit), pp$b.s,
              log(robustlmm:::.sigma(fit)), getME(fit, "theta"))
    names(par0) <- c(paste0("beta", seq_len(p)), paste0("u", seq_len(q)),
                     "log_sigma", paste0("theta", seq_len(L)))
    list(par0 = par0, p = p, q = q, L = L,
         theta_idx = p + q + 1L + seq_len(L))
}

check_fit <- function(name, fit) {
    pp <- fit@pp
    info <- build_par0(fit)
    par0 <- info$par0; y0 <- fit@resp$y
    saved_theta <- getME(fit, "theta")
    on.exit(pp$setTheta(saved_theta), add = TRUE)

    ## A sizeable but feasible theta-excursion (theta >= 0 since theta
    ## parameterises a variance scale).
    parP <- par0
    parP[info$theta_idx] <- pmax(parP[info$theta_idx] * 1.5 + 0.05, 1e-3)

    ## (1) round-trip determinism (forward and reverse).
    s0a <- .scoreVec(par0, fit, y0)
    invisible(.scoreVec(parP, fit, y0))     # excursion
    s0b <- .scoreVec(par0, fit, y0)         # return
    d_round_0 <- max(abs(s0a - s0b))

    sPa <- .scoreVec(parP, fit, y0)
    invisible(.scoreVec(par0, fit, y0))     # excursion the other way
    sPb <- .scoreVec(parP, fit, y0)
    d_round_P <- max(abs(sPa - sPb))

    ## (2) cache-poison resistance. Establish caches at par0, then
    ## corrupt the stored tau and assert the next .scoreVec call
    ## ignores the poison and recomputes from a clean DASvar start.
    invisible(.scoreVec(par0, fit, y0))
    if (length(pp$.tau_e) > 0L) {
        pp$.tau_e    <- pp$.tau_e + 1e3
        pp$.setTau_e <- TRUE
    }
    if (length(pp$.Tbk) > 0L) {
        pp$.Tbk    <- lapply(pp$.Tbk, function(M) M + 1e3)
        pp$.setTbk <- TRUE
    }
    s0_poison <- .scoreVec(par0, fit, y0)
    d_poison  <- max(abs(s0a - s0_poison))

    stopifnot(d_round_0 < TOL_DET)
    stopifnot(d_round_P < TOL_DET)
    stopifnot(d_poison  < TOL_POISON)
    invisible(TRUE)
}

## Diagonal V_b: Dyestuff with DAStau (the only path that uses .tau_e).
fit_dye <- rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, method = "DAStau",
                 rho.e = smoothPsi, rho.b = smoothPsi)
check_fit("Dyestuff (diagonal V_b)", fit_dye)

## Block-diagonal V_b (size 2): sleepstudy with DAStau (uses both
## .tau_e and .Tbk).
fit_sleep <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                   method = "DAStau",
                   rho.e = smoothPsi, rho.b = smoothPsi)
check_fit("Sleepstudy (block-diagonal V_b)", fit_sleep)
