## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## WS8 regression tests for the DAStau tau_b machinery.
##
## 1. Theta-score root consistency for DAStau, block-diagonal V_b:
##    .scoreVec must evaluate the *same* theta-equation the estimator
##    solved, i.e. with T_b from calcTau.nondiag, not the closed-form
##    DASvar pp$Tb(). At the converged fit the theta-score must then
##    be (numerically) zero. Before WS8 the score used pp$Tb()
##    unconditionally; on contaminated data the two T_b's differ and
##    the fitted theta is *not* a root of the mis-specified score.
##    (On clean, balanced data the two nearly coincide, which is how
##    the original port's validation passed.)
##
## 2. Blocks of dimension > 2: rlmer() downgrades DAStau to DASvar
##    up front (rlmer.R, "Falling back to method 'DASvar'"), so the
##    fit must come out as a working DASvar fit, and the full
##    influence function must go through on it. (calcTau.nondiag's
##    own dimension->2 limitation used to be a placeholder stop("yes,
##    do as promised") reachable by direct callers; it now falls back
##    to the DASvar T_b per block type.)

suppressMessages(require(robustlmm))

## ----------------------------------------------------------------
## 1. Theta-score root consistency on contaminated sleepstudy
## ----------------------------------------------------------------
ss <- sleepstudy
set.seed(3)
ii <- sample(nrow(ss), 12L)
ss$Reaction[ii] <- ss$Reaction[ii] + 200

fitT <- rlmer(Reaction ~ Days + (Days | Subject), ss, method = "DAStau")

p <- fitT@pp$p; q <- fitT@pp$q
L <- length(getME(fitT, "theta"))
par0 <- c(robustlmm:::.fixef(fitT), fitT@pp$b.s,
          log(robustlmm:::.sigma(fitT)), getME(fitT, "theta"))

s0       <- robustlmm:::.scoreVec(par0, fitT)
F_theta0 <- s0[p + q + 1L + seq_len(L)]

## reference scale: the theta-score response to a 5% theta perturbation
parP <- par0
parP[p + q + 2L] <- parP[p + q + 2L] * 1.05
sP       <- robustlmm:::.scoreVec(parP, fitT)
F_thetaP <- sP[p + q + 1L + seq_len(L)]

ratio <- max(abs(F_theta0)) / max(abs(F_thetaP - F_theta0))
cat("DAStau block theta-score at fit:", format(max(abs(F_theta0))),
    "; at 5% theta perturbation:", format(max(abs(F_thetaP - F_theta0))),
    "; ratio:", format(ratio), "\n")
stopifnot(ratio < 0.05)

## the (beta, u, sigma) components must vanish at the fit as well
F_head <- s0[seq_len(p + q + 1L)]
sH     <- robustlmm:::.scoreVec(par0 * 1.05, fitT)
stopifnot(max(abs(F_head)) < 0.05 * max(abs(sH[seq_len(p + q + 1L)] - F_head)))

## ----------------------------------------------------------------
## 1b. Sigma-score root consistency on a WEIGHTED DAStau fit
## ----------------------------------------------------------------
## tau_e() performs one alternating (tau | s(tau)) step per call; the
## score must evaluate tau at the joint fixed point the estimator
## solved. Before the fix, .scoreVec restarted tau cold from the
## DASvar initialization, which for weighted fits is ~20% off the
## fixed point: the sigma-score at the fit was -3.65 (ratio 8.09 vs a
## 5% perturbation response) and the full-IF FD errors were flat in h
## at 0.27-0.37. Unweighted fits start ~at the fixed point, so
## section 1 alone cannot detect this.
set.seed(101)
wW <- runif(nrow(Dyestuff), 0.5, 3)
fitW <- suppressWarnings(
    rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = wW,
          method = "DAStau", init = lmerNoFit))

pW <- fitW@pp$p; qW <- fitW@pp$q
LW <- length(getME(fitW, "theta"))
parW <- c(robustlmm:::.fixef(fitW), fitW@pp$b.s,
          log(robustlmm:::.sigma(fitW)), getME(fitW, "theta"))
sW0 <- suppressWarnings(robustlmm:::.scoreVec(parW, fitW))
parWP <- parW
parWP[pW + qW + 2L] <- parWP[pW + qW + 2L] * 1.05
sWP <- suppressWarnings(robustlmm:::.scoreVec(parWP, fitW))

F_sig0 <- sW0[pW + qW + 1L]
ratioW <- abs(F_sig0) / abs(sWP[pW + qW + 1L] - F_sig0)
cat("weighted DAStau sigma-score at fit:", format(F_sig0),
    "; ratio vs 5% theta perturbation:", format(ratioW), "\n")
stopifnot(ratioW < 0.05)

F_thW0 <- sW0[pW + qW + 1L + seq_len(LW)]
ratioWt <- max(abs(F_thW0)) /
    max(abs(sWP[pW + qW + 1L + seq_len(LW)] - F_thW0))
cat("weighted DAStau theta-score at fit ratio:", format(ratioWt), "\n")
stopifnot(ratioWt < 0.05)

## ----------------------------------------------------------------
## 2. DAStau with a block of dimension 3: fallback instead of crash
## ----------------------------------------------------------------
set.seed(11)
nG <- 12L; nPer <- 6L; n <- nG * nPer
g  <- factor(rep(seq_len(nG), each = nPer))
x1 <- rnorm(n); x2 <- rnorm(n)
b  <- matrix(rnorm(3L * nG, sd = c(1, 0.8, 0.6)), nrow = nG, byrow = TRUE)
y  <- 1 + x1 + 0.5 * x2 +
    b[as.integer(g), 1L] + b[as.integer(g), 2L] * x1 +
    b[as.integer(g), 3L] * x2 + rnorm(n, sd = 0.5)
d3 <- data.frame(y, x1, x2, g)

warned <- character(0)
fit3 <- withCallingHandlers(
    rlmer(y ~ x1 + x2 + (x1 + x2 | g), d3, method = "DAStau"),
    warning = function(w) {
        warned <<- c(warned, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
stopifnot(any(grepl("falling back to method 'DASvar'", warned,
                    ignore.case = TRUE)),
          identical(fit3@method, "DASvar"),
          all(is.finite(fixef(fit3))),
          is.finite(robustlmm:::.sigma(fit3)))

## calcTau.nondiag's own fallback (the former placeholder stop), hit
## directly as a direct caller would:
Tb3 <- withCallingHandlers(
    robustlmm:::calcTau.nondiag(fit3,
                                skbs = robustlmm:::.S(fit3),
                                kappas = robustlmm:::.kappa_b(fit3),
                                max.iter = 50L),
    warning = function(w) {
        stopifnot(grepl("falling back to DASvar", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
stopifnot(all(is.finite(as.matrix(Tb3))))

IF3 <- withCallingHandlers(
    implicitIF_full(fit3),
    warning = function(w) invokeRestart("muffleWarning"))
stopifnot(all(is.finite(IF3$IF_beta)),
          all(is.finite(IF3$IF_sigma)),
          all(is.finite(IF3$IF_theta)))

cat("dastau-fallback: all tests passed\n")
