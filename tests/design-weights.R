## WS11 steps 1-2 regression tests: Mallows design weights eta.
##
## Scope of these steps (PLAN-WS11-mallows.md §3.1-3.2): the
## design.weights argument + storage + getME accessor, and the
## fitEffects hook, so the IRLS fixed point solves the eta-weighted
## beta/u estimating equations
##     (U_e X)' H psi_e(r~)                     = 0
##     sigma (U_e Z U_b)' H psi_e(r~) - Lambda_b w_b u = 0
## with H = diag(eta). (The sigma/theta DAS equations, the model vcov,
## the Satterthwaite df and the influence diagnostics also carry eta.)
## The feature has graduated from experimental: it no longer emits an
## experimental warning, and active weights require a single grouping
## factor (it errors otherwise).
##
## Tests:
## 1. design.weights = NULL and eta = rep(1, n) give identical fits
##    (exactness regression, PLAN §5.1).
## 2. smoothPsi + prior weights + eta: the eta-weighted beta/u normal
##    equations hold at the converged fit (machine precision); no
##    experimental warning is emitted.
## 3. cPsi + prior weights + eta: closed-form check -- beta and u
##    equal the solution of the eta-weighted Henderson system at the
##    fitted (sigma, theta) (the "eta-weighted GLS" pin, PLAN §5.2).
## 4. design.weights = "mcd" downweights a planted leverage point and
##    leaves the rest at eta = 1.
## 5. "mcd" on an intercept-only design is an exact no-op (eta == 1,
##    fit identical to NULL, no warning).
## 6. Input validation errors.
## 7. Active design weights on more than one grouping factor error
##    loudly (single-grouping-factor scope, C3); a no-op request on >1
##    factor is still allowed.

suppressMessages(require(robustlmm))

n <- nrow(Dyestuff)
set.seed(42)
w <- runif(n, 0.5, 3)

## ----------------------------------------------------------------
## 1. eta = 1 identical to NULL
## ----------------------------------------------------------------
f0 <- rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = w,
            method = "DASvar", init = lmerNoFit)
f1 <- rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = w,
            method = "DASvar", init = lmerNoFit,
            design.weights = rep(1, n))
stopifnot(identical(fixef(f0), fixef(f1)),
          identical(getME(f0, "theta"), getME(f1, "theta")),
          identical(sigma(f0), sigma(f1)),
          identical(getME(f1, "design.weights"), rep(1, n)),
          identical(getME(f0, "design.weights"), rep(1, n)))
cat("test 1 (eta = 1 == NULL): ok\n")

## ----------------------------------------------------------------
## 2. eta-weighted normal equations at a smoothPsi fit
## ----------------------------------------------------------------
set.seed(7)
eta <- runif(n, 0.4, 1)
warned <- character(0)
fE <- withCallingHandlers(
    rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = w,
          method = "DASvar", init = lmerNoFit, design.weights = eta),
    warning = function(wc) {
        warned <<- c(warned, conditionMessage(wc))
        invokeRestart("muffleWarning")
    })
## the graduated feature no longer emits an experimental warning
stopifnot(!any(grepl("experimental", warned)),
          identical(getME(fE, "design.weights"), eta))

sig   <- sigma(fE)
r_std <- fE@resp$wtres / sig
psi_e <- fE@rho.e@psi(r_std)
U_eX  <- as.matrix(robustlmm:::.U_e(fE) %*% fE@pp$X)
ZL    <- as.matrix(robustlmm:::.U_eZU_b(fE))
u     <- fE@pp$b.s
w_b   <- getME(fE, "w_b_vector")
Lam_b <- diag(as.matrix(robustlmm:::.Lambda_b(fE)))

F_beta <- drop(t(U_eX) %*% (eta * psi_e))
sc_beta <- apply(abs(U_eX * (eta * psi_e)), 2, sum)
relB <- max(abs(F_beta) / sc_beta)
F_u    <- drop(sig * t(ZL) %*% (eta * psi_e)) - Lam_b * w_b * u
sc_u   <- pmax(abs(sig * t(ZL) %*% (eta * psi_e)), abs(Lam_b * w_b * u))
relU <- max(abs(F_u) / pmax(sc_u, 1e-10))
cat("eta-weighted score at fit: beta rel", format(max(relB)),
    " u rel", format(max(relU)), "\n")
## u threshold calibrated 1e-6 -> 5e-4 (2026-06-12): the outer
## alternation stops on parameter change (rel.tol 1e-8), leaving the
## u-score residual at ~1e-5 relative; the eta-free control fit shows
## the same order (6.9e-6), so this is convergence noise, not an
## eta-convention error (which would be O(eta spread), i.e. >= 1e-1).
stopifnot(relB < 1e-6, relU < 5e-4)
cat("test 2 (eta-weighted normal equations): ok\n")

## ----------------------------------------------------------------
## 3. cPsi closed form: eta-weighted Henderson system
## ----------------------------------------------------------------
fC <- suppressWarnings(
    rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = w,
          rho.e = cPsi, rho.b = cPsi, method = "DASvar",
          init = lmerNoFit, design.weights = eta))
U_eXc <- as.matrix(robustlmm:::.U_e(fC) %*% fC@pp$X)
ZLc   <- as.matrix(robustlmm:::.U_eZU_b(fC))
yoff  <- drop(robustlmm:::.U_e(fC) %*%
              (fC@resp$y - fC@resp$offset))
Lamc  <- diag(as.matrix(robustlmm:::.Lambda_b(fC)))
H     <- diag(eta)
p     <- fC@pp$p; q <- fC@pp$q
A11 <- t(U_eXc) %*% H %*% U_eXc
A12 <- t(U_eXc) %*% H %*% ZLc
A22 <- t(ZLc) %*% H %*% ZLc + diag(Lamc, q)
Amat <- rbind(cbind(A11, A12), cbind(t(A12), A22))
rhs  <- c(t(U_eXc) %*% (eta * yoff), t(ZLc) %*% (eta * yoff))
sol  <- drop(solve(Amat, rhs))
relC <- max(abs(sol - c(fixef(fC), fC@pp$b.s))) / max(abs(sol))
cat("cPsi eta-weighted Henderson closed form: max rel err",
    format(relC), "\n")
stopifnot(relC < 1e-7)
cat("test 3 (cPsi eta-weighted closed form): ok\n")

## ----------------------------------------------------------------
## 4. "mcd" downweights a leverage point
## ----------------------------------------------------------------
set.seed(11)
nL <- 60L
gL <- factor(rep(1:12, each = 5L))
xL <- rnorm(nL)
xL[1] <- 8                      # leverage point
yL <- 1 + 0.5 * xL + rnorm(12, sd = 0.5)[as.integer(gL)] + rnorm(nL)
dL <- data.frame(yL, xL, gL)
fM <- suppressWarnings(
    rlmer(yL ~ xL + (1 | gL), dL, method = "DASvar",
          init = lmerNoFit, design.weights = "mcd"))
etaM <- getME(fM, "design.weights")
stopifnot(all(etaM > 0), all(etaM <= 1),
          etaM[1] < 0.5,            # the leverage point is downweighted
          mean(etaM == 1) > 0.8)    # the bulk is untouched
cat("mcd eta at leverage point:", format(etaM[1]),
    "; share at 1:", format(mean(etaM == 1)), "\n")
cat("test 4 (mcd downweights leverage): ok\n")

## ----------------------------------------------------------------
## 5. "mcd" on intercept-only design: exact no-op, no warning
## ----------------------------------------------------------------
warned5 <- character(0)
fNoop <- withCallingHandlers(
    rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = w,
          method = "DASvar", init = lmerNoFit,
          design.weights = "mcd"),
    warning = function(wc) {
        warned5 <<- c(warned5, conditionMessage(wc))
        invokeRestart("muffleWarning")
    })
stopifnot(!any(grepl("experimental", warned5)),
          identical(getME(fNoop, "design.weights"), rep(1, n)),
          identical(fixef(f0), fixef(fNoop)),
          identical(getME(f0, "theta"), getME(fNoop, "theta")))
cat("test 5 (mcd no-op on intercept-only design): ok\n")

## ----------------------------------------------------------------
## 6. input validation
## ----------------------------------------------------------------
stopifnot(
    inherits(tryCatch(
        rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, method = "DASvar",
              init = lmerNoFit, design.weights = rep(1, 5)),
        error = function(e) e), "error"),
    inherits(tryCatch(
        rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, method = "DASvar",
              init = lmerNoFit, design.weights = rep(1.5, n)),
        error = function(e) e), "error"),
    inherits(tryCatch(
        rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, method = "DASvar",
              init = lmerNoFit, design.weights = rep(0, n)),
        error = function(e) e), "error"))
cat("test 6 (input validation): ok\n")

## ----------------------------------------------------------------
## 7. active weights require a single grouping factor (C3)
## ----------------------------------------------------------------
set.seed(13)
dX    <- expand.grid(g1 = factor(1:8), g2 = factor(1:6))
dX$x  <- rnorm(nrow(dX)); dX$x[1] <- 8          # leverage point -> eta < 1
dX$y  <- 1 + 0.5 * dX$x + rnorm(8)[dX$g1] +
             rnorm(6)[dX$g2] + rnorm(nrow(dX))
e7 <- tryCatch(
    rlmer(y ~ x + (1 | g1) + (1 | g2), dX, method = "DASvar",
          init = lmerNoFit, design.weights = "mcd"),
    error = function(e) e)
stopifnot(inherits(e7, "error"),
          grepl("single grouping factor", conditionMessage(e7),
                fixed = TRUE))
## ... but a no-op request (eta == 1 throughout) on > 1 factor is allowed:
## an intercept-only fixed-effects design has no continuous covariate, so
## "mcd" reduces to eta == 1 and never reaches the guard.
fOk <- suppressWarnings(
    rlmer(y ~ 1 + (1 | g1) + (1 | g2), dX, method = "DASvar",
          init = lmerNoFit, design.weights = "mcd"))
stopifnot(all(getME(fOk, "design.weights") == 1))
cat("test 7 (single-grouping-factor enforcement): ok\n")

cat("design-weights: all tests passed\n")
