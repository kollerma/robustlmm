## WS7 regression tests: prior weights and offsets in the influence /
## sandwich machinery.
##
## Convention (matches the estimator, see std.e() in helpers.R and
## fitEffects()): with prior weights w (lme4 convention,
## Var(e_i) = sigma^2 / w_i), the psi functions are evaluated at the
## prior-weight-whitened residuals U_e r / sigma, U_e = diag(sqrt(w)),
## and the whitened design matrices are U_e X, U_e Z U_b. Offsets are
## subtracted from the response before whitening.
##
## Before WS7 the influence machinery divided by U_e instead of
## multiplying (silently wrong for any non-unit weights; exact for
## w == 1), and offsets were dropped from mu entirely.
##
## Tests (all closed-form, no refitting, except test 5):
## 1. Classical closed form (cPsi): the partial local-shift IF equals
##    the GLS sensitivity (X' Omega^-1 X)^-1 X' Omega^-1 with
##    Omega = Z U_b U_b' Z' + W^-1 at the fitted theta.
## 2. caseweightIF (cPsi): column i equals the local-shift column i
##    times the raw residual (y - mu)_i.
## 3. .vcov_model_delta (cPsi) equals the model-based vcov().
## 4. Strict offset equality with lmer: offset of the same order as
##    the response (the old tests/offset.R offset is ~0.1% of the
##    response scale, too small to detect a dropped offset at its
##    tolerances). Plus offset and weights combined.
## 5. Full IF (smoothPsi, weighted fit): finite-difference validation
##    of IF_beta, IF_sigma and IF_theta by refitting at perturbed y.

suppressMessages(require(robustlmm))

## ----------------------------------------------------------------
## Weighted classical fit (cPsi == classical scoring equations)
## ----------------------------------------------------------------
set.seed(42)
w <- runif(nrow(Dyestuff), 0.5, 3)
fitW <- rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = w,
              rho.e = cPsi, rho.b = cPsi, init = lmerNoFit,
              method = "DASvar")

parts <- robustlmm:::implicitIF(fitW)
pp    <- fitW@pp
X     <- as.matrix(pp$X)
Z     <- as.matrix(t(pp$Zt))
U_b   <- as.matrix(pp$U_b)

## 1. Local-shift IF == GLS sensitivity, exact at the fitted theta.
Omega <- Z %*% U_b %*% t(U_b) %*% t(Z) + diag(1 / w)
M_gls <- solve(t(X) %*% solve(Omega, X)) %*% t(X) %*% solve(Omega)
relErr1 <- max(abs(parts$IF_beta - M_gls)) / max(abs(M_gls))
cat("weighted cPsi: max rel err IF_beta vs GLS sensitivity:",
    format(relErr1), "\n")
stopifnot(relErr1 < 1e-7)

## 2. Case-weight IF column i == local-shift column i * raw residual.
cw <- caseweightIF(fitW)
resid_raw <- Dyestuff$Yield - fitted(fitW)
relErr2 <- max(abs(cw$IF_beta - sweep(parts$IF_beta, 2, resid_raw, "*"))) /
    max(abs(cw$IF_beta))
cat("weighted cPsi: max rel err caseweightIF identity:",
    format(relErr2), "\n")
stopifnot(relErr2 < 1e-7)

## 3. Model-delta vcov == model-based vcov for the classical fit.
vd <- robustlmm:::.vcov_model_delta(fitW, parts)
vc <- as.matrix(vcov(fitW))
relErr3 <- max(abs(vd - vc)) / max(abs(vc))
cat("weighted cPsi: max rel err model-delta vs vcov:",
    format(relErr3), "\n")
stopifnot(relErr3 < 1e-4)

## ----------------------------------------------------------------
## 4. Offsets: strict equality with lmer at a realistic offset scale
## ----------------------------------------------------------------
## The offset is added to the response too, so the offset-corrected
## model is the standard, interior Dyestuff fit (adding an offset
## column alone would inject N(0, sd(Yield)^2) noise into Yield - off
## and push theta-hat to the boundary -- singular fits on both sides,
## where ranef comparisons are meaningless). The offset is still of
## the order of the response, so a dropped offset shifts the intercept
## by ~ sd(Yield)/sqrt(n) and every fitted value by ~ sd(Yield):
## detectable at the tolerances below, unlike in tests/offset.R.
set.seed(7)
DyeOff <- within(Dyestuff, {
    off   <- rnorm(length(Yield), sd = sd(Yield))
    Yield <- Yield + off
})

cO <- lmer(Yield ~ 1 + (1 | Batch) + offset(off), DyeOff)
rO <- rlmer(Yield ~ 1 + (1 | Batch) + offset(off), DyeOff,
            rho.e = cPsi, rho.b = cPsi, init = lmerNoFit)
stopifnot(
    all.equal(fixef(cO),  fixef(rO),  tolerance = 1e-5,
              check.attributes = FALSE),
    all.equal(ranef(cO),  ranef(rO),  tolerance = 1e-4,
              check.attributes = FALSE),
    all.equal(fitted(cO), fitted(rO), tolerance = 1e-5,
              check.attributes = FALSE))

## offset + weights combined
cOW <- lmer(Yield ~ 1 + (1 | Batch) + offset(off), DyeOff, weights = w)
rOW <- rlmer(Yield ~ 1 + (1 | Batch) + offset(off), DyeOff, weights = w,
             rho.e = cPsi, rho.b = cPsi, init = lmerNoFit)
stopifnot(
    all.equal(fixef(cOW),  fixef(rOW),  tolerance = 1e-5,
              check.attributes = FALSE),
    all.equal(fitted(cOW), fitted(rOW), tolerance = 1e-5,
              check.attributes = FALSE))

## the influence machinery sees the offset-corrected residuals: on the
## classical fit the GLS-sensitivity identity must hold as in test 1
partsO <- robustlmm:::implicitIF(rOW)
ppO    <- rOW@pp
XO     <- as.matrix(ppO$X)
ZO     <- as.matrix(t(ppO$Zt))
U_bO   <- as.matrix(ppO$U_b)
OmegaO <- ZO %*% U_bO %*% t(U_bO) %*% t(ZO) + diag(1 / w)
M_glsO <- solve(t(XO) %*% solve(OmegaO, XO)) %*% t(XO) %*% solve(OmegaO)
relErr4 <- max(abs(partsO$IF_beta - M_glsO)) / max(abs(M_glsO))
cat("offset+weights cPsi: max rel err IF_beta vs GLS sensitivity:",
    format(relErr4), "\n")
stopifnot(relErr4 < 1e-7)

## ----------------------------------------------------------------
## 5. Weighted robust fit (smoothPsi): finite-difference validation
##    of the full IF (beta, sigma, theta) by refitting at perturbed y.
## ----------------------------------------------------------------
fitR <- rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = w,
              method = "DASvar")
IF <- implicitIF_full(fitR)

refitAt <- function(y) {
    d <- Dyestuff
    d$Yield <- y
    rlmer(Yield ~ 1 + (1 | Batch), d, weights = w, method = "DASvar")
}

h   <- 0.5  ## Yield has sd ~ 50; relative truncation ~ (h / sigma)^2
y0  <- Dyestuff$Yield
idx <- c(1L, 18L)  ## obs 18: the well-known outlying batch
for (i in idx) {
    yp <- y0; yp[i] <- yp[i] + h
    ym <- y0; ym[i] <- ym[i] - h
    fp <- refitAt(yp)
    fm <- refitAt(ym)
    fd_beta  <- (fixef(fp) - fixef(fm)) / (2 * h)
    fd_sigma <- (robustlmm:::.sigma(fp) - robustlmm:::.sigma(fm)) / (2 * h)
    fd_theta <- (getME(fp, "theta") - getME(fm, "theta")) / (2 * h)
    eb <- abs(IF$IF_beta[, i]  - fd_beta)  / max(abs(IF$IF_beta))
    es <- abs(IF$IF_sigma[, i] - fd_sigma) / max(abs(IF$IF_sigma))
    et <- abs(IF$IF_theta[, i] - fd_theta) / max(abs(IF$IF_theta))
    cat(sprintf("weighted smoothPsi FD obs %d: beta %.2e sigma %.2e theta %.2e\n",
                i, max(eb), max(es), max(et)))
    stopifnot(max(eb) < 0.02, max(es) < 0.02, max(et) < 0.05)
}

## ----------------------------------------------------------------
## 6. WS11: eta-weighted (Mallows design weights) full IF -- finite-
##    difference validation. The analytical IF must carry eta on its
##    e-side rows (implicitIF / .scoreVec); refitting at perturbed y
##    confirms it, including at a downweighted observation.
## ----------------------------------------------------------------
set.seed(3)
ne  <- 48L
ge  <- factor(rep(seq_len(8L), each = 6L))
xe  <- rnorm(ne)
etaE <- pmin(1, 1.3 / abs(xe)); etaE[!is.finite(etaE)] <- 1
ye0 <- 1 + 0.5 * xe + rnorm(8L, 0, 1.2)[as.integer(ge)] + rnorm(ne)
de  <- data.frame(y = ye0, x = xe, g = ge)
fitE5 <- suppressWarnings(
    rlmer(y ~ x + (1 | g), de, method = "DASvar", design.weights = etaE))
IFe <- implicitIF_full(fitE5)
refitE <- function(yy) {
    dd <- de; dd$y <- yy
    suppressWarnings(rlmer(y ~ x + (1 | g), dd, method = "DASvar",
                           design.weights = etaE))
}
he <- 0.02
for (i in c(which.min(etaE), 1L)) {      # a downweighted obs + a unit-eta obs
    yp <- ye0; yp[i] <- yp[i] + he
    ym <- ye0; ym[i] <- ym[i] - he
    fp <- refitE(yp); fm <- refitE(ym)
    fb <- (fixef(fp) - fixef(fm)) / (2 * he)
    fs <- (robustlmm:::.sigma(fp) - robustlmm:::.sigma(fm)) / (2 * he)
    ft <- (getME(fp, "theta") - getME(fm, "theta")) / (2 * he)
    eb <- max(abs(IFe$IF_beta[, i]  - fb)) / max(abs(IFe$IF_beta))
    es <- max(abs(IFe$IF_sigma[, i] - fs)) / max(abs(IFe$IF_sigma))
    et <- max(abs(IFe$IF_theta[, i] - ft)) / max(abs(IFe$IF_theta))
    cat(sprintf("eta smoothPsi FD obs %d (eta=%.2f): beta %.2e sigma %.2e theta %.2e\n",
                i, etaE[i], eb, es, et))
    stopifnot(eb < 0.02, es < 0.02, et < 0.05)
}

cat("influence-weights-offset: all tests passed\n")
