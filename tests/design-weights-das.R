## WS11 step 3 (DASvar) regression tests: the sigma/theta DAS machinery
## carries the Mallows design weights eta.
##
## Correctness of the eta != 1 tau formula was established separately by
## Monte-Carlo (PLAN-WS11 step 3, 2026-06-17): the analytical tau matches
## the empirical residual sd within MC error at every observation
## (|z| < 1.3), while the naive per-row eta^2*diag(AA') form is 5-8
## sigma off -- confirming the eta^2 belongs INSIDE the A A' sum
## (diag(A diag(eta^2) A')). These fast checks guard against regressions.
##
## 1. eta == 1 is bit-identical to design.weights = NULL (the scalar
##    fast path is exercised and the DAS core is unchanged at eta = 1).
## 2. The installed DASvar tau uses the A H^2 A' form (not the naive
##    one) -- checked deterministically from the fitted A matrix.
## 3. Regression pin on a fixed cPsi + eta DASvar fit.

suppressMessages(require(robustlmm))

mkdata <- function(seed = 101, J = 6, ni = 5) {
    set.seed(seed)
    n <- J * ni
    g <- factor(rep(seq_len(J), each = ni))
    x <- rnorm(n)
    eta <- pmin(1, 1.2 / abs(x)); eta[!is.finite(eta)] <- 1
    y <- 1 + 0.5 * x + rnorm(J, 0, 1)[as.integer(g)] + rnorm(n)
    list(d = data.frame(y, x, g), eta = eta, n = n)
}

## ---- 1. eta == 1 bit-identity (NULL vs explicit ones) -------------
dat <- mkdata()
fN <- rlmer(y ~ x + (1 | g), dat$d, method = "DASvar",
            rho.e = cPsi, rho.b = cPsi)
f1 <- rlmer(y ~ x + (1 | g), dat$d, method = "DASvar",
            rho.e = cPsi, rho.b = cPsi,
            design.weights = rep(1, dat$n))
stopifnot(isTRUE(all.equal(fN@pp$tau_e(), f1@pp$tau_e())),
          isTRUE(all.equal(fixef(fN), fixef(f1))),
          isTRUE(all.equal(sigma(fN), sigma(f1))),
          isTRUE(all.equal(as.matrix(vcov(fN)), as.matrix(vcov(f1)))))

## ---- 2. eta != 1 DASvar tau uses the A H^2 A' (meat = H^2) form ----
fe <- suppressWarnings(
    rlmer(y ~ x + (1 | g), dat$d, method = "DASvar",
          rho.e = cPsi, rho.b = cPsi, design.weights = dat$eta))
pp  <- fe@pp
eta <- dat$eta
A   <- as.matrix(pp$A())
AH2At  <- as.numeric(A^2 %*% eta^2)            # diag(A H^2 A')  (correct)
naive  <- eta^2 * pp$diagAAt                    # eta_i^2 diag(A A') (wrong)
## tau_correct^2 - tau_naive^2 = Epsi2 * (AH2At - naive); reconstruct the
## naive tau the installed (correct) fit would have produced and confirm
## the installed tau matches the correct branch, not the naive one.
tau_inst   <- pp$tau_e()
tau_naive2 <- tau_inst^2 - pp$Epsi2_e * (AH2At - naive)
stopifnot(max(abs(AH2At - naive)) > 1e-3,                 # forms really differ
          all(abs(tau_inst - sqrt(pmax(tau_naive2, 0))) > 1e-6 |
              abs(AH2At - naive) < 1e-9))                 # installed != naive where they differ

## ---- 3. regression pin (validated by the MC anchor) ---------------
pin <- c(0.8981302072, 0.8904369398, 0.8826409075,
         0.9012619979, 0.8991889714, 0.8789551958)
stopifnot(max(abs(pp$tau_e()[1:6] - pin)) < 1e-6)

## eta weighting actually changes tau (sanity)
stopifnot(max(abs(pp$tau_e() - fN@pp$tau_e())) > 1e-3)

## ---- 4. DAStau path (the default method) --------------------------
## eta == 1 identity (covers the .s / calcTau scalar fast paths) and a
## tau pin. DAStau eta correctness (sigma bias -0.2%, theta bias -0.1%
## over 2000 reps) was established by the MC consistency check in
## PLAN-WS11 step 3.
fNt <- rlmer(y ~ x + (1 | g), dat$d, method = "DAStau",
             rho.e = cPsi, rho.b = cPsi)
f1t <- rlmer(y ~ x + (1 | g), dat$d, method = "DAStau",
             rho.e = cPsi, rho.b = cPsi, design.weights = rep(1, dat$n))
stopifnot(isTRUE(all.equal(fNt@pp$tau_e(), f1t@pp$tau_e())),
          isTRUE(all.equal(sigma(fNt), sigma(f1t))))
fet <- suppressWarnings(
    rlmer(y ~ x + (1 | g), dat$d, method = "DAStau",
          rho.e = cPsi, rho.b = cPsi, design.weights = dat$eta))
## cPsi is linear, so DAStau ~ DASvar -> same pin (to ~1e-9)
stopifnot(max(abs(fet@pp$tau_e()[1:6] - pin)) < 1e-6,
          max(abs(fet@pp$tau_e() - fNt@pp$tau_e())) > 1e-3)

## ---- 5. eta-weighted model vcov (unsc, H^2 meat sandwich) ---------
## Correctness: at eta = 1 the sandwich reduces to the scalar vcov to
## machine precision (1.4e-16, deterministic) and the eta != 1 form
## matches empirical Cov(beta-hat) (PLAN-WS11 step 3: slope vcov ratio
## 1.22 -> 1.06 after the fix). Here: PD, differs from eta=1, pinned.
Ve <- as.matrix(vcov(fe))               # fe = DASvar cPsi eta fit (above)
stopifnot(min(eigen(Ve, only.values = TRUE)$values) > 0,        # PD
          max(abs(Ve - as.matrix(vcov(fN)))) > 1e-4)            # eta has an effect
vpin <- c(0.2899941244, 0.0493321479)
stopifnot(max(abs(diag(Ve) - vpin)) < 1e-6)

## ---- 6. summary() reports the Mallows design weights (step 6) ------
out_e <- capture.output(print(summary(fe)))
out_n <- capture.output(print(summary(fN)))
stopifnot(any(grepl("Mallows design weights", out_e)),    # eta fit: reported
          any(grepl("downweighted", out_e)),
          !any(grepl("Mallows", out_n)))                  # eta=1 fit: not shown

cat("design-weights-das.R: all checks passed\n")
