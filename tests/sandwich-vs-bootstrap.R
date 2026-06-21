## Disable test
quit()

## Cross-check: the cluster-robust sandwich SE for the fixed effects
## (vcov(fit, type = "sandwich")) should agree with a Monte-Carlo SE
## from parametric simulation + refit on a nested design, to within a
## factor of 2. The two estimate different objects -- the sandwich is
## the leading-order partial sandwich at fixed sigma_hat / theta_hat,
## while the MC SE additionally varies (sigma, theta) -- so the
## threshold is loose; the point is to catch order-of-magnitude wiring
## errors, not to demand numerical equality.
##
## Originally listed in the WS4 spec for tests/influence.R as the
## "sandwich vs Maritz-Jarrett" check. That sub-item became obsolete
## when the package's bootstrap was delegated to confintROB
## (confintROB::parametric.rlmerMod refits *lmer*, not rlmer, on the
## simulated data, so its bootstrap_estimates SD measures lmer's
## sampling variability under the rlmer-fitted truth -- apples to
## oranges versus a sandwich SE for rlmer). This file performs the
## proper rlmer-refit Monte-Carlo cross-check instead. ~25 s.

require(robustlmm)
suppressMessages(require(Matrix))

## Nested design (single grouping factor) where the cluster sandwich
## is exact. Sleepstudy: J = 18 subjects, p = 2 fixed effects, n = 180.
fit <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
             method = "DASvar")
nms     <- names(fixef(fit))
beta_h  <- fixef(fit)
sig_h   <- robustlmm:::.sigma(fit)
theta_h <- getME(fit, "theta")
pp      <- fit@pp
X       <- pp$X
Z       <- t(pp$Zt)
U_b     <- pp$U_b
n <- pp$n; p <- pp$p; q <- pp$q
data_   <- fit@frame
yname   <- as.character(formula(fit)[[2L]])
Xb_h    <- as.numeric(X %*% beta_h)

NREP <- 80L
beta_sim <- matrix(NA_real_, NREP, p)
n_fail   <- 0L

t0 <- proc.time()
set.seed(20260603L)
for (r in seq_len(NREP)) {
    b_sim <- sig_h * as.numeric(U_b %*% rnorm(q))
    eps   <- rnorm(n, sd = sig_h)
    y_sim <- Xb_h + as.numeric(Z %*% b_sim) + eps
    d_sim <- data_; d_sim[[yname]] <- y_sim
    f_sim <- tryCatch(
        suppressMessages(suppressWarnings(
            update(fit, data = d_sim, start = theta_h))),
        error = function(e) NULL)
    if (is.null(f_sim)) { n_fail <- n_fail + 1L; next }
    beta_sim[r, ] <- fixef(f_sim)
}
mc_time <- as.numeric((proc.time() - t0)["elapsed"])
ok <- NREP - n_fail
stopifnot(ok >= 0.8 * NREP)

se_mc   <- apply(beta_sim, 2L, sd, na.rm = TRUE)
se_sand <- sqrt(diag(as.matrix(vcov(fit, type = "sandwich"))))
ratio   <- se_sand / se_mc

cat(sprintf("sandwich-vs-bootstrap (nrep=%d, fit failures=%d, %.0fs):\n",
            NREP, n_fail, mc_time))
cat("  coef            sandwich      MC_SE      ratio\n")
for (k in seq_along(nms))
    cat(sprintf("  %-14s %9.4f  %8.4f  %7.3f\n",
                nms[k], se_sand[k], se_mc[k], ratio[k]))

stopifnot(all(is.finite(ratio)))
stopifnot(all(ratio > 0.5 & ratio < 2.0))
