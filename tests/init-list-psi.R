## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## rlmer(init = list(...)) re-solves (beta, u) at the supplied sigma and
## theta with the fit's own psi-functions. It used to run that solve on the
## object coerced from the lmer start, whose rho.e is the classical cPsi, so
## it returned the GLS fit and discarded the supplied fixed effects (and so
## every RANSAC start's).
require(robustlmm)

set.seed(1)
J <- 20; m <- 10
g <- factor(rep(seq_len(J), each = m))
x <- rnorm(J * m)
y <- 2 + x + rnorm(J)[g] + rnorm(J * m)
out <- sample.int(J * m, 40)
y[out] <- y[out] + 8                       # 20% vertical outliers
d <- data.frame(g = g, x = x, y = y)
truth <- list(fixef = c(2, 1), sigma = 1, theta = 1)
fit_b <- function(init, ...)
    suppressWarnings(rlmer(y ~ x + (1 | g), d, method = "DASvar",
                           rho.e = bisquarePsi,
                           rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                           rho.b = chgDefaults(smoothPsi, k = 5.14),
                           rho.sigma.b = chgDefaults(smoothPsi, k = 5.14),
                           init = init, ...))

## 1. With max.iter = 0 the result is the bisquare solve at the supplied
##    sigma and theta, started from the supplied fixed effects: it keeps
##    the outliers rejected (the classical solve gave the lmer intercept).
f0 <- fit_b(truth, max.iter = 0)
lm0 <- lme4::lmer(y ~ x + (1 | g), d)
stopifnot(abs(fixef(f0)[1] - 2) < 0.4,
          fixef(lm0)[1] > 3,
          sigma(f0) == 1, getME(f0, "theta") == 1,
          all(getME(f0, "w_e")[out] < 1e-8))

## 2. The same values applied by hand to a configured fit give the same
##    solve (the start is applied once the psi-functions are in place).
fh <- f0
robustlmm:::setFixef(fh, c(2, 1)); robustlmm:::setSigma(fh, 1)
robustlmm:::setTheta(fh, 1, fit.effects = TRUE, update.sigma = FALSE)
stopifnot(max(abs(robustlmm:::.fixef(fh) - fixef(f0))) < 1e-6)

## 3. Classical psi-functions: the solve is the GLS fit at the supplied
##    sigma and theta whatever fixed effects are supplied.
fc1 <- rlmer(y ~ x + (1 | g), d, method = "DASvar", rho.e = cPsi, rho.b = cPsi,
             rho.sigma.e = cPsi, rho.sigma.b = cPsi,
             init = list(fixef = c(0, 0), sigma = 1, theta = 1), max.iter = 0)
fc2 <- rlmer(y ~ x + (1 | g), d, method = "DASvar", rho.e = cPsi, rho.b = cPsi,
             rho.sigma.e = cPsi, rho.sigma.b = cPsi,
             init = list(fixef = c(5, -3), sigma = 1, theta = 1), max.iter = 0)
stopifnot(max(abs(fixef(fc1) - fixef(fc2))) < 1e-6)

## 4. A supplied u is used as given (no re-solve).
u0 <- rep(0.1, J)
fu <- fit_b(list(fixef = c(2, 1), sigma = 1, theta = 1, u = u0), max.iter = 0)
stopifnot(max(abs(getME(fu, "u") - u0)) < 1e-12,
          max(abs(fixef(fu) - c(2, 1))) < 1e-12)

## 5. ... and the fitted values are those of the supplied values (setU()
##    alone left them at the lmer start's u)
stopifnot(max(abs(fu@resp$mu - robustlmm:::.mu(fu))) < 1e-10)

cat("init-list-psi.R: all checks passed\n")
