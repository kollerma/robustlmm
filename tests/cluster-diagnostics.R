## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## WS15 regression tests: cluster-level Cook's distance and robust
## leverage (hatvalues).
## 1. Per-observation cooks.distance unchanged (backward compatible).
## 2. Cluster Cook's distance flags a planted whole-group outlier.
## 3. groups = "<name>" equals groups = TRUE (the sole grouping factor).
## 4. hatvalues: per-obs leverage; sum is preserved when aggregated per
##    cluster; reduces to the classical lme4 leverage at cPsi.
## 5. Crossed designs error for the cluster-level diagnostics.

suppressMessages({require(robustlmm); require(lme4)})

set.seed(1)
J <- 20L; m <- 6L; n <- J * m
g <- factor(rep(seq_len(J), each = m)); x <- rnorm(n)
y <- 1 + 0.5 * x + rnorm(J, 0, 1)[as.integer(g)] + rnorm(n)
y[g == 7] <- y[g == 7] + 6          # whole-cluster outlier
d <- data.frame(y, x, g)
fit <- suppressWarnings(rlmer(y ~ x + (1 | g), d, method = "DAStau"))

## ---- 1. per-observation cooks.distance unchanged ------------------
cdo <- cooks.distance(fit)
stopifnot(length(cdo) == n, identical(names(cdo)[1], "obs1"),
          all(is.finite(cdo)), all(cdo >= 0))

## ---- 2. cluster Cook's distance flags the planted outlier ----------
cdc <- cooks.distance(fit, groups = TRUE)
stopifnot(length(cdc) == J,
          identical(names(which.max(cdc)), "7"),
          all(is.finite(cdc)))

## ---- 3. groups by name == groups = TRUE ---------------------------
cdn <- cooks.distance(fit, groups = "g")
stopifnot(isTRUE(all.equal(cdc, cdn)))

## ---- 4. hatvalues -------------------------------------------------
hv <- hatvalues(fit)
stopifnot(length(hv) == n, all(hv > 0), all(hv < 1 + 1e-8))
hvc <- hatvalues(fit, groups = TRUE)
stopifnot(length(hvc) == J,
          isTRUE(all.equal(sum(hv), sum(hvc))))      # leverage conserved
## classical limit: at rho = cPsi the robust leverage is the same
## functional form as the classical lme4 leverage. rlmer and lme4 use
## different optimisers, so their theta differ by ~1e-2 and the leverages
## are not identical to FD accuracy; but the shapes coincide (correlation
## ~1) and the effective df (sum) agree closely. Use clean data so cPsi
## and lmer are near the same fit.
set.seed(99)
yc <- 1 + 0.5 * x + rnorm(J, 0, 1.1)[as.integer(g)] + rnorm(n)
dc2 <- data.frame(y = yc, x = x, g = g)
fc  <- suppressWarnings(rlmer(y ~ x + (1 | g), dc2, rho.e = cPsi,
                             rho.b = cPsi, method = "DAStau"))
lmf <- lmer(y ~ x + (1 | g), dc2, REML = FALSE)
stopifnot(cor(hatvalues(fc), hatvalues(lmf)) > 0.999,
          abs(sum(hatvalues(fc)) - sum(hatvalues(lmf))) /
              sum(hatvalues(lmf)) < 0.02)
cat("hatvalues vs lme4 (cPsi): ok (cor ~1, sum within 2%)\n")

## ---- 5. crossed designs error for cluster-level diagnostics -------
dc <- expand.grid(s = factor(1:8), i = factor(1:6)); dc$y <- rnorm(48)
fx <- suppressWarnings(rlmer(y ~ 1 + (1 | s) + (1 | i), dc,
                             method = "DASvar"))
stopifnot(inherits(tryCatch(cooks.distance(fx, groups = TRUE),
                            error = function(e) e), "error"),
          inherits(tryCatch(hatvalues(fx, groups = TRUE),
                            error = function(e) e), "error"))
## but per-observation hatvalues still work on the crossed fit
stopifnot(length(hatvalues(fx)) == nrow(dc))

cat("cluster-diagnostics.R: all checks passed\n")
