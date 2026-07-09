## Disable test
quit()

## Tests for the EXPERIMENTAL robust variance-component score test,
## anova(fit0, fit1, test = "score"). The test is one-sided, for a
## single added independent scalar variance component, computed from
## the robust null fit only and calibrated by a score-only parametric
## bootstrap (null refits only). See ?anova.rlmerMod.
##
## Checks:
## 1. Scope errors (fail loudly): two grouping factors; correlated-slope
##    alternative; non-nested RE pair; fixed-effects-only pair.
## 2. Determinism: same seed => identical p-value and bootstrap sample.
## 3. Cross-validation: the package statistic equals an in-test
##    reimplementation of the validation study's analytic exchangeable
##    formula to 1e-8 (the package uses a per-cluster eigendecomposition
##    of V_j; in the balanced exchangeable case both are the same
##    symmetric square root).
## 4. Null behaviour smoke: on clean H0 data p is in (0, 1], no NAs,
##    >= 10 effective replicates; s_j attached and named per cluster.
## 5. Power smoke: strong random slope (tau = 1.5, J = 30) rejects at
##    alpha = 0.10 (seeded).
## 6. Contaminated smoke: (a) a planted mean-shifted cluster: runs with
##    no warnings besides the small-nsim one, valid p; (b) a cluster
##    contaminated in the tested (slope) direction: its s_j is the top
##    per-cluster contribution (the diagnosability product).

suppressMessages(require(robustlmm))

mkData <- function(J = 12L, ni = 6L, tau_slope = 0, seed = 1L) {
    set.seed(seed)
    g <- factor(rep(seq_len(J), each = ni))
    x <- rep(as.numeric(scale(seq_len(ni))), J)
    gi <- as.integer(g)
    s <- if (tau_slope > 0) rnorm(J, 0, tau_slope) else rep(0, J)
    y <- 1 + 0.5 * x + rnorm(J)[gi] + s[gi] * x + rnorm(J * ni)
    data.frame(y, x, g)
}

## Fits are built at the top level so that the bootstrap's
## update(fit, data = ...) re-evaluates the call cleanly -- the same
## convention as tests/anova-robust-null.R.

## ---- 1. scope errors --------------------------------------------------
## (a) two grouping factors
set.seed(2L)
dcr <- expand.grid(g = factor(1:8), h = factor(letters[1:4]),
                   rep = 1:2)
dcr$x <- rep(as.numeric(scale(1:8)), 8)
dcr$y <- 1 + 0.5 * dcr$x + rnorm(8)[as.integer(dcr$g)] +
    0.7 * rnorm(4)[as.integer(dcr$h)] + rnorm(nrow(dcr))
f0cr <- suppressWarnings(rlmer(y ~ x + (1 | g) + (1 | h), dcr,
                               method = "DASvar"))
f1cr <- suppressWarnings(rlmer(y ~ x + (1 | g) + (1 | h) + (0 + x | g),
                               dcr, method = "DASvar"))
e1a <- tryCatch(anova(f0cr, f1cr, test = "score"), error = identity)
stopifnot(inherits(e1a, "error"),
          grepl("single grouping factor", conditionMessage(e1a)),
          grepl("test = \"boot\"", conditionMessage(e1a)))
cat("test 1a (two grouping factors: stop): ok\n")

## (b) correlated-slope alternative (theta jumps by 2, not 1)
d1  <- mkData(seed = 1L)
f0  <- suppressWarnings(rlmer(y ~ x + (1 | g), d1, method = "DASvar"))
f1  <- suppressWarnings(rlmer(y ~ x + (1 | g) + (0 + x | g), d1,
                              method = "DASvar"))
f1c <- suppressWarnings(rlmer(y ~ x + (1 + x | g), d1,
                              method = "DASvar"))
e1b <- tryCatch(anova(f0, f1c, test = "score"), error = identity)
stopifnot(inherits(e1b, "error"),
          grepl("exactly one variance parameter", conditionMessage(e1b)),
          grepl("test = \"boot\"", conditionMessage(e1b)))
cat("test 1b (correlated-slope alternative: stop): ok\n")

## (c) non-nested RE pair: (1|g) + (0+x|g) vs (1+x|g) -- theta lengths
## differ by 1 but no single scalar column is added
e1c <- tryCatch(anova(f1, f1c, test = "score"), error = identity)
stopifnot(inherits(e1c, "error"),
          grepl("could not identify a single added scalar",
                conditionMessage(e1c)))
cat("test 1c (non-nested RE pair: stop): ok\n")

## (d) fixed-effects-only pair
f0f <- suppressWarnings(rlmer(y ~ 1 + (1 | g), d1, method = "DASvar"))
e1d <- tryCatch(anova(f0f, f0, test = "score"), error = identity)
stopifnot(inherits(e1d, "error"),
          grepl("do not differ in their\\s+random-effects structure",
                conditionMessage(e1d)))
cat("test 1d (fixed-effects-only pair: stop): ok\n")

## ---- 2. determinism under a fixed seed --------------------------------
a2a <- suppressWarnings(anova(f0, f1, test = "score", nsim = 19L,
                              seed = 42L))
a2b <- suppressWarnings(anova(f0, f1, test = "score", nsim = 19L,
                              seed = 42L))
stopifnot(identical(a2a[["Pr(>=Score)"]][2L], a2b[["Pr(>=Score)"]][2L]),
          identical(attr(a2a, "boot")$S_boot, attr(a2b, "boot")$S_boot))
cat("test 2 (determinism under fixed seed): ok\n")

## ---- 3. cross-validation vs the study formula -------------------------
## In-test reimplementation of the validation study's score_stat
## (analytic symmetric square root for the balanced exchangeable case;
## robust-score-test.R in the study record).
study_stat <- function(fit, xv, gv) {
    sig  <- sigma(fit)
    th   <- as.numeric(robustlmm::getME(fit, "theta"))
    beta <- robustlmm:::.fixef(fit)
    yv   <- fit@resp$y
    r    <- as.numeric(yv - fit@pp$X %*% beta)
    psi  <- fit@rho.e@psi
    kap1 <- fit@rho.e@Epsi2()
    tau2 <- (sig * th)^2
    idx  <- split(seq_along(r), gv)
    sj <- vapply(idx, function(ii) {
        ni <- length(ii)
        ri <- r[ii]; xi <- xv[ii]
        a1 <- 1 / sig
        a2 <- 1 / sqrt(sig^2 + ni * tau2)
        mr <- mean(ri); mx <- mean(xi)
        rt <- a1 * (ri - mr) + a2 * mr
        v  <- a1 * (xi - mx) + a2 * mx
        tj <- sum(v * psi(rt))
        tj^2 - kap1 * sum(v^2)
    }, numeric(1))
    den <- sqrt(sum((sj - mean(sj))^2))
    list(S = sum(sj) / den, sj = sj)
}
st  <- study_stat(f0, d1$x, d1$g)
sc  <- robustlmm:::.anova_score_scope(f0, f1)
pk  <- robustlmm:::.anova_score_stat(f0, sc$z_add)
stopifnot(max(abs(sc$z_add - d1$x)) < 1e-12,
          abs(pk$S - st$S) < 1e-8,
          max(abs(pk$s_j - st$sj)) < 1e-8)
## ... and the observed statistic reported by anova() is the same
stopifnot(abs(a2a$Score[2L] - st$S) < 1e-8)
cat("test 3 (package statistic == study formula to 1e-8): ok\n")

## ---- 4. null behaviour smoke -------------------------------------------
pv4 <- a2a[["Pr(>=Score)"]][2L]
b4  <- attr(a2a, "boot")
stopifnot(is.finite(pv4), pv4 > 0, pv4 <= 1,
          !anyNA(b4$S_boot), length(b4$S_boot) >= 10L,
          b4$nsim == 19L,
          identical(names(b4$s_j), levels(d1$g)),
          all(is.finite(b4$s_j)),
          isTRUE(all.equal(b4$S_obs, a2a$Score[2L])))
## table shape consistent with the boot path's (npar, Df, stat, p)
stopifnot(identical(rownames(a2a), c("fit0", "fit1")),
          identical(colnames(a2a),
                    c("npar", "Df", "Score", "Pr(>=Score)")),
          any(grepl("EXPERIMENTAL", attr(a2a, "heading"))))
cat("test 4 (clean-null smoke: valid p, s_j named per cluster): ok\n")

## ---- 5. power smoke -----------------------------------------------------
d5 <- mkData(J = 30L, tau_slope = 1.5, seed = 5L)
f05 <- suppressWarnings(rlmer(y ~ x + (1 | g), d5, method = "DASvar"))
f15 <- suppressWarnings(rlmer(y ~ x + (1 | g) + (0 + x | g), d5,
                              method = "DASvar"))
a5 <- suppressWarnings(anova(f05, f15, test = "score", nsim = 39L,
                             seed = 6L))
stopifnot(a5[["Pr(>=Score)"]][2L] <= 0.10)
cat("test 5 (power smoke, tau = 1.5, J = 30: rejects at 0.10): ok\n")

## ---- 6. contaminated smoke ----------------------------------------------
## (a) planted mean-shifted cluster: the score path must NOT emit the
## boot path's downweighted-group warning (no null cleaning by design);
## the only tolerated warning is the small-nsim one.
d6 <- mkData(seed = 7L)
d6$y[d6$g == 4] <- d6$y[d6$g == 4] + 8
f06 <- suppressWarnings(rlmer(y ~ x + (1 | g), d6, method = "DASvar"))
f16 <- suppressWarnings(rlmer(y ~ x + (1 | g) + (0 + x | g), d6,
                              method = "DASvar"))
warns6 <- character(0)
a6 <- withCallingHandlers(
    anova(f06, f16, test = "score", nsim = 19L, seed = 8L),
    warning = function(w) {
        warns6 <<- c(warns6, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
stopifnot(all(grepl("effective replicates", warns6)),
          !any(grepl("downweight", warns6)))
pv6 <- a6[["Pr(>=Score)"]][2L]
stopifnot(is.finite(pv6), pv6 > 0, pv6 <= 1)
cat("test 6a (mean-shifted cluster: no warning-spam, valid p): ok\n")

## (b) cluster contaminated in the tested (slope) direction: its
## bounded contribution s_j saturates and tops the attached vector --
## the diagnosability product documented in ?anova.rlmerMod.
d6b <- mkData(seed = 9L)
d6b$y[d6b$g == 5] <- d6b$y[d6b$g == 5] + 4 * d6b$x[d6b$g == 5]
f06b <- suppressWarnings(rlmer(y ~ x + (1 | g), d6b, method = "DASvar"))
f16b <- suppressWarnings(rlmer(y ~ x + (1 | g) + (0 + x | g), d6b,
                               method = "DASvar"))
a6b <- suppressWarnings(anova(f06b, f16b, test = "score", nsim = 19L,
                              seed = 10L))
sj6 <- attr(a6b, "boot")$s_j
stopifnot(identical(names(which.max(sj6)), "5"))
cat("test 6b (slope-direction cluster tops s_j): ok\n")

cat("anova-score: all tests passed\n")
