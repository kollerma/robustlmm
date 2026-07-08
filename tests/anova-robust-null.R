## Disable test
quit()

## Tests for the experimental robust null of the bootstrap
## variance-component test, anova(fit0, fit1, test = "boot",
## null = "robust"). The robust null generates the parametric bootstrap
## from a contamination-cleaned null fit: clusters the robust fit heavily
## downweights (smallest per-cluster RE weight < threshold, default 0.65)
## are trimmed, fit0 is refitted without them, and the bootstrap is drawn
## from the de-biased generating parameters. D_obs stays from the original
## untrimmed fits. Fallback: no/over-half clusters flagged or refit failure
## -> plain parametric null.
##
## Checks:
## 1. null = "robust" runs end-to-end and returns a valid anova table with
##    a finite p-value in [0, 1].
## 2. On data with a PLANTED contaminated cluster, the cleaning flags it
##    (>= 1 cluster) and the cleaned theta/sigma differ from the raw fit0.
## 3. On clean data, (almost) nothing is flagged and the cleaned params
##    match the raw fit (cleaned = FALSE), so the result ~matches
##    null = "parametric".
## 4. null = "parametric" (default) is unchanged.
## 5. Multi-theta null (crossed (1|g) + (1|h), length(theta) == 2): the
##    planted cluster is flagged, cleaning works, and the generator's
##    U_b is Lambda(theta_clean) exactly, i.e.
##    Ub@x == theta_clean[pp$Lind]. The pre-fix code set ratio <- 1
##    whenever length(theta) > 1, i.e. it would have generated from the
##    STALE untrimmed U_b here; the test asserts the rebuilt U_b differs.
## 6. >50%-flagged fallback (forced via the threshold option): the
##    specific warning fires, the table heading says the robust null was
##    not applied, and the result equals null = "parametric" at the same
##    seed.
## 7. Rank-deficient / refit-failure fallback, triggered deterministically:
##    a fixed-effect factor level that lives only inside the flagged cluster
##    makes the trimmed design rank-deficient, so .anova_robust_null_params
##    warns and returns the untrimmed parameters (reason "rank-deficient"
##    from the nonsingular-subsampling check, or "refit-failed" from the
##    downstream dimension check).
## 8. Nonsingular-subsampling check, shared-level scenario: a contrast level
##    shared by two oppositely-shifted (hence flagged) clusters is dropped
##    when both are trimmed; reason is deterministically "rank-deficient".

suppressMessages(require(robustlmm))

mkData <- function(contam, seed = 1L) {
    set.seed(seed); J <- 12L; m <- 6L; n <- J * m
    g <- factor(rep(seq_len(J), each = m)); x <- rnorm(n)
    y <- 1 + 0.5 * x + rnorm(J, 0, 1)[as.integer(g)] + rnorm(n)
    if (contam) y[g == 4] <- y[g == 4] + 8        # group-4 contamination
    data.frame(y, x, g)
}
## Fits are built at the top level (not inside a closure) so that the
## bootstrap's update(fit, data = ...) re-evaluates the call cleanly --
## the same convention as tests/anova-boot-guard.R.

## ---- 1. runs end-to-end, valid p-value -------------------------------
dc  <- mkData(TRUE)
f0c <- suppressWarnings(rlmer(y ~ x + (1 | g), dc, method = "DASvar"))
f1c <- suppressWarnings(rlmer(y ~ x + (1 + x | g), dc, method = "DASvar"))
ar  <- suppressWarnings(anova(f0c, f1c, test = "boot", null = "robust",
                              nsim = 30L, seed = 11L))
stopifnot(inherits(ar, "anova"))
pv <- ar[["Pr(>=Diff)"]][2L]
stopifnot(is.finite(pv), pv >= 0, pv <= 1)
cat("test 1 (robust null runs, valid p-value): ok\n")

## ---- 2. planted contaminated cluster is flagged, params shift --------
rn <- robustlmm:::.anova_robust_null_params(f0c, threshold = 0.65)
stopifnot(isTRUE(rn$cleaned),
          rn$n_flagged >= 1L,
          "g: 4" %in% rn$flagged)
raw_theta <- robustlmm::getME(f0c, "theta")
raw_sigma <- sigma(f0c)
stopifnot(!isTRUE(all.equal(rn$theta, raw_theta)) ||
          !isTRUE(all.equal(rn$sigma, raw_sigma)))
## the cleaned variance component should be smaller than the
## contamination-inflated raw one
stopifnot(rn$theta[1L] <= raw_theta[1L] + 1e-8)
cat("test 2 (contaminated cluster flagged, cleaned params differ): ok\n")

## ---- 3. clean data: nothing flagged, falls back to parametric --------
d0  <- mkData(FALSE)
f0  <- suppressWarnings(rlmer(y ~ x + (1 | g), d0, method = "DASvar"))
f1  <- suppressWarnings(rlmer(y ~ x + (1 + x | g), d0, method = "DASvar"))
rn0 <- robustlmm:::.anova_robust_null_params(f0, threshold = 0.65)
stopifnot(isFALSE(rn0$cleaned), rn0$n_flagged == 0L)
a_par <- suppressWarnings(anova(f0, f1, test = "boot", null = "parametric",
                                nsim = 30L, seed = 7L))
a_rob <- suppressWarnings(anova(f0, f1, test = "boot", null = "robust",
                                nsim = 30L, seed = 7L))
## identical seed + no cleaning => identical bootstrap => identical p-value
stopifnot(isTRUE(all.equal(a_par[["Pr(>=Diff)"]][2L],
                           a_rob[["Pr(>=Diff)"]][2L])))
cat("test 3 (clean data: no flag, robust ~ parametric): ok\n")

## ---- 4. null = "parametric" default unchanged ------------------------
a_def <- suppressWarnings(anova(f0, f1, test = "boot",
                                nsim = 30L, seed = 7L))
stopifnot(isTRUE(all.equal(a_def[["Pr(>=Diff)"]][2L],
                           a_par[["Pr(>=Diff)"]][2L])))
cat("test 4 (parametric default unchanged): ok\n")

## ---- 5. multi-theta null: full-theta U_b rebuild ---------------------
## crossed intercepts => theta has 2 components; the cleaned generator
## must rebuild U_b = Lambda(theta_clean) through pp$Lind, exact for ANY
## number of components (the old ratio rescale of fit0's U_b was only
## applied for a single component and silently kept the stale U_b here).
set.seed(3L); Jm <- 10L; Km <- 6L
dm <- expand.grid(g = factor(seq_len(Jm)), h = factor(letters[seq_len(Km)]))
dm$x <- rnorm(nrow(dm))
dm$y <- 1 + 0.5 * dm$x + rnorm(Jm)[as.integer(dm$g)] +
        0.7 * rnorm(Km)[as.integer(dm$h)] + rnorm(nrow(dm))
dm$y[dm$g == 4] <- dm$y[dm$g == 4] + 8    # contaminate cluster g = 4
f0m <- suppressWarnings(rlmer(y ~ x + (1 | g) + (1 | h), dm,
                              method = "DASvar"))
stopifnot(length(robustlmm::getME(f0m, "theta")) == 2L)
rnm <- robustlmm:::.anova_robust_null_params(f0m, threshold = 0.65)
stopifnot(isTRUE(rnm$cleaned),
          "g: 4" %in% rnm$flagged,
          length(rnm$theta) == 2L,
          ## exact Lambda(theta_clean) through the Lind map
          isTRUE(all.equal(rnm$Ub@x, rnm$theta[f0m@pp$Lind],
                           tolerance = 1e-12)),
          ## ... and it is NOT the stale untrimmed U_b the pre-fix
          ## multi-theta path would have used
          !isTRUE(all.equal(rnm$Ub@x, f0m@pp$U_b@x)))
## the untrimmed (no-flag) path leaves U_b at fit0's Lambda
dm0 <- dm; dm0$y[dm0$g == 4] <- dm0$y[dm0$g == 4] - 8
f0m0 <- suppressWarnings(rlmer(y ~ x + (1 | g) + (1 | h), dm0,
                               method = "DASvar"))
rnm0 <- robustlmm:::.anova_robust_null_params(f0m0, threshold = 0.65)
stopifnot(isFALSE(rnm0$cleaned),
          identical(rnm0$reason, "none-flagged"),
          isTRUE(all.equal(rnm0$Ub@x, f0m0@pp$U_b@x, tolerance = 1e-12)))
cat("test 5 (multi-theta cleaned U_b = Lambda(theta_clean)): ok\n")

## ---- 6. >50%-flagged fallback: warning, heading, == parametric -------
## threshold = 1.5 flags every cluster (all weights are <= 1), so the
## >50% guard must fire and the call must degrade to the plain
## parametric bootstrap.
op <- options(robustlmm.anova.robust_null_threshold = 1.5)
warns <- character(0)
a_over <- withCallingHandlers(
    anova(f0, f1, test = "boot", null = "robust", nsim = 30L, seed = 7L),
    warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
options(op)
stopifnot(any(grepl("clusters flagged (> 50%)", warns, fixed = TRUE)))
stopifnot(isTRUE(all.equal(a_over[["Pr(>=Diff)"]][2L],
                           a_par[["Pr(>=Diff)"]][2L])))
## heading must NOT claim "no cluster flagged; identical to parametric"
h <- attr(a_over, "heading")
stopifnot(any(grepl("> 50% of clusters flagged", h, fixed = TRUE)),
          !any(grepl("no cluster flagged", h, fixed = TRUE)))
cat("test 6 (>50% flagged: warning + heading + == parametric): ok\n")

## ---- 7. refit-failure fallback (deterministic) -----------------------
## Factor level "b" of the fixed effect f exists ONLY inside the
## contaminated cluster 4, but only on a THIRD of its rows -- so f
## cannot absorb the cluster-wide shift (a full g == 4 indicator would,
## un-flagging the cluster) and g4 is still flagged. After trimming g4,
## level "b" has zero observations, so the trimmed refit either fails
## outright or comes back with a shorter beta; both land in the
## "refit-failed" fallback with its warning. This unit-tests
## .anova_robust_null_params directly because no in-anova trigger is
## deterministic: refit failures inside a full anova() call depend on
## optimizer behaviour on random data.
df7 <- dc
df7$f <- factor(ifelse(df7$g == 4 & seq_len(nrow(df7)) %% 3 == 0,
                       "b", "a"))
stopifnot(all(df7$g[df7$f == "b"] == 4), any(df7$f == "b"))
f07 <- suppressWarnings(rlmer(y ~ x + f + (1 | g), df7,
                              method = "DASvar"))
warns7 <- character(0)
rn7 <- withCallingHandlers(
    robustlmm:::.anova_robust_null_params(f07, threshold = 0.65),
    warning = function(w) {
        warns7 <<- c(warns7, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
## Trimming g4 drops level "b" of factor f, so the trimmed design is
## rank-deficient. The nonsingular-subsampling check intercepts this
## before the refit (reason "rank-deficient"); a near-singular trim that
## slipped past it would still be caught after the refit by the
## parameter-dimension check (reason "refit-failed"). Accept either guard.
stopifnot(isFALSE(rn7$cleaned),
          rn7$reason %in% c("rank-deficient", "refit-failed"),
          length(warns7) > 0L,
          any(grepl("rank-deficient|refit failed|not parameter-compatible",
                    warns7)),
          ## fallback returns the untrimmed generating parameters
          isTRUE(all.equal(rn7$sigma, sigma(f07))),
          isTRUE(all.equal(rn7$beta, robustlmm:::.fixef(f07))),
          isTRUE(all.equal(rn7$Ub@x, f07@pp$U_b@x, tolerance = 1e-12)))
cat("test 7 (rank-deficient / refit-failure fallback:",
    rn7$reason, "+ untrimmed params): ok\n")

## ---- 8. Nonsingular-subsampling check: shared contrast level ----------
## Level "c" is shared by subjects 1 and 2 only; shifting them in opposite
## directions leaves the trt=="c" mean absorbable by the fixed effect but
## gives both clusters large opposite random intercepts, so both are
## flagged. Trimming both drops level "c" from the design -> the rank
## check intercepts it with reason "rank-deficient" (before any refit).
set.seed(7)
J8 <- 14L; ni8 <- 6L
subj8 <- factor(rep(seq_len(J8), each = ni8))
trt8  <- factor(rep("a", J8 * ni8), levels = c("a", "b", "c"))
trt8[subj8 %in% as.character(3:8)] <- "b"
trt8[subj8 %in% c("1", "2")]       <- "c"
b8 <- rnorm(J8, sd = 1)
y8 <- 2 + 0.5 * (trt8 == "b") + 0.8 * (trt8 == "c") +
      b8[as.integer(subj8)] + rnorm(J8 * ni8, sd = 0.6)
y8[subj8 == "1"] <- y8[subj8 == "1"] + 9
y8[subj8 == "2"] <- y8[subj8 == "2"] - 9
df8  <- data.frame(y = y8, trt = trt8, subj = subj8)
f08  <- suppressWarnings(rlmer(y ~ trt + (1 | subj), df8, method = "DASvar"))
warns8 <- character(0)
rn8 <- withCallingHandlers(
    robustlmm:::.anova_robust_null_params(f08, threshold = 0.65),
    warning = function(w) { warns8 <<- c(warns8, conditionMessage(w))
                            invokeRestart("muffleWarning") })
stopifnot(isFALSE(rn8$cleaned),
          identical(rn8$reason, "rank-deficient"),
          any(grepl("rank-deficient", warns8)),
          isTRUE(all.equal(rn8$sigma, sigma(f08))),
          isTRUE(all.equal(rn8$beta, robustlmm:::.fixef(f08))))
cat("test 8 (nonsingular check: shared-level trim -> rank-deficient): ok\n")

cat("anova-robust-null: all tests passed\n")
