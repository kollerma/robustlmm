## Disable test
quit()

## github-transition: structured-covariance test (a few structured DASvar
## fits) -- disable on the github-transition (CRAN-release) branch with a
## top-of-file quit() per that branch's convention if it is too slow there.
##
## Stabilization tests for structured cs/ar1 covariances (PLAN item 6):
## 1. the inference surfaces (vcov, Satterthwaite df, sandwich vcov,
##    anova) behave GRACEFULLY on a structured single-factor fit -- vcov
##    and the sandwich are finite, and summary(df = "satterthwaite") never
##    errors: it either reports a finite df or falls back to the t-table
##    with an explanatory note. (The df machinery differentiates over the
##    unconstrained theta, so it works for cs but cannot handle some
##    structured blocks, e.g. a larger ar1 -- there it degrades to the
##    t-table rather than crashing. Certifying the df for structured fits
##    is part of why the feature stays experimental.)
## 2. the scope guard warns when cs/ar1 is combined with more than one
##    grouping factor (the untested/unvalidated case) and stays silent for
##    a single grouping factor and for diag (fitted directly, exempt).
## The feature stays experimental: cs/ar1 use a manifold projection with
## no consistency proof, validated only for single-factor balanced designs.

suppressMessages({ require(robustlmm); require(lme4) })
if (packageVersion("lme4") < "2.0.0") {
    cat("structured-inference: skipped (needs lme4 >= 2.0-0)\n")
    quit(save = "no")
}

## evaluate expr for its side effects, returning the warning messages
catchW <- function(expr) {
    w <- character(0)
    withCallingHandlers(expr, warning = function(c) {
        w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
    w
}
hasGuard <- function(w) any(grepl("single grouping factor", w, fixed = TRUE))

## ---- cs data: single grouping factor, balanced repeated measures ----
set.seed(1); nsubj <- 30L; truesd <- c(2, 1.5, 1.2); rho <- 0.4
R <- matrix(rho, 3, 3); diag(R) <- 1
Lt <- t(chol(diag(truesd) %*% R %*% diag(truesd)))
dat <- expand.grid(f = factor(c("a", "b", "c")), rep = 1:5,
                   subj = factor(seq_len(nsubj)))
bvals <- t(Lt %*% matrix(rnorm(3 * nsubj), 3, nsubj))
dat$y <- 10 + bvals[cbind(as.integer(dat$subj), as.integer(dat$f))] +
    rnorm(nrow(dat), sd = 0.8)

## graceful-inference check: vcov + sandwich finite; summary(df = "sat")
## never errors and is internally consistent -- a finite df column, or a
## fallback to the t-table with an explanatory note.
checkInference <- function(fit) {
    stopifnot(all(is.finite(as.matrix(vcov(fit)))),
              nrow(vcov(fit)) == fit@pp$p,
              all(is.finite(as.matrix(vcov(fit, type = "sandwich")))))
    s  <- summary(fit, df = "satterthwaite")
    cm <- s$coefficients
    if ("df" %in% colnames(cm))
        stopifnot(all(is.finite(cm[, "df"])), all(cm[, "df"] >= 1),
                  "Pr(>|t|)" %in% colnames(cm))
    else                          # graceful fallback
        stopifnot(identical(colnames(cm),
                            c("Estimate", "Std. Error", "t value")),
                  !is.null(s$dfNote))
    av <- anova(fit)
    stopifnot(inherits(av, "anova") || is.data.frame(av))
    invisible(TRUE)
}

## ---- 1. cs single-factor fit: inference is graceful (df works here) -
fcs <- suppressWarnings(rlmer(y ~ 1 + f + cs(0 + f | subj), dat,
                              method = "DASvar"))
checkInference(fcs)
stopifnot("df" %in% colnames(summary(fcs, df = "satterthwaite")$coefficients))
cat("test 1 (cs inference graceful, df available): ok\n")

## ---- 1b. ar1 fit: inference graceful (df may fall back to t-table) --
set.seed(2); nl <- 4L; arsd <- 1.5; arrho <- 0.6
Rar <- arrho^abs(outer(seq_len(nl), seq_len(nl), "-"))
Lta <- t(chol((arsd^2) * Rar))
da <- expand.grid(f = factor(seq_len(nl)), rep = 1:4,
                  subj = factor(seq_len(25L)))
bv <- t(Lta %*% matrix(rnorm(nl * 25L), nl, 25L))
da$y <- 5 + bv[cbind(as.integer(da$subj), as.integer(da$f))] +
    rnorm(nrow(da), sd = 0.8)
far <- suppressWarnings(rlmer(y ~ 1 + f + ar1(0 + f | subj), da,
                             method = "DASvar"))
checkInference(far)
cat("test 1b (ar1 inference graceful): ok\n")

## ---- 2. scope guard: cs/ar1 + >1 grouping factor warns --------------
dat$site <- factor(rep(seq_len(6L), length.out = nrow(dat)))
m2 <- suppressWarnings(lmer(y ~ 1 + f + cs(0 + f | subj) + (1 | site),
                            dat, REML = TRUE))
stopifnot(hasGuard(catchW(robustlmm:::.checkStructuredCovariance(m2))))
## single grouping factor cs -> no guard warning
m1 <- suppressWarnings(lmer(y ~ 1 + f + cs(0 + f | subj), dat, REML = TRUE))
stopifnot(!hasGuard(catchW(robustlmm:::.checkStructuredCovariance(m1))))
## diag + >1 grouping factor -> exempt (fitted directly), no guard warning
md <- suppressWarnings(lmer(y ~ 1 + f + diag(0 + f | subj) + (1 | site),
                            dat, REML = TRUE))
stopifnot(!hasGuard(catchW(robustlmm:::.checkStructuredCovariance(md))))
cat("test 2 (cs/ar1 >1-factor scope guard): ok\n")

## ---- 3. the guard propagates through rlmer() ------------------------
stopifnot(hasGuard(catchW(rlmer(y ~ 1 + f + cs(0 + f | subj) + (1 | site),
                                dat, method = "DASvar", max.iter = 2L))))
cat("test 3 (guard fires through rlmer): ok\n")

cat("structured-inference: all tests passed\n")
