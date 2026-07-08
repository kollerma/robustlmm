## Disable test
quit()

## github-transition: MONTE-CARLO TEST (a few seconds per design) -- when
## merging to the github-transition (CRAN-release) branch, disable it there
## by adding `quit()` at the top, per that branch's convention.
##
## Fast, seeded coverage smoke test for the robust Satterthwaite df
## (PLAN-3.5.0-stabilization.md, item 1.2). For each design type the
## Satterthwaite t-CI is meant to attain ~nominal coverage of the fixed
## effects; the full validation (large nsim) lives in the build-ignored
## inference study, this is a regression guard against gross breakage.
##
## It checks, for a single-factor, a nested, a crossed (CGM) and a
## variance-component-boundary design:
##   (a) every fit yields a FINITE Satterthwaite df >= 1, and the MEAN df
##       lands in a per-design range -- the sharp guard, independent of
##       nsim: a regression to the z / Inf-df fallback (or NaN) blows past
##       the upper bound, a collapse to df ~ 1 falls below the lower, and
##   (b) the empirical 95% t-CI coverage of a known fixed-effect
##       coefficient clears a generous floor (gross under-coverage guard).
## The bounds are wide on purpose (small nsim) and carry margin around the
## fixed seeded values, so they are not flaky yet flag a catastrophically
## wrong df or variance. Precise coverage validation lives in the
## build-ignored inference study.

suppressMessages(require(robustlmm))

## small nsim: enough to catch gross breakage, fast enough to run serially
## on win-builder. Seeds are per-replicate (and rlmer is deterministic given
## the data) so the result is fully reproducible regardless of how the work
## is scheduled across cores -- the coverage bands below are therefore not
## flaky, they are fixed numbers with margin around the seeded values.
nsim    <- 20L
## 2 cores (CRAN's check ceiling) off Windows, serial on Windows/win-builder
mc.cores <- if (.Platform$OS.type == "windows") 1L else 2L

## 95% t-CI from the default-vcov SE and the robust Satterthwaite df, for
## one fixed-effect coefficient; returns the df used and whether the CI
## covered the truth. Uses the same internal path as summary(df = "sat").
## .satterthwaite_df returns one value per fixed effect in fixef() order
## (unnamed), so locate the coefficient positionally, as summary() does.
coverOne <- function(fit, coef, truth) {
    res <- tryCatch(robustlmm:::.satterthwaite_df(fit),
                    error = function(e) NULL)
    if (is.null(res)) return(c(df = NA_real_, covered = NA_real_,
                               boundary = NA_real_))
    j   <- match(coef, names(fixef(fit)))
    df  <- pmax(as.numeric(res)[j], 1)
    se  <- attr(res, "se")[j]
    est <- fixef(fit)[[j]]
    hw  <- qt(0.975, df) * se
    c(df = df,
      covered = as.numeric(truth >= est - hw && truth <= est + hw),
      boundary = as.numeric(isTRUE(attr(res, "boundary"))))
}

## run one design: simData(seed) -> data.frame, fit it, return coverOne()
runDesign <- function(simData, fitOne, coef, truth, base) {
    rows <- parallel::mclapply(seq_len(nsim), function(i) {
        set.seed(base + i)
        dat <- simData()
        fit <- tryCatch(suppressWarnings(suppressMessages(fitOne(dat))),
                        error = function(e) NULL)
        if (is.null(fit)) return(c(df = NA_real_, covered = NA_real_,
                                   boundary = NA_real_))
        coverOne(fit, coef, truth)
    }, mc.cores = mc.cores)
    do.call(rbind, rows)
}

## Assertions (all deterministic for the fixed per-replicate seeds, so the
## numeric bounds are not flaky -- they carry margin around the seeded
## values):
##   * almost all replicates fit and produce a finite df,
##   * the mean df lands in [dfLo, dfHi] -- the sharp guard: a regression to
##     the z / Inf-df fallback blows past dfHi, a collapse to df ~ 1 falls
##     below dfLo,
##   * the empirical 95% t-CI coverage is at least `covLo` (gross
##     under-coverage guard; CGM crossed fits legitimately over-cover, so
##     there is no upper coverage bound).
## `expectBoundary` additionally requires the boundary path to fire.
checkDesign <- function(name, M, covLo, dfLo, dfHi, expectBoundary = FALSE) {
    ok  <- is.finite(M[, "df"]) & !is.na(M[, "covered"])
    cov <- mean(M[ok, "covered"]); mdf <- mean(M[ok, "df"])
    cat(sprintf("%-14s  fits=%d/%d  mean.df=%.1f  coverage=%.3f%s\n",
                name, sum(ok), nrow(M), mdf, cov,
                if (expectBoundary)
                    sprintf("  boundary=%.2f", mean(M[ok, "boundary"]))
                else ""))
    stopifnot(sum(ok) >= 0.9 * nrow(M),          # almost all fits usable
              all(M[ok, "df"] >= 1),             # finite, >= 1 (not z/Inf)
              mdf >= dfLo, mdf <= dfHi,           # df not collapsed / Inf
              cov >= covLo)                       # ~nominal coverage
    if (expectBoundary)
        stopifnot(mean(M[ok, "boundary"]) > 0)   # boundary path exercised
}

## ---- 1. single grouping factor: y ~ time + (1 | id) ----------------
beta1 <- c(1, 0.5)
M1 <- runDesign(
    simData = function() {
        J <- 18L; n <- 4L
        id <- factor(rep(seq_len(J), each = n)); time <- rep(seq_len(n), J)
        b  <- rnorm(J, 0, 1)[as.integer(id)]
        data.frame(y = beta1[1] + beta1[2] * time + b + rnorm(J * n),
                   time, id)
    },
    fitOne = function(d) rlmer(y ~ time + (1 | id), d, method = "DASvar"),
    coef = "time", truth = beta1[2], base = 1000L)
## seeded: coverage 0.90, mean df ~37
checkDesign("single-factor", M1, covLo = 0.70, dfLo = 10, dfHi = 150)

## ---- 2. nested: y ~ 1 + (1 | school/class), strong components ------
M2 <- runDesign(
    simData = function() {
        d <- expand.grid(class = factor(1:3), school = factor(1:12),
                         rep = 1:2)
        cs <- interaction(d$school, d$class, drop = TRUE)
        d$y <- 1 + rnorm(12, 0, 1.4)[d$school] +
                   rnorm(nlevels(cs), 0, 1.1)[cs] + rnorm(nrow(d))
        d
    },
    fitOne = function(d) rlmer(y ~ 1 + (1 | school/class), d,
                              method = "DASvar"),
    coef = "(Intercept)", truth = 1, base = 2000L)
## seeded: coverage 1.00, mean df ~11
checkDesign("nested", M2, covLo = 0.75, dfLo = 3, dfHi = 40)

## ---- 3. crossed (CGM multiway): y ~ 1 + (1 | g1) + (1 | g2) ---------
M3 <- runDesign(
    simData = function() {
        d <- expand.grid(g1 = factor(1:10), g2 = factor(1:7), rep = 1L)
        d$y <- 1 + rnorm(10, 0, 1.2)[d$g1] + rnorm(7, 0, 1.0)[d$g2] +
                   rnorm(nrow(d))
        d
    },
    fitOne = function(d) rlmer(y ~ 1 + (1 | g1) + (1 | g2), d,
                              method = "DASvar"),
    coef = "(Intercept)", truth = 1, base = 3000L)
## seeded: coverage 1.00 (CGM is conservative at this size), mean df ~10
checkDesign("crossed", M3, covLo = 0.70, dfLo = 3, dfHi = 40)

## ---- 4. boundary: nested with the INNER component truly 0 -----------
## strong school variance, NO class-within-school variance. The
## class:school theta is estimated at the 0 boundary; the df is then
## reported conditional on it (reducible, since the strong school
## component remains), giving a small finite df where the t-correction
## genuinely matters (cf. WS14). Coverage of the intercept should stay
## ~nominal.
M4 <- runDesign(
    simData = function() {
        d <- expand.grid(class = factor(1:3), school = factor(1:12),
                         rep = 1:2)
        d$y <- 1 + rnorm(12, 0, 1.4)[d$school] + rnorm(nrow(d))  # no class RE
        d
    },
    fitOne = function(d) rlmer(y ~ 1 + (1 | school/class), d,
                              method = "DASvar"),
    coef = "(Intercept)", truth = 1, base = 4000L)
## seeded: coverage 1.00, mean df ~9, boundary path fires on every replicate
checkDesign("boundary", M4, covLo = 0.75, dfLo = 2, dfHi = 40,
            expectBoundary = TRUE)

cat("satterthwaite-coverage.R: all checks passed\n")
