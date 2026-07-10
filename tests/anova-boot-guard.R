## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Regression test for the bootstrap-VC-test influential-group guard
## (PLAN-3.5.0-stabilization.md item 3.1, C4). anova(test = "boot")
## simulates the null from fit0's estimates, so a contaminated group that
## biases fit0 makes the test anti-conservative. The guard screens fit0's
## own random-effects robustness weights (a heavily downweighted group is
## one the robust fit treated as an outlier) and warns, naming the groups.
##
## NB: cluster Cook's distance does NOT work as this screen -- the robust
## fit absorbs a group shift into a downweighted random effect, so the
## *influence* the distance measures is removed along with the bias. The
## RE weights are the direct signal; that is what this guard uses.
##
## Checks:
## 1. .re_min_weight_by_group: labeled per-group minima, one per level,
##    for an intercept and a random-slope model.
## 2. a heavily downweighted group fires the guard, naming it; clean data
##    does not.
## 3. anova(test = "boot") emits the guard warning end-to-end on
##    contaminated data and is silent (re: the guard) on clean data.

suppressMessages(require(robustlmm))

mkData <- function(contam, seed = 1L) {
    set.seed(seed); J <- 12L; m <- 6L; n <- J * m
    g <- factor(rep(seq_len(J), each = m)); x <- rnorm(n)
    y <- 1 + 0.5 * x + rnorm(J, 0, 1)[as.integer(g)] + rnorm(n)
    if (contam) y[g == 4] <- y[g == 4] + 8        # group-4 contamination
    data.frame(y, x, g)
}
catchW <- function(expr) {
    w <- character(0)
    withCallingHandlers(expr, warning = function(c) {
        w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
    w
}
guardMsg <- function(w) w[grepl("heavily downweights", w)]

## ---- 1. .re_min_weight_by_group labels per group --------------------
fc <- suppressWarnings(rlmer(y ~ x + (1 | g), mkData(TRUE),
                             method = "DASvar"))
mw <- robustlmm:::.re_min_weight_by_group(fc)
stopifnot(isTRUE(attr(mw, "labeled")),
          length(mw) == 12L,            # one per group level
          "g: 4" %in% names(mw),
          mw[["g: 4"]] < 0.5)           # the contaminated group
## random-slope model: still one entry per level (min over the 2 terms)
set.seed(2); J <- 8L
ds <- data.frame(y = NA, x = rnorm(80), g = factor(rep(1:J, each = 10)))
ds$y <- 1 + ds$x + rnorm(J)[as.integer(ds$g)] + rnorm(80)
fs <- suppressWarnings(rlmer(y ~ x + (1 + x | g), ds, method = "DASvar"))
mws <- robustlmm:::.re_min_weight_by_group(fs)
stopifnot(isTRUE(attr(mws, "labeled")), length(mws) == J)
cat("test 1 (.re_min_weight_by_group labels groups): ok\n")

## ---- 2. guard fires on contamination, silent on clean ---------------
stopifnot(length(catchW(robustlmm:::.bootstrap_group_guard(
    suppressWarnings(rlmer(y ~ x + (1 | g), mkData(FALSE),
                           method = "DASvar"))))) == 0)
wd <- catchW(robustlmm:::.bootstrap_group_guard(fc))
stopifnot(length(guardMsg(wd)) == 1L,
          grepl("g: 4", guardMsg(wd), fixed = TRUE))
cat("test 2 (guard fires only under contamination): ok\n")

## ---- 3. anova(test = "boot") emits the guard end-to-end -------------
f0c <- suppressWarnings(rlmer(y ~ x + (1 | g), mkData(TRUE),
                             method = "DASvar"))
f1c <- suppressWarnings(rlmer(y ~ x + (1 + x | g), mkData(TRUE),
                             method = "DASvar"))
wboot <- catchW(anova(f0c, f1c, test = "boot", nsim = 10L, seed = 1L))
stopifnot(length(guardMsg(wboot)) == 1L)
## clean data: no guard warning from the same call
f0 <- suppressWarnings(rlmer(y ~ x + (1 | g), mkData(FALSE),
                            method = "DASvar"))
f1 <- suppressWarnings(rlmer(y ~ x + (1 + x | g), mkData(FALSE),
                            method = "DASvar"))
wclean <- catchW(anova(f0, f1, test = "boot", nsim = 10L, seed = 1L))
stopifnot(length(guardMsg(wclean)) == 0L)
cat("test 3 (anova boot emits guard end-to-end): ok\n")

cat("anova-boot-guard: all tests passed\n")
