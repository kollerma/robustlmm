## Stabilization regression tests for the RANSAC initial estimator
## (PLAN-3.5.0-stabilization.md item 4): two graduation claims that must
## not silently regress --
##   (a) determinism: the start is fully reproducible from a `seed`
##       (and, for the `init = "ransac"` string form, from an outer
##       set.seed), and
##   (b) scope: single, crossed and >2-level nested grouping structures
##       all yield a usable high-breakdown start (subsampling is
##       stratified by the finest grouping factor).
## The numbers themselves (phony-rate recovery etc.) are validated in the
## build-ignored inst/simulationStudy/ransacBasin.R and ransacConsensus.R.

suppressMessages(require(robustlmm))

## ---- (a) determinism ------------------------------------------------
set.seed(1)
d <- data.frame(y = rnorm(80), x = rnorm(80), g = factor(rep(1:16, 5)))

## same seed -> identical subsample and fit
r1 <- ransac_lme4(y ~ x + (1 | g), d, K = 12L, seed = 99L)
r2 <- ransac_lme4(y ~ x + (1 | g), d, K = 12L, seed = 99L)
stopifnot(!is.null(r1$fit),
          identical(r1$subset, r2$subset),
          isTRUE(all.equal(fixef(r1$fit), fixef(r2$fit))))

## rlmer_ransac reproducible from a seed
fa <- suppressWarnings(rlmer_ransac(y ~ x + (1 | g), d, K = 12L, seed = 7L))
fb <- suppressWarnings(rlmer_ransac(y ~ x + (1 | g), d, K = 12L, seed = 7L))
stopifnot(isTRUE(all.equal(fixef(fa), fixef(fb))))

## init = "ransac" reproducible under an outer set.seed (the string form
## draws from the ambient RNG, by design)
set.seed(3); fi1 <- rlmer(y ~ x + (1 | g), d, method = "DASvar",
                          init = "ransac")
set.seed(3); fi2 <- rlmer(y ~ x + (1 | g), d, method = "DASvar",
                          init = "ransac")
stopifnot(isTRUE(all.equal(fixef(fi1), fixef(fi2))))
cat("determinism (seed-reproducible start): ok\n")

## ---- (b) scope: crossed ---------------------------------------------
set.seed(5)
dc <- expand.grid(g1 = factor(1:10), g2 = factor(1:6))
dc$x <- rnorm(nrow(dc)); dc$y <- 1 + dc$x + rnorm(nrow(dc))
rc <- ransac_lme4(y ~ x + (1 | g1) + (1 | g2), dc, K = 12L, seed = 1L)
stopifnot(!is.null(rc$fit), is.finite(rc$scale))
fcx <- suppressWarnings(
    rlmer_ransac(y ~ x + (1 | g1) + (1 | g2), dc, K = 12L, seed = 1L))
stopifnot(inherits(fcx, "rlmerMod"))
cat("scope (crossed): ok\n")

## ---- (b) scope: >2-level nested -------------------------------------
set.seed(6)
dn <- expand.grid(rep = 1:4, c = factor(1:2), b = factor(1:2),
                  a = factor(1:5))
dn$x <- rnorm(nrow(dn)); dn$y <- 1 + dn$x + rnorm(nrow(dn))
rn <- suppressWarnings(ransac_lme4(y ~ x + (1 | a/b/c), dn, K = 12L,
                                   seed = 1L))
stopifnot(!is.null(rn$fit), is.finite(rn$scale))
cat("scope (3-level nested): ok\n")

cat("ransac-scope: all tests passed\n")
