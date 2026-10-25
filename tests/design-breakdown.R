## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## design_breakdown(): breakdown ceiling of the fixed-effects design
## (Mili and Coakley 1996) and Yohai's MM guarantee at b = 1/2.
require(robustlmm)
set.seed(1)
d <- data.frame(y = rnorm(40), g = gl(2, 20), x = rnorm(40), h = gl(4, 10))

## balanced dummy: k = n/2, ceiling 1/4, MM guarantee empty
r1 <- design_breakdown(model.matrix(~ g, d))
stopifnot(r1$k == 20, r1$max_bp == 0.25, r1$mm_bp == 0, r1$b_opt == 0.25,
          r1$exact)
## general position: k = p - 1
r2 <- design_breakdown(model.matrix(~ x, d))
stopifnot(r2$k == 1, r2$max_bp == 0.5, isTRUE(all.equal(r2$mm_bp, 0.475)))
## a continuous covariate does not lift the dummy's ceiling
r3 <- design_breakdown(model.matrix(~ g + x, d))
stopifnot(r3$k == 20, r3$max_bp == 0.25)
## four levels of 10: three levels on one hyperplane
r4 <- design_breakdown(model.matrix(~ h, d))
stopifnot(r4$k == 30, r4$max_bp == 0.125, r4$b_opt == 0.125)
## fits are accepted; sleepstudy's Days = 0 rows form the largest set
fit <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
             method = "DASvar")
r5 <- design_breakdown(fit)
stopifnot(r5$k == 18, r5$n == 180, isTRUE(all.equal(r5$mm_bp, 0.4)))
## the random-search fallback returns a lower bound on k
set.seed(42); u0 <- runif(1); set.seed(42)
r6 <- design_breakdown(model.matrix(~ h, d), max_subsets = 1)
## the random search leaves the caller's RNG stream untouched
stopifnot(identical(runif(1), u0))
stopifnot(!r6$exact, r6$k <= 30)
cat("design-breakdown.R: all checks passed\n")
