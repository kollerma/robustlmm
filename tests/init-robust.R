## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## rlmer(init = "robust"): high-breakdown start from a marginal S-estimate
## screen and lmer on the retained observations.
require(robustlmm)

set.seed(11)
J <- 20; m <- 10
g <- factor(rep(seq_len(J), each = m))
x <- rnorm(J * m); b <- rnorm(J); e <- rnorm(J * m)
bad <- sample.int(J, 6)
b[bad] <- b[bad] + 8                      # 30% outlying clusters
y <- 2 + x + b[g] + e
d <- data.frame(g = g, x = x, y = y)

## 1. The start is the lmer fit, on the observations the S-estimate of
##    the model's own fixed-effects part keeps (fixed resampling seed).
linit <- lme4::lmer(y ~ x + (1 | g), d, REML = TRUE)
st <- robustlmm:::.robustInitList(
    linit, quote(lme4::lmer(y ~ x + (1 | g), data = d, REML = TRUE)),
    environment())
set.seed(1)
S <- robustbase::lmrob(y ~ x, d, method = "S")
keep <- abs(residuals(S)) <= qnorm(0.995) * S$scale
lk <- lme4::lmer(y ~ x + (1 | g), d[keep, ])
stopifnot(sum(!keep) > 0,
          max(abs(st$fixef - unname(lme4::fixef(lk)))) < 1e-8,
          abs(st$sigma - sigma(lk)) < 1e-8,
          abs(st$theta - unname(lme4::getME(lk, "theta"))) < 1e-8)

## 1b. It is computed for the model as specified (review M1): contrasts=,
##     subset= and variables outside `data` behave as in the fit itself, and
##     the start leaves the random number stream untouched.
set.seed(2)
d$f <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
d$y2 <- d$y + c(a = 0, b = 5, c = 10)[as.character(d$f)]
fk1 <- rlmer(y2 ~ x + f + (1 | g), d, setting = "RSEhb",
             contrasts = list(f = "contr.sum"))
d2 <- d; contrasts(d2$f) <- contr.sum(3)
fk2 <- rlmer(y2 ~ x + f + (1 | g), d2, setting = "RSEhb")
stopifnot(isTRUE(all.equal(fixef(fk1), fixef(fk2))))
rows <- which(as.integer(d$g) <= 15)
ws <- character()
fs1 <- withCallingHandlers(
    rlmer(y ~ x + (1 | g), d, setting = "RSEhb", subset = rows),
    warning = function(w) {
        ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
    })
fs2 <- rlmer(y ~ x + (1 | g), d[rows, ], setting = "RSEhb")
## (subset= used to trigger the warning meant for an init fitted on a
## subsample, although the internal start uses the same observations)
stopifnot(isTRUE(all.equal(fixef(fs1), fixef(fs2))),
          !any(grepl("was fitted on", ws)))
yy <- d$y; xx <- d$x; gg <- d$g
fe <- rlmer(yy ~ xx + (1 | gg), setting = "RSEhb")
fd <- rlmer(y ~ x + (1 | g), d, setting = "RSEhb")
stopifnot(isTRUE(all.equal(unname(fixef(fe)), unname(fixef(fd)))))
set.seed(5); u0 <- runif(1); set.seed(5)
fd2 <- rlmer(y ~ x + (1 | g), d, setting = "RSEhb")
stopifnot(identical(runif(1), u0), identical(fixef(fd2), fixef(fd)))

## 1c. Prior weights and offsets are supported: unit weights change
##     nothing, and an offset is the same as subtracting it from y.
d$w <- 1
fw <- rlmer(y ~ x + (1 | g), d, weights = w, setting = "RSEhb")
stopifnot(isTRUE(all.equal(fixef(fw), fixef(fd))))
d$o <- 0.5 * d$x
d$yo <- d$y + d$o
fo <- rlmer(yo ~ x + offset(o) + (1 | g), d, setting = "RSEhb")
stopifnot(max(abs(fixef(fo) - fixef(fd))) < 1e-5)

## 2. With the high-breakdown psi- and scale-functions the robust start
##    recovers the intercept where the default lmer start does not.
hb <- function(init)
    rlmer(y ~ x + (1 | g), d, method = "DASvar",
          rho.e = makeBisquarePsi(3.443689), rho.b = bisquarePsi,
          rho.sigma.e = sScalePsi, rho.sigma.b = sScalePsi, init = init)
f_rob <- hb("robust")
f_def <- hb(NULL)
stopifnot(abs(fixef(f_rob)[1] - 2) < 0.5, abs(fixef(f_def)[1] - 2) > 1)

## 3. Works with a correlated random-effects block, and with monotone
##    psi-functions reaches the same fit as the default start.
fs_rob <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy, init = "robust")
fs_def <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy)
stopifnot(max(abs(fixef(fs_rob) - fixef(fs_def))) < 1e-4)

## 5. Redescending rho.b: a random effect whose own weight and whose
##    observations' weights are all zero leaves the effects equations
##    singular; the solve now sets it to zero instead of failing (this data
##    set failed with "LU factorization ... near-singular").
set.seed(4030040L)
J <- 50; m <- 10
g <- factor(rep(seq_len(J), each = m))
x <- rnorm(J * m); b <- rnorm(J); e <- rnorm(J * m)
y <- 2 + x + b[g] + e
bb <- sample.int(J * m, round(0.1 * J * m)); y[bb] <- y[bb] + 8
d5 <- data.frame(g = g, x = x, y = y)
f5 <- rlmer(y ~ x + (1 | g), d5, method = "DASvar", setting = "RSEhb")
stopifnot(all(is.finite(fixef(f5))), abs(fixef(f5)[1] - 2) < 0.5)

cat("init-robust.R: all checks passed\n")
