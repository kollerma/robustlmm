## init = "ransac" must fit the FULL data, not the RANSAC subsample.
##
## Regression test for the defect found 2026-07-27 while running the
## thread5-sugasawa study: `.rlmerInit()` builds the working object from
## any merMod passed as `init` (`lobj <- as(init, "rlmerMod")`), so the
## model frame came from that object and `formula`/`data` were ignored.
## `ransac_lme4()$fit` lives on the winning subsample, so
## `rlmer(init = "ransac")` estimated the model on roughly half the data:
## 510 of 1000 observations, 27 of 50 clusters, on the design below.
require(robustlmm)

set.seed(20260728L)
m <- 30L; n_i <- rep(c(8L, 12L), each = m / 2L); n <- sum(n_i)
id <- factor(rep(seq_len(m), times = n_i))
x1 <- rnorm(n); x2 <- rnorm(n)
b  <- matrix(rnorm(2 * m, sd = c(1, 0.5)), m, 2, byrow = TRUE)
y  <- 0.5 + 0.3 * x1 + 0.5 * x2 + b[as.integer(id), 1] +
      b[as.integer(id), 2] * x2 + rnorm(n, sd = 1.2)
dat <- data.frame(y = y, x1 = x1, x2 = x2, id = id)
form <- y ~ x1 + x2 + (1 + x2 | id)

fit <- suppressWarnings(rlmer(form, dat, method = "DASvar", init = "ransac",
                              rho.e = smoothPsi, rho.b = smoothPsi))

## the whole point: every observation and every cluster is in the fit
stopifnot(nobs(fit) == n,
          nrow(ranef(fit)[[1]]) == m,
          length(getME(fit, "u")) == 2L * m)

## rlmer_ransac() routes through the same branch and must agree
fit2 <- suppressWarnings(rlmer_ransac(form, dat, K = 10L, method = "DASvar",
                                      rho.e = smoothPsi, rho.b = smoothPsi))
stopifnot(nobs(fit2) == n, nrow(ranef(fit2)[[1]]) == m)

## The list form of init accepts partial input; anything omitted keeps
## the value from the internal lmer initialisation.  Omitting `u` is the
## supported way to start from a fit on other observations -- setting it
## to zero instead drives theta onto the boundary, which is why
## .ransacInitList() leaves it out.
base <- suppressWarnings(rlmer(form, dat, method = "DASvar",
                               rho.e = smoothPsi, rho.b = smoothPsi))
part <- suppressWarnings(
    rlmer(form, dat, method = "DASvar", rho.e = smoothPsi, rho.b = smoothPsi,
          init = list(fixef = unname(fixef(base)))))
## Starting beta at base's own converged value must land back on
## base: same model, same data, same solution.  (theta[2] is an
## off-diagonal and may legitimately be <= 0, so do not test its sign.)
stopifnot(nobs(part) == n,
          all.equal(unname(fixef(part)), unname(fixef(base)), tolerance = 1e-5),
          all.equal(unname(getME(part, "theta")), unname(getME(base, "theta")),
                    tolerance = 1e-5))

## An init fitted on fewer observations than `data` has rows is the
## silent-data-loss case; it must warn.
sub <- droplevels(dat[1:(n %/% 2L), ])
lsub <- lme4::lmer(form, sub, REML = TRUE)
ww <- tryCatch({
    suppressMessages(rlmer(form, dat, method = "DASvar", init = lsub,
                           rho.e = smoothPsi, rho.b = smoothPsi))
    character(0)
}, warning = function(w) conditionMessage(w))
stopifnot(length(ww) > 0, grepl("was fitted on", ww))

## fitDatasets_rlmer_ransac / _ransac_bisq carry their own inline copy of
## the RANSAC-start call -- they do not route through rlmer_ransac() --
## so they need checking separately.  This is where the first pass at the
## fix missed them.
ds <- generateAnovaDatasets(2, 1, 8, 3, lmeFormula = y ~ 1,
                            heavyTailedMethod = 0)
for (fn in list(fitDatasets_rlmer_ransac, fitDatasets_rlmer_ransac_bisq)) {
    fits <- suppressWarnings(fn(ds, K = 10L))
    for (f in fits)
        stopifnot(nobs(f) == ds$numberOfRows,
                  nrow(ranef(f)[[1]]) == ds$numberOfSubjects)
}

cat("ransac-init-fullwdata.R: OK (incl. fitDatasets wrappers)\n")
