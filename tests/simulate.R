## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## simulate.rlmerMod: parametric and wild bootstrap samples.
require(robustlmm)

num <- function(x) as.numeric(if (is.data.frame(x)) unlist(x) else x)
chk <- function(a, b, tol = 1e-10)
    stopifnot(isTRUE(all.equal(num(a), num(b), tolerance = tol)))
err <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")

fm <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy)
n <- nobs(fm)

## ---- shape, names, seed ----------------------------------------------------
s1 <- simulate(fm, nsim = 3, seed = 11)
stopifnot(is.data.frame(s1), identical(dim(s1), c(180L, 3L)),
          identical(names(s1), paste0("sim_", 1:3)),
          identical(rownames(s1), rownames(model.frame(fm))),
          !is.null(attr(s1, "seed")))
## reproducible under a seed, which does not disturb the RNG stream
set.seed(1); u1 <- runif(1)
set.seed(1); s1b <- simulate(fm, nsim = 3, seed = 11); u2 <- runif(1)
stopifnot(identical(s1, s1b), identical(u1, u2))
## without a seed the stream advances and the state is recorded
set.seed(2); a <- simulate(fm); b <- simulate(fm)
stopifnot(!identical(a$sim_1, b$sim_1))
set.seed(2); a2 <- simulate(fm)
stopifnot(identical(a$sim_1, a2$sim_1))
for (tp in c("parametric", "wild"))
    stopifnot(identical(simulate(fm, 2, seed = 5, type = tp),
                        simulate(fm, 2, seed = 5, type = tp)))

## ---- parametric: moments match the fitted model -----------------------------
X <- getME(fm, "X"); Z <- as.matrix(getME(fm, "Z"))
Lam <- as.matrix(getME(fm, "Lambda")); s <- sigma(fm)
V <- s^2 * (Z %*% tcrossprod(Lam) %*% t(Z) + diag(n))
mu0 <- drop(X %*% fixef(fm))
nsim <- 4000
Y <- as.matrix(simulate(fm, nsim = nsim, seed = 3))
## mean: standardised deviations behave like N(0, 1)
zm <- (rowMeans(Y) - mu0) / sqrt(diag(V) / nsim)
stopifnot(abs(mean(zm)) < 4 / sqrt(n), max(abs(zm)) < 4.5)
## marginal variances and the within-subject covariances
S <- cov(t(Y))
stopifnot(max(abs(diag(S) / diag(V) - 1)) < 0.15)
i1 <- which(sleepstudy$Subject == "308")
stopifnot(max(abs(S[i1, i1] - V[i1, i1])) / max(V[i1, i1]) < 0.12)
## different subjects are uncorrelated
i2 <- which(sleepstudy$Subject == "309")
stopifnot(max(abs(cov2cor(S)[i1, i2])) < 0.08)

## use.u = TRUE: conditional on the fitted random effects
Yc <- as.matrix(simulate(fm, nsim = nsim, seed = 4, use.u = TRUE))
zc <- (rowMeans(Yc) - fitted(fm)) / (s / sqrt(nsim))
stopifnot(max(abs(zc)) < 4.5,
          max(abs(apply(Yc, 1, var) / s^2 - 1)) < 0.15,
          max(abs(cor(t(Yc))[i1, i1][upper.tri(diag(10))])) < 0.08)

## prior weights and an offset
set.seed(7)
ss <- sleepstudy
ss$w <- runif(n, 0.5, 2); ss$off <- rnorm(n, sd = 20)
fw <- rlmer(Reaction ~ Days + offset(off) + (1 | Subject), ss, weights = w)
Yw <- as.matrix(simulate(fw, nsim = nsim, seed = 8, use.u = TRUE))
stopifnot(max(abs((rowMeans(Yw) - fitted(fw)) /
                  (sigma(fw) / sqrt(ss$w * nsim)))) < 4.5,
          max(abs(apply(Yw, 1, var) * ss$w / sigma(fw)^2 - 1)) < 0.15)
Yw2 <- as.matrix(simulate(fw, nsim = nsim, seed = 9))
mw <- drop(getME(fw, "X") %*% fixef(fw)) + ss$off
tau2 <- as.numeric(VarCorr(fw)$Subject)
stopifnot(max(abs(apply(Yw2, 1, var) /
                  (tau2 + sigma(fw)^2 / ss$w) - 1)) < 0.15,
          abs(mean(rowMeans(Yw2) - mw)) < 4 * sqrt(tau2 / 18 / nsim) + 0.5)

## ---- wild: equal to confintROB's generator in the simple case --------------
if (requireNamespace("confintROB", quietly = TRUE) &&
    requireNamespace("MASS", quietly = TRUE)) {
    gen <- confintROB:::createWildSampleFunction(
        y = getME(fm, "y"), X = as.matrix(getME(fm, "X")),
        id = getME(fm, "flist")[[1]], bet = unname(fixef(fm)))
    set.seed(21); ref <- gen(4)
    mine <- simulate(fm, nsim = 4, seed = 21, type = "wild")
    for (i in 1:4) stopifnot(identical(unname(mine[[i]]), unname(ref[[i]])))
    f1 <- rlmer(Reaction ~ Days + (1 | Subject), sleepstudy)
    gen1 <- confintROB:::createWildSampleFunction(
        y = getME(f1, "y"), X = as.matrix(getME(f1, "X")),
        id = getME(f1, "flist")[[1]], bet = unname(fixef(f1)))
    set.seed(22); ref1 <- gen1(2)
    mine1 <- simulate(f1, nsim = 2, seed = 22, type = "wild")
    for (i in 1:2) stopifnot(identical(unname(mine1[[i]]), unname(ref1[[i]])))
}

## wild: structure of the samples. Per cluster, (y* - X beta) is a multiple
## of the rescaled residual, the multiplier being one of the two weights
yw <- simulate(fm, nsim = 50, seed = 31, type = "wild")
r0 <- getME(fm, "y") - mu0
ratio <- (as.matrix(yw) - mu0) / r0
mammen <- c(-(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2)
## constant sign pattern within a subject, and two distinct multipliers
sgn <- sign(ratio[, 1]); grp <- sleepstudy$Subject
for (g in levels(grp)) {
    q <- ratio[grp == g, , drop = FALSE]
    ## each column: all observations of the subject share one weight
    wt <- q / (q[, 1] / ifelse(q[1, 1] > 0, mammen[2], mammen[1]))
    stopifnot(all(abs(wt - rep(wt[1, ], each = nrow(wt))) < 1e-8),
              all(apply(abs(outer(wt[1, ], mammen, "-")), 1, min) < 1e-8))
}
## Rademacher weights: the residuals only change sign
yr <- simulate(fm, nsim = 20, seed = 32, type = "wild",
               wild.weights = "rademacher")
ym <- simulate(fm, nsim = 1, seed = 33, type = "wild")
rr <- abs(as.matrix(yr) - mu0)
stopifnot(max(abs(rr - rr[, 1])) < 1e-8)

## the offset is not rescaled: it is added back unchanged
fo <- rlmer(Reaction ~ Days + offset(off) + (1 | Subject), ss)
fo0 <- rlmer(I(Reaction - off) ~ Days + (1 | Subject), ss)
chk(simulate(fo, 3, seed = 41, type = "wild"),
    simulate(fo0, 3, seed = 41, type = "wild") + ss$off, tol = 1e-6)
## prior weights: the leverage is that of the weighted design
yww <- simulate(fw, 2, seed = 42, type = "wild")
Xw <- getME(fw, "X") * sqrt(ss$w)
rw <- sqrt(ss$w) * (ss$Reaction - ss$off - drop(getME(fw, "X") %*% fixef(fw)))
i1 <- which(ss$Subject == "308")
P1 <- Xw[i1, ] %*% solve(crossprod(Xw), t(Xw[i1, ]))
sc <- sqrt(diag(MASS::ginv(diag(10) - P1))) * rw[i1] / sqrt(ss$w[i1])
dev <- yww$sim_1[i1] - (drop(getME(fw, "X") %*% fixef(fw)) + ss$off)[i1]
stopifnot(any(abs(dev / sc - mammen[1]) < 1e-8) ||
          any(abs(dev / sc - mammen[2]) < 1e-8),
          diff(range(dev / sc)) < 1e-8)

## ---- clusters of the wild bootstrap -------------------------------------------
## nested: the default is the outermost factor
set.seed(5)
d <- expand.grid(rep = 1:4, sub = 1:3, top = 1:8)
d$sub <- interaction(d$top, d$sub)
d$y <- rnorm(8)[d$top] + rnorm(24, sd = 0.7)[as.integer(d$sub)] +
    rnorm(nrow(d))
d$top <- factor(d$top)
fn <- rlmer(y ~ 1 + (1 | top) + (1 | sub), d)
chk(simulate(fn, 3, seed = 51, type = "wild"),
    simulate(fn, 3, seed = 51, type = "wild", cluster = "top"))
ytop <- simulate(fn, 1, seed = 52, type = "wild")$sim_1
ysub <- simulate(fn, 1, seed = 52, type = "wild", cluster = "sub")$sim_1
stopifnot(!isTRUE(all.equal(ytop, ysub)))
## a vector of clusters
chk(simulate(fn, 2, seed = 53, type = "wild", cluster = as.integer(d$top)),
    simulate(fn, 2, seed = 53, type = "wild", cluster = "top"))
## crossed: refused unless cluster is given
fc <- rlmer(diameter ~ 1 + (1 | plate) + (1 | sample), Penicillin)
stopifnot(err(simulate(fc, type = "wild")),
          is.data.frame(simulate(fc, 2, seed = 1, type = "wild",
                                 cluster = "plate")),
          is.data.frame(simulate(fc, 2, seed = 1)),       # parametric is fine
          err(simulate(fc, type = "wild", cluster = "nofactor")),
          err(simulate(fc, type = "wild", cluster = 1:3)),
          err(simulate(fm, cluster = "Subject")),          # parametric + cluster
          err(simulate(fm, type = "wild", use.u = TRUE)),
          err(simulate(fm, nsim = 0)), err(simulate(fm, nsim = 1.5)))

## ---- the bootstrap of anova() is unchanged ------------------------------------
## its generator is .sim_parametric(); the first replicate of a seeded run
## equals simulate()'s first sample
f0 <- rlmer(Reaction ~ Days + (1 | Subject), sleepstudy)
set.seed(61)
dsim <- robustlmm:::.boot_sim_data(sleepstudy, "Reaction",
                                   drop(getME(f0, "X") %*% fixef(f0)),
                                   t(getME(f0, "Zt")), f0@pp$U_b, sigma(f0),
                                   nobs(f0), f0@pp$q)
chk(dsim$Reaction, simulate(f0, 1, seed = 61)$sim_1, tol = 1e-14)
## and the same stream as before the generators were shared
set.seed(61)
q0 <- f0@pp$q
b_sim <- sigma(f0) * as.numeric(f0@pp$U_b %*% rnorm(q0))
eps <- rnorm(nobs(f0), sd = sigma(f0))
stopifnot(identical(dsim$Reaction,
                    unname(drop(getME(f0, "X") %*% fixef(f0)) +
                           as.numeric(t(getME(f0, "Zt")) %*% b_sim) + eps)))

## ---- prior weights and offset in the bootstrap of anova() ----------------------
## the generator adds the offset and uses the error variance sigma^2 / w,
## like simulate(); for the weighted model with offset `fw` its replicate
## equals simulate()'s under the same seed
set.seed(62)
dw <- robustlmm:::.boot_sim_data(ss, "Reaction",
                                 drop(getME(fw, "X") %*% fixef(fw)) + ss$off,
                                 t(getME(fw, "Zt")), fw@pp$U_b, sigma(fw),
                                 nobs(fw), fw@pp$q,
                                 robustlmm:::.boot_prior_weights(fw))
chk(dw$Reaction, simulate(fw, 1, seed = 62)$sim_1, tol = 1e-14)
stopifnot(is.null(robustlmm:::.boot_prior_weights(f0)))
## the refits find the weights and the offset, also for models fitted
## inside a function; the replicates are then centred at the fitted null
## model (offset included)
mkw <- function() {
    dl <- ss
    list(rlmer(Reaction ~ Days + offset(off) + (1 | Subject), dl,
               weights = w),
         rlmer(Reaction ~ Days + offset(off) + (Days || Subject), dl,
               weights = w))
}
mw <- mkw()
for (tst in c("boot", "score")) {
    aw <- suppressWarnings(anova(mw[[1]], mw[[2]], test = tst, nsim = 12,
                                 seed = 63))
    stopifnot(attr(aw, "boot")$n_fail == 0,
              is.finite(aw[2, ncol(aw)]))
}
bd <- robustlmm:::.boot_data_frame(mw[[1]], mw[[1]])
stopifnot(all(c("w", "off") %in% names(bd)), nrow(bd) == n)

## ---- refit(): equal to update() with the response replaced --------------------
par <- function(f) c(fixef(f), sigma(f), getME(f, "theta"))
for (form in c(Reaction ~ Days + (1 | Subject),
               Reaction ~ Days + offset(off) + (Days | Subject))) {
    ## non-default psi-function, prior weights and offset are kept
    fr <- rlmer(form, ss, weights = w, rho.e = chgDefaults(smoothPsi, k = 2))
    ynew <- simulate(fr, seed = 71)
    r <- refit(fr, ynew)                  # a one-column data frame
    dnew <- ss; dnew$Reaction <- ynew$sim_1
    u <- update(fr, data = dnew)
    chk(par(r), par(u), tol = 1e-5)
    stopifnot(identical(getME(r, "y"), ynew$sim_1),
              identical(r@rho.e@tDefs, fr@rho.e@tDefs),
              identical(weights(r), weights(fr)),
              r@optinfo$conv$opt == 0,
              identical(par(refit(fr, ynew$sim_1)), par(r)),
              ## the original fit is not modified
              identical(getME(fr, "y"), ss$Reaction))
}
stopifnot(err(refit(fm, 1:3)), err(refit(fm, rep(NA_real_, n))))
## nothing of the original call is evaluated again: a fit made inside a
## function, whose data no longer exists, and a fit on a subset
mk <- function() {
    dd <- sleepstudy
    rlmer(Reaction ~ Days + (1 | Subject), dd, subset = Days >= 2)
}
fl <- mk()
w <- NULL
rl <- withCallingHandlers(refit(fl, simulate(fl, seed = 73)$sim_1),
                          warning = function(c) w <<- c(w, conditionMessage(c)))
stopifnot(nobs(rl) == 144, is.null(w), rl@optinfo$conv$opt == 0)
bl <- lme4::bootMer(fl, fixef, nsim = 3, seed = 74)
stopifnot(all(is.finite(bl$t)))
## lme4::bootMer() runs on an rlmerMod through simulate() and refit()
bm <- suppressWarnings(lme4::bootMer(f0, fixef, nsim = 4, seed = 72))
stopifnot(identical(dim(bm$t), c(4L, 2L)), all(is.finite(bm$t)))
chk(bm$t0, fixef(f0))

## ---- seeded bootstrap results are those of the version before simulate() -----
## (reference values computed with robustlmm 3.5.1 before this method was
## added; confint delegates to confintROB and must not be affected)
f1 <- rlmer(Reaction ~ Days + (Days || Subject), sleepstudy)
sc <- suppressWarnings(anova(f0, f1, test = "score", nsim = 20, seed = 102))
chk(sc[2, "Pr(>=Score)"], 0.0476190476, tol = 1e-8)
chk(mean(attr(sc, "boot")$S_boot), -0.7194372310, tol = 1e-4)
if (requireNamespace("confintROB", quietly = TRUE)) {
    ciw <- suppressWarnings(confint(f0, method = "boot", boot.type = "wild",
                                    nsim = 30, seed = 103))
    chk(ciw, c(237.6933382383, 8.4716515804, 259.6761156410, 13.1279700009),
        tol = 1e-4)
    cip <- suppressWarnings(confint(f0, method = "boot",
                                    boot.type = "parametric",
                                    nsim = 30, seed = 103))
    chk(cip, c(243.8365178482, 9.4172454430, 257.3980871406, 11.8986410839),
        tol = 1e-4)
}
cat("simulate.R: all checks passed\n")
