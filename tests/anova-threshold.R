## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Thresholded bootstrap generator for anova(test = "boot")
## (.boot_threshold_theta; Guedon, Baey and Kuhn 2024).
require(robustlmm)

thr_fun <- robustlmm:::.boot_threshold_theta

set.seed(101)
J <- 30; m <- 8
d <- data.frame(g = factor(rep(seq_len(J), each = m)),
                x = rnorm(J * m), z = rnorm(J * m))
b0 <- rnorm(J)
d$y <- 1 + 0.5 * d$x + b0[d$g] + rnorm(J * m)

## 1. A nuisance slope with true variance zero is thresholded; the
##    clearly positive intercept variance is kept.
fit0 <- rlmer(y ~ x + (1 | g) + (0 + x | g), d, method = "DASvar")
th0 <- getME(fit0, "theta")
r <- thr_fun(fit0, th0)
stopifnot(length(r$zeroed) == 1L,
          grepl("x", r$zeroed),
          r$theta[grepl("Intercept", names(th0))] ==
              th0[grepl("Intercept", names(th0))],
          r$theta[!grepl("Intercept", names(th0))] == 0)

## 2. Nothing is thresholded when every component is clearly positive.
fitS <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
              method = "DASvar")
rS <- thr_fun(fitS, getME(fitS, "theta"))
stopifnot(length(rS$zeroed) == 0L,
          identical(rS$theta, getME(fitS, "theta")))

## 3. Zeroing the middle component of a 3 x 3 block keeps the retained
##    covariance and a lower-triangular factor in the block's pattern.
fit3 <- rlmer(y ~ x + z + (1 + x + z | g), d, method = "DASvar")
th3 <- getME(fit3, "theta")
B <- as.matrix(fit3@blocks[[1]])
L <- matrix(0, 3, 3); L[B != 0] <- th3[B[B != 0]]
L[2, ] <- c(1e-4, 1e-4, 0)         # component 2: tiny variance
L[3, ] <- c(0.3, 0.2, 0.8)
thmod <- th3; thmod[B[B != 0]] <- L[B != 0]
r3 <- thr_fun(fit3, thmod)
Ln <- matrix(0, 3, 3); Ln[B != 0] <- r3$theta[B[B != 0]]
S <- tcrossprod(L); Sn <- tcrossprod(Ln)
stopifnot(length(r3$zeroed) == 1L,
          all(Ln[upper.tri(Ln)] == 0),
          all(diag(Ln) >= 0),
          all(Sn[2, ] == 0), all(Sn[, 2] == 0),
          isTRUE(all.equal(Sn[-2, -2], S[-2, -2], tolerance = 1e-12)))

## 3b. The same with a singular retained covariance whose dependent
##     component is not the last retained one (components 1 and 2
##     perfectly correlated, component 4 zeroed): the rebuilt factor must
##     not permute the retained components.
d$w <- rnorm(J * m)
fit4 <- suppressWarnings(
    rlmer(y ~ x + z + w + (1 + x + z + w | g), d, method = "DASvar"))
th4 <- getME(fit4, "theta")
B4 <- as.matrix(fit4@blocks[[1]])
L4 <- rbind(c(1, 0, 0, 0), c(2, 0, 0, 0), c(0.5, 0.7, 0.3, 0),
            c(1e-4, 0, 0, 1e-4))
thmod4 <- th4; thmod4[B4[B4 != 0]] <- L4[B4 != 0]
r4 <- thr_fun(fit4, thmod4)
Ln4 <- matrix(0, 4, 4); Ln4[B4 != 0] <- r4$theta[B4[B4 != 0]]
S4 <- tcrossprod(L4); Sn4 <- tcrossprod(Ln4)
stopifnot(length(r4$zeroed) == 1L,
          all(Ln4[upper.tri(Ln4)] == 0),
          all(Sn4[4, ] == 0),
          isTRUE(all.equal(Sn4[-4, -4], S4[-4, -4], tolerance = 1e-12)))

## 4. anova(test = "boot") reports the thresholded generator in its
##    heading, nuisance = "raw" restores the plain generator, and the
##    seeded p-values are reproducible.
fit1 <- rlmer(y ~ x + (1 | g) + (0 + x | g) + (0 + z | g), d,
              method = "DASvar")
a1 <- suppressWarnings(anova(fit0, fit1, test = "boot", nsim = 12,
                             seed = 7))
a1b <- suppressWarnings(anova(fit0, fit1, test = "boot", nsim = 12,
                              seed = 7))
a2 <- suppressWarnings(anova(fit0, fit1, test = "boot", nsim = 12,
                             seed = 7, nuisance = "raw"))
h1 <- attr(a1, "heading"); h2 <- attr(a2, "heading")
stopifnot(any(grepl("generator: nuisance variance component", h1)),
          !any(grepl("generator:", h2)),
          identical(a1[["Pr(>=Diff)"]], a1b[["Pr(>=Diff)"]]),
          length(attr(a1, "boot")$threshold$zeroed) == 1L,
          is.null(attr(a2, "boot")$threshold))

## 5. A quasi-deviance difference at its atom (both fits put the tested
##    variance on the boundary) is rounding noise; it is treated as 0, so
##    the p-value is 1 rather than an arbitrary noise comparison.
##    (With this seed fit1 puts the z-slope variance at 0: D_obs ~ -1e-6.)
Dobs <- attr(a1, "boot")$D_obs
stopifnot(abs(Dobs) <= 0.01,
          a1[["Pr(>=Diff)"]][2L] == 1, a2[["Pr(>=Diff)"]][2L] == 1,
          any(grepl("at its atom", attr(a1, "heading"))),
          any(grepl("consider test = \"score\"", attr(a1, "heading"))))

## 6. A covariate that fit1 uses only transformed ((0 + log(zp) | g): the
##    model frame stores "log(zp)", not zp) is taken from the data, so the
##    bootstrap refits of fit1 succeed (they all failed before).
d$zp <- exp(d$z)
fitL <- rlmer(y ~ x + (1 | g) + (0 + log(zp) | g), d, method = "DASvar")
fit00 <- rlmer(y ~ x + (1 | g), d, method = "DASvar")
aL <- suppressWarnings(anova(fit00, fitL, test = "boot", nsim = 12,
                             seed = 3))
stopifnot(attr(aL, "boot")$n_fail == 0L)

cat("anova-threshold.R: all checks passed\n")
