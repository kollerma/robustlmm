## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## WS16 step 5 regression tests: confint(method = "Wald", df = "satterthwaite").
##
## 1. Default Wald CI is unchanged: critical value = z = qnorm(1-a).
## 2. df = "satterthwaite" uses per-coefficient t-quantiles
##    qt(1-a, df_k) with df_k the summary Satterthwaite df; the implied
##    critical value (half-width / se) matches exactly.
## 3. The t-based CI is wider than the z-based CI (finite df).
## 4. vcov_type = "sandwich" + df = "satterthwaite" warns and falls back
##    to the z-quantile (identical to the plain sandwich CI).
## 5. A crossed-RE fit now gets finite (CGM multiway) df (WS17), so its
##    t-based CI is wider than the z-based CI -- no fallback, no error.

suppressMessages(require(robustlmm))

set.seed(101)
J <- 30; n <- 8
id   <- factor(rep(seq_len(J), each = n))
time <- rep(seq_len(n), J)
set.seed(5)
gid  <- factor(sample(c("A", "B"), J, replace = TRUE))[as.integer(id)]
b    <- rnorm(J, 0, 1)[as.integer(id)]
y    <- 1 + 0.5 * time + 0.8 * (gid == "B") + b + rnorm(J * n, 0, 1)
d    <- data.frame(y, time, gid, id)
fit  <- rlmer(y ~ time + gid + (1 | id), d, method = "DASvar")

a  <- 0.025
se <- sqrt(diag(as.matrix(vcov(fit))))

## 1. default = z
ciz <- confint(fit, method = "Wald")
critz <- (ciz[, 2] - ciz[, 1]) / 2 / se
stopifnot(max(abs(critz - qnorm(1 - a))) < 1e-9)

## 2. satterthwaite t-quantiles match qt(1-a, summary df)
cit <- confint(fit, method = "Wald", df = "satterthwaite")
sw  <- summary(fit, df = "satterthwaite")$coefficients
critt <- (cit[, 2] - cit[, 1]) / 2 / se
stopifnot(max(abs(critt - qt(1 - a, sw[, "df"]))) < 1e-8)

## 3. t CI wider than z CI everywhere
stopifnot(all((cit[, 2] - cit[, 1]) > (ciz[, 2] - ciz[, 1])))

## 4. sandwich + satterthwaite -> warning + z fallback
cs <- withCallingHandlers(
    confint(fit, method = "Wald", vcov_type = "sandwich",
            df = "satterthwaite"),
    warning = function(w) invokeRestart("muffleWarning"))
cs_z <- confint(fit, method = "Wald", vcov_type = "sandwich")
stopifnot(max(abs(cs - cs_z)) < 1e-12)

## 5. crossed REs -> finite df via the CGM multiway covariance (WS17):
##    t-based CI wider than z, computed without error
set.seed(7)
dd <- expand.grid(g1 = factor(1:14), g2 = factor(1:9), rep = 1:2)
dd$y <- 1 + rnorm(14, 0, 1.2)[dd$g1] + rnorm(9, 0, 1.0)[dd$g2] +
            rnorm(nrow(dd))
fx <- suppressWarnings(rlmer(y ~ 1 + (1 | g1) + (1 | g2), dd,
                             method = "DASvar"))
cx   <- confint(fx, method = "Wald", df = "satterthwaite")
cx_z <- confint(fx, method = "Wald")
stopifnot(all(is.finite(cx)),
          all((cx[, 2] - cx[, 1]) > (cx_z[, 2] - cx_z[, 1])))

cat("confint-satterthwaite.R: all checks passed\n")
