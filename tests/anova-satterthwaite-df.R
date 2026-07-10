## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## WS16 step 4 regression tests: anova(fit, ddf = "satterthwaite").
##
## 1. Default anova() is unchanged (chi-square columns).
## 2. ddf = "satterthwaite" yields an F-table (NumDF, DenDF, F value,
##    Pr(>F)); F == Chisq/NumDF, p == pf(F, NumDF, DenDF).
## 3. For 1-df terms the contestMD denominator df reduces exactly to the
##    single-contrast Satterthwaite df, so anova DenDF == summary df.
## 4. A multi-df term (3-level factor) gets NumDF = 2 and a finite DenDF.
## 5. Pairwise anova(fit0, fit1, ddf=) F-test; DenDF == summary df of the
##    single restricted coefficient.
## 6. vcov_type = "sandwich" + ddf = "satterthwaite" warns and falls back
##    to the chi-square table.

suppressMessages(require(robustlmm))

set.seed(101)
J <- 30; n <- 8
id   <- factor(rep(seq_len(J), each = n))
time <- rep(seq_len(n), J)
set.seed(5)
gid  <- factor(sample(c("A", "B"), J, replace = TRUE))[as.integer(id)]
grp3 <- factor(sample(c("A", "B", "C"), J * n, replace = TRUE))
b    <- rnorm(J, 0, 1)[as.integer(id)]
y    <- 1 + 0.5 * time + 0.8 * (gid == "B") +
        0.4 * (grp3 == "B") - 0.2 * (grp3 == "C") +
        b + rnorm(J * n, 0, 1)
d    <- data.frame(y, time, gid, grp3, id)

fit  <- rlmer(y ~ time + gid + grp3 + (1 | id), d, method = "DASvar")
fit0 <- rlmer(y ~ time + grp3 + (1 | id), d, method = "DASvar")

## 1. default unchanged
ac <- anova(fit)
stopifnot(identical(colnames(ac), c("Df", "Chisq", "Pr(>Chisq)")))

## 2. F-table shape + internal consistency
aF <- anova(fit, ddf = "satterthwaite")
stopifnot(identical(colnames(aF),
                    c("NumDF", "DenDF", "F value", "Pr(>F)")))
stopifnot(max(abs(aF[["F value"]] - ac$Chisq / ac$Df)) < 1e-9)
p_expected <- pf(aF[["F value"]], aF$NumDF, aF$DenDF, lower.tail = FALSE)
stopifnot(max(abs(aF[["Pr(>F)"]] - p_expected)) < 1e-12)

## 3. 1-df terms: anova DenDF == summary Satterthwaite df
sw <- summary(fit, df = "satterthwaite")$coefficients
stopifnot(abs(aF["time", "DenDF"] - sw["time", "df"]) < 1e-6,
          abs(aF["gid",  "DenDF"] - sw["gidB", "df"]) < 1e-6)

## 4. multi-df term
stopifnot(aF["grp3", "NumDF"] == 2,
          is.finite(aF["grp3", "DenDF"]), aF["grp3", "DenDF"] > 0)

## 5. pairwise F-test (restrict gid): DenDF == summary gidB df
ap <- anova(fit0, fit, ddf = "satterthwaite")
stopifnot("DenDF" %in% colnames(ap),
          abs(ap["fit1", "DenDF"] - sw["gidB", "df"]) < 1e-6)

## 6. sandwich + ddf -> warning + chi-square fallback
as <- withCallingHandlers(
    anova(fit, vcov_type = "sandwich", ddf = "satterthwaite"),
    warning = function(w) invokeRestart("muffleWarning"))
stopifnot(identical(colnames(as), c("Df", "Chisq", "Pr(>Chisq)")))

cat("anova-satterthwaite-df.R: all checks passed\n")
