## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## contrast_test(): agreement with anova(), summary(), confint() and
## emmeans for the matching contrasts, on the model-based covariance
## (Satterthwaite df) and on the cluster sandwich.
require(robustlmm)

chk <- function(a, b, tol = 1e-8)
    stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = tol)))

ss <- sleepstudy
ss$grp <- factor(rep(c("a", "b", "c"), length.out = 18))[ss$Subject]
fm <- rlmer(Reaction ~ Days + grp + (Days | Subject), ss)
p <- length(fixef(fm))
nms <- names(fixef(fm))

## ---- one contrast: the coefficient table of summary() --------------------
co <- coef(summary(fm, df = "satterthwaite"))
for (j in seq_len(p)) {
    L <- numeric(p); L[j] <- 1
    ct <- contrast_test(fm, L)
    chk(c(ct$Estimate, ct$`Std. Error`, ct$df, ct$`t value`, ct$`Pr(>|t|)`),
        co[j, c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")],
        tol = 1e-6)
}

## named columns, in any order and a subset
ct1 <- contrast_test(fm, c(Days = 1))
ct2 <- contrast_test(fm, matrix(c(0, 1, 0, 0), 1))
chk(unlist(ct1), unlist(ct2), tol = 1e-12)
ct3 <- contrast_test(fm, c(Days = 5, "(Intercept)" = 1))
chk(ct3$Estimate, sum(fixef(fm)[1:2] * c(1, 5)))
## rhs shifts the estimate, not the standard error
ct4 <- contrast_test(fm, c(Days = 5, "(Intercept)" = 1), rhs = 300)
chk(ct4$Estimate, ct3$Estimate - 300)
chk(ct4$`Std. Error`, ct3$`Std. Error`)

## ---- several contrasts: anova()'s per-term F test ------------------------
av <- anova(fm, ddf = "satterthwaite")
Lg <- diag(p)[grep("^grp", nms), , drop = FALSE]
ctg <- contrast_test(fm, Lg)
chk(c(ctg$NumDF, ctg$DenDF, ctg$`F value`, ctg$`Pr(>F)`),
    unlist(av["grp", c("NumDF", "DenDF", "F value", "Pr(>F)")]), tol = 1e-6)
## a one-row matrix gives the t version of anova's 1-df F test
ctd <- contrast_test(fm, diag(p)[2, , drop = FALSE])
chk(ctd$`t value`^2, av["Days", "F value"], tol = 1e-6)
chk(ctd$df, av["Days", "DenDF"], tol = 1e-6)

## ddf = "none": the chi-square Wald table
avn <- anova(fm, ddf = "none")
ctn <- contrast_test(fm, Lg, ddf = "none")
chk(ctn$`F value` * ctn$NumDF, avn["grp", "Chisq"])
chk(ctn$`Pr(>F)`, avn["grp", "Pr(>Chisq)"])
stopifnot(is.infinite(ctn$DenDF),
          is.infinite(contrast_test(fm, c(Days = 1), ddf = "none")$df))

## ---- the sandwich path: anova(vcov_type = "sandwich") and confint --------
avs <- anova(fm, vcov_type = "sandwich")
cts <- contrast_test(fm, Lg, vcov_type = "sandwich")
chk(c(cts$NumDF, cts$DenDF, cts$`F value`, cts$`Pr(>F)`),
    unlist(avs["grp", c("NumDF", "DenDF", "F value", "Pr(>F)")]))
ct1s <- contrast_test(fm, c(Days = 1), vcov_type = "sandwich")
chk(ct1s$`t value`^2, avs["Days", "F value"])
chk(ct1s$df, avs["Days", "DenDF"])
## its interval is confint()'s
ci <- confint(fm, "Days", vcov_type = "sandwich")
chk(ct1s$Estimate + c(-1, 1) * qt(0.975, ct1s$df) * ct1s$`Std. Error`,
    as.numeric(ci))
## G1 sandwich with the normal reference
cts0 <- contrast_test(fm, c(Days = 1), vcov_type = "sandwich", ddf = "none")
chk(cts0$`Std. Error`, sqrt(vcov(fm, type = "sandwich")["Days", "Days"]))

## ---- emmeans: df and SE of a pairwise contrast ---------------------------
if (requireNamespace("emmeans", quietly = TRUE)) {
    em <- emmeans::emmeans(fm, "grp")
    pw <- as.data.frame(emmeans::contrast(em, "pairwise", adjust = "none"))
    ## a - b = -grpb; a - c = -grpc; b - c = grpb - grpc
    Lp <- rbind(c(-1, 0), c(0, -1), c(1, -1))
    colnames(Lp) <- c("grpb", "grpc")
    for (i in 1:3) {
        ct <- contrast_test(fm, Lp[i, ])
        chk(c(ct$Estimate, ct$`Std. Error`, ct$df),
            c(pw$estimate[i], pw$SE[i], pw$df[i]), tol = 1e-5)
    }
}

## ---- input validation ----------------------------------------------------
err <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
stopifnot(err(contrast_test(fm, c(1, 2))),                 # wrong length
          err(contrast_test(fm, c(Dayz = 1))),             # unknown name
          err(contrast_test(fm, rbind(c(Days = 1), c(Days = 2)))),  # dependent
          err(contrast_test(fm, c(Days = 1), rhs = c(1, 2))),
          err(contrast_test(fm, c(Days = 0))),
          err(contrast_test(fm, c(Days = 1), ddf = "satterthwaite",
                            vcov_type = "sandwich")),
          err(contrast_test(fm, c(Days = 1), ddf = "cluster")),
          err(contrast_test(lmer(Reaction ~ Days + (1 | Subject), ss),
                            c(Days = 1))))
cat("contrast-test.R: all checks passed\n")
