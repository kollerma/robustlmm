## Disable test
quit()

## github-transition: SLOW TEST (~66s) -- when merging to the
## github-transition (CRAN-release) branch, disable it there by adding
## `quit()` at the top, per that branch's convention.
## Tests for anova.rlmerMod (single-fit Wald, pairwise Wald, and
## bootstrap quasi-deviance for variance-component tests).

require(robustlmm)

fit_full <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                  method = "DASvar")
fit_int  <- rlmer(Reaction ~ 1    + (Days | Subject), sleepstudy,
                  method = "DASvar")
fit_int_re <- rlmer(Reaction ~ Days + (1 | Subject), sleepstudy,
                    method = "DASvar")

## ----------------------------------------------------------------
## 1. Single-fit anova(fit): per-term Wald chi-sq table.
## ----------------------------------------------------------------
a1 <- anova(fit_full)
stopifnot(inherits(a1, "anova"))
stopifnot(is.data.frame(a1))
stopifnot(identical(rownames(a1), c("(Intercept)", "Days")))
stopifnot(identical(colnames(a1), c("Df", "Chisq", "Pr(>Chisq)")))
stopifnot(all(a1$Df == 1L))
stopifnot(all(is.finite(a1$Chisq)) && all(a1$Chisq > 0))
stopifnot(all(a1[["Pr(>Chisq)"]] >= 0 & a1[["Pr(>Chisq)"]] <= 1))

## vcov_type = "sandwich" yields a different table but same structure.
a1s <- anova(fit_full, vcov_type = "sandwich")
stopifnot(identical(dim(a1), dim(a1s)))
stopifnot(any(abs(a1$Chisq - a1s$Chisq) > 0.5))    # genuinely different
stopifnot(grepl("sandwich", attr(a1s, "heading")[2]))

## ----------------------------------------------------------------
## 2. Pairwise Wald for fixed-effects-only nested models.
## ----------------------------------------------------------------
a2 <- anova(fit_int, fit_full)              # H0: Days = 0
stopifnot(inherits(a2, "anova"))
stopifnot(identical(rownames(a2), c("fit0", "fit1")))
stopifnot(is.na(a2$Df[1L]) && a2$Df[2L] == 1L)
stopifnot(is.na(a2$Chisq[1L]) && a2$Chisq[2L] > 0)
## The pair-Wald Chisq should match the single-fit Wald for the Days
## term in the larger model (they test the same restriction in this
## simple case).
stopifnot(abs(a2$Chisq[2L] - a1["Days", "Chisq"]) < 1e-10)

## Swap order: anova(fit_full, fit_int) should canonicalise to the
## same result.
a2r <- anova(fit_full, fit_int)
stopifnot(identical(a2[, "Chisq"], a2r[, "Chisq"]))

## Non-nested: error.
fit_other <- rlmer(Reaction ~ I(Days^2) + (Days | Subject),
                   sleepstudy, method = "DASvar")
err <- tryCatch(anova(fit_full, fit_other), error = identity)
stopifnot(inherits(err, "error"))
stopifnot(grepl("nested", conditionMessage(err)))

## Same fixed-effect design: error.
err2 <- tryCatch(anova(fit_full, fit_full), error = identity)
stopifnot(inherits(err2, "error"))
stopifnot(grepl("same fixed-effect design", conditionMessage(err2)))

## ----------------------------------------------------------------
## 3. Variance-component test: forced bootstrap when test = "Wald"
##    is requested (with a warning), and bootstrap path works.
##    nsim = 40 is enough for the smoke test; the package warns about
##    < 200 effective reps.
## ----------------------------------------------------------------
w <- tryCatch(anova(fit_int_re, fit_full, nsim = 40L, seed = 1L,
                    test = "Wald"),
              warning = function(w) w)
stopifnot(inherits(w, "warning"))
stopifnot(grepl("random-effects structure", conditionMessage(w)))

set.seed(123)
a3 <- suppressWarnings(
    anova(fit_int_re, fit_full, nsim = 40L, seed = 1L))
stopifnot(inherits(a3, "anova"))
stopifnot(identical(rownames(a3), c("fit0", "fit1")))
stopifnot(identical(colnames(a3),
                    c("npar", "QuasiDev", "Df", "Diff", "Pr(>=Diff)")))
stopifnot(is.na(a3$Diff[1L]) && a3$Diff[2L] > 0)
## On sleepstudy the random-slope is well-supported, so p should be
## small.
stopifnot(a3[["Pr(>=Diff)"]][2L] < 0.1)
## Bootstrap attribute records the empirical null sample.
boot_attr <- attr(a3, "boot")
stopifnot(!is.null(boot_attr))
stopifnot(length(boot_attr$D_boot) >= 30L)
stopifnot(boot_attr$nsim == 40L)
stopifnot(boot_attr$D_obs == a3$Diff[2L])

## Reproducibility under fixed seed.
a3b <- suppressWarnings(
    anova(fit_int_re, fit_full, nsim = 40L, seed = 1L))
stopifnot(identical(attr(a3, "boot")$D_boot,
                    attr(a3b, "boot")$D_boot))

## ----------------------------------------------------------------
## 4. test = "boot" works on fixed-effects-only nested models too
##    (explicit opt-in; default is Wald). Verify it returns a
##    bootstrap-style table.
## ----------------------------------------------------------------
a4 <- suppressWarnings(
    anova(fit_int, fit_full, test = "boot", nsim = 40L, seed = 1L))
stopifnot(identical(colnames(a4),
                    c("npar", "QuasiDev", "Df", "Diff", "Pr(>=Diff)")))
stopifnot(!is.null(attr(a4, "boot")))
