## Test ransac_lme4 and rlmer_ransac.

require(robustlmm)

## --- 1. ransac_lme4 returns a valid lmer fit ----------------------
set.seed(42)
res <- ransac_lme4(Reaction ~ Days + (Days | Subject),
                   data = sleepstudy, K = 30L, sub_frac = 0.5)
stopifnot(methods::is(res$fit, "lmerMod"))
stopifnot(is.finite(res$scale))
stopifnot(length(res$subset) == ceiling(0.5 * nrow(sleepstudy)))
stopifnot(length(res$scales) == 30L)
stopifnot(res$K == 30L)
stopifnot(res$n_sub == ceiling(0.5 * nrow(sleepstudy)))

## --- 2. Best subset really minimises scale ------------------------
stopifnot(res$scale <= min(res$scales, na.rm = TRUE) + 1e-12)

## --- 3. rlmer_ransac yields a usable rlmer fit --------------------
set.seed(99)
fit <- rlmer_ransac(Reaction ~ Days + (Days | Subject),
                     data = sleepstudy, K = 30L, method = "DASvar")
stopifnot(methods::is(fit, "rlmerMod"))
stopifnot(all(is.finite(robustlmm:::.fixef(fit))))
stopifnot(robustlmm:::.sigma(fit) > 0)

## --- 4. Reproducibility via seed argument -------------------------
res1 <- ransac_lme4(Reaction ~ Days + (Days | Subject),
                    data = sleepstudy, K = 10L, seed = 123)
res2 <- ransac_lme4(Reaction ~ Days + (Days | Subject),
                    data = sleepstudy, K = 10L, seed = 123)
stopifnot(identical(res1$subset, res2$subset))
stopifnot(all.equal(res1$scales, res2$scales))

## --- 5. Argument validation ---------------------------------------
err <- tryCatch(
    ransac_lme4(Reaction ~ Days + (Days | Subject),
                data = sleepstudy, K = 0L), error = identity)
stopifnot(inherits(err, "error"))
err <- tryCatch(
    ransac_lme4(Reaction ~ Days + (Days | Subject),
                data = sleepstudy, sub_frac = 1.5), error = identity)
stopifnot(inherits(err, "error"))

## --- 6. init = "ransac" string form ------------------------------
## Plumbs the string form through .rlmerInit -> ransac_lme4 -> rlmerMod path.
set.seed(11)
fitR <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
              method = "DASvar", init = "ransac")
stopifnot(methods::is(fitR, "rlmerMod"))
stopifnot(all(is.finite(robustlmm:::.fixef(fitR))))
stopifnot(robustlmm:::.sigma(fitR) > 0)

## On clean sleepstudy the RANSAC start should not move the fit
## materially compared with the default init: the fixed-effects
## differences should be within a few standard errors of the default fit.
set.seed(12)
fitD <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
              method = "DASvar")
seD <- sqrt(diag(as.matrix(vcov(fitD))))
stopifnot(all(abs(robustlmm:::.fixef(fitR) - robustlmm:::.fixef(fitD)) <
              5 * seD))

## --- 7. init = <unknown string> falls through to the existing error ---
err <- tryCatch(
    rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
          method = "DASvar", init = "ranzac"),
    error = identity)
stopifnot(inherits(err, "error"))
stopifnot(grepl("Unsuitable init", conditionMessage(err)))
