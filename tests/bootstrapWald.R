## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Tests for confint.rlmerMod (Wald path + dispatch to confintROB).

require(robustlmm)

fit <- rlmer(Yield ~ (1 | Batch), Dyestuff, method = "DASvar")
beta <- fixef(fit)
V    <- as.matrix(vcov(fit))
se   <- sqrt(diag(V))

## ----------------------------------------------------------------
## 1. method = "Wald" (closed form): shape, dimnames, ordering,
##    and a byte-for-byte check against beta_hat +/- z * SE.
## ----------------------------------------------------------------
ci_w <- confint(fit)
stopifnot(is.matrix(ci_w))
stopifnot(identical(dim(ci_w), c(length(beta), 2L)))
stopifnot(identical(colnames(ci_w), c("2.5 %", "97.5 %")))
stopifnot(identical(rownames(ci_w), names(beta)))
z95   <- qnorm(0.975)
ref   <- cbind(beta - z95 * se, beta + z95 * se)
stopifnot(max(abs(unname(ci_w) - unname(ref))) < 1e-10)
stopifnot(all(ci_w[, 1] < beta & beta < ci_w[, 2]))
stopifnot(attr(ci_w, "method") == "Wald")
stopifnot(attr(ci_w, "vcov_type") == "default")

## level = 0.99 -> 0.5 % / 99.5 % columns and a wider interval.
ci_99 <- confint(fit, level = 0.99)
stopifnot(identical(colnames(ci_99), c("0.5 %", "99.5 %")))
stopifnot(all((ci_99[, 2] - ci_99[, 1]) > (ci_w[, 2] - ci_w[, 1])))

## parm selection via character and integer paths.
ci_int  <- confint(fit, parm = 1L)
ci_name <- confint(fit, parm = "(Intercept)")
stopifnot(identical(unname(ci_int), unname(ci_name)))
err <- tryCatch(confint(fit, parm = "nope"), error = identity)
stopifnot(inherits(err, "error"))

## ----------------------------------------------------------------
## 2. method = "Wald" with vcov_type = "sandwich": uses the robust
##    cluster sandwich. Same closed form, different V.
## ----------------------------------------------------------------
##    The default for the sandwich (df = "auto") is the CR2 sandwich with a
##    Bell-McCaffrey t df per coefficient (theory-methods F3); df = "none"
##    gives the G1 normal-quantile interval of robustlmm 3.5.0.
ci_sw  <- confint(fit, vcov_type = "sandwich")
V_cr2  <- vcov(fit, type = "sandwich", correction = "CR2")
se_cr2 <- sqrt(diag(as.matrix(V_cr2)))
df_bm  <- vapply(seq_along(beta), function(i) {
    L <- matrix(0, 1L, length(beta)); L[1L, i] <- 1
    robustlmm:::.cr2_df(V_cr2, L)
}, numeric(1))
ref_sw <- cbind(beta - qt(0.975, df_bm) * se_cr2,
                beta + qt(0.975, df_bm) * se_cr2)
stopifnot(max(abs(unname(ci_sw) - unname(ref_sw))) < 1e-10)
stopifnot(attr(ci_sw, "vcov_type") == "sandwich",
          attr(ci_sw, "correction") == "CR2",
          max(abs(attr(ci_sw, "df") - df_bm)) < 1e-12)
## intercept-only balanced design: CR2 is close to G1 and the df to J - 1
V_swr  <- suppressWarnings(vcov(fit, type = "sandwich"))
J_sw   <- attr(V_swr, "n.clusters")
se_sw  <- sqrt(diag(as.matrix(V_swr)))
stopifnot(abs(se_cr2 / se_sw - 1) < 0.05, abs(df_bm - (J_sw - 1)) < 0.5)
ci_swz  <- suppressWarnings(confint(fit, vcov_type = "sandwich", df = "none"))
ref_swz <- cbind(beta - z95 * se_sw, beta + z95 * se_sw)
stopifnot(max(abs(unname(ci_swz) - unname(ref_swz))) < 1e-10)
ci_swc  <- confint(fit, vcov_type = "sandwich", df = "cluster")
stopifnot(max(abs(unname(ci_swc) - unname(ci_sw))) < 1e-12)

## ----------------------------------------------------------------
## 3. method = "boot" / "BCa" delegate to confintROB. Skipped if
##    confintROB is not installed (it is in Suggests).
## ----------------------------------------------------------------
have_conf <- requireNamespace("confintROB", quietly = TRUE)
if (have_conf) {
    nsim_test <- 50L
    ci_b <- suppressWarnings(
        confint(fit, method = "boot", nsim = nsim_test, seed = 20260602L))
    stopifnot(is.matrix(ci_b))
    stopifnot(identical(dim(ci_b), c(length(beta), 2L)))
    stopifnot(identical(rownames(ci_b), names(beta)))
    stopifnot(identical(colnames(ci_b), c("2.5 %", "97.5 %")))
    stopifnot(all(ci_b[, 1] < beta & beta < ci_b[, 2]))
    stopifnot(attr(ci_b, "method")    == "boot")
    ## default boot.type is "wild" (confintROB's own recommendation)
    stopifnot(attr(ci_b, "boot.type") == "wild")
    ## explicit boot.type = "parametric" is still honoured
    ci_bp <- suppressWarnings(
        confint(fit, method = "boot", boot.type = "parametric",
                nsim = nsim_test, seed = 20260602L))
    stopifnot(attr(ci_bp, "boot.type") == "parametric")

    ## Reproducibility under the same seed.
    ci_b2 <- suppressWarnings(
        confint(fit, method = "boot", nsim = nsim_test, seed = 20260602L))
    stopifnot(identical(ci_b, ci_b2))

    ## parm selection (subsetting done in our wrapper, not by
    ## confintROB; the wrapper omits parm in its delegation precisely
    ## because confintROB treats an explicit parm = NULL as "no
    ## parameters").
    ci_b_p <- suppressWarnings(
        confint(fit, parm = "(Intercept)", method = "boot",
                nsim = nsim_test, seed = 20260602L))
    stopifnot(identical(rownames(ci_b_p), "(Intercept)"))
    stopifnot(max(abs(unname(ci_b_p) - unname(ci_b))) < 1e-10)

    ## BCa returns a proper 2-col matrix; this is the regression that
    ## motivated dropping parm from the delegated call.
    ci_bca <- suppressWarnings(
        confint(fit, method = "BCa", nsim = nsim_test, seed = 20260602L))
    stopifnot(is.matrix(ci_bca))
    stopifnot(identical(dim(ci_bca), c(length(beta), 2L)))
    stopifnot(all(is.finite(ci_bca)))
    stopifnot(attr(ci_bca, "method") == "BCa")

    ## vcov_type = "sandwich" cannot be combined with method = "boot"
    ## (confintROB does not use a covariance of the fixed effects): an
    ## error, where earlier versions warned and ignored vcov_type.
    w <- tryCatch(confint(fit, method = "boot", nsim = nsim_test,
                          seed = 20260602L, vcov_type = "sandwich"),
                  error = function(e) e)
    stopifnot(inherits(w, "error"))
    stopifnot(grepl("vcov_type", conditionMessage(w)))
} else {
    cat("confintROB not installed; boot / BCa tests skipped.\n")
    err <- tryCatch(confint(fit, method = "boot", nsim = 30L),
                    error = identity)
    stopifnot(inherits(err, "error"))
    stopifnot(grepl("confintROB", conditionMessage(err)))
}
