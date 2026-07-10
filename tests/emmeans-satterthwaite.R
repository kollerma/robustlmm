## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## WS16 step 3 regression tests: emm_basis.rlmerMod reports a finite
## Satterthwaite df (not Inf) for the default vcov, consistent with
## summary(fit, df = "satterthwaite"), and falls back to Inf for a
## user-supplied vcov.. emmeans is Suggests, so skip when absent.

suppressMessages(require(robustlmm))

if (!requireNamespace("emmeans", quietly = TRUE)) {
    cat("emmeans not installed; skipping emmeans-satterthwaite.R\n")
} else {
    suppressMessages({require(emmeans); require(lme4)})

    set.seed(101)
    J <- 30; n <- 8
    id   <- factor(rep(seq_len(J), each = n))
    time <- rep(seq_len(n), J)
    set.seed(5)
    gid  <- factor(sample(c("A", "B"), J, replace = TRUE))[as.integer(id)]
    b    <- rnorm(J, 0, 1)[as.integer(id)]
    y    <- 1 + 0.5 * time + 0.8 * (gid == "B") + b + rnorm(J * n, 0, 1)
    fit  <- rlmer(y ~ time + gid + (1 | id),
                  data.frame(y, time, gid, id), method = "DASvar")

    sw <- summary(fit, df = "satterthwaite")$coefficients

    ## 1. default vcov: finite df everywhere in a ref_grid
    s_gid <- summary(emmeans(fit, ~ gid))
    stopifnot(all(is.finite(s_gid$df)), all(s_gid$df > 0))

    ## 2. the pure "time" trend contrast == the "time" fixef contrast,
    ##    so its emmeans df must equal summary()'s time df.
    et <- summary(emtrends(fit, ~ 1, var = "time"))
    stopifnot(abs(et$df - sw["time", "df"]) < 1e-6)

    ## 3. the A - B pairwise contrast == the gidB fixef contrast.
    pw <- summary(contrast(emmeans(fit, ~ gid), "pairwise"))
    stopifnot(abs(pw$df[1] - sw["gidB", "df"]) < 1e-6)

    ## 4. user-supplied vcov. -> z-based (Inf) fallback, unchanged
    s_v <- summary(emmeans(fit, ~ gid, vcov. = vcov(fit)))
    stopifnot(all(!is.finite(s_v$df)))

    cat("emmeans-satterthwaite.R: all checks passed\n")
}
