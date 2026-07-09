########################################################
## inferenceStudy.R                                     #
##                                                      #
## Phases A, C, D of the WS6 simulation study           #
## evaluating the new inference machinery in robustlmm: #
##   A: confidence-interval coverage for beta           #
##      (vcov default / sandwich; confint Wald)         #
##   C: single-fit anova(fit) Type I + power            #
##   D: pairwise anova(f0, f1) Type I + power           #
##      (fixed-effects-only nested models; Wald)        #
##                                                      #
## Phase B (bootstrap CI via confintROB) and Phase E    #
## (VC bootstrap anova) are out of scope here -- they   #
## need substantially more compute and are tracked      #
## separately.                                          #
##                                                      #
## Parallelised via parallel::mclapply. Sequential      #
## fallback when n_cores = 1L or on Windows.            #
########################################################

suppressPackageStartupMessages({
    require(robustlmm)
    require(lme4)
})
src_path <- system.file("simulationStudy", package = "robustlmm")
source(file.path(src_path, "randomNumberGenerators.R"))
source(file.path(src_path, "contaminationStrategies.R"))

## --------------------------------------------------------------
## Configuration defaults
## --------------------------------------------------------------
DEFAULTS <- list(
    J_values     = c(8L, 18L, 50L, 100L),
    n_per_subj   = 10L,            # observations per cluster (within-cluster size)
    n_levels     = 3L,             # fixed-effects levels for the ANOVA design
    n_reps       = 500L,
    contamination = c("clean", "t3", "CN_e", "CN_re", "shift_subj"),
    seed         = 20260603L,
    n_cores      = max(1L, parallel::detectCores() - 1L),
    method       = "DASvar",       # rlmer fitting method
    eff_size     = 0.5             # power-cell effect size in sigma units
)

## --------------------------------------------------------------
## Generators                                                    #
##                                                               #
## A single (J, n_per_subj, contamination, effect_size) cell     #
## yields n_reps datasets. We use generateAnovaDatasets() for    #
## the diagonal-V_b workhorse design                             #
##     y = beta_0 + beta_1 * Var1 + b_{Var2} + epsilon           #
## (Var1: 3-level factor, Var2: subject grouping factor).        #
## --------------------------------------------------------------

.make_generator <- function(contamination, groupname = "Var2") {
    switch(contamination,
           clean      = list(err = srnorm,  re = srnorm,  contam_fn = NULL),
           t3         = list(err = srt3,    re = srt3,    contam_fn = NULL),
           CN_e       = list(err = srcnorm, re = srnorm,  contam_fn = NULL),
           CN_re      = list(err = srnorm,  re = srcnorm, contam_fn = NULL),
           shift_subj = list(
               err = srnorm, re = srnorm,
               contam_fn = local({
                   gn <- groupname
                   function(data, eps = 0.10, sigma_e = 1.0)
                       strategy_within(data, eps = eps,
                                        yname = "y", groupname = gn,
                                        magnitude_sigma = 5,
                                        sigma_e = sigma_e)
               })),
           stop("unknown contamination type: ", contamination))
}

##' Build n_reps datasets for one cell.
##' @return list(generator = ds, contam_fn = optional contamination
##'   transformer; trueBeta = the data-generating beta)
.cell_datasets <- function(J, n_per_subj, n_levels, n_reps,
                           contamination, eff_size, cell_seed) {
    g <- .make_generator(contamination)
    ## Build trueBeta: (Intercept) = 1, plus (n_levels - 1) Var1-level
    ## coefficients each = eff_size * seq_along (spreads the levels).
    ## eff_size = 0 -> Var1 has no effect (Phase D Type-I); > 0 ->
    ## Var1 has a graded effect (power cell).
    n_coef <- n_levels                                   # 1 + (n_levels - 1)
    trueBeta <- c(1, eff_size * seq_len(n_levels - 1L))
    ds <- generateAnovaDatasets(
        numberOfDatasetsToGenerate  = n_reps,
        numberOfLevelsInFixedFactor = n_levels,
        numberOfSubjects            = J,
        numberOfReplicates          = n_per_subj,
        errorGenerator       = g$err,
        randomEffectGenerator = g$re,
        lmeFormula = y ~ Var1,
        trueBeta = trueBeta,
        lower = 0, arrange = TRUE)
    list(ds = ds, contam_fn = g$contam_fn,
         trueBeta = trueBeta,
         eff_size = eff_size)
}

##' Longitudinal design for Phase E (VC bootstrap test). Random
##' intercept + slope structure controllable via the trueTheta
##' parameter -- eff_size = 0 means no slope variance (Type-I cell),
##' eff_size > 0 means positive slope SD (power cell). The lme4
##' Cholesky convention used by generateLongitudinalDatasets:
##'   theta = c(sd_intercept,        # = sigma * theta_1
##'             cov_int_slope_via_chol,
##'             sd_slope_orthogonal) # = sigma * theta_3
##' so theta = (1, 0, 0) is "intercept only"; theta = (1, 0, s) gives
##' an uncorrelated intercept/slope with slope SD = s.
.cell_datasets_long <- function(J, n_timepoints, n_reps,
                                 contamination, eff_size, cell_seed) {
    g <- .make_generator(contamination, groupname = "id")
    trueTheta <- c(1, 0, eff_size)
    ds <- generateLongitudinalDatasets(
        numberOfDatasetsToGenerate = n_reps,
        numberOfSubjects           = J,
        numberOfTimepoints         = n_timepoints,
        numberOfTreatmentLevels    = 1L,
        errorGenerator       = g$err,
        randomEffectGenerator = g$re,
        trueBeta  = c(0, 1),      # intercept, time
        trueTheta = trueTheta)
    list(ds = ds, contam_fn = g$contam_fn,
         trueBeta = c(0, 1),
         eff_size = eff_size)
}

## --------------------------------------------------------------
## Method evaluators                                              #
##                                                                #
## Each returns a data.frame row per fixed-effect coefficient    #
## with columns (method, coef, lower, upper, length, covered)    #
## given the truth `beta_true`.                                   #
## --------------------------------------------------------------

.eval_ci_wald <- function(fit, beta_true, vcov_type) {
    ci <- suppressWarnings(
        confint(fit, method = "Wald", vcov_type = vcov_type))
    data.frame(method = paste0("Wald.", vcov_type),
               coef   = rownames(ci),
               lower  = ci[, 1L], upper = ci[, 2L],
               length = ci[, 2L] - ci[, 1L],
               covered = (ci[, 1L] <= beta_true) &
                         (beta_true <= ci[, 2L]),
               stringsAsFactors = FALSE)
}

## Same shape as .eval_ci_wald but for an lmerMod: lme4's
## confint(., method = "Wald") returns CIs for variance components
## too, so we filter to the fixed-effect rows by name.
.eval_ci_wald_lmer <- function(fit_lmer, beta_true, fixef_names) {
    ci <- suppressWarnings(
        confint(fit_lmer, method = "Wald"))
    keep <- rownames(ci) %in% fixef_names
    ci   <- ci[keep, , drop = FALSE]
    ci   <- ci[match(fixef_names, rownames(ci)), , drop = FALSE]
    data.frame(method = "Wald.lmer",
               coef   = rownames(ci),
               lower  = ci[, 1L], upper = ci[, 2L],
               length = ci[, 2L] - ci[, 1L],
               covered = (ci[, 1L] <= beta_true) &
                         (beta_true <= ci[, 2L]),
               stringsAsFactors = FALSE)
}

## Per-term Wald chi-sq for a lmerMod, applied identically to
## .anova_single_wald's logic for rlmer (extract V-block per term,
## form bk' V_k^-1 bk). Returns a 3-col data.frame: term, Chisq, p.
.anova_lmer_single <- function(fit_lmer) {
    mm <- model.matrix(fit_lmer)
    assign <- attr(mm, "assign")
    tlabs  <- attr(terms(fit_lmer), "term.labels")
    hasInt <- attr(terms(fit_lmer), "intercept") == 1L
    V  <- as.matrix(vcov(fit_lmer))
    bh <- fixef(fit_lmer)
    unique_aids <- if (hasInt) c(0L, seq_along(tlabs)) else seq_along(tlabs)
    term_names  <- if (hasInt) c("(Intercept)", tlabs) else tlabs
    Df <- integer(length(unique_aids))
    Chisq <- numeric(length(unique_aids))
    for (k in seq_along(unique_aids)) {
        idx <- which(assign == unique_aids[k])
        Df[k] <- length(idx)
        Vk <- V[idx, idx, drop = FALSE]
        bk <- bh[idx]
        Chisq[k] <- as.numeric(t(bk) %*% solve(Vk, bk))
    }
    data.frame(term = term_names,
               Df = Df, Chisq = Chisq,
               p = stats::pchisq(Chisq, Df, lower.tail = FALSE),
               stringsAsFactors = FALSE)
}

##' One rep -> one data.frame of Phase-A coverage records (one row
##' per (method, coef)). Fits rlmer + lmer side by side on the same
##' simulated dataset so the rlmer-vs-lmer comparison is paired.
.phaseA_rep <- function(rep_idx, cell, formula) {
    dat <- cell$ds$generateData(rep_idx)
    if (!is.null(cell$contam_fn))
        dat <- cell$contam_fn(dat)
    fit <- tryCatch(
        suppressMessages(suppressWarnings(
            rlmer(formula, data = dat, method = "DASvar"))),
        error = function(e) NULL)
    fit_lm <- tryCatch(
        suppressMessages(suppressWarnings(
            lmer(formula, data = dat, REML = TRUE))),
        error = function(e) NULL)
    if (is.null(fit) && is.null(fit_lm)) return(NULL)
    beta_true <- cell$trueBeta
    out <- list()
    if (!is.null(fit)) {
        out[[length(out) + 1L]] <- .eval_ci_wald(fit, beta_true, "default")
        out[[length(out) + 1L]] <- .eval_ci_wald(fit, beta_true, "sandwich")
    }
    if (!is.null(fit_lm))
        out[[length(out) + 1L]] <- .eval_ci_wald_lmer(
            fit_lm, beta_true, names(fixef(fit_lm)))
    do.call(rbind, out)
}

##' One rep -> Phase-C records (per-term anova(fit) p-values), with
##' rlmer + lmer paired on the same dataset.
.phaseC_rep <- function(rep_idx, cell, formula) {
    dat <- cell$ds$generateData(rep_idx)
    if (!is.null(cell$contam_fn))
        dat <- cell$contam_fn(dat)
    fit <- tryCatch(
        suppressMessages(suppressWarnings(
            rlmer(formula, data = dat, method = "DASvar"))),
        error = function(e) NULL)
    fit_lm <- tryCatch(
        suppressMessages(suppressWarnings(
            lmer(formula, data = dat, REML = TRUE))),
        error = function(e) NULL)
    if (is.null(fit) && is.null(fit_lm)) return(NULL)
    out <- list()
    if (!is.null(fit)) {
        a_d <- anova(fit, vcov_type = "default")
        a_s <- anova(fit, vcov_type = "sandwich")
        out[[length(out) + 1L]] <- data.frame(
            method = "anova.Wald.default",
            term = rownames(a_d), Df = a_d$Df, Chisq = a_d$Chisq,
            p = a_d[["Pr(>Chisq)"]], stringsAsFactors = FALSE)
        out[[length(out) + 1L]] <- data.frame(
            method = "anova.Wald.sandwich",
            term = rownames(a_s), Df = a_s$Df, Chisq = a_s$Chisq,
            p = a_s[["Pr(>Chisq)"]], stringsAsFactors = FALSE)
    }
    if (!is.null(fit_lm)) {
        a_l <- .anova_lmer_single(fit_lm)
        out[[length(out) + 1L]] <- data.frame(
            method = "anova.lmer", term = a_l$term, Df = a_l$Df,
            Chisq = a_l$Chisq, p = a_l$p, stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
}

##' One rep -> Phase-D records (pairwise Wald anova: fit_null vs
##' fit_full, where fit_full adds Var1 as a fixed effect).
##'
##' Type-I cell  : fit_null is the truth (eff_size = 0) -> p should
##'   be Uniform(0, 1).
##' Power cell   : eff_size > 0 -> p should be small.
##' One rep -> Phase-E records (VC bootstrap quasi-deviance test on
##' the PSD-cone boundary): fit_null = random intercept only,
##' fit_full = intercept + slope. Type-I cell: true theta_slope = 0;
##' Power cell: theta_slope > 0. Runs anova(test = "boot") with
##' nsim_boot inner replicates.
## Custom bootstrap for Phase F: *unweighted* residual bootstrap
## instead of the Gaussian parametric simulation. Resamples the
## rlmer-conditional residuals of fit0 with replacement (the `prob`
## argument is set from `weights(fit0)`, which under rlmer returns
## the lme4 model weights -- all 1s for an unweighted fit. So this
## is effectively uniform resampling of the residuals, NOT
## robustness-weighted resampling).
##
## The idea behind residual bootstrap: the empirical residual
## distribution carries the contamination signature that Phase E's
## Gaussian step drops, so under shift_subj the bootstrap null may
## become wider and the Type-I inflation should shrink.
##
## The follow-up variant (Phase G, `weighted = TRUE`) replaces
## `weights(fit0)` with `getME(fit0, "w_e")` -- the actual rlmer
## robustness weights -- to downweight (not drop) the contaminated
## residuals at sampling time.
.anova_pair_boot_resid <- function(fit0, fit1, nsim, seed = NULL,
                                   weighted = FALSE) {
    if (!is.null(seed)) set.seed(seed)
    sig0_obs <- robustlmm:::.sigma(fit0)
    quasi_dev <- function(fit, scale = robustlmm:::.sigma(fit)) {
        r_std <- residuals(fit) / scale
        2 * sum(fit@rho.e@rho(r_std))
    }
    D_obs <- quasi_dev(fit0, scale = sig0_obs) -
             quasi_dev(fit1, scale = sig0_obs)

    pp0    <- fit0@pp
    sig0   <- sig0_obs
    th0    <- getME(fit0, "theta")
    th1    <- getME(fit1, "theta")
    beta0  <- robustlmm:::.fixef(fit0)
    X      <- pp0$X
    Z0     <- t(pp0$Zt)
    Ub0    <- pp0$U_b
    n      <- pp0$n; q0 <- pp0$q
    yname  <- as.character(formula(fit0)[[2L]])
    data_obs <- fit0@frame
    Xb0    <- as.numeric(X %*% beta0)

    ## Residuals + sampling weights from fit0. weighted = FALSE
    ## (Phase F) uses weights(fit0) = lme4 model weights = all 1s,
    ## i.e. uniform resampling. weighted = TRUE (Phase G) uses the
    ## rlmer robustness weights, downweighting contaminated residuals.
    r_obs <- residuals(fit0)
    w_obs <- if (weighted)
                 pmax(getME(fit0, "w_e"), .Machine$double.eps)
             else
                 pmax(weights(fit0), .Machine$double.eps)

    D_boot <- rep(NA_real_, nsim)
    n_fail <- 0L
    for (r in seq_len(nsim)) {
        b_sim <- sig0 * as.numeric(Ub0 %*% rnorm(q0))
        eps   <- sample(r_obs, n, replace = TRUE, prob = w_obs)
        y_sim <- Xb0 + as.numeric(Z0 %*% b_sim) + eps
        d_sim <- data_obs; d_sim[[yname]] <- y_sim
        f0_s <- tryCatch(suppressMessages(suppressWarnings(
                    update(fit0, data = d_sim, start = th0))),
                    error = function(e) NULL)
        f1_s <- tryCatch(suppressMessages(suppressWarnings(
                    update(fit1, data = d_sim, start = th1))),
                    error = function(e) NULL)
        if (is.null(f0_s) || is.null(f1_s)) {
            n_fail <- n_fail + 1L; next
        }
        sig0_s <- robustlmm:::.sigma(f0_s)
        D_boot[r] <- quasi_dev(f0_s, scale = sig0_s) -
                     quasi_dev(f1_s, scale = sig0_s)
    }
    D_clean <- D_boot[!is.na(D_boot)]
    if (length(D_clean) < 10L)
        return(NULL)
    list(D_obs = D_obs, D_boot = D_clean,
         p.value = mean(D_clean >= D_obs),
         n_fail = n_fail, nsim = nsim)
}

##' One rep -> Phase-F record (VC bootstrap quasi-deviance test with
##' residual bootstrap instead of Gaussian parametric simulation).
##' Same dataset structure as Phase E.
.phaseF_rep <- function(rep_idx, cell, formula_full, formula_null,
                        nsim_boot) {
    dat <- cell$ds$generateData(rep_idx)
    if (!is.null(cell$contam_fn))
        dat <- cell$contam_fn(dat)
    fit_full <- tryCatch(
        suppressMessages(suppressWarnings(
            do.call(rlmer, list(formula = formula_full, data = dat,
                                 method = "DASvar")))),
        error = function(e) NULL)
    fit_null <- tryCatch(
        suppressMessages(suppressWarnings(
            do.call(rlmer, list(formula = formula_null, data = dat,
                                 method = "DASvar")))),
        error = function(e) NULL)
    if (is.null(fit_full) || is.null(fit_null)) return(NULL)
    boot_out <- tryCatch(
        .anova_pair_boot_resid(fit_null, fit_full, nsim = nsim_boot,
                               seed = 31L * rep_idx + 1L),
        error = function(e) NULL)
    if (is.null(boot_out)) return(NULL)
    data.frame(method = "pair.boot.QD.resid",
               Df       = length(getME(fit_full, "theta")) -
                          length(getME(fit_null, "theta")),
               Diff     = boot_out$D_obs,
               p        = boot_out$p.value,
               nsim_eff = length(boot_out$D_boot),
               stringsAsFactors = FALSE)
}

##' One rep -> Phase-G record: same as Phase F but with the rlmer
##' robustness weights (`getME(fit0, "w_e")`) as resampling
##' probabilities -- the variant Phase F was originally meant to be.
.phaseG_rep <- function(rep_idx, cell, formula_full, formula_null,
                        nsim_boot) {
    dat <- cell$ds$generateData(rep_idx)
    if (!is.null(cell$contam_fn))
        dat <- cell$contam_fn(dat)
    fit_full <- tryCatch(
        suppressMessages(suppressWarnings(
            do.call(rlmer, list(formula = formula_full, data = dat,
                                 method = "DASvar")))),
        error = function(e) NULL)
    fit_null <- tryCatch(
        suppressMessages(suppressWarnings(
            do.call(rlmer, list(formula = formula_null, data = dat,
                                 method = "DASvar")))),
        error = function(e) NULL)
    if (is.null(fit_full) || is.null(fit_null)) return(NULL)
    boot_out <- tryCatch(
        .anova_pair_boot_resid(fit_null, fit_full, nsim = nsim_boot,
                               seed = 31L * rep_idx + 2L,
                               weighted = TRUE),
        error = function(e) NULL)
    if (is.null(boot_out)) return(NULL)
    data.frame(method = "pair.boot.QD.residw",
               Df       = length(getME(fit_full, "theta")) -
                          length(getME(fit_null, "theta")),
               Diff     = boot_out$D_obs,
               p        = boot_out$p.value,
               nsim_eff = length(boot_out$D_boot),
               stringsAsFactors = FALSE)
}

## Custom bootstrap for Phase H: *block* (subject-level) residual
## bootstrap. The Phase F mechanism analysis showed why per-obs
## residual resampling can't fix shift_subj: the bootstrap's
## random-effects draw `b_sim ~ N(0, sigma^2 Ub(theta_hat_0) Ub^T)`
## inherits the inflated theta_intercept_hat that fit0 absorbed from
## the shifted subjects, so the bootstrap null is biased in the SAME
## direction as the data. Phase H breaks that chain by dropping the
## Gaussian b_sim draw entirely: it resamples whole subjects' TOTAL
## residual vectors r_tot_j = y_j - X_j beta_hat_0 (random effect +
## conditional residual together, preserving within-subject
## correlation and the empirical, possibly contaminated, distribution
## of subject effects). Shifted subjects then enter the bootstrap
## null as discrete outlying blocks -- exactly how they enter the
## data -- rather than as inflated Gaussian variance.
##
## Anticipated cost: in power cells the unmodeled slope variance also
## lives in r_tot and is preserved by block transplantation, so the
## bootstrap null is no longer slope-free and power should drop
## towards alpha. Quantifying that Type-I/power trade-off is the
## point of running both cell types.
##
## Requires a balanced design (equal block sizes, common time
## ordering) -- true for .cell_datasets_long.
.anova_pair_boot_block <- function(fit0, fit1, nsim, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    sig0_obs <- robustlmm:::.sigma(fit0)
    quasi_dev <- function(fit, scale = robustlmm:::.sigma(fit)) {
        r_std <- residuals(fit) / scale
        2 * sum(fit@rho.e@rho(r_std))
    }
    D_obs <- quasi_dev(fit0, scale = sig0_obs) -
             quasi_dev(fit1, scale = sig0_obs)

    pp0    <- fit0@pp
    th0    <- getME(fit0, "theta")
    th1    <- getME(fit1, "theta")
    beta0  <- robustlmm:::.fixef(fit0)
    X      <- pp0$X
    yname  <- as.character(formula(fit0)[[2L]])
    data_obs <- fit0@frame
    Xb0    <- as.numeric(X %*% beta0)

    grp <- getME(fit0, "flist")[[1L]]
    idx_by_subj <- split(seq_along(grp), grp)
    block_sizes <- lengths(idx_by_subj)
    if (length(unique(block_sizes)) != 1L)
        stop("block bootstrap requires a balanced design")
    n_subj <- length(idx_by_subj)

    y_obs <- data_obs[[yname]]
    r_tot <- y_obs - Xb0

    D_boot <- rep(NA_real_, nsim)
    n_fail <- 0L
    for (r in seq_len(nsim)) {
        donor <- sample.int(n_subj, n_subj, replace = TRUE)
        y_sim <- Xb0
        for (j in seq_len(n_subj))
            y_sim[idx_by_subj[[j]]] <-
                Xb0[idx_by_subj[[j]]] + r_tot[idx_by_subj[[donor[j]]]]
        d_sim <- data_obs; d_sim[[yname]] <- y_sim
        f0_s <- tryCatch(suppressMessages(suppressWarnings(
                    update(fit0, data = d_sim, start = th0))),
                    error = function(e) NULL)
        f1_s <- tryCatch(suppressMessages(suppressWarnings(
                    update(fit1, data = d_sim, start = th1))),
                    error = function(e) NULL)
        if (is.null(f0_s) || is.null(f1_s)) {
            n_fail <- n_fail + 1L; next
        }
        sig0_s <- robustlmm:::.sigma(f0_s)
        D_boot[r] <- quasi_dev(f0_s, scale = sig0_s) -
                     quasi_dev(f1_s, scale = sig0_s)
    }
    D_clean <- D_boot[!is.na(D_boot)]
    if (length(D_clean) < 10L)
        return(NULL)
    list(D_obs = D_obs, D_boot = D_clean,
         p.value = mean(D_clean >= D_obs),
         n_fail = n_fail, nsim = nsim)
}

##' One rep -> Phase-H record (VC bootstrap quasi-deviance test with
##' subject-block residual bootstrap). Same dataset structure as
##' Phase E.
.phaseH_rep <- function(rep_idx, cell, formula_full, formula_null,
                        nsim_boot) {
    dat <- cell$ds$generateData(rep_idx)
    if (!is.null(cell$contam_fn))
        dat <- cell$contam_fn(dat)
    fit_full <- tryCatch(
        suppressMessages(suppressWarnings(
            do.call(rlmer, list(formula = formula_full, data = dat,
                                 method = "DASvar")))),
        error = function(e) NULL)
    fit_null <- tryCatch(
        suppressMessages(suppressWarnings(
            do.call(rlmer, list(formula = formula_null, data = dat,
                                 method = "DASvar")))),
        error = function(e) NULL)
    if (is.null(fit_full) || is.null(fit_null)) return(NULL)
    boot_out <- tryCatch(
        .anova_pair_boot_block(fit_null, fit_full, nsim = nsim_boot,
                               seed = 31L * rep_idx + 3L),
        error = function(e) NULL)
    if (is.null(boot_out)) return(NULL)
    data.frame(method = "pair.boot.QD.block",
               Df       = length(getME(fit_full, "theta")) -
                          length(getME(fit_null, "theta")),
               Diff     = boot_out$D_obs,
               p        = boot_out$p.value,
               nsim_eff = length(boot_out$D_boot),
               stringsAsFactors = FALSE)
}

.phaseE_rep <- function(rep_idx, cell, formula_full, formula_null,
                        nsim_boot) {
    dat <- cell$ds$generateData(rep_idx)
    if (!is.null(cell$contam_fn))
        dat <- cell$contam_fn(dat)
    ## Use do.call so the formula object is captured by value in the
    ## fit's $call; bare `rlmer(formula_full, ...)` would capture the
    ## symbol `formula_full`, which is out of scope when the bootstrap
    ## later does update(fit, data = d_sim).
    fit_full <- tryCatch(
        suppressMessages(suppressWarnings(
            do.call(rlmer, list(formula = formula_full, data = dat,
                                 method = "DASvar")))),
        error = function(e) NULL)
    fit_null <- tryCatch(
        suppressMessages(suppressWarnings(
            do.call(rlmer, list(formula = formula_null, data = dat,
                                 method = "DASvar")))),
        error = function(e) NULL)
    if (is.null(fit_full) || is.null(fit_null)) return(NULL)
    boot_out <- tryCatch(
        suppressMessages(suppressWarnings(
            anova(fit_null, fit_full, test = "boot",
                  nsim = nsim_boot,
                  seed = 31L * rep_idx))),
        error = function(e) NULL)
    if (is.null(boot_out)) return(NULL)
    data.frame(method = "pair.boot.QD",
               Df     = boot_out$Df[2L],
               Diff   = boot_out$Diff[2L],
               p      = boot_out[["Pr(>=Diff)"]][2L],
               nsim_eff = length(attr(boot_out, "boot")$D_boot),
               stringsAsFactors = FALSE)
}

.phaseD_rep <- function(rep_idx, cell, formula_full, formula_null) {
    dat <- cell$ds$generateData(rep_idx)
    if (!is.null(cell$contam_fn))
        dat <- cell$contam_fn(dat)
    fit_full <- tryCatch(
        suppressMessages(suppressWarnings(
            rlmer(formula_full, data = dat, method = "DASvar"))),
        error = function(e) NULL)
    fit_null <- tryCatch(
        suppressMessages(suppressWarnings(
            rlmer(formula_null, data = dat, method = "DASvar"))),
        error = function(e) NULL)
    fit_full_lm <- tryCatch(
        suppressMessages(suppressWarnings(
            lmer(formula_full, data = dat, REML = FALSE))),
        error = function(e) NULL)
    fit_null_lm <- tryCatch(
        suppressMessages(suppressWarnings(
            lmer(formula_null, data = dat, REML = FALSE))),
        error = function(e) NULL)
    out <- list()
    if (!is.null(fit_full) && !is.null(fit_null)) {
        a_d <- suppressWarnings(
            anova(fit_null, fit_full, vcov_type = "default"))
        a_s <- suppressWarnings(
            anova(fit_null, fit_full, vcov_type = "sandwich"))
        out[[length(out) + 1L]] <- data.frame(
            method = c("pair.Wald.default", "pair.Wald.sandwich"),
            Df     = c(a_d$Df[2L],     a_s$Df[2L]),
            Chisq  = c(a_d$Chisq[2L],  a_s$Chisq[2L]),
            p      = c(a_d[["Pr(>Chisq)"]][2L],
                       a_s[["Pr(>Chisq)"]][2L]),
            stringsAsFactors = FALSE)
    }
    if (!is.null(fit_full_lm) && !is.null(fit_null_lm)) {
        a_l <- suppressMessages(suppressWarnings(
            anova(fit_null_lm, fit_full_lm)))
        ## lme4's anova.merMod LRT puts the larger model on the
        ## second row; extract its Chisq/Df/p.
        out[[length(out) + 1L]] <- data.frame(
            method = "pair.LRT.lmer",
            Df     = a_l[["Df"]][2L],
            Chisq  = a_l[["Chisq"]][2L],
            p      = a_l[["Pr(>Chisq)"]][2L],
            stringsAsFactors = FALSE)
    }
    if (length(out) == 0L) return(NULL)
    do.call(rbind, out)
}

## --------------------------------------------------------------
## Cell drivers                                                   #
## --------------------------------------------------------------

.with_seed <- function(seed, body) {
    if (!is.null(seed)) set.seed(seed)
    body
}

.run_phase_cell <- function(phase, J, contamination, eff_size, n_reps,
                            n_per_subj, n_levels, cell_seed, n_cores,
                            nsim_boot = 200L, verbose = TRUE) {
    if (verbose)
        cat(sprintf("  phase=%s  J=%d  contam=%-12s  eff_size=%.2f  n_reps=%d\n",
                    phase, J, contamination, eff_size, n_reps))
    cell <- .with_seed(cell_seed,
                       if (phase %in% c("E", "F", "G", "H"))
                           .cell_datasets_long(J,
                                                n_timepoints = n_per_subj,
                                                n_reps,
                                                contamination,
                                                eff_size, cell_seed)
                       else
                           .cell_datasets(J, n_per_subj, n_levels, n_reps,
                                          contamination, eff_size, cell_seed))

    rep_fun <- switch(phase,
        "A" = function(r) .phaseA_rep(r, cell, formula = y ~ Var1 + (1 | Var2)),
        "C" = function(r) .phaseC_rep(r, cell, formula = y ~ Var1 + (1 | Var2)),
        "E" = function(r) .phaseE_rep(r, cell,
                                       formula_full = y ~ time + (time | id),
                                       formula_null = y ~ time + (1   | id),
                                       nsim_boot = nsim_boot),
        "F" = function(r) .phaseF_rep(r, cell,
                                       formula_full = y ~ time + (time | id),
                                       formula_null = y ~ time + (1   | id),
                                       nsim_boot = nsim_boot),
        "G" = function(r) .phaseG_rep(r, cell,
                                       formula_full = y ~ time + (time | id),
                                       formula_null = y ~ time + (1   | id),
                                       nsim_boot = nsim_boot),
        "H" = function(r) .phaseH_rep(r, cell,
                                       formula_full = y ~ time + (time | id),
                                       formula_null = y ~ time + (1   | id),
                                       nsim_boot = nsim_boot),
        "D" = function(r) .phaseD_rep(r, cell,
                                       formula_full = y ~ Var1 + (1 | Var2),
                                       formula_null = y ~ 1 + (1 | Var2)),
        stop("unknown phase: ", phase))

    if (n_cores > 1L && phase %in% c("E", "F", "G", "H")) {
        ## Phase E does ~400 rlmer refits per outer rep (200 boot x 2
        ## fits). mclapply (with either mc.preschedule choice) silently
        ## fails at this scale with "sendMaster: ignoring SIGPIPE
        ## signal" -- the fork-after-rlmer state interacts badly with
        ## something native (likely BLAS thread pool or an Rcpp module
        ## slot in the rlmer predictor). PSOCK clusters use persistent
        ## R worker processes spawned via Rscript over sockets, so
        ## they sidestep the fork issue entirely. Setup cost is ~5s
        ## per worker; per-rep work is ~100s, so the relative overhead
        ## is well under 1%.
        cl <- parallel::makeCluster(n_cores, type = "PSOCK")
        on.exit(parallel::stopCluster(cl), add = TRUE)
        ## The worker R processes need (a) robustlmm + lme4 loaded
        ## and (b) the unexported helpers (.phaseE_rep, .quasi_deviance,
        ## ...) sourced. inferenceStudy.R itself carries those helpers,
        ## so each worker sources it from system.file().
        script_dir <- system.file("simulationStudy",
                                   package = "robustlmm")
        if (!nzchar(script_dir))
            script_dir <- normalizePath("inst/simulationStudy",
                                          mustWork = FALSE)
        parallel::clusterExport(cl, varlist = "script_dir",
                                 envir = environment())
        parallel::clusterEvalQ(cl, {
            suppressMessages({
                library(robustlmm)
                library(lme4)
            })
            source(file.path(script_dir, "randomNumberGenerators.R"))
            source(file.path(script_dir, "contaminationStrategies.R"))
            source(file.path(script_dir, "inferenceStudy.R"))
        })
        ## rep_fun, plus the free variables in its body that the
        ## worker can't resolve via library() + source(): `cell` (the
        ## per-cell dataset list) and `nsim_boot` (the bootstrap
        ## depth). When rep_fun's enclosing env is GlobalEnv (it is
        ## here, because the function was constructed inside
        ## .run_phase_cell which itself was sourced at script top
        ## level), parLapply doesn't ship the closure environment; the
        ## worker has to look up free vars in its own GlobalEnv. Going
        ## via a tiny passthrough fixes the lookup: rep_fun lives in
        ## the worker's GlobalEnv (via clusterExport), and the wrapper
        ## just calls it.
        parallel::clusterExport(cl,
            varlist = c("cell", "nsim_boot", "rep_fun"),
            envir = environment())
        recs <- parallel::parLapply(cl, seq_len(n_reps),
                                     function(r) rep_fun(r))
    } else if (n_cores > 1L && .Platform$OS.type != "windows") {
        recs <- parallel::mclapply(seq_len(n_reps), rep_fun,
                                    mc.cores = n_cores,
                                    mc.preschedule = TRUE)
    } else {
        recs <- lapply(seq_len(n_reps), rep_fun)
    }
    n_fail <- sum(vapply(recs, is.null, logical(1)))
    df <- do.call(rbind, recs[!vapply(recs, is.null, logical(1))])
    if (!is.null(df) && nrow(df) > 0L) {
        df$phase         <- phase
        df$J             <- J
        df$contamination <- contamination
        df$eff_size      <- eff_size
    }
    list(records = df, n_fail = n_fail)
}

## --------------------------------------------------------------
## Aggregators                                                    #
## --------------------------------------------------------------

.aggregate <- function(records_by_phase) {
    cov_A <- type1_C <- power_C <- type1_D <- power_D <-
        type1_E <- power_E <- NULL

    a <- records_by_phase$A
    if (!is.null(a) && nrow(a) > 0L)
        cov_A <- aggregate(
            cbind(covered = a$covered, length = a$length) ~
                J + contamination + method + coef,
            data = a, FUN = mean, na.rm = TRUE)

    c_all <- records_by_phase$C
    if (!is.null(c_all) && nrow(c_all) > 0L) {
        ## Restrict the Type-I cells to the Var1 terms (the ones the
        ## generator zeroed) -- the intercept is always nonzero in the
        ## truth, so testing it tells us nothing about Type I.
        c_t1 <- c_all[c_all$eff_size == 0 & c_all$term != "(Intercept)",
                       , drop = FALSE]
        if (nrow(c_t1) > 0L)
            type1_C <- aggregate(
                cbind(rej05 = c_t1$p < 0.05, rej01 = c_t1$p < 0.01) ~
                    J + contamination + method + term,
                data = c_t1, FUN = mean, na.rm = TRUE)
        c_pow <- c_all[c_all$eff_size > 0 & c_all$term != "(Intercept)",
                        , drop = FALSE]
        if (nrow(c_pow) > 0L)
            power_C <- aggregate(
                cbind(rej05 = c_pow$p < 0.05, rej01 = c_pow$p < 0.01) ~
                    J + contamination + eff_size + method + term,
                data = c_pow, FUN = mean, na.rm = TRUE)
    }

    d_all <- records_by_phase$D
    if (!is.null(d_all) && nrow(d_all) > 0L) {
        d_t1 <- d_all[d_all$eff_size == 0, , drop = FALSE]
        if (nrow(d_t1) > 0L)
            type1_D <- aggregate(
                cbind(rej05 = d_t1$p < 0.05, rej01 = d_t1$p < 0.01) ~
                    J + contamination + method,
                data = d_t1, FUN = mean, na.rm = TRUE)
        d_pow <- d_all[d_all$eff_size > 0, , drop = FALSE]
        if (nrow(d_pow) > 0L)
            power_D <- aggregate(
                cbind(rej05 = d_pow$p < 0.05, rej01 = d_pow$p < 0.01) ~
                    J + contamination + eff_size + method,
                data = d_pow, FUN = mean, na.rm = TRUE)
    }

    e_all <- records_by_phase$E
    if (!is.null(e_all) && nrow(e_all) > 0L) {
        e_t1 <- e_all[e_all$eff_size == 0, , drop = FALSE]
        if (nrow(e_t1) > 0L)
            type1_E <- aggregate(
                cbind(rej05 = e_t1$p < 0.05, rej01 = e_t1$p < 0.01) ~
                    J + contamination + method,
                data = e_t1, FUN = mean, na.rm = TRUE)
        e_pow <- e_all[e_all$eff_size > 0, , drop = FALSE]
        if (nrow(e_pow) > 0L)
            power_E <- aggregate(
                cbind(rej05 = e_pow$p < 0.05, rej01 = e_pow$p < 0.01) ~
                    J + contamination + eff_size + method,
                data = e_pow, FUN = mean, na.rm = TRUE)
    }

    type1_F <- power_F <- NULL
    f_all <- records_by_phase$F
    if (!is.null(f_all) && nrow(f_all) > 0L) {
        f_t1 <- f_all[f_all$eff_size == 0, , drop = FALSE]
        if (nrow(f_t1) > 0L)
            type1_F <- aggregate(
                cbind(rej05 = f_t1$p < 0.05, rej01 = f_t1$p < 0.01) ~
                    J + contamination + method,
                data = f_t1, FUN = mean, na.rm = TRUE)
        f_pow <- f_all[f_all$eff_size > 0, , drop = FALSE]
        if (nrow(f_pow) > 0L)
            power_F <- aggregate(
                cbind(rej05 = f_pow$p < 0.05, rej01 = f_pow$p < 0.01) ~
                    J + contamination + eff_size + method,
                data = f_pow, FUN = mean, na.rm = TRUE)
    }

    ## Phases G and H share Phase F's record layout, so one helper
    ## covers both.
    type1_G <- power_G <- type1_H <- power_H <- NULL
    .agg_boot_phase <- function(all) {
        if (is.null(all) || nrow(all) == 0L) return(list(NULL, NULL))
        t1 <- all[all$eff_size == 0, , drop = FALSE]
        pw <- all[all$eff_size > 0, , drop = FALSE]
        list(
            if (nrow(t1) > 0L)
                aggregate(
                    cbind(rej05 = t1$p < 0.05, rej01 = t1$p < 0.01) ~
                        J + contamination + method,
                    data = t1, FUN = mean, na.rm = TRUE),
            if (nrow(pw) > 0L)
                aggregate(
                    cbind(rej05 = pw$p < 0.05, rej01 = pw$p < 0.01) ~
                        J + contamination + eff_size + method,
                    data = pw, FUN = mean, na.rm = TRUE))
    }
    g_agg <- .agg_boot_phase(records_by_phase$G)
    type1_G <- g_agg[[1L]]; power_G <- g_agg[[2L]]
    h_agg <- .agg_boot_phase(records_by_phase$H)
    type1_H <- h_agg[[1L]]; power_H <- h_agg[[2L]]

    list(coverage_A = cov_A,
         type1_C    = type1_C, power_C = power_C,
         type1_D    = type1_D, power_D = power_D,
         type1_E    = type1_E, power_E = power_E,
         type1_F    = type1_F, power_F = power_F,
         type1_G    = type1_G, power_G = power_G,
         type1_H    = type1_H, power_H = power_H)
}

## --------------------------------------------------------------
## Top-level driver                                              #
## --------------------------------------------------------------

##' Run phases A + C + D of the WS6 inference study.
##'
##' @param phases character subset of {"A", "C", "D", "E", "F", "G", "H"}
##' @param J_values integer vector of cluster counts to sweep
##' @param n_reps Monte-Carlo replicates per cell
##' @param contamination subset of DEFAULTS$contamination
##' @param eff_sizes numeric vector; 0 = Type-I cell, >0 = power cell
##'   (only used by phases C/D)
##' @param n_cores parallel::mclapply mc.cores
##' @param master_seed deterministic seed for the whole sweep
##' @param outfile optional .rds path to save results
##' @return a list with $results (per-rep data.frame), $aggregates
##'   (see .aggregate), $provenance.
runInferenceStudy <- function(phases     = c("A", "C", "D"),
                              J_values   = DEFAULTS$J_values,
                              n_reps     = DEFAULTS$n_reps,
                              contamination = DEFAULTS$contamination,
                              eff_sizes  = c(0, DEFAULTS$eff_size),
                              n_cores    = DEFAULTS$n_cores,
                              n_per_subj = DEFAULTS$n_per_subj,
                              n_levels   = DEFAULTS$n_levels,
                              nsim_boot  = 200L,
                              master_seed = DEFAULTS$seed,
                              outfile    = NULL,
                              verbose    = TRUE) {

    t_start <- Sys.time()
    provenance <- list(
        R_version         = R.version.string,
        robustlmm_version = as.character(packageVersion("robustlmm")),
        robustlmm_sha     = tryCatch(packageDescription("robustlmm")$GithubSHA1,
                                       error = function(e) NA_character_),
        phases            = phases,
        J_values          = J_values,
        n_reps            = n_reps,
        contamination     = contamination,
        eff_sizes         = eff_sizes,
        n_per_subj        = n_per_subj,
        n_levels          = n_levels,
        nsim_boot         = nsim_boot,
        master_seed       = master_seed,
        n_cores           = n_cores,
        started           = format(t_start))
    if (verbose) {
        cat("=== robustlmm inference simulation ===\n")
        cat(sprintf("R: %s | robustlmm %s\n",
                    provenance$R_version, provenance$robustlmm_version))
        cat(sprintf("phases: %s | J: %s | n_reps: %d | n_cores: %d\n",
                    paste(phases, collapse = ","),
                    paste(J_values, collapse = ","),
                    n_reps, n_cores))
    }

    set.seed(master_seed)
    records_by_phase <- list()
    cell_id <- 0L

    for (phase in phases) {
        ## Phase E / F use a longitudinal generator (.cell_datasets_long);
        ## phase E uses .cell_datasets (one-way ANOVA). Either way the
        ## cell builder is dispatched inside .run_phase_cell.
        ##
        ## Phase A is a Type-I-of-coverage problem -> eff_size irrelevant.
        es_grid <- if (phase == "A") 0 else eff_sizes
        phase_records <- list()
        for (J in J_values) {
            for (contam in contamination) {
                for (es in es_grid) {
                    cell_id <- cell_id + 1L
                    cell_seed <- master_seed + cell_id
                    cell_res <- .run_phase_cell(
                        phase = phase, J = J, contamination = contam,
                        eff_size = es, n_reps = n_reps,
                        n_per_subj = n_per_subj, n_levels = n_levels,
                        nsim_boot = nsim_boot,
                        cell_seed = cell_seed, n_cores = n_cores,
                        verbose = verbose)
                    if (!is.null(cell_res$records))
                        phase_records[[length(phase_records) + 1L]] <- cell_res$records
                    if (verbose && cell_res$n_fail > 0L)
                        cat(sprintf("    fit failures: %d\n", cell_res$n_fail))
                }
            }
        }
        if (length(phase_records) > 0L)
            records_by_phase[[phase]] <- do.call(rbind, phase_records)
    }

    aggregates <- .aggregate(records_by_phase)
    results <- records_by_phase
    runtime <- as.numeric(Sys.time() - t_start, units = "secs")
    provenance$finished <- format(Sys.time())
    provenance$runtime_secs <- runtime

    out <- list(results = results,
                aggregates = aggregates,
                provenance = provenance)
    if (!is.null(outfile)) {
        saveRDS(out, outfile)
        if (verbose) cat(sprintf("saved %s (runtime %.0fs)\n", outfile, runtime))
    }
    out
}

## --------------------------------------------------------------
## Top-level Rscript invocation: small-config smoke run by default.
## --------------------------------------------------------------
if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    cfg  <- list(phases   = c("A", "C", "D"),
                 J_values = c(8L, 18L),
                 n_reps   = 25L,
                 outfile  = sprintf("inferenceStudy_smoke_%s.rds",
                                    format(Sys.Date())))
    for (a in args) {
        kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
        if (length(kv) == 2L) {
            key <- kv[1L]; val <- kv[2L]
            cfg[[key]] <- if (key %in% c("J_values"))
                              as.integer(strsplit(val, ",")[[1L]])
                          else if (key %in% c("n_reps", "n_cores",
                                              "n_per_subj", "n_levels",
                                              "nsim_boot",
                                              "master_seed"))
                              as.integer(val)
                          else if (key %in% c("eff_sizes"))
                              as.numeric(strsplit(val, ",")[[1L]])
                          else strsplit(val, ",")[[1L]]
        }
    }
    do.call(runInferenceStudy, cfg)
}
