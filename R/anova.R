## anova.rlmerMod: per-term Wald table (single fit), pairwise Wald
## restriction test (fixed-effects-only nested models), and a
## parametric-bootstrap quasi-deviance test (variance-component
## comparisons, where Wald/score asymptotics are invalid on the PSD-cone
## boundary -- see paper2-boundary-singularity).

## Robust quasi-deviance: 2 sum rho_e(r_i / scale) at the conditional
## residuals of `fit`, evaluated at an externally supplied `scale`.
## With rho_e(z) = z^2 / 2 (the classical case) this reduces to sum
## r_i^2 / scale^2.
##
## For the nested-model discrepancy QD(fit0; sigma) - QD(fit1; sigma)
## both fits MUST use the same scale, otherwise the difference is
## contaminated by sigma_0 != sigma_1 and can flip sign even when the
## larger model genuinely fits better (Cantoni and Ronchetti 2001
## evaluate the quasi-deviance at a common scale for exactly this
## reason). The convention here is to use sigma_hat from fit0 (the
## null) throughout.
.quasi_deviance <- function(fit, scale = .sigma(fit)) {
    r_std <- residuals(fit) / scale
    2 * sum(fit@rho.e@rho(r_std))
}

## Two rlmer fits differ in their random-effects structure iff their
## spherical-RE vectors have different lengths, their theta vectors
## have different lengths, or the random-only sub-formulas are not
## identical after normalisation.
.differ_in_re <- function(fit0, fit1) {
    if (fit0@pp$q != fit1@pp$q) return(TRUE)
    if (length(getME(fit0, "theta")) != length(getME(fit1, "theta")))
        return(TRUE)
    rf0 <- formula(fit0, random.only = TRUE)
    rf1 <- formula(fit1, random.only = TRUE)
    !isTRUE(all.equal(rf0, rf1))
}

## Are the fixed-effect column-name sets nested? Returns the extra
## column names in fit1 (if fit0 nests in fit1).
.extra_fixed_cols <- function(fit0, fit1) {
    nms0 <- colnames(fit0@pp$X)
    nms1 <- colnames(fit1@pp$X)
    list(extra = setdiff(nms1, nms0),
         nested = all(nms0 %in% nms1))
}

## Per-term robust Wald test from a single rlmer fit. ddf = "none"
## (default) gives the chi-square table; ddf = "satterthwaite" gives an
## F-test with a Satterthwaite denominator df per term (default vcov).
.anova_single_wald <- function(fit, vcov_type = "default",
                               ddf = "none", IF = NULL) {
    mm     <- model.matrix(fit)
    assign <- attr(mm, "assign")
    tlabs  <- attr(terms(fit), "term.labels")
    hasInt <- attr(terms(fit), "intercept") == 1L

    V  <- as.matrix(vcov(fit, type = vcov_type))
    bh <- .fixef(fit)
    names(bh) <- colnames(mm)

    if (hasInt) {
        unique_aids <- c(0L, seq_along(tlabs))
        term_names  <- c("(Intercept)", tlabs)
    } else {
        unique_aids <- seq_along(tlabs)
        term_names  <- tlabs
    }
    p     <- ncol(mm)
    Df    <- integer(length(unique_aids))
    Chisq <- numeric(length(unique_aids))
    DenDF <- rep(NA_real_, length(unique_aids))
    for (k in seq_along(unique_aids)) {
        idx <- which(assign == unique_aids[k])
        Df[k] <- length(idx)
        Vk <- V[idx, idx, drop = FALSE]
        bk <- bh[idx]
        Chisq[k] <- as.numeric(crossprod(bk, solve(Vk, bk)))
        if (ddf == "satterthwaite") {
            L <- matrix(0, length(idx), p)
            L[cbind(seq_along(idx), idx)] <- 1
            fd <- tryCatch(.satterthwaite_Fdenom(fit, L, IF = IF),
                           error = function(e) list(df2 = NA_real_))
            DenDF[k] <- fd$df2
        }
    }
    if (ddf == "satterthwaite") {
        Fval <- Chisq / Df
        Pr   <- stats::pf(Fval, Df, DenDF, lower.tail = FALSE)
        out  <- data.frame(NumDF = Df, DenDF = DenDF,
                           `F value` = Fval, `Pr(>F)` = Pr,
                           check.names = FALSE)
        rownames(out) <- term_names
        return(structure(out,
            heading = c(
                "Robust Wald F-test for fixed effects (rlmer)",
                "  denominator df: Satterthwaite (IF-based), vcov_type = \"default\""),
            class = c("anova", "data.frame")))
    }
    Pr <- stats::pchisq(Chisq, Df, lower.tail = FALSE)
    out <- data.frame(Df = Df, Chisq = Chisq,
                      `Pr(>Chisq)` = Pr, check.names = FALSE)
    rownames(out) <- term_names
    structure(out,
              heading = c("Robust Wald test for fixed effects (rlmer)",
                          sprintf("  vcov_type = \"%s\"", vcov_type)),
              class = c("anova", "data.frame"))
}

## Robust Wald restriction test for two nested rlmer fits that differ
## only in the fixed-effects design.
.anova_pair_wald <- function(fit0, fit1, vcov_type = "default",
                             ddf = "none", IF = NULL) {
    extras <- .extra_fixed_cols(fit0, fit1)
    if (!extras$nested)
        stop("Models are not nested: fit0's fixed-effect columns ",
             "must be a subset of fit1's.", call. = FALSE)
    if (length(extras$extra) == 0L)
        stop("Models share the same fixed-effect design; nothing to ",
             "test.", call. = FALSE)
    V1   <- as.matrix(vcov(fit1, type = vcov_type))
    bh1  <- .fixef(fit1)
    nms1 <- colnames(fit1@pp$X)
    idx  <- match(extras$extra, nms1)
    Vk   <- V1[idx, idx, drop = FALSE]
    bk   <- bh1[idx]
    Chi  <- as.numeric(crossprod(bk, solve(Vk, bk)))
    Df   <- length(idx)
    p0 <- fit0@pp$p; p1 <- fit1@pp$p
    if (ddf == "satterthwaite") {
        L <- matrix(0, length(idx), p1)
        L[cbind(seq_along(idx), idx)] <- 1
        fd <- tryCatch(.satterthwaite_Fdenom(fit1, L, IF = IF),
                       error = function(e) list(df2 = NA_real_))
        Fval <- Chi / Df
        Pr   <- stats::pf(Fval, Df, fd$df2, lower.tail = FALSE)
        out  <- data.frame(npar = c(p0, p1),
                           NumDF = c(NA_integer_, Df),
                           DenDF = c(NA_real_, fd$df2),
                           `F value` = c(NA_real_, Fval),
                           `Pr(>F)` = c(NA_real_, Pr),
                           check.names = FALSE)
        rownames(out) <- c("fit0", "fit1")
        return(structure(out,
            heading = c(
                "Robust Wald restriction F-test (rlmer)",
                sprintf("  H0: %s = 0",
                        paste(extras$extra, collapse = " = ")),
                "  denominator df: Satterthwaite (IF-based), vcov_type = \"default\""),
            class = c("anova", "data.frame")))
    }
    Pr   <- stats::pchisq(Chi, Df, lower.tail = FALSE)
    out <- data.frame(npar  = c(p0, p1),
                      Df    = c(NA_integer_, Df),
                      Chisq = c(NA_real_, Chi),
                      `Pr(>Chisq)` = c(NA_real_, Pr),
                      check.names = FALSE)
    rownames(out) <- c("fit0", "fit1")
    structure(out,
              heading = c(
                  "Robust Wald restriction test (rlmer)",
                  sprintf("  H0: %s = 0",
                          paste(extras$extra, collapse = " = ")),
                  sprintf("  vcov_type = \"%s\"", vcov_type)),
              class = c("anova", "data.frame"))
}

## Parametric bootstrap of the quasi-deviance difference D = QD(fit0)
## - QD(fit1) under H_0: fit0 is the data-generating model. Simulates
## from the fitted central LMM (Z U_b epsilon_b + sigma epsilon) at
## fit0's parameter estimates, refits BOTH rlmer models on each
## simulated dataset, and returns the empirical p-value for the
## observed discrepancy.
##
## The bootstrap is required for variance-component tests (Wald/score
## asymptotics are invalid on the PSD-cone boundary; see Self-Liang
## 1987 and the paper2-boundary-singularity discussion). It also runs
## for fixed-effects-only nested comparisons when test = "boot" is
## explicit.
## Smallest robustness weight per random-effects group. The robust fit
## downweights the spherical random effects through rho.b, so a group
## whose random effect is heavily downweighted is one the fit treated as
## a group outlier. getME(., "w_b_vector") is aligned with the spherical
## random effects (b-index order, see wgt.b in helpers.R); each b-index
## is mapped to its "factor: level" through the Gp/cnms/flist layout
## (term i owns b-indices Gp[i]+1 .. Gp[i+1], level-major: the
## length(cnms[[i]]) coordinates of one level are consecutive). Minima
## are aggregated over all b's carrying the same label, so several terms
## sharing a grouping factor contribute to one entry per level. Returns
## a named numeric ("factor: level" -> smallest weight over that level's
## random effects); attr "labeled" is FALSE (and names dropped) if the
## layout cannot be resolved, in which case the raw per-coefficient
## weights are returned so the caller can still screen on the minimum.
.re_min_weight_by_group <- function(fit) {
    w    <- getME(fit, "w_b_vector")
    cnms <- fit@cnms
    fl   <- fit@flist
    asg  <- attr(fl, "assign")
    Gp   <- fit@Gp
    unlabeled <- structure(as.numeric(w), labeled = FALSE)
    if (is.null(cnms) || is.null(fl) || is.null(asg) ||
        length(asg) != length(cnms) ||
        length(Gp) != length(cnms) + 1L ||
        Gp[length(Gp)] != length(w))
        return(unlabeled)
    ## explicit b-index -> "factor: level" labels, term by term
    lab <- character(length(w))
    for (i in seq_along(cnms)) {
        nc <- length(cnms[[i]])
        fi <- asg[i]
        levs <- levels(fl[[fi]]); nl <- length(levs)
        if (Gp[i + 1L] - Gp[i] != nc * nl)
            return(unlabeled)
        lab[Gp[i] + seq_len(nc * nl)] <-
            paste0(names(fl)[fi], ": ", rep(levs, each = nc))
    }
    if (any(!nzchar(lab)))
        return(unlabeled)
    grp <- factor(lab, levels = unique(lab))
    res <- vapply(split(as.numeric(w), grp), min, numeric(1))
    structure(res, labeled = TRUE)
}

## Auto-guard for the bootstrap variance-component test (C4): the
## parametric bootstrap simulates from the central LMM at fit0's
## estimates, so it is anti-conservative when a group is contaminated and
## biases those estimates. The robust fit's own RE weights detect this
## directly -- a heavily downweighted group is a group outlier -- whereas
## cluster Cook's distance does not (the downweighting that removes the
## bias also removes the influence the distance would measure). Warns,
## naming the offending groups, when any group's minimum RE weight falls
## below `threshold` (default 0.5).
.bootstrap_group_guard <- function(fit0, threshold = 0.5) {
    mw <- tryCatch(.re_min_weight_by_group(fit0), error = function(e) NULL)
    if (is.null(mw) || !length(mw)) return(invisible(NULL))
    flagged <- mw[is.finite(mw) & mw < threshold]
    if (!length(flagged)) return(invisible(NULL))
    flagged <- flagged[order(flagged)]
    if (isTRUE(attr(mw, "labeled"))) {
        shown <- utils::head(flagged, 5L)
        detail <- paste0(names(shown), " (min w = ",
                         sprintf("%.2f", shown), ")", collapse = "; ")
        if (length(flagged) > length(shown))
            detail <- paste0(detail, "; and ",
                             length(flagged) - length(shown), " more")
        lead <- paste0("the robust fit heavily downweights random-effect ",
                       "group(s): ", detail)
    } else {
        lead <- sprintf(paste0("the robust fit heavily downweights %d ",
                               "random effect(s) (smallest weight %.2f)"),
                        length(flagged), min(flagged))
    }
    warning("anova(test = \"boot\"): ", lead, ". The parametric-bootstrap ",
            "variance-component test simulates from the central model at ",
            "the null fit's estimates, so it can be anti-conservative when ",
            "groups are contaminated; treat the p-value with caution and ",
            "inspect these groups (see ?anova.rlmerMod).", call. = FALSE)
    invisible(flagged)
}

## Default cleaning threshold for the experimental robust null (see
## .anova_robust_null_params). A cluster is flagged when its smallest
## random-effect robustness weight falls below this value; the value 0.65
## was tuned by simulation (power preserved, contamination-induced Type-I
## inflation reduced at larger J). Overridable, but not part of the public
## interface.
.anova_robust_null_threshold <- function()
    getOption("robustlmm.anova.robust_null_threshold", 0.65)

## EXPERIMENTAL contamination-cleaned null for the bootstrap
## variance-component test. The plain parametric bootstrap simulates from
## fit0's robust estimates, so any bias the contaminated clusters induced
## in (sigma0, theta0, beta0) propagates into a too-narrow bootstrap null
## and the test becomes anti-conservative. This helper de-biases the
## *generating* parameters: it flags the clusters the robust fit heavily
## downweights (smallest per-cluster random-effect weight < threshold,
## reusing .re_min_weight_by_group), refits fit0 with those clusters
## removed, and returns the cleaned generating parameters together with
## the generator's U_b = Lambda(theta) rebuilt at the cleaned theta.
## D_obs is still computed from the original untrimmed fits, so the test
## uses all the data; only the bootstrap null GENERATION is cleaned.
##
## Guard / fallback: if no cluster is flagged, the weight layout cannot
## be resolved, more than 50% of clusters are flagged, or the trimmed
## refit fails, the untrimmed fit0 parameters are returned (cleaned =
## FALSE, with the fallback cause in $reason) so the call degrades to
## the plain parametric bootstrap. A single warning is emitted on the
## >50%/failure fallback.
.anova_robust_null_params <- function(fit0,
                                      threshold = .anova_robust_null_threshold()) {
    sig0    <- .sigma(fit0)
    th_orig <- getME(fit0, "theta")
    beta0   <- .fixef(fit0)
    Ub0     <- fit0@pp$U_b
    untrimmed <- function(reason)
        list(sigma = sig0, theta = th_orig, beta = beta0,
             Ub = Ub0, flagged = character(0), n_flagged = 0L,
             cleaned = FALSE, reason = reason)

    mw <- tryCatch(.re_min_weight_by_group(fit0), error = function(e) NULL)
    if (is.null(mw) || !length(mw) || !isTRUE(attr(mw, "labeled")))
        return(untrimmed("layout-unresolved"))
    flagged_names <- names(mw)[is.finite(mw) & mw < threshold]
    n_flagged <- length(flagged_names)
    if (n_flagged == 0L) return(untrimmed("none-flagged"))

    ## Map the labelled "<factor>: <level>" entries back to (factor, level)
    ## and build a subset that keeps every observation NOT in a flagged
    ## cluster. The labels are produced by .re_min_weight_by_group as
    ## paste0(factor, ": ", level).
    fl <- fit0@flist
    keep <- rep(TRUE, nrow(fit0@frame))
    n_clusters <- 0L
    for (fi in seq_along(fl)) {
        fac  <- names(fl)[fi]
        levs <- levels(fl[[fi]])
        n_clusters <- n_clusters + length(levs)
        prefix <- paste0(fac, ": ")
        hit <- flagged_names[startsWith(flagged_names, prefix)]
        if (!length(hit)) next
        bad_levs <- substring(hit, nchar(prefix) + 1L)
        keep <- keep & !(as.character(fl[[fi]]) %in% bad_levs)
    }
    ## >50% of clusters flagged -> do not clean.
    if (n_clusters == 0L || n_flagged > floor(0.5 * n_clusters)) {
        warning("anova(test = \"boot\", null = \"robust\"): ", n_flagged,
                " of ", n_clusters, " clusters flagged (> 50%); falling ",
                "back to the plain parametric null (no cleaning).",
                call. = FALSE)
        return(untrimmed("over-half-flagged"))
    }
    if (!any(keep)) return(untrimmed("no-observations-left"))

    ## Nonsingular-subsampling check on the trimmed design (Koller and
    ## Stahel 2017; the same test the RANSAC subsampler uses). Removing the
    ## flagged (contaminated) clusters must not leave the fixed-effects
    ## design rank-deficient, drop a fixed-factor contrast level, collapse a
    ## grouping factor below two levels, or make a random-slope covariate
    ## constant. Unlike the RANSAC repair we cannot add units back -- the
    ## excluded clusters are exactly the ones we are deliberately removing --
    ## so a deficient trim falls back to the plain parametric null. This
    ## catches near-singular trims (e.g. a contrast column that becomes
    ## near-constant) that the post-refit parameter-dimension check below
    ## would miss.
    rank_info <- tryCatch(.ransac_design_info(formula(fit0), fit0@frame),
                          error = function(e) NULL)
    if (!is.null(rank_info)) {
        chk <- tryCatch(.ransac_check_subsample(which(keep), fit0@frame,
                                                rank_info),
                        error = function(e) list(ok = TRUE))
        if (!isTRUE(chk$ok)) {
            warning("anova(test = \"boot\", null = \"robust\"): trimming the ",
                    "flagged clusters leaves the design rank-deficient ",
                    "(dropped contrast level, collapsed grouping factor, or ",
                    "constant random-slope covariate); falling back to the ",
                    "plain parametric null (no cleaning).", call. = FALSE)
            return(untrimmed("rank-deficient"))
        }
    }

    data_keep <- fit0@frame[keep, , drop = FALSE]
    for (fi in seq_along(fl)) {
        fac <- names(fl)[fi]
        if (is.factor(data_keep[[fac]]))
            data_keep[[fac]] <- droplevels(data_keep[[fac]])
    }
    fc <- .boot_refit(fit0, data_keep, th_orig)
    if (is.null(fc)) {
        warning("anova(test = \"boot\", null = \"robust\"): the trimmed ",
                "null refit failed; falling back to the plain parametric ",
                "null (no cleaning).", call. = FALSE)
        return(untrimmed("refit-failed"))
    }

    sig_c  <- .sigma(fc)
    th_c   <- getME(fc, "theta")
    beta_c <- .fixef(fc)
    ## The trimmed refit must be parameter-compatible with fit0 (same
    ## theta / beta dimensions); a level dropped from a fixed-effect
    ## factor, say, makes the cleaned parameters unusable as generating
    ## parameters for the full design.
    if (length(th_c) != length(th_orig) || anyNA(th_c) ||
        length(beta_c) != length(beta0) || anyNA(beta_c)) {
        warning("anova(test = \"boot\", null = \"robust\"): the trimmed ",
                "null refit is not parameter-compatible with the full ",
                "fit; falling back to the plain parametric null (no ",
                "cleaning).", call. = FALSE)
        return(untrimmed("refit-failed"))
    }
    ## Build the generator's U_b = Lambda(theta) directly at the cleaned
    ## theta: robustlmm's internal (Cholesky-space) theta maps onto the
    ## non-zeros of U_b through pp$Lind (exactly what rlmerPredD$setTheta
    ## does), so this is exact for any number of variance components.
    ## For a single component it coincides with scaling fit0's U_b by
    ## theta_clean / theta_orig.
    Ub_c   <- Ub0
    Ub_c@x <- th_c[fit0@pp$Lind]

    list(sigma = sig_c, theta = th_c, beta = beta_c, Ub = Ub_c,
         flagged = flagged_names, n_flagged = n_flagged, cleaned = TRUE,
         reason = "cleaned")
}

## Shared machinery of the parametric-bootstrap paths (test = "boot"
## and test = "score"): simulate one dataset from the null generating
## parameters, y* = X beta0 + Z Lambda(theta0) (sigma0 z_b) + sigma0
## eps. Consumes exactly rnorm(q0) followed by rnorm(n), so extracting
## this helper leaves the callers' RNG streams (and hence seeded
## p-values) bit-identical.
.boot_sim_data <- function(data_obs, yname, Xb0, Z0, Ub0, sig0, n, q0) {
    b_sim <- sig0 * as.numeric(Ub0 %*% rnorm(q0))
    eps   <- rnorm(n, sd = sig0)
    d_sim <- data_obs
    d_sim[[yname]] <- Xb0 + as.numeric(Z0 %*% b_sim) + eps
    d_sim
}

## Quiet bootstrap refit at a warm start; NULL on error.
.boot_refit <- function(fit, d_sim, start)
    tryCatch(suppressMessages(suppressWarnings(
        update(fit, data = d_sim, start = start))),
        error = function(e) NULL)

.anova_pair_boot <- function(fit0, fit1, nsim = 1000L, seed = NULL,
                             verbose = FALSE, robust.null = FALSE) {
    if (!is.null(seed)) set.seed(seed)
    sig0_obs <- .sigma(fit0)
    D_obs <- .quasi_deviance(fit0, scale = sig0_obs) -
             .quasi_deviance(fit1, scale = sig0_obs)

    pp0    <- fit0@pp
    th0    <- getME(fit0, "theta")
    X      <- pp0$X
    Z0     <- t(pp0$Zt)
    n      <- pp0$n; q0 <- pp0$q
    yname  <- as.character(formula(fit0)[[2L]])
    data_obs <- fit0@frame
    th1    <- getME(fit1, "theta")

    ## Generating parameters for the bootstrap null. Plain parametric:
    ## fit0's raw robust estimates. Robust null: a contamination-cleaned
    ## null fit (D_obs above is unchanged -- the test still uses all data).
    rn_info <- NULL
    if (isTRUE(robust.null)) {
        rn_info <- .anova_robust_null_params(fit0)
        sig0  <- rn_info$sigma
        beta0 <- rn_info$beta
        Ub0   <- rn_info$Ub
    } else {
        sig0  <- sig0_obs
        beta0 <- .fixef(fit0)
        Ub0   <- pp0$U_b
    }
    Xb0    <- as.numeric(X %*% beta0)

    D_boot <- rep(NA_real_, nsim)
    n_fail <- 0L
    t0     <- Sys.time()
    for (r in seq_len(nsim)) {
        d_sim <- .boot_sim_data(data_obs, yname, Xb0, Z0, Ub0, sig0, n, q0)
        f0_s <- .boot_refit(fit0, d_sim, th0)
        f1_s <- .boot_refit(fit1, d_sim, th1)
        if (is.null(f0_s) || is.null(f1_s)) {
            n_fail <- n_fail + 1L; next
        }
        sig0_s <- .sigma(f0_s)
        D_boot[r] <- .quasi_deviance(f0_s, scale = sig0_s) -
                     .quasi_deviance(f1_s, scale = sig0_s)
        if (verbose && r %% 50L == 0L) {
            elapsed <- as.numeric(Sys.time() - t0, units = "secs")
            message(sprintf("  boot rep %d/%d, elapsed %.1fs",
                            r, nsim, elapsed))
        }
    }
    D_clean <- D_boot[!is.na(D_boot)]
    if (length(D_clean) < 10L)
        stop("Bootstrap: only ", length(D_clean), " of ", nsim,
             " replicates succeeded; need >= 10 to compute a p-value.",
             call. = FALSE)
    if (length(D_clean) < 200L)
        warning("Bootstrap p-value based on only ", length(D_clean),
                " effective replicates; consider raising nsim for a ",
                "stable tail estimate.", call. = FALSE)
    pval <- mean(D_clean >= D_obs)
    list(D_obs = D_obs, D_boot = D_clean,
         p.value = pval, n_fail = n_fail, nsim = nsim,
         robust.null = rn_info)
}

## ------------------------------------------------------------------
## EXPERIMENTAL robust variance-component SCORE test (test = "score").
##
## One-sided score test for a single added (scalar, uncorrelated)
## variance component, computed from the robust NULL fit only and
## calibrated by a parametric SCORE bootstrap that refits ONLY the
## null model per replicate (~4x cheaper than the deviance bootstrap's
## double refits). Statistic (per cluster j, from the null fit):
##
##   V_j        = sigma0^2 (I + Z0_j Lambda0 Lambda0' Z0_j')
##   rtil_j     = V_j^{-1/2} (y_j - X_j beta0)   (whitened marginal resid,
##                                                symmetric square root)
##   v_j        = V_j^{-1/2} z_add_j             (whitened tested direction;
##                                                z_add = the added RE term's
##                                                design column)
##   s_j        = (v_j' psi(rtil_j))^2 - kappa_1 ||v_j||^2,
##                psi = fit0's rho.e psi, kappa_1 = rho.e@Epsi2()
##   S          = sum_j s_j / sqrt(sum_j (s_j - s_bar)^2)
##
## kappa_1 ||v_j||^2 is the exact E_Phi[(v' psi(W))^2], W ~ N(0, I)
## (cross terms vanish by independence + odd psi), so each s_j is
## Fisher-consistency centered; the self-normalisation makes any
## j-constant factor drop out. Reject for S large positive (extra
## variance in the tested direction inflates the v-weighted quadratic
## form). Validated by simulation (robust-score-test.R phase 1 +
## robust-score-test-phase2.R, robustlmm-open-problems): all
## pre-registered gates passed; see ?anova.rlmerMod for the numbers.
##
## NO null cleaning and NO .bootstrap_group_guard call here, by
## design: the per-cluster contributions are psi-bounded, so the
## bootstrap null law is insensitive to the O(eps)-biased generating
## parameters -- in 600+ paired replicates (including the
## uncentered-x honesty cells) the trim65-cleaned generating
## parameters changed at most 1-2 decisions vs the raw ones, so the
## cleaning composition is validated unnecessary and the C4 guard's
## "poisoned null" warning would be misleading noise for this path.
## ------------------------------------------------------------------

## Scope check for test = "score" (fail loudly, matching the package's
## C3 discipline): both fits rlmerMod on the same data, a single
## grouping factor shared by both, identical fixed-effects design, and
## the alternative adds exactly ONE independent scalar random-effect
## column relative to the null: length(theta1) == length(theta0) + 1
## and the added column is a 1-dimensional block, structurally
## uncorrelated with every other random effect (checked on the sparsity
## pattern of fit1's Lambda). Covers both `(1|g)` vs `(1|g) + (0+x|g)`
## and diagonal structures adding one component, e.g. `(1|g)` vs
## `diag(x|g)`.
##
## Returns list(z_add, label) where z_add is the added term's design
## column per observation, extracted from fit1's Z: the added block
## owns one Z column per cluster and each observation loads on its own
## cluster's column only, so the per-observation design value is the
## row sum over the added block's columns.
.anova_score_scope <- function(fit0, fit1) {
    .score_stop <- function(...)
        stop(..., " -- not supported by test = \"score\"; use ",
             "test = \"boot\".", call. = FALSE)
    if (length(fit0@flist) != 1L || length(fit1@flist) != 1L)
        .score_stop("test = \"score\" requires a single grouping ",
                    "factor in both models; crossed or nested ",
                    "grouping factors are out of scope")
    if (!identical(names(fit0@flist)[1L], names(fit1@flist)[1L]) ||
        !identical(as.character(fit0@flist[[1L]]),
                   as.character(fit1@flist[[1L]])))
        .score_stop("both models must share the same grouping factor ",
                    "(same variable, same observation order)")
    if (nrow(fit0@frame) != nrow(fit1@frame) ||
        !isTRUE(all.equal(as.numeric(fit0@resp$y),
                          as.numeric(fit1@resp$y))))
        .score_stop("both models must be fitted to the same data")
    if (!identical(colnames(fit0@pp$X), colnames(fit1@pp$X)))
        .score_stop("the fixed-effects designs differ; test = ",
                    "\"score\" tests a single added variance ",
                    "component only")
    th0 <- getME(fit0, "theta")
    th1 <- getME(fit1, "theta")
    if (length(th1) != length(th0) + 1L)
        .score_stop("the alternative must add exactly one variance ",
                    "parameter (scalar component) relative to the ",
                    "null; got length(theta) ", length(th0), " vs ",
                    length(th1), ". Multi-component or ",
                    "correlated-slope alternatives are out of scope")
    cn0 <- fit0@cnms
    cn1 <- fit1@cnms
    Gp1 <- fit1@Gp
    nl  <- length(levels(fit1@flist[[1L]]))
    gname <- names(fit1@flist)[1L]

    ## Which of fit1's b-columns belong to the added scalar block?
    add_cols <- NULL
    add_name <- NULL
    if (length(cn1) == length(cn0) + 1L) {
        ## case (a): the alternative has one extra 1-column term,
        ## e.g. (1|g) vs (1|g) + (0 + x|g)
        for (k in seq_along(cn1)) {
            if (length(cn1[[k]]) != 1L) next
            if (identical(unname(cn1[-k]), unname(cn0))) {
                add_cols <- (Gp1[k] + 1L):Gp1[k + 1L]
                add_name <- cn1[[k]]
                break
            }
        }
    } else if (length(cn1) == length(cn0)) {
        ## case (b): one term gains one column under a diagonal
        ## structure, e.g. (1|g) vs diag(x|g). b-layout within a term
        ## is level-major: the length(cnms[[k]]) coordinates of one
        ## level are consecutive.
        deleted_pos <- function(small, big) {
            if (length(big) != length(small) + 1L) return(NA_integer_)
            for (pos in seq_along(big))
                if (identical(big[-pos], small)) return(pos)
            NA_integer_
        }
        differs <- which(!vapply(seq_along(cn1), function(k)
            identical(cn1[[k]], cn0[[k]]), logical(1)))
        if (length(differs) == 1L) {
            k   <- differs
            pos <- deleted_pos(cn0[[k]], cn1[[k]])
            nc1 <- length(cn1[[k]])
            if (!is.na(pos) && Gp1[k + 1L] - Gp1[k] == nc1 * nl) {
                add_cols <- Gp1[k] + (seq_len(nl) - 1L) * nc1 + pos
                add_name <- cn1[[k]][pos]
            }
        }
    }
    if (is.null(add_cols))
        .score_stop("could not identify a single added scalar ",
                    "random-effect column in the alternative model ",
                    "(models must be nested with exactly one added ",
                    "independent component)")

    ## The added block must be structurally uncorrelated with every
    ## other random effect: no structural nonzero of fit1's Lambda
    ## pattern may couple the added b-coordinates to any other
    ## b-coordinate (a correlated-slope alternative would).
    Ub1 <- methods::as(methods::as(fit1@pp$U_b, "generalMatrix"),
                       "TsparseMatrix")
    in_add_i <- (Ub1@i + 1L) %in% add_cols
    in_add_j <- (Ub1@j + 1L) %in% add_cols
    if (any(in_add_i != in_add_j))
        .score_stop("the added random-effect column is correlated ",
                    "with other random effects in the alternative ",
                    "model; only an independent (uncorrelated) added ",
                    "component is in scope")

    Z1 <- t(fit1@pp$Zt)
    z_add <- Matrix::rowSums(Z1[, add_cols, drop = FALSE])
    list(z_add = as.numeric(z_add),
         label = paste0("(", add_name, " | ", gname, ")"))
}

## The score statistic S and the per-cluster contributions s_j, from
## the null fit and the added term's design column. The per-cluster
## marginal covariance is built from the null fit's own estimates:
## with M = Z Lambda(theta0) (= Z %*% pp$U_b), V_j / sigma0^2 =
## I + M_j M_j' restricted to cluster j's rows (exact for a single
## grouping factor, where V is block-diagonal by cluster). The
## symmetric inverse square root is computed by a per-cluster
## eigendecomposition (cluster sizes are small); in the balanced
## exchangeable special case this reproduces the analytic root used in
## the validation study to machine precision.
.anova_score_stat <- function(fit0, z_add) {
    sig  <- .sigma(fit0)
    beta <- .fixef(fit0)
    r    <- as.numeric(fit0@resp$y - fit0@resp$offset -
                       fit0@pp$X %*% beta)
    psi  <- fit0@rho.e@psi
    kap1 <- fit0@rho.e@Epsi2()
    M    <- Matrix::crossprod(fit0@pp$Zt, fit0@pp$U_b)
    idx  <- split(seq_along(r), fit0@flist[[1L]])
    s_j  <- vapply(idx, function(ii) {
        Mi <- as.matrix(M[ii, , drop = FALSE])
        Mi <- Mi[, colSums(abs(Mi)) > 0, drop = FALSE]
        Vi <- diag(length(ii)) + tcrossprod(Mi)   # V_j / sigma^2
        ee <- eigen(Vi, symmetric = TRUE)
        ## symmetric V_j^{-1/2} (eigenvalues >= 1 by construction)
        W  <- ee$vectors %*% (t(ee$vectors) / sqrt(ee$values)) / sig
        rt <- as.numeric(W %*% r[ii])
        v  <- as.numeric(W %*% z_add[ii])
        tj <- sum(v * psi(rt))
        tj^2 - kap1 * sum(v^2)
    }, numeric(1))
    den <- sqrt(sum((s_j - mean(s_j))^2))
    S   <- if (den > 0) sum(s_j) / den else NA_real_
    list(S = S, s_j = s_j)
}

## Parametric SCORE bootstrap: simulate nsim datasets from fit0's
## estimates with the same generation machinery as .anova_pair_boot,
## refit ONLY the null model per replicate, and recompute S*. p-value
## = (1 + #{S* >= S_obs}) / (n_effective + 1), one-sided.
.anova_pair_score <- function(fit0, fit1, nsim = 199L, seed = NULL,
                              verbose = FALSE) {
    scope <- .anova_score_scope(fit0, fit1)
    if (!is.null(seed)) set.seed(seed)
    obs <- .anova_score_stat(fit0, scope$z_add)
    if (!is.finite(obs$S))
        stop("test = \"score\": the observed statistic is degenerate ",
             "(zero variance of the per-cluster contributions).",
             call. = FALSE)

    pp0   <- fit0@pp
    th0   <- getME(fit0, "theta")
    X     <- pp0$X
    Z0    <- t(pp0$Zt)
    Ub0   <- pp0$U_b
    n     <- pp0$n; q0 <- pp0$q
    sig0  <- .sigma(fit0)
    beta0 <- .fixef(fit0)
    yname <- as.character(formula(fit0)[[2L]])
    data_obs <- fit0@frame
    Xb0   <- as.numeric(X %*% beta0)

    S_boot <- rep(NA_real_, nsim)
    t0     <- Sys.time()
    for (r in seq_len(nsim)) {
        d_sim <- .boot_sim_data(data_obs, yname, Xb0, Z0, Ub0, sig0, n, q0)
        f0_s <- .boot_refit(fit0, d_sim, th0)
        if (is.null(f0_s)) next
        S_boot[r] <- tryCatch(.anova_score_stat(f0_s, scope$z_add)$S,
                              error = function(e) NA_real_)
        if (verbose && r %% 50L == 0L) {
            elapsed <- as.numeric(Sys.time() - t0, units = "secs")
            message(sprintf("  score boot rep %d/%d, elapsed %.1fs",
                            r, nsim, elapsed))
        }
    }
    S_clean <- S_boot[!is.na(S_boot)]
    n_fail  <- nsim - length(S_clean)
    if (length(S_clean) < 10L)
        stop("Score bootstrap: only ", length(S_clean), " of ", nsim,
             " replicates succeeded; need >= 10 to compute a p-value.",
             call. = FALSE)
    if (length(S_clean) < 100L)
        warning("Score-bootstrap p-value based on only ",
                length(S_clean), " effective replicates; consider ",
                "raising nsim for a stable tail estimate.",
                call. = FALSE)
    pval <- (1 + sum(S_clean >= obs$S)) / (length(S_clean) + 1)
    list(S_obs = obs$S, s_j = obs$s_j, S_boot = S_clean,
         p.value = pval, n_fail = n_fail, nsim = nsim,
         label = scope$label)
}

##' Analysis of variance for an rlmer fit.
##'
##' Three modes.
##' \describe{
##'   \item{Single fit \code{anova(fit)}}{returns a per-term robust Wald
##'     chi-square table from \code{vcov(fit, type = vcov_type)}. Each
##'     term is tested marginally as
##'     \eqn{T_t = \hat\beta_t^\top \hat V_t^{-1} \hat\beta_t \sim
##'       \chi^2_{k_t}} under \eqn{H_0: \beta_t = 0}.}
##'   \item{Nested fits differing only in fixed effects (default
##'     \code{test = "Wald"})}{robust Wald restriction test on the
##'     extra coefficients. The same \eqn{V} used by the single-fit
##'     table is used here; \code{vcov_type = "sandwich"} carries through.}
##'   \item{Nested fits differing in random-effects structure}{the
##'     Wald and score asymptotics are invalid (boundary problem; see
##'     Self-Liang 1987 and Koller 2026 paper 2). \code{test = "Wald"}
##'     warns and switches to the parametric bootstrap automatically;
##'     \code{test = "boot"} is the default valid path. The bootstrap
##'     simulates \code{nsim} datasets from the fitted central LMM at
##'     \code{fit0}'s estimates, refits both rlmer models per
##'     replicate, and uses the quasi-deviance difference
##'     \eqn{D = 2 \sum [\rho_e(r_i^{(0)}/\hat\sigma_0) -
##'       \rho_e(r_i^{(1)}/\hat\sigma_1)]} as the discrepancy
##'     statistic (Heritier and Ronchetti 1994; Cantoni and Ronchetti
##'     2001; Heritier, Cantoni, Copt and Victoria-Feser 2009 sec. 5).
##'     For the special case of a \emph{single added variance
##'     component} the experimental \code{test = "score"} offers a
##'     contamination-robust alternative (see Details).}
##' }
##'
##' Only pairwise comparison is implemented; chains \code{anova(fit0,
##' fit1, fit2, ...)} are not yet supported.
##'
##' \strong{Small-J caveat for the sandwich.} \code{vcov_type =
##'   "sandwich"} is markedly anti-conservative at small \eqn{J}: in
##' a simulation study the pairwise Wald Type-I rate reached
##' 0.15-0.20 (vs. nominal 0.05) at \eqn{J = 8} and was still elevated
##' at \eqn{J = 18}, returning to nominal only by \eqn{J \gtrsim 50}.
##' \code{vcov_sandwich} emits a warning for \eqn{J < 20}; prefer
##' \code{vcov_type = "default"} for hypothesis tests at small
##' \eqn{J}.
##'
##' \strong{Subject-contamination caveat for the bootstrap
##' (experimental).} The parametric-bootstrap quasi-deviance path
##' (\code{test = "boot"}) is experimental --- it is currently the only
##' exposed variance-component path. It calibrates correctly under clean
##' Gaussian and heavy-tailed errors but is markedly anti-conservative
##' when a subset of subjects is contaminated. In a simulation study with 10\% of subjects
##' shifted by 5 standard errors (\code{shift_subj}), Type-I climbed
##' to 0.13-0.24 (2.5-5x nominal) across \eqn{J \in \{8, 18, 50\}}.
##' The mechanism: the bootstrap simulates from the fitted central
##' LMM at \code{fit0}'s estimates, so any bias \code{fit0}'s
##' \eqn{\hat\sigma_0, \hat\theta_0} picked up from the contaminated
##' subjects propagates into a too-narrow bootstrap null. As an automatic
##' guard, \code{anova(test = "boot")} now inspects the null fit's own
##' random-effects robustness weights and \strong{warns} when a group is
##' heavily downweighted (smallest weight below \eqn{0.5}) --- a direct
##' signal that the fit treated that group as an outlier, so the bootstrap
##' null may be poisoned. (The robust fit absorbs a contaminated group
##' into a downweighted random effect, which is why
##' \code{\link[=cooks.distance.rlmerMod]{cooks.distance}} --- an
##' \emph{influence} measure --- does \emph{not} reliably flag it here:
##' the downweighting that removes the bias also removes the influence.)
##' When the warning fires, treat the p-value with caution.
##'
##' \strong{Experimental robust null (\code{null = "robust"}).} As a
##' mitigation for the subject-contamination anti-conservativeness above,
##' the bootstrap path accepts \code{null = "robust"}: instead of
##' generating the parametric bootstrap from \code{fit0}'s raw estimates,
##' it generates from a \emph{contamination-cleaned} null fit. Clusters the
##' robust fit heavily downweights (smallest per-cluster random-effect
##' weight below a tuned threshold, the same signal the automatic guard
##' uses) are flagged, \code{fit0} is refitted with those clusters removed,
##' and the bootstrap is generated from the resulting de-biased
##' \eqn{(\hat\sigma, \hat\theta, \hat\beta)} (with the generator's
##' \eqn{U_b = \Lambda(\theta)} rebuilt at the full cleaned \eqn{\hat\theta}
##' vector). The observed discrepancy \eqn{D} is still
##' taken from the original untrimmed fits, so the test uses all the data;
##' only the bootstrap null generation is cleaned. If no cluster is
##' flagged, more than half the clusters are flagged, trimming would leave
##' the design rank-deficient (a dropped contrast level, a collapsed
##' grouping factor, or a constant random-slope covariate --- checked by
##' the same nonsingular-subsampling test used by the RANSAC initial
##' estimator), or the trimmed refit fails, the call falls back to the
##' plain parametric null (so on clean data \code{null = "robust"} closely
##' matches \code{null = "parametric"}); the table heading states which
##' case applied. In simulation this reduces the contamination-induced
##' Type-I inflation at larger \eqn{J} while preserving power, at the cost
##' of mild conservatism on clean data. It is \strong{experimental},
##' validated by simulation rather than a finite-sample theorem; the table
##' heading notes when the robust null was applied.
##'
##' \strong{Experimental robust score test (\code{test = "score"}).}
##' For the common special case where the alternative adds exactly
##' \emph{one independent scalar variance component} relative to the
##' null --- e.g. \code{(1|g)} vs \code{(1|g) + (0 + x|g)}, or diagonal
##' structures adding one component --- \code{test = "score"} runs a
##' one-sided robust score test computed from the robust \emph{null}
##' fit only and calibrated by a score-only parametric bootstrap that
##' refits just the null model per replicate (roughly 4x cheaper than
##' \code{test = "boot"}'s double refits; hence the smaller default
##' \code{nsim = 199}). Per cluster \eqn{j}, the whitened marginal
##' residuals \eqn{\tilde r_j = V_j^{-1/2}(y_j - X_j \hat\beta_0)} and
##' the whitened tested direction \eqn{v_j = V_j^{-1/2} z_j} (with
##' \eqn{V_j} the null fit's marginal covariance and \eqn{z_j} the
##' added term's design column) give the bounded contribution
##' \eqn{s_j = (v_j^\top \psi(\tilde r_j))^2 - \kappa_1 \|v_j\|^2}
##' (\eqn{\psi} the fit's \code{rho.e} psi-function, \eqn{\kappa_1 =
##' E[\psi(Z)^2]}); the statistic is the self-normalised cluster sum
##' \eqn{S = \sum_j s_j / \sqrt{\sum_j (s_j - \bar s)^2}}, and the
##' one-sided p-value is \eqn{(1 + \#\{S^* \ge S\})/(n_{\mathrm{eff}} +
##' 1)}. Because each \eqn{s_j} is psi-bounded, whole-cluster
##' contamination shifts the statistic and its bootstrap reference law
##' together instead of inflating the test: in simulation (Gaussian
##' balanced designs, one scalar tested component; not a theorem)
##' contaminated-null Type-I was 0.035 vs 0.115 for \code{test =
##' "boot"} with 10\% of clusters shifted at \eqn{J = 50}, clean-null
##' Type-I 0.045, and power 0.920 vs 0.900 --- with no null cleaning
##' needed (raw vs cleaned generating parameters changed at most 1-2
##' decisions in 600+ paired replicates, so \code{null = "robust"} is
##' ignored and no downweighted-group warning is issued for this
##' path). Across an adversarial sweep the pattern held: \eqn{J = 18}
##' clean/contaminated 0.055/0.015, uncentered-\eqn{x} designs
##' clean/contaminated 0.045/0.075, 5\% single-observation outliers
##' 0.060. \emph{Scope}: both fits on the same data with identical
##' fixed effects and a single shared grouping factor, and the
##' alternative adds exactly one uncorrelated scalar component
##' (\code{length(theta)} differs by 1); anything else ---
##' multi-component or correlated-slope alternatives, crossed or
##' nested factors --- stops with an error, use \code{test = "boot"}
##' there. \emph{Caveat}: contamination aligned with the tested
##' direction is indistinguishable from the alternative for any test
##' with power; inspect the attached per-cluster contributions
##' (\code{attr(., "boot")$s_j}; the largest values identify the
##' clusters driving the statistic) and cluster diagnostics
##' (\code{cooks.distance(fit, groups = )},
##' \code{\link[=hatvalues.rlmerMod]{hatvalues}}) when a rejection is
##' suspect. \strong{Experimental}: simulation-validated, not proven;
##' the deviance bootstrap \code{test = "boot"} remains the default
##' variance-component path.
##'
##' @param object An \code{rlmerMod} object.
##' @param ... A second \code{rlmerMod} object for pairwise comparison.
##' @param test One of \code{"Wald"} (default; closed form for nested
##'   fixed-effects-only tests), \code{"boot"} (parametric bootstrap
##'   quasi-deviance; the default valid path for variance-component
##'   tests), or \code{"score"} (experimental one-sided robust score
##'   test for a single added variance component, calibrated by a
##'   score-only parametric bootstrap; see Details for its scope and
##'   validation record).
##' @param null Bootstrap null generation for \code{test = "boot"}.
##'   \code{"parametric"} (default) generates from \code{fit0}'s robust
##'   estimates (the exact current behaviour). \code{"robust"} is an
##'   experimental contamination-robust null that generates from a cleaned
##'   null fit (see Details); it only affects the variance-component
##'   bootstrap path and is ignored for the Wald paths.
##' @param vcov_type Forwarded to \code{\link[=vcov.rlmerMod]{vcov}}
##'   for the Wald paths; ignored for the bootstrap.
##' @param ddf Denominator degrees of freedom for the Wald paths.
##'   \code{"none"} (default) reports the chi-square table as before;
##'   \code{"satterthwaite"} reports an F-test with a Satterthwaite
##'   denominator df (the multivariate \code{lmerTest::contestMD}
##'   generalisation, built on the robust IF-based covariance of the
##'   variance parameters). It requires \code{vcov_type =
##'   "default"} and a single grouping factor; otherwise it warns and
##'   falls back to the chi-square table.
##' @param nsim Bootstrap replicates when \code{test = "boot"} or
##'   \code{test = "score"}; default 1000 for \code{"boot"} and 199
##'   for \code{"score"} (whose null-only refits are ~4x cheaper).
##' @param seed Optional RNG seed for reproducibility of the bootstrap.
##' @param verbose Bootstrap progress messages.
##' @return An \code{"anova"} data.frame; the bootstrap path attaches
##'   \code{attr(., "boot") = list(D_boot, n_fail, nsim, D_obs)}. The
##'   score path attaches \code{attr(., "boot") = list(S_boot, n_fail,
##'   nsim, S_obs, s_j)} with \code{s_j} the named per-cluster
##'   contributions to the observed statistic.
##' @references Heritier S, Ronchetti E (1994). \emph{Robust
##'   bounded-influence tests in general parametric models}. JASA
##'   89(427): 897--904.
##' @references Cantoni E, Ronchetti E (2001). \emph{Robust inference
##'   for generalized linear models}. JASA 96(455): 1022--1030.
##' @references Heritier S, Cantoni E, Copt S, Victoria-Feser MP
##'   (2009). \emph{Robust Methods in Biostatistics}. Wiley.
##' @seealso \code{\link[=vcov.rlmerMod]{vcov}},
##'   \code{\link[=confint.rlmerMod]{confint}}
##' @importFrom stats anova pchisq pf
##' @method anova rlmerMod
##' @export
anova.rlmerMod <- function(object, ..., test = c("Wald", "boot", "score"),
                           null = c("parametric", "robust"),
                           vcov_type = c("default", "sandwich"),
                           ddf = c("none", "satterthwaite"),
                           nsim = 1000L, seed = NULL,
                           verbose = FALSE) {
    stopifnot(is(object, "rlmerMod"))
    vcov_type <- match.arg(vcov_type)
    test_arg  <- match.arg(test)
    null      <- match.arg(null)
    ddf       <- match.arg(ddf)
    ## Score refits only the null model per replicate (~4x cheaper),
    ## so its validated default is nsim = 199, not the boot's 1000.
    if (test_arg == "score" && missing(nsim))
        nsim <- 199L
    if (ddf == "satterthwaite" && vcov_type != "default") {
        warning("ddf = \"satterthwaite\" is only available with ",
                "vcov_type = \"default\"; reverting to the chi-square ",
                "Wald table.", call. = FALSE)
        ddf <- "none"
    }
    ## Compute the (cached) influence function once for the fit that
    ## carries the tested coefficients; fall back to the chi-square
    ## table if it is unavailable (unsupported design, etc.).
    .satt_IF <- function(fit) {
        IF <- tryCatch(implicitIF_full(fit), error = function(e) NULL)
        if (is.null(IF))
            warning("Satterthwaite df unavailable for this fit; ",
                    "reverting to the chi-square Wald table.",
                    call. = FALSE)
        IF
    }

    dots <- list(...)
    if (length(dots) == 0L) {
        IF <- NULL
        if (ddf == "satterthwaite") {
            IF <- .satt_IF(object)
            if (is.null(IF)) ddf <- "none"
        }
        return(.anova_single_wald(object, vcov_type = vcov_type,
                                  ddf = ddf, IF = IF))
    }

    if (length(dots) > 1L)
        stop("anova.rlmerMod currently supports pairwise comparison ",
             "only (anova(fit0, fit1)).", call. = FALSE)
    other <- dots[[1L]]
    if (!is(other, "rlmerMod"))
        stop("Both arguments must be rlmerMod objects.", call. = FALSE)

    ## Canonicalise so fit0 = smaller, fit1 = larger.
    n_par <- function(f) f@pp$p + length(getME(f, "theta")) + 1L
    if (n_par(other) < n_par(object)) {
        fit0 <- other; fit1 <- object
    } else {
        fit0 <- object; fit1 <- other
    }

    is_vc <- .differ_in_re(fit0, fit1)
    if (is_vc && test_arg == "Wald") {
        warning("Models differ in random-effects structure; the Wald ",
                "test on the variance-component boundary is invalid. ",
                "Switching to test = \"boot\" (parametric bootstrap).",
                call. = FALSE)
        test_arg <- "boot"
    }

    if (test_arg == "score") {
        if (!is_vc)
            stop("test = \"score\" tests a single added variance ",
                 "component; the models do not differ in their ",
                 "random-effects structure. Use test = \"Wald\" for ",
                 "fixed-effects comparisons.", call. = FALSE)
        ## No .bootstrap_group_guard call here, deliberately: the score
        ## statistic's per-cluster contributions are psi-bounded, so
        ## its bootstrap null law is insensitive to the O(eps)-biased
        ## generating parameters (validated: cleaned vs raw generating
        ## parameters changed at most 1-2 decisions in 600+ paired
        ## replicates, including the uncentered-x cells). The guard's
        ## "poisoned null" warning is a deviance-bootstrap concern and
        ## would be misleading noise for this path.
        sc <- .anova_pair_score(fit0, fit1, nsim = nsim, seed = seed,
                                verbose = verbose)
        out <- data.frame(npar  = c(n_par(fit0), n_par(fit1)),
                          Df    = c(NA_integer_, 1L),
                          Score = c(NA_real_, sc$S_obs),
                          `Pr(>=Score)` = c(NA_real_, sc$p.value),
                          check.names = FALSE)
        rownames(out) <- c("fit0", "fit1")
        attr(out, "boot") <- list(S_boot = sc$S_boot,
                                  n_fail = sc$n_fail,
                                  nsim   = sc$nsim,
                                  S_obs  = sc$S_obs,
                                  s_j    = sc$s_j)
        heading <- c(
            "Robust score test for one added variance component (rlmer)",
            "  EXPERIMENTAL: simulation-validated, not proven; see ?anova.rlmerMod",
            sprintf("  tested component: %s, one-sided", sc$label),
            sprintf(paste0("  parametric score bootstrap (null refits ",
                           "only): nsim = %d (effective = %d, fit ",
                           "failures = %d)"),
                    sc$nsim, length(sc$S_boot), sc$n_fail),
            paste0("  largest per-cluster contributions attr(., ",
                   "\"boot\")$s_j identify the driving clusters"))
        return(structure(out, heading = heading,
                         class = c("anova", "data.frame")))
    }

    if (test_arg == "Wald") {
        IF <- NULL
        if (ddf == "satterthwaite") {
            IF <- .satt_IF(fit1)
            if (is.null(IF)) ddf <- "none"
        }
        return(.anova_pair_wald(fit0, fit1, vcov_type = vcov_type,
                                ddf = ddf, IF = IF))
    }

    ## C4 auto-guard: warn if the robust null fit downweights a group
    ## (the bootstrap null is built from fit0, so group contamination
    ## there makes the test anti-conservative).
    .bootstrap_group_guard(fit0)
    robust_null <- identical(null, "robust")
    boot  <- .anova_pair_boot(fit0, fit1, nsim = nsim, seed = seed,
                              verbose = verbose, robust.null = robust_null)
    p0 <- fit0@pp$p; p1 <- fit1@pp$p
    L0 <- length(getME(fit0, "theta")); L1 <- length(getME(fit1, "theta"))
    df_diff <- (p1 - p0) + (L1 - L0)
    sig0 <- .sigma(fit0)
    out <- data.frame(npar  = c(n_par(fit0), n_par(fit1)),
                      QuasiDev = c(.quasi_deviance(fit0, scale = sig0),
                                   .quasi_deviance(fit1, scale = sig0)),
                      Df    = c(NA_integer_, df_diff),
                      Diff  = c(NA_real_, boot$D_obs),
                      `Pr(>=Diff)` = c(NA_real_, boot$p.value),
                      check.names = FALSE)
    rownames(out) <- c("fit0", "fit1")
    attr(out, "boot") <- list(D_boot = boot$D_boot, n_fail = boot$n_fail,
                              nsim = nsim, D_obs = boot$D_obs,
                              robust.null = boot$robust.null)
    heading <- c(
        "Parametric bootstrap quasi-deviance test (rlmer)",
        sprintf("  nsim = %d (effective = %d, fit failures = %d)",
                nsim, length(boot$D_boot), boot$n_fail))
    if (robust_null) {
        rn <- boot$robust.null
        cleaned <- !is.null(rn) && isTRUE(rn$cleaned)
        reason <- if (is.null(rn)) NULL else rn$reason
        heading <- c(heading,
            if (cleaned)
                sprintf(paste0("  experimental robust null: bootstrap ",
                               "generated from a contamination-cleaned ",
                               "fit (%d cluster(s) trimmed)"),
                        rn$n_flagged)
            else if (identical(reason, "layout-unresolved"))
                paste0("  experimental robust null: random-effect weight ",
                       "layout unresolved; robust null not applied ",
                       "(plain parametric null used)")
            else if (identical(reason, "over-half-flagged"))
                paste0("  experimental robust null: > 50% of clusters ",
                       "flagged; robust null not applied (plain ",
                       "parametric null used)")
            else if (identical(reason, "refit-failed") ||
                     identical(reason, "no-observations-left"))
                paste0("  experimental robust null: trimmed refit not ",
                       "usable; robust null not applied (plain ",
                       "parametric null used)")
            else if (identical(reason, "rank-deficient"))
                paste0("  experimental robust null: trimming the flagged ",
                       "clusters would make the design rank-deficient; ",
                       "robust null not applied (plain parametric null used)")
            else
                paste0("  experimental robust null: no cluster flagged; ",
                       "identical to null = \"parametric\""))
    }
    structure(out, heading = heading,
              class = c("anova", "data.frame"))
}
