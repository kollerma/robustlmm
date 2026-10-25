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
    r_std <- fit@resp$wtres / scale
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

## Small-J cluster-sandwich reference (theory-methods F3): the CR2
## (bias-reduced) sandwich with the approximate-Hotelling F test of
## Pustejovsky and Tipton (2018) -- Bell-McCaffrey t for one df. When CR2 is
## not defined (random effects not nested within clusters) the G1 sandwich
## with F(k, J - 1) is used instead.
.cluster_vcov <- function(fit, cluster = NULL) {
    ## only the "random effects not nested within the clusters" condition
    ## falls back to G1; any other error surfaces. Warnings of the CR2
    ## attempt are held back until it succeeds (the G1 fallback re-issues
    ## the cluster warnings).
    held <- list()
    V <- tryCatch(withCallingHandlers(
                      .vcov_sandwich(fit, cluster = cluster, correction = "CR2"),
                      warning = function(w) {
                          held[[length(held) + 1L]] <<- w
                          invokeRestart("muffleWarning")
                      }),
                  rlmer_cr_not_nested = function(e) NULL)
    if (is.null(V))
        return(.vcov_sandwich(fit, cluster = cluster, correction = "G1"))
    for (w in held) warning(w)
    V
}
## Reference used by ddf = "auto" for vcov_type = "sandwich" (F3/F3b
## decision in research/theory-methods).
.ddf_auto_sandwich <- "cluster"

## F statistic and denominator df for the Wald statistic W of L beta = 0
## (L: k x p) with a .cluster_vcov() covariance. A contrast without an
## estimable sandwich variance (eta = NA; .vcov_sandwich has warned) or
## with eta <= k - 1 gets NA; fewer than one denominator df warns.
.cluster_Fref <- function(W, L, V, label = "") {
    k <- nrow(rbind(L))
    if (identical(attr(V, "correction"), "CR2")) {
        eta <- .cr2_df(V, L)
        df2 <- eta - k + 1
        if (is.finite(eta) && df2 <= 0)
            warning(sprintf(paste0("%s: the cluster sandwich has %.2g ",
                                   "effective degrees of freedom for %d ",
                                   "numerator df; no F reference (NA)."),
                            label, eta, k), call. = FALSE)
        if (!is.finite(df2) || df2 <= 0)
            return(list(F = NA_real_, df2 = NA_real_))
        if (df2 < 1)
            warning(sprintf(paste0("%s: %.2g denominator df -- the ",
                                   "clusters carry almost no information ",
                                   "on this term."), label, df2),
                    call. = FALSE)
        return(list(F = W * df2 / (eta * k), df2 = df2))
    }
    list(F = W / k, df2 = attr(V, "n.clusters") - 1)
}
.cluster_heading <- function(V) {
    if (identical(attr(V, "correction"), "CR2"))
        sprintf(paste0("  vcov_type = \"sandwich\" (CR2, %d clusters); ",
                       "denominator df: approximate Hotelling ",
                       "(Bell-McCaffrey for 1 df)"),
                attr(V, "n.clusters"))
    else
        sprintf(paste0("  vcov_type = \"sandwich\" (G1: CR2 needs random ",
                       "effects nested within clusters); denominator df: ",
                       "J - 1 = %d"), attr(V, "n.clusters") - 1L)
}

## Per-term robust Wald test from a single rlmer fit. ddf = "none" gives the
## chi-square table; ddf = "satterthwaite" an F-test with a Satterthwaite
## denominator df per term (default vcov); ddf = "cluster" the small-J
## cluster-sandwich F-test (.cluster_vcov / .cluster_Fref).
.anova_single_wald <- function(fit, vcov_type = "default",
                               ddf = "none", IF = NULL, cluster = NULL) {
    mm     <- model.matrix(fit)
    assign <- attr(mm, "assign")
    tlabs  <- attr(terms(fit), "term.labels")
    hasInt <- attr(terms(fit), "intercept") == 1L

    V  <- if (ddf == "cluster") .cluster_vcov(fit, cluster)
          else as.matrix(vcov(fit, type = vcov_type, cluster = cluster))
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
    Fval  <- Chisq / pmax(Df, 1L)
    for (k in seq_along(unique_aids)) {
        idx <- which(assign == unique_aids[k])
        Df[k] <- length(idx)
        Vk <- V[idx, idx, drop = FALSE]
        bk <- bh[idx]
        Chisq[k] <- as.numeric(crossprod(bk, solve(Vk, bk)))
        Fval[k]  <- Chisq[k] / Df[k]
        if (ddf %in% c("satterthwaite", "cluster")) {
            L <- matrix(0, length(idx), p)
            L[cbind(seq_along(idx), idx)] <- 1
        }
        if (ddf == "satterthwaite") {
            fd <- tryCatch(.satterthwaite_Fdenom(fit, L, IF = IF),
                           error = function(e) list(df2 = NA_real_))
            DenDF[k] <- fd$df2
        } else if (ddf == "cluster") {
            fr <- .cluster_Fref(Chisq[k], L, V, label = term_names[k])
            Fval[k]  <- fr$F
            DenDF[k] <- fr$df2
        }
    }
    if (ddf %in% c("satterthwaite", "cluster")) {
        Pr   <- stats::pf(Fval, Df, DenDF, lower.tail = FALSE)
        out  <- data.frame(NumDF = Df, DenDF = DenDF,
                           `F value` = Fval, `Pr(>F)` = Pr,
                           check.names = FALSE)
        rownames(out) <- term_names
        return(structure(out,
            heading = c(
                "Robust Wald F-test for fixed effects (rlmer)",
                if (ddf == "cluster") .cluster_heading(V) else
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
                             ddf = "none", IF = NULL, cluster = NULL) {
    extras <- .extra_fixed_cols(fit0, fit1)
    if (!extras$nested)
        stop("Models are not nested: fit0's fixed-effect columns ",
             "must be a subset of fit1's.", call. = FALSE)
    if (length(extras$extra) == 0L)
        stop("Models share the same fixed-effect design; nothing to ",
             "test.", call. = FALSE)
    V1   <- if (ddf == "cluster") .cluster_vcov(fit1, cluster)
            else as.matrix(vcov(fit1, type = vcov_type, cluster = cluster))
    bh1  <- .fixef(fit1)
    nms1 <- colnames(fit1@pp$X)
    idx  <- match(extras$extra, nms1)
    Vk   <- V1[idx, idx, drop = FALSE]
    bk   <- bh1[idx]
    Chi  <- as.numeric(crossprod(bk, solve(Vk, bk)))
    Df   <- length(idx)
    p0 <- fit0@pp$p; p1 <- fit1@pp$p
    if (ddf %in% c("satterthwaite", "cluster")) {
        L <- matrix(0, length(idx), p1)
        L[cbind(seq_along(idx), idx)] <- 1
        if (ddf == "cluster") {
            fr <- .cluster_Fref(Chi, L, V1,
                                label = paste(extras$extra, collapse = ", "))
            fd <- list(df2 = fr$df2)
            Fval <- fr$F
        } else {
            fd <- tryCatch(.satterthwaite_Fdenom(fit1, L, IF = IF),
                           error = function(e) list(df2 = NA_real_))
            Fval <- Chi / Df
        }
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
                if (ddf == "cluster") .cluster_heading(V1) else
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
## p-values) bit-identical. The draw itself is .sim_parametric(), the
## generator of simulate(type = "parametric"): the callers add the offset
## of the null model to Xb0 and pass its prior weights w, so the errors
## have variance sigma0^2 / w, as the refits assume.
.boot_sim_data <- function(data_obs, yname, Xb0, Z0, Ub0, sig0, n, q0,
                           w = NULL) {
    d_sim <- data_obs
    d_sim[[yname]] <- .sim_parametric(Xb0, Z0, Ub0, sig0, n, q0, w)
    d_sim
}

## ------------------------------------------------------------------
## Thresholded bootstrap generator (Guedon, Baey and Kuhn 2024).
##
## When a nuisance variance component of the null model is estimated
## near zero, generating the bootstrap from the raw estimate makes the
## quasi-deviance bootstrap inconsistent: the null law of the statistic
## depends on whether the nuisance component lies on the boundary, and
## a small positive estimate puts the generator off the boundary. The
## remedy is to set such components to exactly zero before generating,
## with a threshold c_J -> 0 such that sqrt(J) c_J -> infinity (J the
## number of groups). Machine-checked in the robustlmm theory notes
## (Topic4.FaceSelection: zero_selected, positive_retained;
## Topic4.ShrinkBootstrap: gbk_deviance_size). The weaker condition
## J c_J -> infinity is not sufficient
## (Topic4.FaceSelection.note_condition_insufficient).
##
## The threshold acts on the variance scale, where a zero component is
## sqrt(J)-tight (the standard-deviation scale theta is only
## J^(1/4)-tight at zero). Component c of a block is set to zero when
## its relative variance tau_c = (L L')_cc is below lambda_J * se_c,
## where se_c = sqrt(2 / sum_j ||z_jc||^4) is the Gaussian standard
## error of tau_c at zero (z_jc = the design column of component c in
## group j) and lambda_J = sqrt(2 log J). The retained components keep
## their covariance; their Cholesky factor is rebuilt by a QR of the
## retained rows of L (valid also when that covariance is singular).
## cs / ar1 blocks are left unchanged: their structure ties the
## components together.
## ------------------------------------------------------------------
.boot_threshold_theta <- function(fit0, theta,
                                  lambda = function(J) sqrt(2 * log(J))) {
    theta_new <- theta
    thnames <- names(getME(fit0, "theta"))
    Z <- t(fit0@pp$Zt)
    struct <- tryCatch(.structuredCovInfo(fit0),
                       error = function(e) vector("list", length(fit0@blocks)))
    zeroed <- character(0)
    info <- vector("list", length(fit0@blocks))
    for (b in seq_along(fit0@blocks)) {
        if (!is.null(struct[[b]])) next
        B <- as.matrix(fit0@blocks[[b]])
        d <- NROW(B)
        idx <- fit0@idx[[b]]
        J <- NCOL(idx)
        if (J < 2L) next
        nz <- B != 0
        L <- matrix(0, d, d)
        L[nz] <- theta[B[nz]]
        tau <- rowSums(L^2)
        z4 <- vapply(seq_len(d), function(cc) {
            zc <- Z[, idx[cc, ], drop = FALSE]
            sum(Matrix::colSums(zc^2)^2)
        }, numeric(1))
        se <- sqrt(2 / z4)
        lam <- lambda(J)
        zero <- is.finite(se) & tau < lam * se
        info[[b]] <- list(tau = tau, se = se, lambda = lam, J = J,
                          zero = zero, names = thnames[diag(B)])
        if (!any(zero)) next
        keep <- !zero
        Lnew <- matrix(0, d, d)
        if (any(keep)) {
            k <- sum(keep)
            ## tol = 0: no column pivoting, which would permute the
            ## retained components when their covariance is singular
            R <- qr.R(qr(t(L[keep, , drop = FALSE]), tol = 0))
            R <- R[seq_len(k), seq_len(k), drop = FALSE]
            sgn <- sign(diag(R))
            sgn[sgn == 0] <- 1
            Lnew[keep, keep] <- t(sgn * R)
        }
        ## the rebuilt factor must fit the block's theta pattern
        if (any(Lnew[!nz] != 0)) next
        theta_new[B[nz]] <- Lnew[nz]
        zeroed <- c(zeroed, thnames[diag(B)][zero])
    }
    list(theta = theta_new, zeroed = zeroed, info = info)
}

## Data frame for the bootstrap refits: fit0's model frame plus every
## column that only fit1 uses (e.g. the covariate of an added random
## slope). Before this, refits of fit1 failed whenever fit1 used a
## variable absent from fit0. A variable that enters a model only
## transformed (the frame of (0 + log(z) | g) stores "log(z)", not z) is
## taken from the fit's data, row-matched to the model frame.
.boot_data_frame <- function(fit0, fit1) {
    f0 <- fit0@frame
    f1 <- fit1@frame
    if (nrow(f0) != nrow(f1))
        stop("anova(test = \"boot\"): the two fits use different ",
             "observations (", nrow(f0), " vs ", nrow(f1), " rows); ",
             "fit both models to the same data.", call. = FALSE)
    for (v in setdiff(names(f1), names(f0)))
        f0[[v]] <- f1[[v]]
    ## variables of the formulas, and of the `weights` and `offset`
    ## arguments, which the refits evaluate in the bootstrap data as well
    argVars <- function(fit)
        c(all.vars(fit@call$weights), all.vars(fit@call$offset))
    need <- setdiff(unique(c(all.vars(formula(fit0)), all.vars(formula(fit1)),
                             argVars(fit0), argVars(fit1))),
                    names(f0))
    if (!length(need)) return(f0)
    env <- environment(formula(fit1))
    d <- tryCatch(eval(fit1@call$data, env), error = function(e) NULL)
    if (is.data.frame(d)) {
        idx <- match(rownames(f0), rownames(d))
        if (anyNA(idx)) {
            idx <- if (nrow(d) == nrow(f0)) seq_len(nrow(d)) else NULL
        }
        if (!is.null(idx))
            for (v in intersect(need, names(d))) f0[[v]] <- d[[v]][idx]
    }
    ## anything else must be found by the refit's formula environment
    miss <- need[!need %in% names(f0) &
                 !vapply(need, exists, logical(1), envir = env)]
    if (length(miss))
        stop("anova(test = \"boot\"): cannot find ",
             paste(miss, collapse = ", "), ", which the models use only ",
             "transformed; fit both models with `data` containing it.",
             call. = FALSE)
    f0
}

## Message before a long bootstrap (see .slowMessage): the workload is
## nsim * refits * the workload measure of the fit, `refits` being the
## number of robust refits per replicate.
.bootMessage <- function(fit, nsim, refits, what, avoid)
    .slowMessage(nsim * refits * .satterthwaite_df_workload(fit),
                 "robustlmm.message.boot.max", 1e7,
                 sprintf("%s with %d replicates", what, as.integer(nsim)),
                 avoid)

## Prior weights of the null model for the generator; NULL if all one
.boot_prior_weights <- function(fit) {
    w <- fit@resp$weights
    if (all(w == 1)) NULL else w
}

## Quiet bootstrap refit at a warm start; NULL on error.
.boot_refit <- function(fit, d_sim, start)
    tryCatch(suppressMessages(suppressWarnings(
        update(fit, data = d_sim, start = start))),
        error = function(e) NULL)

.anova_pair_boot <- function(fit0, fit1, nsim = 1000L, seed = NULL,
                             verbose = FALSE, robust.null = FALSE,
                             nuisance = "threshold") {
    .bootMessage(fit1, nsim, 2, "anova(test = \"boot\")",
                 paste0("Reduce 'nsim', or use test = \"score\" where it ",
                        "applies (one refit per replicate)."))
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
    data_obs <- .boot_data_frame(fit0, fit1)
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
    ## Nuisance variance components estimated near zero are set to zero
    ## in the generator (see .boot_threshold_theta).
    thr <- NULL
    if (identical(nuisance, "threshold")) {
        th_gen <- if (!is.null(rn_info) && !is.null(rn_info$theta))
                      rn_info$theta else th0
        thr <- .boot_threshold_theta(fit0, th_gen)
        if (length(thr$zeroed))
            Ub0@x <- thr$theta[pp0$Lind]
    }
    Xb0    <- as.numeric(X %*% beta0) + fit0@resp$offset
    w0     <- .boot_prior_weights(fit0)

    D_boot <- rep(NA_real_, nsim)
    n_fail <- 0L
    t0     <- Sys.time()
    for (r in seq_len(nsim)) {
        d_sim <- .boot_sim_data(data_obs, yname, Xb0, Z0, Ub0, sig0, n, q0,
                                w0)
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
    ## D is a likelihood-ratio-type difference with an atom at 0 (both fits
    ## put the tested variance on the boundary). There its computed value is
    ## rounding noise of either sign (|D| < 0.01 in practice, while genuine
    ## positive values are O(1)), and comparing noise with noise made the
    ## p-value arbitrary. Treat D <= D_TIE as exactly 0: an observed D at the
    ## atom gets p = 1. This also maps the occasional clearly negative D (the
    ## larger robust fit in a worse local solution, about 1 in 2000
    ## replicates in T2) to 0: it carries no evidence against the null
    ## either.
    D_TIE <- 0.01
    D_eff <- function(D) ifelse(D <= D_TIE, 0, D)
    pval <- mean(D_eff(D_clean) >= D_eff(D_obs))
    list(D_obs = D_obs, D_boot = D_clean,
         p.value = pval, n_fail = n_fail, nsim = nsim,
         robust.null = rn_info, threshold = thr)
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
                    "\"score\" tests added variance components only")
    th0 <- getME(fit0, "theta")
    th1 <- getME(fit1, "theta")
    q_add <- length(th1) - length(th0)
    if (q_add < 1L)
        .score_stop("the alternative must add at least one variance ",
                    "parameter relative to the null; got length(theta) ",
                    length(th0), " vs ", length(th1))
    cn0 <- fit0@cnms
    cn1 <- fit1@cnms
    Gp1 <- fit1@Gp
    nl  <- length(levels(fit1@flist[[1L]]))
    gname <- names(fit1@flist)[1L]

    ## The added components are fit1's component names that fit0 does not
    ## have (one grouping factor, so names identify components). Removing
    ## them from fit1's terms must give fit0's terms exactly, and each
    ## added component must carry exactly one variance parameter (its own
    ## scalar block): q_add components for q_add extra theta entries. A
    ## correlated slope adds two parameters (variance and covariance) for
    ## one component and is rejected here.
    added <- setdiff(unlist(cn1), unlist(cn0))
    reduced <- Filter(length, lapply(cn1, setdiff, added))
    if (!length(added) || !identical(unname(reduced), unname(cn0)))
        .score_stop("could not identify a single added scalar ",
                    "random-effect column, or a set of them, in the ",
                    "alternative model (models must be nested, the ",
                    "alternative adding independent components)")
    if (length(added) != q_add)
        .score_stop("each added component must carry exactly one variance ",
                    "parameter (an independent scalar block); got ",
                    length(added), " added component(s) for ", q_add,
                    " added parameters (length(theta) ", length(th0),
                    " vs ", length(th1), "). Correlated-slope ",
                    "alternatives are out of scope")

    ## b-columns of each added component; b-layout within a term is
    ## level-major: the length(cnms[[k]]) coordinates of one level are
    ## consecutive.
    add_cols <- lapply(added, function(nm) {
        for (k in seq_along(cn1)) {
            pos <- match(nm, cn1[[k]])
            if (is.na(pos)) next
            nc1 <- length(cn1[[k]])
            if (Gp1[k + 1L] - Gp1[k] != nc1 * nl) return(NULL)
            return(Gp1[k] + (seq_len(nl) - 1L) * nc1 + pos)
        }
        NULL
    })
    if (any(vapply(add_cols, is.null, logical(1))))
        .score_stop("could not identify the added random-effect columns ",
                    "in the alternative model")

    ## Each added component must be structurally uncorrelated with every
    ## other random effect: no structural nonzero of fit1's Lambda pattern
    ## may couple its b-coordinates to any other b-coordinate.
    Ub1 <- methods::as(methods::as(fit1@pp$U_b, "generalMatrix"),
                       "TsparseMatrix")
    for (cols in add_cols) {
        in_i <- (Ub1@i + 1L) %in% cols
        in_j <- (Ub1@j + 1L) %in% cols
        if (any(in_i != in_j))
            .score_stop("an added random-effect column is correlated ",
                        "with other random effects in the alternative ",
                        "model; only independent (uncorrelated) added ",
                        "components are in scope")
    }

    Z1 <- t(fit1@pp$Zt)
    labels <- paste0("(", added, " | ", gname, ")")
    if (q_add == 1L) {
        z_add <- Matrix::rowSums(Z1[, add_cols[[1L]], drop = FALSE])
        return(list(z_add = as.numeric(z_add), label = labels))
    }
    z_add <- vapply(add_cols, function(cols)
        as.numeric(Matrix::rowSums(Z1[, cols, drop = FALSE])),
        numeric(nrow(Z1)))
    colnames(z_add) <- added
    list(z_add = z_add, label = labels)
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

## ------------------------------------------------------------------
## Several added components (q >= 2): cone (chi-bar) statistic of the
## whitened robust score.
##
## Per cluster j and tested direction k, s_jk = (v_jk' psi(rt_j))^2 -
## kappa_1 ||v_jk||^2 exactly as in the scalar statistic. With
## U = sum_j s_j, W = sum_j (s_j - sbar)(s_j - sbar)' and the sensitivity
## M = d E[U] / d lambda (lambda = the added variances), the local model
## is U ~ N(M lambda, W) with lambda >= 0, and the statistic is
##     T = ||b||^2 - min_{lambda >= 0} ||b - A lambda||^2
##       = ||Pi_K b||^2,   b = W^{-1/2} U,  A = W^{-1/2} M,
## the squared norm of the projection of b onto the cone K = A R^q_+
## (the cone statistic of the robustlmm theory notes, Topic 4; its
## locally uniform limit and the consistency of the parametric
## bootstrap at an invertible whitening are machine-checked there,
## Topic4.GbkChiBar). Under r ~ N(0, I + lambda v_l v_l'),
##     M_kl = sum_j [ a (v_jk'v_jl)^2 + c sum_i v_jki^2 v_jli^2 ],
## a = (E psi')^2 and c = E[Z psi(Z) psi'(Z)] - (E psi')^2 (Z ~ N(0,1);
## E[psi psi''] = -E[psi'^2] + E[Z psi psi'] by parts). For psi(z) = z,
## c = 0 and M is the Gaussian score information.
## ------------------------------------------------------------------
.score_multi_consts <- function(rho) {
    psi <- rho@psi
    dpsi <- rho@Dpsi
    Ed <- rho@EDpsi()
    Ezpp <- stats::integrate(function(z) z * psi(z) * dpsi(z) *
                                 stats::dnorm(z), -Inf, Inf)$value
    c(a = Ed^2, c = Ezpp - Ed^2)
}

## ||Pi_K b||^2 for K = A R^q_+ by exact enumeration of the active sets.
.cone_proj_sq <- function(A, b) {
    q <- ncol(A)
    best <- sum(b^2)                       # residual at lambda = 0
    for (m in seq_len(2^q - 1L)) {
        S <- which(bitwAnd(m, 2^(seq_len(q) - 1L)) > 0)
        AS <- A[, S, drop = FALSE]
        lam <- tryCatch(solve(crossprod(AS), crossprod(AS, b)),
                        error = function(e) NULL)
        if (is.null(lam) || any(lam < 0)) next
        res <- sum((b - AS %*% lam)^2)
        if (res < best) best <- res
    }
    sum(b^2) - best
}

.anova_score_stat_multi <- function(fit0, Zadd, consts = NULL) {
    if (is.null(consts)) consts <- .score_multi_consts(fit0@rho.e)
    sig  <- .sigma(fit0)
    beta <- .fixef(fit0)
    r    <- as.numeric(fit0@resp$y - fit0@resp$offset -
                       fit0@pp$X %*% beta)
    psi  <- fit0@rho.e@psi
    kap1 <- fit0@rho.e@Epsi2()
    Mz   <- Matrix::crossprod(fit0@pp$Zt, fit0@pp$U_b)
    idx  <- split(seq_along(r), fit0@flist[[1L]])
    q    <- ncol(Zadd)
    Msens <- matrix(0, q, q)
    s <- t(vapply(idx, function(ii) {
        Mi <- as.matrix(Mz[ii, , drop = FALSE])
        Mi <- Mi[, colSums(abs(Mi)) > 0, drop = FALSE]
        Vi <- diag(length(ii)) + tcrossprod(Mi)   # V_j / sigma^2
        ee <- eigen(Vi, symmetric = TRUE)
        W  <- ee$vectors %*% (t(ee$vectors) / sqrt(ee$values)) / sig
        rt <- as.numeric(W %*% r[ii])
        V  <- W %*% Zadd[ii, , drop = FALSE]      # columns v_jk
        pr <- psi(rt)
        G  <- crossprod(V)                         # v_jk' v_jl
        H  <- crossprod(V^2)                       # sum_i v_jki^2 v_jli^2
        Msens <<- Msens + consts[["a"]] * G^2 + consts[["c"]] * H
        as.numeric(crossprod(V, pr))^2 - kap1 * colSums(V^2)
    }, numeric(q)))
    if (q == 1L) s <- t(s)
    U  <- colSums(s)
    Sc <- sweep(s, 2L, colMeans(s))
    Wc <- crossprod(Sc)
    ev <- eigen(Wc, symmetric = TRUE)
    if (any(ev$values <= 0))
        return(list(T = NA_real_, s_j = s, U = U, W = Wc, M = Msens,
                    min_eig_cor = 0))
    Wmh <- ev$vectors %*% (t(ev$vectors) / sqrt(ev$values))   # W^{-1/2}
    d <- sqrt(diag(Wc))
    min_eig_cor <- min(eigen(Wc / tcrossprod(d), symmetric = TRUE,
                             only.values = TRUE)$values)
    T <- .cone_proj_sq(Wmh %*% Msens, as.numeric(Wmh %*% U))
    list(T = T, s_j = s, U = U, W = Wc, M = Msens,
         min_eig_cor = min_eig_cor)
}

## Parametric SCORE bootstrap: simulate nsim datasets from fit0's
## estimates with the same generation machinery as .anova_pair_boot,
## refit ONLY the null model per replicate, and recompute S*. p-value
## = (1 + #{S* >= S_obs}) / (n_effective + 1), one-sided.
.anova_pair_score <- function(fit0, fit1, nsim = 199L, seed = NULL,
                              verbose = FALSE) {
    scope <- .anova_score_scope(fit0, fit1)
    if (is.matrix(scope$z_add))
        return(.anova_pair_score_multi(fit0, scope, nsim = nsim,
                                       seed = seed, verbose = verbose))
    .bootMessage(fit0, nsim, 1, "anova(test = \"score\")", "Reduce 'nsim'.")
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
    data_obs <- .boot_data_frame(fit0, fit0)
    Xb0   <- as.numeric(X %*% beta0) + fit0@resp$offset
    w0    <- .boot_prior_weights(fit0)

    S_boot <- rep(NA_real_, nsim)
    t0     <- Sys.time()
    for (r in seq_len(nsim)) {
        d_sim <- .boot_sim_data(data_obs, yname, Xb0, Z0, Ub0, sig0, n, q0,
                                w0)
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

## Parametric score bootstrap for q >= 2 added components: as in the
## scalar case, refit only the null model per replicate and recompute T*.
.anova_pair_score_multi <- function(fit0, scope, nsim = 199L, seed = NULL,
                                    verbose = FALSE) {
    .bootMessage(fit0, nsim, 1, "anova(test = \"score\")", "Reduce 'nsim'.")
    if (!is.null(seed)) set.seed(seed)
    Zadd   <- scope$z_add
    q      <- ncol(Zadd)
    if (q > 6L)
        stop("test = \"score\" supports at most 6 added components; got ",
             q, ".", call. = FALSE)
    consts <- .score_multi_consts(fit0@rho.e)
    obs <- .anova_score_stat_multi(fit0, Zadd, consts)
    if (!is.finite(obs$T))
        stop("test = \"score\": the covariance of the per-cluster score ",
             "contributions is singular; the tested directions cannot be ",
             "separated. Use test = \"boot\".", call. = FALSE)
    if (obs$min_eig_cor < 0.01)
        warning("test = \"score\": the score contributions of the tested ",
                "components are nearly collinear (smallest eigenvalue of ",
                "their correlation matrix ", signif(obs$min_eig_cor, 2),
                "); near a singular whitening the cone statistic's ",
                "calibration is unreliable. Consider test = \"boot\".",
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
    data_obs <- .boot_data_frame(fit0, fit0)
    Xb0   <- as.numeric(X %*% beta0) + fit0@resp$offset
    w0    <- .boot_prior_weights(fit0)

    T_boot <- rep(NA_real_, nsim)
    t0     <- Sys.time()
    for (r in seq_len(nsim)) {
        d_sim <- .boot_sim_data(data_obs, yname, Xb0, Z0, Ub0, sig0, n, q0,
                                w0)
        f0_s <- .boot_refit(fit0, d_sim, th0)
        if (is.null(f0_s)) next
        T_boot[r] <- tryCatch(.anova_score_stat_multi(f0_s, Zadd, consts)$T,
                              error = function(e) NA_real_)
        if (verbose && r %% 50L == 0L) {
            elapsed <- as.numeric(Sys.time() - t0, units = "secs")
            message(sprintf("  score boot rep %d/%d, elapsed %.1fs",
                            r, nsim, elapsed))
        }
    }
    T_clean <- T_boot[!is.na(T_boot)]
    n_fail  <- nsim - length(T_clean)
    if (length(T_clean) < 10L)
        stop("Score bootstrap: only ", length(T_clean), " of ", nsim,
             " replicates succeeded; need >= 10 to compute a p-value.",
             call. = FALSE)
    if (length(T_clean) < 100L)
        warning("Score-bootstrap p-value based on only ",
                length(T_clean), " effective replicates; consider ",
                "raising nsim for a stable tail estimate.",
                call. = FALSE)
    pval <- (1 + sum(T_clean >= obs$T)) / (length(T_clean) + 1)
    list(S_obs = obs$T, s_j = obs$s_j, S_boot = T_clean,
         p.value = pval, n_fail = n_fail, nsim = nsim,
         label = scope$label, q = q, min_eig_cor = obs$min_eig_cor)
}

##' Analysis of variance for an rlmer fit.
##'
##' Three modes.
##' \describe{
##'   \item{Single fit \code{anova(fit)}}{a per-term robust Wald table:
##'     each term is tested with \eqn{\hat\beta_t^\top \hat V_t^{-1}
##'     \hat\beta_t}, \eqn{\hat V} from \code{vcov(fit, type = vcov_type)},
##'     against \eqn{\chi^2} (default covariance) or a small-sample
##'     \eqn{F} reference (sandwich; see \code{ddf}).}
##'   \item{Nested fits differing only in the fixed effects}{a robust Wald
##'     test of the extra coefficients with the same covariance.}
##'   \item{Nested fits differing in the random effects}{Wald asymptotics
##'     fail on the boundary of the parameter space (Self and Liang 1987),
##'     so \code{test = "Wald"} switches to \code{test = "boot"}: a
##'     parametric bootstrap from \code{fit0}'s estimates of the
##'     quasi-deviance difference \eqn{D = 2 \sum_i
##'     [\rho_e(r_i^{(0)}/\hat\sigma_0) - \rho_e(r_i^{(1)}/\hat\sigma_0)]}
##'     (Heritier and Ronchetti 1994; Cantoni and Ronchetti 2001; Heritier
##'     et al. 2009).
##'     \code{test = "score"} is an experimental robust score test for
##'     added independent variance components.}
##' }
##' Only pairwise comparisons are implemented.
##'
##' \strong{Sandwich covariance.} With few clusters the G1 cluster sandwich
##' is biased down for between-cluster terms, and \eqn{\chi^2} tests based
##' on it are anti-conservative. For \code{vcov_type = "sandwich"} the
##' default reference is therefore the bias-reduced CR2 sandwich
##' (\code{\link{vcov_sandwich}}; Bell and McCaffrey 2002) with the
##' approximate Hotelling \eqn{F}-test of Pustejovsky and Tipton (2018),
##' \eqn{F = (\eta - k + 1) W / (\eta k)} on \eqn{(k, \eta - k + 1)}
##' degrees of freedom. It is much closer to its nominal level but can stay
##' somewhat anti-conservative for between-cluster terms; see the section
##' \dQuote{Small-sample references for the cluster sandwich} of
##' \code{vignette("simulationStudies")}. CR2 needs every random effect
##' nested within one cluster; otherwise the G1 sandwich with
##' \eqn{F(k, J - 1)} is used, and the table heading says which. A term
##' that a single cluster informs on its own (e.g. a factor level present
##' in one cluster) gets no reference distribution (\code{NA}), with a
##' warning.
##'
##' \strong{Bootstrap test (experimental).} \code{test = "boot"} holds its
##' level under error contamination but is anti-conservative when whole
##' groups are contaminated, because the bootstrap is generated from
##' estimates those groups have biased. It warns when the null fit heavily
##' downweights a group (smallest random-effect robustness weight below
##' 0.5). \code{null = "robust"} generates from a refit without the heavily
##' downweighted groups instead (the observed statistic still uses all the
##' data); it falls back to the plain generator when no group or more than
##' half of the groups are flagged, when removing them would leave the
##' design rank-deficient, or when the refit fails, and the table heading
##' says which applied. When both fits put the tested variance on the
##' boundary, \eqn{D} sits at its atom at zero and the p-value is 1. See the
##' ANOVA section of \code{vignette("robustlmm-3.5.0")}.
##'
##' \strong{Score test (experimental).} \code{test = "score"} uses the
##' robust null fit only. Per cluster \eqn{j} it forms the bounded
##' contribution \eqn{s_j = (v_j^\top \psi(\tilde r_j))^2 - \kappa_1
##' \|v_j\|^2}, with \eqn{\tilde r_j} the whitened marginal residuals,
##' \eqn{v_j} the whitened design column of a tested component and
##' \eqn{\kappa_1 = E[\psi(Z)^2]}{kappa_1 = E[psi(Z)^2]}, and calibrates
##' the statistic by a parametric bootstrap that refits only the null
##' model. One added
##' component is tested with the self-normalised sum of the \eqn{s_j};
##' several are tested jointly with the cone (chi-bar) statistic of the
##' whitened score vector, with a warning when the tested directions are
##' nearly collinear. Scope: both fits on the same data, with the same fixed
##' effects and a single grouping factor, and an alternative that adds
##' uncorrelated scalar components; otherwise it stops with an error.
##' Contamination aligned with a tested direction cannot be distinguished
##' from the alternative; the per-cluster contributions
##' (\code{attr(., "boot")$s_j}) show which clusters drive a rejection.
##' Simulation results: the ANOVA section of
##' \code{vignette("robustlmm-3.5.0")} (one component) and the section
##' \dQuote{Several added variance components} of
##' \code{vignette("simulationStudies")}.
##'
##' @param object An \code{rlmerMod} object.
##' @param ... A second \code{rlmerMod} object for pairwise comparison.
##' @param test \code{"Wald"} (default; nested fixed-effects models),
##'   \code{"boot"} (parametric bootstrap; the default for
##'   variance-component tests) or \code{"score"} (experimental; see
##'   Details).
##' @param null Bootstrap generator for \code{test = "boot"}:
##'   \code{"parametric"} (default) generates from \code{fit0}'s estimates,
##'   \code{"robust"} (experimental) from a refit without heavily
##'   downweighted groups (see Details).
##' @param nuisance How the \code{test = "boot"} generator treats variance
##'   components of the null model estimated close to zero:
##'   \code{"threshold"} (default) sets a component to zero when its
##'   relative variance is below \eqn{\sqrt{2\log J}} times its standard
##'   error at zero (\eqn{J} the number of groups), which keeps the
##'   bootstrap consistent when a nuisance component is zero (Guedon, Baey
##'   and Kuhn 2024); \code{"raw"} uses the estimates as they are (the
##'   behaviour before version 3.5.1). Compound-symmetric and AR(1) blocks are
##'   never thresholded. See the section \dQuote{Near-zero nuisance variance
##'   components in the bootstrap test} of \code{vignette("simulationStudies")}.
##' @param vcov_type Covariance for the Wald paths, passed to
##'   \code{\link[=vcov.rlmerMod]{vcov}}.
##' @param ddf Reference distribution for the Wald paths. \code{"auto"}
##'   (default): \code{"none"} for \code{vcov_type = "default"} and
##'   \code{"cluster"} for \code{vcov_type = "sandwich"}. \code{"none"}:
##'   \eqn{\chi^2}. \code{"satterthwaite"}: \eqn{F} with robust
##'   Satterthwaite denominator degrees of freedom (default covariance and a
##'   single grouping factor only; otherwise the default reference, with a
##'   warning). \code{"cluster"}: the CR2 sandwich with approximate Hotelling
##'   degrees of freedom (sandwich only; see Details).
##' @param cluster For \code{vcov_type = "sandwich"}, the clusters, passed
##'   to \code{\link{vcov_sandwich}}; needed for models with several
##'   grouping factors, e.g. \code{cluster = "school"} for
##'   \code{(1 | school/class)}.
##' @param nsim Bootstrap replicates; default 1000 for \code{test = "boot"}
##'   and 199 for \code{test = "score"}.
##' @param seed Optional RNG seed for the bootstrap.
##' @param verbose Bootstrap progress messages.
##' @return An \code{"anova"} data.frame; the bootstrap path attaches
##'   \code{attr(., "boot") = list(D_boot, n_fail, nsim, D_obs)}. The
##'   score path attaches \code{attr(., "boot") = list(S_boot, n_fail,
##'   nsim, S_obs, s_j)} with \code{s_j} the named per-cluster
##'   contributions to the observed statistic.
##' @references Self SG, Liang KY (1987). Asymptotic properties of maximum
##'   likelihood estimators and likelihood ratio tests under nonstandard
##'   conditions. \emph{Journal of the American Statistical Association}
##'   82(398), 605--610.
##' @references Guedon T, Baey C, Kuhn E (2024). \emph{Bootstrap test
##'   procedure for variance components in nonlinear mixed effects models
##'   in the presence of nuisance parameters and a singular Fisher
##'   information matrix}. arXiv:2306.10779.
##' @references Bell RM, McCaffrey DF (2002). Bias reduction in standard
##'   errors for linear regression with multi-stage samples. \emph{Survey
##'   Methodology} 28(2), 169--181.
##' @references Pustejovsky JE, Tipton E (2018). Small-sample methods for
##'   cluster-robust variance estimation and hypothesis testing in fixed
##'   effects models. \emph{Journal of Business & Economic Statistics}
##'   36(4), 672--683.
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
##' @section Lifecycle:
##'   The Wald tests are stable. Experimental, and subject to change:
##'   \code{test = "score"}, and for \code{test = "boot"} the options
##'   \code{null = "robust"} and \code{nuisance = "threshold"} (its
##'   default).
##' @export
anova.rlmerMod <- function(object, ..., test = c("Wald", "boot", "score"),
                           null = c("parametric", "robust"),
                           nuisance = c("threshold", "raw"),
                           vcov_type = c("default", "sandwich"),
                           ddf = c("auto", "none", "satterthwaite",
                                   "cluster"),
                           nsim = 1000L, seed = NULL,
                           verbose = FALSE, cluster = NULL) {
    stopifnot(is(object, "rlmerMod"))
    bootArgGiven <- !missing(null) || !missing(nuisance)
    vcov_type <- match.arg(vcov_type)
    test_arg  <- match.arg(test)
    null      <- match.arg(null)
    nuisance  <- match.arg(nuisance)
    ddf       <- match.arg(ddf)
    if (!is.null(cluster) && vcov_type != "sandwich")
        stop("'cluster' is used only with vcov_type = \"sandwich\"",
             call. = FALSE)
    ## (with the default test, a variance-component comparison switches
    ## to the bootstrap below, so these arguments may apply)
    if (bootArgGiven && (test_arg == "score" ||
                         (test_arg == "Wald" && !missing(test))))
        stop("'null' and 'nuisance' are used only with test = \"boot\"",
             call. = FALSE)
    .check_number(nsim, "nsim", lower = 1, integer = TRUE)
    ## Score refits only the null model per replicate (~4x cheaper),
    ## so its validated default is nsim = 199, not the boot's 1000.
    if (test_arg == "score" && missing(nsim))
        nsim <- 199L
    ## Each vcov has a default reference (ddf = "auto"): the chi-square
    ## table for the model-based vcov, the small-J cluster F-test for the
    ## sandwich (theory-methods F3). A reference unavailable for the chosen
    ## vcov falls back to that default.
    if (ddf == "satterthwaite" && vcov_type != "default") {
        warning("ddf = \"satterthwaite\" is only available with ",
                "vcov_type = \"default\"; using the sandwich's default ",
                "reference (ddf = \"auto\").", call. = FALSE)
        ddf <- "auto"
    }
    if (ddf == "cluster" && vcov_type != "sandwich") {
        warning("ddf = \"cluster\" applies to vcov_type = \"sandwich\"; ",
                "using the chi-square Wald table.", call. = FALSE)
        ddf <- "none"
    }
    if (ddf == "auto")
        ddf <- if (vcov_type == "sandwich") .ddf_auto_sandwich else "none"
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
                                  ddf = ddf, IF = IF, cluster = cluster))
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
            stop("test = \"score\" tests added variance ",
                 "components; the models do not differ in their ",
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
        q_sc <- if (is.null(sc$q)) 1L else sc$q
        out <- data.frame(npar  = c(n_par(fit0), n_par(fit1)),
                          Df    = c(NA_integer_, q_sc),
                          Score = c(NA_real_, sc$S_obs),
                          `Pr(>=Score)` = c(NA_real_, sc$p.value),
                          check.names = FALSE)
        rownames(out) <- c("fit0", "fit1")
        attr(out, "boot") <- list(S_boot = sc$S_boot,
                                  n_fail = sc$n_fail,
                                  nsim   = sc$nsim,
                                  S_obs  = sc$S_obs,
                                  s_j    = sc$s_j)
        heading <- if (q_sc == 1L) c(
            "Robust score test for one added variance component (rlmer)",
            "  EXPERIMENTAL; see ?anova.rlmerMod",
            sprintf("  tested component: %s, one-sided", sc$label)) else c(
            sprintf(paste0("Robust score test for %d added variance ",
                           "components (rlmer)"), q_sc),
            paste0("  EXPERIMENTAL: cone (chi-bar) statistic of the ",
                   "whitened robust score; see ?anova.rlmerMod"),
            sprintf("  tested components: %s, all one-sided",
                    paste(sc$label, collapse = ", ")))
        heading <- c(heading,
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
                                ddf = ddf, IF = IF, cluster = cluster))
    }

    ## C4 auto-guard: warn if the robust null fit downweights a group
    ## (the bootstrap null is built from fit0, so group contamination
    ## there makes the test anti-conservative).
    .bootstrap_group_guard(fit0)
    robust_null <- identical(null, "robust")
    boot  <- .anova_pair_boot(fit0, fit1, nsim = nsim, seed = seed,
                              verbose = verbose, robust.null = robust_null,
                              nuisance = nuisance)
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
                              robust.null = boot$robust.null,
                              threshold = boot$threshold)
    heading <- c(
        "Parametric bootstrap quasi-deviance test (rlmer)",
        sprintf("  nsim = %d (effective = %d, fit failures = %d)",
                nsim, length(boot$D_boot), boot$n_fail))
    if (abs(boot$D_obs) <= 0.01) {
        in_score_scope <- !inherits(tryCatch(.anova_score_scope(fit0, fit1),
                                             error = function(e) e), "error")
        heading <- c(heading, paste0(
            "  note: the alternative fit puts the tested variance at 0, ",
            "so the deviance difference is at its atom (p = 1); the ",
            "deviance bootstrap has little information near the null",
            if (in_score_scope) "; consider test = \"score\"" else ""))
    } else if (boot$D_obs < 0) {
        heading <- c(heading, paste0(
            "  note: the alternative fit has the larger quasi-deviance ",
            "(Diff < 0: the robust fits reached different local ",
            "solutions), so there is no evidence against the null (p = 1)"))
    }
    if (!is.null(boot$threshold) && length(boot$threshold$zeroed))
        heading <- c(heading,
            sprintf(paste0("  generator: nuisance variance component(s) ",
                           "%s estimated near zero set to 0 ",
                           "(threshold; nuisance = \"raw\" to disable)"),
                    paste(boot$threshold$zeroed, collapse = ", ")))
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
