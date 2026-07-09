## Validation study for the EXPERIMENTAL robust variance-component score
## test, anova(fit0, fit1, test = "score") -- see ?anova.rlmerMod. The
## script calls the PACKAGE interface end-to-end (statistic + score-only
## parametric bootstrap), so it doubles as a usage demo. Design (one
## added independent random slope): y ~ x + (1|g) with J = 50 clusters
## of size 6 and x = scale(1:6) within cluster, vs the alternative
## adding (0 + x | g). Headline cells, comparing test = "score" against
## the deviance bootstrap test = "boot" on identical data and seeds:
##   A  contaminated Type-I: 10% of clusters shifted by +5 (tau_slope = 0)
##   B  clean Type-I
##   C  clean power: tau_slope = 0.5
##
## ---- FINDING (nrep = 200/cell, score nsim = 199, boot nsim = 59; -----
## 2026-07-03/04; shipped in scoreTest_results.rds) -- rejection rates at
## alpha = 0.05:
##   A (contaminated null): score 0.035  vs boot 0.115
##   B (clean null):        score 0.045  vs boot 0.040
##   C (clean power):       score 0.920  vs boot 0.900
## The psi-bounded per-cluster contributions keep the score test's size
## under whole-cluster shift contamination where the deviance bootstrap
## (whose null generation is poisoned by the contaminated estimates) is
## anti-conservative; these are the numbers cited in the 3.5.0 feature
## vignette. The shipped rds carries the per-replicate p-values and
## observed statistics (the convention of the other results files here:
## per-replicate estimates and decisions, not raw bootstrap draws) plus
## the summary rates of the adversarial sweep (small J, uncentered x,
## slope-direction and observation-level contamination) cited in
## ?anova.rlmerMod. Original prototype study: robust-score-test.R +
## robust-score-test-phase2.R in the robustlmm-open-problems archive
## (statistic verified identical to the package implementation in
## robust-score-test-package-crosscheck.R).
## ----------------------------------------------------------------------

suppressMessages({ library(robustlmm); library(parallel) })

## NREP = 200 for the shipped results; the default here is modest so a
## test run finishes quickly (~40s per replicate on an M1 Max, dominated
## by the two bootstraps). Note: the score p-value uses the (1 + #)/(n +
## 1) convention, the boot p-value the mean(D* >= D) convention.
NREP       <- 20L
NCORES     <- 5L
NSIM_SCORE <- 199L   # package default for test = "score"
NSIM_BOOT  <- 59L    # as in the original study (boot refits both models)

CELLS <- rbind(
    data.frame(cell = "A", J = 50L, eps = 0.10, tau_slope = 0.0,
               kind = "TypeI-contam", base_seed = 30000L),
    data.frame(cell = "B", J = 50L, eps = 0.00, tau_slope = 0.0,
               kind = "TypeI-clean",  base_seed = 95000L),
    data.frame(cell = "C", J = 50L, eps = 0.00, tau_slope = 0.5,
               kind = "POWER",        base_seed = 60000L))

## data generator (identical to the original study)
gen_data <- function(J, ni = 6, b0 = 1, b1 = 0.5, tau_int = 1,
                     tau_slope = 0, eps = 0, shift_sd = 5) {
    g <- factor(rep(seq_len(J), each = ni))
    x <- rep(as.numeric(scale(seq_len(ni))), J)
    a <- rnorm(J, 0, tau_int)
    s <- if (tau_slope > 0) rnorm(J, 0, tau_slope) else rep(0, J)
    gi <- as.integer(g)
    y <- b0 + b1 * x + a[gi] + s[gi] * x + rnorm(J * ni)
    if (eps > 0) {
        bad <- sample.int(J, max(1L, round(eps * J)))
        y[gi %in% bad] <- y[gi %in% bad] + shift_sd
    }
    data.frame(y = y, x = x, g = g)
}

one <- function(rep_id, J, eps, tau_slope, base_seed) {
    set.seed(base_seed + rep_id)
    d <- gen_data(J = J, eps = eps, tau_slope = tau_slope)
    f0 <- tryCatch(suppressMessages(suppressWarnings(
        rlmer(y ~ x + (1 | g), data = d, method = "DASvar"))),
        error = function(e) NULL)
    f1 <- tryCatch(suppressMessages(suppressWarnings(
        rlmer(y ~ x + (1 | g) + (0 + x | g), data = d, method = "DASvar"))),
        error = function(e) NULL)
    if (is.null(f0) || is.null(f1))
        return(c(p_score = NA, p_boot = NA, S_obs = NA, D_obs = NA))
    bs <- base_seed * 1000L + rep_id   # common bootstrap seed
    aS <- tryCatch(suppressMessages(suppressWarnings(
        anova(f0, f1, test = "score", nsim = NSIM_SCORE, seed = bs))),
        error = function(e) NULL)
    aB <- tryCatch(suppressMessages(suppressWarnings(
        anova(f0, f1, test = "boot", nsim = NSIM_BOOT, seed = bs))),
        error = function(e) NULL)
    c(p_score = if (is.null(aS)) NA else aS[2L, "Pr(>=Score)"],
      p_boot  = if (is.null(aB)) NA else aB[2L, "Pr(>=Diff)"],
      S_obs   = if (is.null(aS)) NA else attr(aS, "boot")$S_obs,
      D_obs   = if (is.null(aB)) NA else attr(aB, "boot")$D_obs)
}

rate_ <- function(pv, alpha = 0.05) {
    pv <- pv[!is.na(pv)]
    c(rate = mean(pv <= alpha),
      se = sqrt(mean(pv <= alpha) * (1 - mean(pv <= alpha)) /
                max(length(pv), 1)),
      n = length(pv))
}

t0 <- Sys.time()
res <- agg <- list()
for (i in seq_len(nrow(CELLS))) {
    cn <- CELLS$cell[i]
    cat(sprintf("cell %s (%s): J=%d eps=%.2f tau_slope=%.1f nrep=%d\n",
                cn, CELLS$kind[i], CELLS$J[i], CELLS$eps[i],
                CELLS$tau_slope[i], NREP))
    M <- do.call(rbind, mclapply(seq_len(NREP), one, J = CELLS$J[i],
                                 eps = CELLS$eps[i],
                                 tau_slope = CELLS$tau_slope[i],
                                 base_seed = CELLS$base_seed[i],
                                 mc.cores = NCORES))
    res[[cn]] <- data.frame(cell = cn, rep = seq_len(NREP), M,
                            row.names = NULL)
    rs <- rate_(M[, "p_score"]); rb <- rate_(M[, "p_boot"])
    cat(sprintf("    score rate = %.3f (se %.3f), boot rate = %.3f (se %.3f), n = %d\n",
                rs["rate"], rs["se"], rb["rate"], rb["se"], rs["n"]))
    agg[[cn]] <- data.frame(CELLS[i, ], score_rate = rs[["rate"]],
                            score_se = rs[["se"]], boot_rate = rb[["rate"]],
                            boot_se = rb[["se"]], n = rs[["n"]],
                            row.names = NULL)
}
results <- do.call(rbind, res); rownames(results) <- NULL
agg <- do.call(rbind, agg); rownames(agg) <- NULL
runtime <- as.numeric(Sys.time() - t0, units = "mins")
cat(sprintf("total %.1f min\n", runtime))
print(agg, row.names = FALSE, digits = 3)

saveRDS(list(results = results, agg = agg,
             config = list(NREP = NREP, NSIM_SCORE = NSIM_SCORE,
                           NSIM_BOOT = NSIM_BOOT, runtime_min = runtime)),
        file.path("inst", "simulationStudy", "scoreTest_results.rds"))
cat("saved inst/simulationStudy/scoreTest_results.rds\n")
