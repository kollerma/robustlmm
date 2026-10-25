## Precompute the numbers quoted in the section "Studies for the
## small-sample and variance-component methods" of
## vignettes/simulationStudies.Rnw, so that the vignette pulls them from
## data via \Sexpr{} rather than hard-coding them. Reads the (build-ignored)
## result files under research/theory-methods and writes the small, shipped
## RDS vignettes/methods-studies-numbers.rds.
##
## Re-run from the package root after re-running any of the studies:
##   Rscript vignettes/precompute-methods-studies.R
##
## Build-ignored (see .Rbuildignore); the studies it reads are not shipped.

rd <- "research/theory-methods"

## --- cluster sandwich, intercept and within-cluster slope (T5) ----------
t5 <- readRDS(file.path(rd, "t5_results/summary.rds"))
sandwichWithin <- list(
    J = t5$J,
    covG1z = cbind(int = t5$a_int, slope = t5$a_slope),     # G1 + normal
    covG1t = cbind(int = t5$b_int, slope = t5$b_slope))     # G1 + t(J - 1)

## --- cluster sandwich with between-cluster terms (F3b, F3c, F3d) --------
## F3d re-ran both designs with fresh seeds, end-to-end through anova() and
## confint(); its table is the one reported.
f3d <- readRDS(file.path(rd, "f3d_results/summary.rds"))$table
f3d$termLabel <- c(w = "between covariate", f = "between factor (2 df)",
                   x = "within covariate", z = "between binary",
                   u = "between log-normal", v = "within covariate",
                   zu = "joint $z = u = 0$ (2 df)")[f3d$term]
f3b <- readRDS(file.path(rd, "f3b_results/summary.rds"))
f3b <- f3b[f3b$scen == "clean", ]
f3c <- readRDS(file.path(rd, "f3c_results/summary.rds"))
sandwichBetween <- list(
    table = f3d[, c("design", "J", "term", "termLabel", "n_ok", "old", "new",
                    "cov_int", "cov_slope")],
    nreps = max(f3d$n),
    oldRange = range(f3d$old),
    newRangeJ20 = range(f3d$new[f3d$J >= 20]),
    newMax = max(f3d$new),
    covRangeJ20 = range(c(f3d$cov_int, f3d$cov_slope)[f3d$J >= 20]),
    covRange = range(c(f3d$cov_int, f3d$cov_slope)),
    ## the first runs (F3b: design A, F3c: design B)
    firstNewRangeJ20 = range(c(f3b$CR2_HTZ[f3b$J >= 20], f3c$new[f3c$J >= 20])),
    firstOldMax = max(c(f3b$G1_chisq, f3c$old)),
    firstG1tIntMin = min(f3b$covG1_int),
    firstG1tIntMax = max(f3b$covG1_int))

## --- nuisance threshold of the bootstrap test (T2) -----------------------
t2 <- readRDS(file.path(rd, "t2_results/recount.rds"))
bootThreshold <- list(
    table = t2,
    atomRange = range(t2$atom_obs),
    nNegative = sum(t2$neg_far),
    nValues = 3 * sum(t2$n),
    nIdenticalTied = sum(t2$raw_tied == t2$thr_tied),
    nCells = nrow(t2),
    sizeJ50 = range(c(t2$raw_tied, t2$thr_tied)[grepl("J50", t2$cell)]))

## --- several added variance components, score test (T1) -----------------
t1 <- readRDS(file.path(rd, "t1_results/summary.rds"))
scoreMulti <- list(
    nreps = max(t1$n_ok),
    sizeClean = t1$reject_own[t1$scen == "null"],
    sizeContam = t1$reject_own[t1$scen == "contam"],
    power = t1$reject_own[t1$scen == "power"],
    J = t1$J[t1$scen == "null"])

## --- RANSAC candidate ranking (T4, F2 re-run, F1) ------------------------
t4 <- readRDS(file.path(rd, "t4_results/summary.rds"))
t4c <- t4[t4$contam, ]; t4n <- t4[!t4$contam, ]
old <- options(f2.source_only = TRUE)
Sys.setenv(F2_DIR = "f2_rerun_results")
source(file.path(rd, "f2_f1_summarise.R"), local = TRUE)
options(old)
b <- f2_basin(); iv <- f2_invar()
f1 <- lapply(list.files(file.path(rd, "f1_results"), "^cell_.*[.]rds$",
                        full.names = TRUE), readRDS)
names(f1) <- vapply(f1, function(o) o$cell, "")
rmse <- function(x) sqrt(mean((x - 2)^2, na.rm = TRUE))
ransac <- list(
    nreps = t4c$n,
    startIntMedQn = t4c$med_int_err_start_Qn,
    startIntMedCentred = t4c$med_int_err_start_centred,
    startSlopeRmseQn = t4c$rmse_slope_start_Qn,
    startSlopeRmseCentred = t4c$rmse_slope_start_centred,
    cleanRmseRatio = t4n$rmse_int_start_centred / t4n$rmse_int_start_Qn,
    sameAsDefault = range(iv$same_as_default),
    nInvar = max(iv$n),
    phonyRange = range(b$pooled$phony),
    phonyMaxDiff = max(abs(b$pooled$phony_minus_A)),
    nBasin = max(b$pooled$n),
    maxDiscordant = max(b$pooled$mcnemar_b + b$pooled$mcnemar_c),
    f1StartB50 = rmse(f1$clus3$res[, "b50_arm_start"]),
    f1StartB25 = rmse(f1$clus3$res[, "b25_arm_start"]),
    f1StartMarg = rmse(f1$clus3$res[, "marg_arm_start"]),
    f1Nreps = nrow(f1$clus3$res))

numbers <- list(sandwichWithin = sandwichWithin,
                sandwichBetween = sandwichBetween,
                bootThreshold = bootThreshold, scoreMulti = scoreMulti,
                ransac = ransac)
saveRDS(numbers, "vignettes/methods-studies-numbers.rds")
cat("wrote vignettes/methods-studies-numbers.rds\n")
str(numbers, max.level = 2)
