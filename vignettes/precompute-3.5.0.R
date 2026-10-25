## Precompute the numbers quoted in vignettes/robustlmm-3.5.0.Rnw (and
## vignettes/rlmer.Rnw, vignettes/simulationStudies.Rnw), so that the
## vignettes pull them from data via \Sexpr{} rather than hard-coding them. Reads the (build-ignored) simulation-study result files and writes
## the small, shipped RDS vignettes/robustlmm-3.5.0-numbers.rds.
##
## Re-run this after re-running the underlying studies:
##   Rscript inst/simulationStudy/structuredCovariances.R
##   Rscript inst/simulationStudy/mallowsEfficiency.R
##   Rscript inst/simulationStudy/scoreTest.R
##   Rscript inst/simulationStudy/dasmcValidation.R
##   Rscript inst/simulationStudy/ransacConsensus.R
##   Rscript vignettes/precompute-3.5.0.R
##
## Build-ignored (see .Rbuildignore); the studies it reads are not shipped.

studyDir <- "inst/simulationStudy"

## --- structured covariances (Section 5.2) -------------------------------
sc <- readRDS(file.path(studyDir, "structuredCovariances_results.rds"))
sm <- sc$summary
rownames(sm) <- sm$scenario
getRow <- function(scn) sm[scn, ]

structured <- list(
    nreps     = sc$config$nreps,
    contamPct = round(100 * sc$config$eps),
    cs = list(
        trueRho    = getRow("cs_clean")$true_rho,
        cleanLmer  = getRow("cs_clean")$lmer_rho,
        cleanRlmer = getRow("cs_clean")$rlmer_rho,
        contLmer   = getRow("cs_cont")$lmer_rho,
        contRlmer  = getRow("cs_cont")$rlmer_rho,
        trueSd     = getRow("cs_cont")$true_sd,
        contSdLmer = getRow("cs_cont")$lmer_sd,
        contSdRlmer = getRow("cs_cont")$rlmer_sd),
    ar1 = list(
        trueRho    = getRow("ar1_clean")$true_rho,
        cleanLmer  = getRow("ar1_clean")$lmer_rho,
        cleanRlmer = getRow("ar1_clean")$rlmer_rho,
        contLmer   = getRow("ar1_cont")$lmer_rho,
        contRlmer  = getRow("ar1_cont")$rlmer_rho))

## --- Mallows design-weight efficiency (Section 5.1) ---------------------
## The shipped default is gamma = 1, cutoff 0.975. Report its clean-model
## efficiency cost relative to the plain RSE (1 - RE_vs_RSE).
me <- readRDS(file.path(studyDir, "mallowsEfficiency_results.rds"))
defaultRow <- me$tab[me$tab$method == "Mallows g1 c0.975", ]
mallows <- list(
    nreps       = me$config$nreps,
    costVsRsePct = round(100 * (1 - defaultRow$RE_vs_RSE), 1),
    pctDown      = defaultRow$pct_down)

## --- robust variance-component score test (Section "ANOVA") ------------
## Headline cells A-C (J = 50) and the adversarial sweep of
## inst/simulationStudy/scoreTest.R; rejection rates at alpha = 0.05.
st <- readRDS(file.path(studyDir, "scoreTest_results.rds"))
agg <- st$agg
rownames(agg) <- agg$cell
sw <- st$sweep
rownames(sw) <- sw$cell
sizeCells <- c("S1a", "S1b", "S2a", "S2b", "S4")   # null hypothesis true
scoreTest <- list(
    nreps       = st$config$NREP,
    J           = agg["A", "J"],
    contamPct   = round(100 * agg["A", "eps"]),
    contamScore = agg["A", "score_rate"],
    contamBoot  = agg["A", "boot_rate"],
    cleanScore  = agg["B", "score_rate"],
    cleanBoot   = agg["B", "boot_rate"],
    powerScore  = agg["C", "score_rate"],
    powerBoot   = agg["C", "boot_rate"],
    sweepSize   = range(sw[sizeCells, "score_rate"]),
    sweepSlopeDirection = sw["S3", "score_rate"])

## --- Monte-Carlo DAS-tau calibration (Section "Structured covariances") -
## Mean paired difference to lmer in the fitted ar1 correlation, largest
## absolute value over the two true correlations.
dm <- readRDS(file.path(studyDir, "dasmcValidation_results.rds"))
dasmc <- list(
    nreps          = dm$config$NREP,
    J              = dm$config$J,
    blockDim       = dm$config$NC,
    dasvarResidual = max(abs(dm$p_var_lmer$d_rho)),
    dasmcResidual  = max(abs(dm$p_mc_lmer$d_rho)))

## --- RANSAC multi-start consensus (Section "RANSAC") ---------------------
## Pooled over the six generators of inst/simulationStudy/ransacConsensus.R:
## phony-correlation rate and slope RMSE of the single best start
## (n_starts = 1) and of the consensus, and the paired discordances.
rc <- readRDS(file.path(studyDir, "ransacConsensus_results.rds"))
r <- rc$results
ok <- is.finite(r$phony1) & is.finite(r$phonyc)
slopeTrue <- unname(lme4::fixef(lme4::lmer(Reaction ~ Days + (Days | Subject),
                                            lme4::sleepstudy))[2])
ransacConsensus <- list(
    nreps           = rc$config$NREPS,
    nGenerators     = length(unique(r$gen)),
    nStarts         = rc$config$NSTART,
    nPairs          = sum(ok),
    phonySingle     = mean(r$phony1[ok]),
    phonyConsensus  = mean(r$phonyc[ok]),
    consensusBetter = sum(r$phony1[ok] == 1 & r$phonyc[ok] == 0),
    singleBetter    = sum(r$phony1[ok] == 0 & r$phonyc[ok] == 1),
    slopeRmseSingle    = sqrt(mean((r$slope1 - slopeTrue)^2, na.rm = TRUE)),
    slopeRmseConsensus = sqrt(mean((r$slopec - slopeTrue)^2, na.rm = TRUE)))

numbers <- list(structured = structured, mallows = mallows,
                scoreTest = scoreTest, dasmc = dasmc,
                ransacConsensus = ransacConsensus)
saveRDS(numbers, "vignettes/robustlmm-3.5.0-numbers.rds")
cat("wrote vignettes/robustlmm-3.5.0-numbers.rds\n")
str(numbers)
