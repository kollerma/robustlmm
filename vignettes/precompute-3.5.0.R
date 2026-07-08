## Precompute the numbers quoted in vignettes/robustlmm-3.5.0.Rnw, so that
## the vignette pulls them from data via \Sexpr{} rather than hard-coding
## them. Reads the (build-ignored) simulation-study result files and writes
## the small, shipped RDS vignettes/robustlmm-3.5.0-numbers.rds.
##
## Re-run this after re-running the underlying studies:
##   Rscript inst/simulationStudy/structuredCovariances.R
##   Rscript inst/simulationStudy/mallowsEfficiency.R
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

numbers <- list(structured = structured, mallows = mallows)
saveRDS(numbers, "vignettes/robustlmm-3.5.0-numbers.rds")
cat("wrote vignettes/robustlmm-3.5.0-numbers.rds\n")
str(numbers)
