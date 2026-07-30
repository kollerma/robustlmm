########################################################
## robustnessBlockDiagonalExtendedRansacRefit.R        ##
##                                                     ##
## Re-fits ONLY the two RANSAC arms of the Phase-C     ##
## Extended study, after the init = "ransac" fix       ##
## (robustlmm 3.5.0-2).                                ##
##                                                     ##
## Why only two arms. `fitDatasets_rlmer_ransac` and   ##
## `_ransac_bisq` built the RANSAC start and passed it ##
## to rlmer() as a fitted model object, which made the ##
## fit use the RANSAC subsample rather than the data:  ##
## on this design 510 of 1000 observations, 27 of 50   ##
## clusters. `DAStau_bisq` and `DAStau_sizeOBR` do not ##
## use RANSAC and are unaffected, as are all seven     ##
## methods of the base robustnessBlockDiagonal.R study ##
## -- including RSEa, whose 23.5% is the baseline the  ##
## 16.6% is compared against. So the baseline does not ##
## move and only the two RANSAC rows are re-fitted.    ##
##                                                     ##
## The datasets are regenerated identically to         ##
## robustnessBlockDiagonalExtended.R (same set.seed,   ##
## same generators, same order), so the re-fitted arms ##
## stay paired with both the retained Extended arms    ##
## and the base study.                                 ##
##                                                     ##
## Writes to its own baseFilename; the pre-fix archive ##
## is left intact as the record of what was published. ##
##                                                     ##
##   Rscript robustnessBlockDiagonalExtendedRansacRefit.R [outdir] [ncores]
########################################################

require(robustlmm)
source(system.file("simulationStudy/randomNumberGenerators.R",
                   package = "robustlmm"))

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else getwd()
## 5 cores on the M1 Max for long jobs, per project convention.
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("robustlmm   :", as.character(packageVersion("robustlmm")), "\n")
cat("R           :", R.version.string, "\n")
cat("outdir      :", normalizePath(outdir), "\n")
cat("ncores      :", ncores, "\n\n")
flush.console()

########################################################
## generate datasets -- verbatim from
## robustnessBlockDiagonalExtended.R, so the same
## set.seed(101) produces the same 6000 datasets
########################################################
z1 <- rep(1, 10)
z2 <- 0:9
K <- list()
K[[1]] <- tcrossprod(z1, z1)
K[[2]] <- tcrossprod(z1, z2) + tcrossprod(z2, z1)
K[[3]] <- tcrossprod(z2, z2)
names(K) <- c("Subject.Intercept.",
              "Subject.Days.Intercept.",
              "Subject.Days")
p <- length(z2)
n <- length(unique(sleepstudy$Subject))
groups <- cbind(rep(1:p, each = n), rep(1:n, p))

preparedDataset <-
    prepareMixedEffectDataset(
        Reaction ~ Days + (Days | Subject),
        data = sleepstudy,
        lmeFormula = Reaction ~ Days,
        heavyLmeRandom = ~ Days,
        heavyLmeGroups = ~ Subject,
        lqmmRandom = ~ Days,
        lqmmGroup = "Subject",
        lqmmCovariance = "pdSymm",
        groups = groups,
        varcov = K,
        lower = c(0, -Inf, 0)
    )

generateDatasets <-
    function(errorGenerator, randomEffectGenerator) {
        generateMixedEffectDatasets(
            numberOfDatasetsToGenerate = 1000,
            preparedDataset,
            errorGenerator,
            randomEffectGenerator)
    }

set.seed(101)
datasetList <- list(
    "N/N"       = generateDatasets(srnorm, srnorm),
    "N/CN"      = generateDatasets(srnorm, srcnorm),
    "CN/N"      = generateDatasets(srcnorm, srnorm),
    "CN/CN"     = generateDatasets(srcnorm, srcnorm),
    "t3/t3"     = generateDatasets(srt3, srt3),
    "skt3/skt3" = generateDatasets(srskt3, srskt3)
)
datasets <- bindDatasets(datasetList = datasetList)
datasetIndexToGeneratorMap <-
    rep(factor(names(datasetList), names(datasetList)), each = 1000)
saveRDS(datasetIndexToGeneratorMap,
        file.path(outdir, "datasetIndexToGeneratorMap.rds"))

########################################################
## the two affected arms only
########################################################
fittingFunctions <- list(
    fitDatasets_rlmer_ransac,
    fitDatasets_rlmer_ransac_bisq
)

########################################################
## Pre-flight: assert the fix is actually in the build
## being used.  This is the whole point of the re-run,
## and it is cheap to check before spending hours: fit
## one dataset through each arm and require that every
## observation and every subject is in the fit.  Pre-fix
## these came back at roughly half.
########################################################
cat("=== pre-flight: RANSAC arms must fit all the data ===\n")
probe <- datasets
probe$numberOfDatasets <- 1L
for (fn in list(fitDatasets_rlmer_ransac, fitDatasets_rlmer_ransac_bisq)) {
    f <- suppressWarnings(suppressMessages(fn(probe, datasetIndices = 1L)))[[1]]
    lab <- attr(f, "label")
    nobs_f <- nobs(f)
    nsub_f <- nrow(ranef(f)[[1]])
    cat(sprintf("  %-32s nobs %d / %d, subjects %d / %d\n",
                lab, nobs_f, nrow(sleepstudy), nsub_f,
                length(unique(sleepstudy$Subject))))
    if (nobs_f != nrow(sleepstudy) ||
        nsub_f != length(unique(sleepstudy$Subject)))
        stop("pre-flight failed for ", lab, ": the installed robustlmm ",
             "still fits the RANSAC subsample. Install a build with the ",
             "init = \"ransac\" fix (3.5.0-2 or later) before re-running.")
}
cat("pre-flight OK\n\n")
flush.console()

baseFilename <- "datasets_robustnessBlockDiagonalExtendedRansacRefit"
t0 <- proc.time()[["elapsed"]]
results <- processDatasetsInParallel(
    datasets, outdir, baseFilename, fittingFunctions,
    chunkSize = 50, checkProcessed = TRUE,
    createMinimalSaveFile = TRUE,
    ncores = ncores, stdErrors = TRUE)
cat(sprintf("\nfitting took %.2f hours\n",
            (proc.time()[["elapsed"]] - t0) / 3600))

########################################################
## the check that matters: every fit must use all 180
## observations.  Pre-fix these arms used ~half.
########################################################
cat("\n=== rows per method ===\n")
print(table(results$label))

saveRDS(list(results = results,
             versions = list(
                 R = R.version.string,
                 robustlmm = as.character(packageVersion("robustlmm")),
                 lme4 = as.character(packageVersion("lme4")),
                 when = Sys.time()),
             seed = 101L,
             note = paste("Post-fix re-fit of the two RANSAC arms of the",
                          "Phase-C Extended study; see",
                          "IF-thread1/thread4/ransac-subsample-provenance.md")),
        file.path(outdir, "ransacRefit_results.rds"))
cat("\nsaved", file.path(outdir, "ransacRefit_results.rds"), "\n")
