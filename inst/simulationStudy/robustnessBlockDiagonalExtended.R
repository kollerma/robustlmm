########################################################
## robustnessBlockDiagonalExtended.R                   ##
##                                                     ##
## Runs the new fitDatasets_rlmer_* wrappers (DAStau_  ##
## bisq, DAStau_sizeOBR, ransac, ransac_bisq) on the   ##
## same 6000 contaminated sleepstudy datasets used by  ##
## robustnessBlockDiagonal.R (KS-2022 Section 4.4).    ##
## The phony-rho table for the new methods is computed ##
## and compared against the stored KS-2022 baseline    ##
## (RSEn 16.5% at CN/N, RSEa 23.5% at CN/N).           ##
########################################################

require(robustlmm)
require(ggplot2)
source(system.file("simulationStudy/randomNumberGenerators.R",
                   package = "robustlmm"))

## path where to find the processed results
path <- system.file("simulationStudy", package = "robustlmm")
ncores <- parallel::detectCores()

########################################################
## generate datasets (identical to robustnessBlockDiagonal.R)
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
stopifnot(all.equal(sleepstudy, dplyr::arrange(sleepstudy, groups[, 1], groups[, 2])))

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

########################################################
## Fitting functions: only the NEW methods.
########################################################
fittingFunctions <- list(
    fitDatasets_rlmer_DAStau_bisq,
    fitDatasets_rlmer_DAStau_sizeOBR,
    fitDatasets_rlmer_ransac,
    fitDatasets_rlmer_ransac_bisq
)

########################################################
## Fit. Use a new baseFilename to avoid overwriting the
## original stored results.
########################################################
baseFilename <- "datasets_robustnessBlockDiagonalExtended"
createMinimalSaveFile <-
    path != system.file("simulationStudy", package = "robustlmm")
results <- processDatasetsInParallel(
    datasets, path, baseFilename, fittingFunctions,
    chunkSize = 50, checkProcessed = TRUE,
    createMinimalSaveFile = createMinimalSaveFile,
    ncores = ncores, stdErrors = TRUE)

########################################################
## Convert (theta_1, theta_2, theta_3) -> (sigma_0, sigma_1, rho).
########################################################
convertTheta <- function(thetas, sigma) {
    if (!is.matrix(thetas))
        thetas <- matrix(thetas, 1, 3, byrow = TRUE)
    v11 <- thetas[, 1] ^ 2
    v21 <- thetas[, 1] * thetas[, 2]
    v22 <- thetas[, 2] ^ 2 + thetas[, 3] ^ 2
    correlation <- v21 / sqrt(v11 * v22)
    cbind(B0.sigma = sigma * sqrt(v11),
          B1.sigma = sigma * sqrt(v22),
          B.corr   = correlation)
}

plotData <- cbind(
    data.frame(Method = as.factor(results$label),
               Generator = datasetIndexToGeneratorMap[results$datasetIndex]),
    results$coefficients,
    results$sigma,
    convertTheta(results$thetas, results$sigma[, 1]),
    results$datasetIndex
)
names(plotData)[3:4] <- c("beta0", "beta1")

## Use the same shortener as the original script. Falls back to
## stripping the "fitDatasets_rlmer_" prefix for methods missing
## from the shortener.
new_short <- function(labels) {
    out <- shortenLabelsKS2022(labels)
    miss <- is.na(out) | out == labels
    out[miss] <- sub("^fitDatasets_rlmer_", "", labels[miss])
    out
}
levels(plotData$Method) <- new_short(levels(plotData$Method))

########################################################
## Phony rho table for the new methods, side by side
## with the stored original rates.
########################################################
correlationData_new <- plotData[, c("Method", "Generator", "B.corr")]
correlationTable_new <-
    aggregate(B.corr ~ Method + Generator, correlationData_new,
              function(x) mean(abs(x) - 1 > -1e-4))
correlationTable_new <-
    stats::reshape(correlationTable_new, direction = "wide",
                   idvar = "Method", timevar = "Generator")
names(correlationTable_new) <- sub("B.corr.", "", names(correlationTable_new))

cat("\n=== New-methods phony-rho rates ===\n")
print(correlationTable_new)

cat("\n=== Stored KS-2022 baseline (from robustnessBlockDiagonal) ===\n")
aggregatedOriginal <- file.path(path,
    "datasets_robustnessBlockDiagonal-aggregated.Rdata")
if (file.exists(aggregatedOriginal)) {
    e <- new.env()
    load(aggregatedOriginal, envir = e)
    if (!is.null(e$correlationTable)) {
        print(e$correlationTable)
    }
}

########################################################
## Plots.
########################################################
plotData_long <- reshape2::melt(plotData[, -ncol(plotData)], 1:2)
plotDataTmp <-
    aggregate(plotData_long[["value"]], plotData_long[1:3],
              function(x) unlist(MASS::hubers(x, k = 1.345)))
plotDataTmp <- cbind(plotDataTmp[1:3], plotDataTmp[[4]])
names(plotDataTmp)[4:5] <- c("location", "scale")
plotDataAggr_new <-
    reshape2::melt(plotDataTmp, 1:3, variable.name = "type")

plotDataTruth <-
    data.frame(
        variable = factor(levels(plotDataAggr_new$variable),
                          levels(plotDataAggr_new$variable)),
        type = factor("location", levels = levels(plotDataAggr_new$type)),
        value = c(datasets$trueBeta, datasets$trueSigma,
                  convertTheta(datasets$trueTheta, datasets$trueSigma))
    )

plot_robustnessBlockDiagonalExtended <-
    ggplot(plotDataAggr_new,
           aes(Generator, value, color = Method)) +
    geom_hline(data = plotDataTruth, aes(yintercept = value)) +
    lemon::geom_pointline(aes(group = Method)) +
    xlab("") +
    ggh4x::facet_grid2(variable ~ type, scales = "free_y",
                        independent = "y")

if (interactive())
    print(plot_robustnessBlockDiagonalExtended)
