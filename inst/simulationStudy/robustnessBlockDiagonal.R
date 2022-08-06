########################################################
## Replicate simulation study in Section 4.4 in       ##
## Koller and Stahel 2022 (figures 4 and 5)           ##
########################################################

require(robustlmm)
require(ggplot2)
source(system.file("simulationStudy/randomNumberGenerators.R",
                   package = "robustlmm"))

## path where to find the processed results
path <- system.file("simulationStudy", package = "robustlmm")
## number of cores to use for fitting datasets
ncores <- parallel::detectCores()

########################################################
## generate datasets and configure fitting functions  ##
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

## verify groups definition
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
        lower = c(0,-Inf, 0)
    )

generateDatasets <-
    function(errorGenerator, randomEffectGenerator) {
        datasets <-
            generateMixedEffectDatasets(
                numberOfDatasetsToGenerate = 1000,
                preparedDataset,
                errorGenerator,
                randomEffectGenerator
            )
        return(datasets)
    }

set.seed(101)
datasetList <- list(
    "N/N" = generateDatasets(srnorm, srnorm),
    "N/CN" = generateDatasets(srnorm, srcnorm),
    "CN/N" = generateDatasets(srcnorm, srnorm),
    "CN/CN" = generateDatasets(srcnorm, srcnorm),
    "t3/t3" = generateDatasets(srt3, srt3),
    "skt3/skt3" = generateDatasets(srskt3, srskt3)
)
datasets <- bindDatasets(datasetList = datasetList)
datasetIndexToGeneratorMap <-
    rep(factor(names(datasetList), names(datasetList)), each = 1000)

fittingFunctions <-
    list(
        fitDatasets_lmer,
        fitDatasets_rlmer_DAStau,
        fitDatasets_rlmer_DAStau_noAdj,
        fitDatasets_heavyLme,
        fitDatasets_lqmm,
        fitDatasets_varComprob_compositeTau_OGK,
        fitDatasets_varComprob_S_OGK
    )

########################################################
## fit all datasets or load and verify results        ##
########################################################

baseFilename <- "datasets_robustnessBlockDiagonal"
createMinimalSaveFile <-
    path != system.file("simulationStudy", package = "robustlmm")
results <- processDatasetsInParallel(
    datasets,
    path,
    baseFilename,
    fittingFunctions,
    chunkSize = 50,
    checkProcessed = TRUE,
    createMinimalSaveFile = createMinimalSaveFile,
    ncores = ncores,
    stdErrors = TRUE
)

########################################################
## prepare data for plotting                          ##
########################################################

convertTheta <- function(thetas, sigma) {
    if (!is.matrix(thetas)) {
        thetas <- matrix(thetas, 1, 3, byrow = TRUE)
    }
    v11 <- thetas[, 1] ^ 2
    v21 <- thetas[, 1] * thetas[, 2]
    v22 <- thetas[, 2] ^ 2 + thetas[, 3] ^ 2
    correlation <- v21 / sqrt(v11 * v22)
    return(cbind(
        B0.sigma = sigma * sqrt(v11),
        B1.sigma = sigma * sqrt(v22),
        B.corr = correlation
    ))
}

plotData <- cbind(
    data.frame(Method = as.factor(results$label),
               Generator = datasetIndexToGeneratorMap[results$datasetIndex]),
    results$coefficients,
    results$sigma,
    convertTheta(results$thetas, results$sigma[, 1]),
    results$datasetIndex
)
plotData <-
    subset(plotData,
           Method != "fitDatasets_heavyLme" &
               Method != "fitDatasets_lqmm")
levels(plotData$Method) <-
    shortenLabelsKS2022(levels(plotData$Method))
names(plotData)[3:4] <- c("beta0", "beta1")
plotDataLong <- reshape2::melt(plotData[-ncol(plotData)], 1:2)
plotDataTmp <-
    aggregate(plotDataLong[["value"]], plotDataLong[1:3], function(x)
        unlist(MASS::hubers(x, k = 1.345)))
plotDataTmp <- cbind(plotDataTmp[1:3],
                     plotDataTmp[[4]])
names(plotDataTmp)[4:5] <- c("location", "scale")
plotDataAggr <-
    reshape2::melt(plotDataTmp, 1:3, variable.name = "type")

plotDataTruth <-
    data.frame(
        variable = factor(
            levels(plotDataAggr$variable),
            levels(plotDataAggr$variable)
        ),
        type = factor("location", levels = levels(plotDataAggr$type)),
        value = c(
            datasets$trueBeta,
            datasets$trueSigma,
            convertTheta(datasets$trueTheta, datasets$trueSigma)
        )
    )

########################################################
## create correlation == 1 or -1 table                ##
########################################################

correlationData <-
    subset(plotDataLong, variable == "B.corr", select = -3)
correlationTmp <-
    aggregate(value ~ Method + Generator, correlationData,
              function(x)
                  mean(abs(x) - 1 > -0.0001))
correlationTable <-
    stats::reshape(
        correlationTmp,
        direction = "wide",
        idvar = "Method",
        timevar = "Generator"
    )
names(correlationTable) <-
    sub("value.", "", names(correlationTable))

########################################################
## load and verify aggregated data from full results  ##
########################################################

aggregatedFile <-
    file.path(path, paste0(baseFilename, "-aggregated.Rdata"))
runningOnMinimalProcessedResults <- max(plotData$datasetIndex) == 3
if (runningOnMinimalProcessedResults) {
    if (file.exists(aggregatedFile)) {
        load(aggregatedFile)
        stopifnot(all.equal(plotData, partialPlotData, check.attributes = FALSE))
    } else {
        warning("Running on minimal processed results, ",
                "but aggregated plot data is missing.")
    }
} else if (!file.exists(aggregatedFile) && createMinimalSaveFile) {
    partialPlotData <- subset(plotData, datasetIndex <= 3)
    save(partialPlotData,
         plotDataLong,
         plotDataAggr,
         plotDataTruth,
         correlationTable,
         file = aggregatedFile)
}

########################################################
## plot results                                       ##
########################################################

plot_robustnessBlockDiagonal <-
    ggplot(plotDataAggr, aes(Generator, value, color = Method)) +
    geom_hline(data = plotDataTruth, aes(yintercept = value)) +
    lemon::geom_pointline(aes(group = Method)) +
    xlab("") +
    ggh4x::facet_grid2(variable ~ type, scales = "free_y", independent = "y")

if (interactive()) {
    print(plot_robustnessBlockDiagonal)
}

plot_violinBlockDiagonal <-
    ggplot(plotDataLong, aes(Method, value)) +
    geom_violin() + xlab("") +
    facet_grid(variable ~ Generator, scales = "free_y") +
    theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
    ))

if (interactive()) {
    print(plot_violinBlockDiagonal)
}
