########################################################
## Replicate simulation study in Section 4.2 in       ##
## Koller and Stahel 2022 (figures 2 and 3)           ##
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

generateDatasets <-
    function(errorGenerator, randomEffectGenerator) {
        datasets <-
            generateAnovaDatasets(
                numberOfDatasetsToGenerate = 1000,
                numberOfLevelsInFixedFactor = 3,
                numberOfSubjects = 5,
                numberOfReplicates = 10,
                errorGenerator,
                randomEffectGenerator,
                lmeFormula = y ~ Var1,
                lower = 0,
                arrange = TRUE
            )
        ## add additional arguments used by varComprob
        data <- datasets$generateData(1)
        numberOfObservationsPerSample <-
            datasets$numberOfLevelsInFixedFactor * datasets$numberOfReplicates
        groups <-
            cbind(
                rep(1:numberOfObservationsPerSample, each = datasets$numberOfSubjects),
                as.numeric(data$Var2)
            )
        ## verify groups definition
        stopifnot(
            length(unique(groups[, 1])) == numberOfObservationsPerSample,
            length(unique(groups[, 2])) == datasets$numberOfSubjects,
            all.equal(data, dplyr::arrange(data, groups[, 1], groups[, 2]))
        )
        varcov <-
            matrix(1,
                   numberOfObservationsPerSample,
                   numberOfObservationsPerSample)
        datasets$groups <- groups
        datasets$varcov <- varcov
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

fitDatasets_varComprob_compositeTau_OGK_with_init_scale <-
    function(datasets, postFit) {
        if (require(robustvarComp)) {
            lcontrol <-
                varComprob.control(
                    lower = datasets[["lower"]],
                    cov.init = "covOGK",
                    init = list(scale = 15)
                )
        } else {
            lcontrol <- NULL
        }
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_compositeTau_OGK_with_init_scale",
                postFit = postFit
            )
        )
    }

fitDatasets_varComprob_S_OGK_with_init_scale <-
    function(datasets, postFit) {
        if (require(robustvarComp)) {
            lcontrol <-
                varComprob.control(
                    lower = datasets[["lower"]],
                    method = "S",
                    psi = "optimal",
                    cov.init = "covOGK",
                    init = list(scale = 15)
                )
        } else {
            lcontrol <- NULL
        }
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_S_OGK_with_init_scale",
                postFit = postFit
            )
        )
    }

fittingFunctions <-
    list(
        fitDatasets_lmer,
        fitDatasets_rlmer_DAStau,
        fitDatasets_rlmer_DAStau_noAdj,
        fitDatasets_varComprob_compositeTau_OGK_with_init_scale,
        fitDatasets_varComprob_S_OGK_with_init_scale
    )

########################################################
## fit all datasets or load and verify results        ##
########################################################

baseFilename <- "datasets_robustnessDiagonal"
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

computeSignificance <- function(results, datasets) {
    trueBeta <-
        matrix(
            datasets$trueBeta,
            nrow = NROW(results$coefficients),
            ncol = datasets$numberOfLevelsInFixedFactor,
            byrow = TRUE
        )
    tValues <-
        (results$coefficients - trueBeta) / results$standardErrors
    df.resid <-
        datasets$numberOfRows - 1 -
        (datasets$numberOfLevelsInFixedFactor - 1) -
        (datasets$numberOfSubjects - 1)
    criticalValue <- qt(0.975, df.resid)
    significant <- abs(tValues) > criticalValue
    colnames(significant) <-
        paste(colnames(significant), "significance", sep = ".")
    return(significant)
}

plotData <- cbind(
    data.frame(Method = as.factor(results$label),
               Generator = datasetIndexToGeneratorMap[results$datasetIndex]),
    results$coefficients,
    results$sigma,
    results$sigma * results$thetas,
    computeSignificance(results, datasets),
    results$datasetIndex
)
levels(plotData$Method) <-
    shortenLabelsKS2022(levels(plotData$Method))
names(plotData)[-(1:2)] <-
    c(
        "beta0",
        "beta1",
        "beta2",
        "sigma",
        "B.sigma",
        "beta0.significance",
        "beta1.significance",
        "beta2.significance",
        "datasetIndex"
    )

plotDataRobustnessTmp <- reshape2::melt(plotData[1:7], 1:2)
plotDataRobustnessAggr <-
    aggregate(plotDataRobustnessTmp[["value"]], plotDataRobustnessTmp[1:3], function(x)
        unlist(MASS::hubers(x, k = 1.345)))
plotDataRobustnessAggr <- cbind(plotDataRobustnessAggr[1:3],
                                plotDataRobustnessAggr[[4]])
names(plotDataRobustnessAggr)[4:5] <- c("location", "scale")
plotDataRobustness <-
    reshape2::melt(plotDataRobustnessAggr, 1:3, variable.name = "type")

plotDataTruth <-
    data.frame(
        variable = factor(
            levels(plotDataRobustness$variable),
            levels(plotDataRobustness$variable)
        ),
        type = factor("location", levels = levels(plotDataRobustness$type)),
        value = c(
            datasets$trueBeta,
            datasets$trueSigma,
            datasets$trueSigma * datasets$trueTheta
        )
    )

plotDataCoverageTmp <-
    reshape2::melt(plotData[-c(3:7, ncol(plotData))], 1:2)
plotDataCoverage <-
    aggregate(plotDataCoverageTmp[["value"]], plotDataCoverageTmp[1:3], function(x)
        1 - mean(x))
names(plotDataCoverage)[4] <- "coverage probability"
levels(plotDataCoverage$variable) <-
    sub(".significance", "", levels(plotDataCoverage$variable))

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
         plotDataRobustness,
         plotDataTruth,
         plotDataCoverage,
         file = aggregatedFile)
}

########################################################
## plot results                                       ##
########################################################

plot_robustnessDiagonal <-
    ggplot(plotDataRobustness, aes(Generator, value, color = Method)) +
    geom_hline(data = plotDataTruth, aes(yintercept = value)) +
    lemon::geom_pointline(aes(group = Method)) +
    xlab("") +
    ggh4x::facet_grid2(variable ~ type, scales = "free_y", independent = "y")

if (interactive()) {
    print(plot_robustnessDiagonal)
}

plot_coverageDiagonal <-
    ggplot(plotDataCoverage,
           aes(Generator, `coverage probability`, color = Method)) +
    geom_hline(yintercept = 0.95,
               color = "gray",
               size = 0.5) +
    lemon::geom_pointline(aes(group = Method)) +
    xlab("") +
    facet_wrap( ~ variable)

if (interactive()) {
    print(plot_coverageDiagonal)
}
