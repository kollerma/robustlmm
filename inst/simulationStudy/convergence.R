########################################################
## Simulation study looking at effect of increasing   ##
## number of observations on the efficiency of the    ##
## resulting estimates. As requested by one reviewer  ##
## of Koller and Stahel 2022.                         ##
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

fittingFunctions <-
    list(fitDatasets_lmer,
         fitDatasets_rlmer_DAStau,
         fitDatasets_rlmer_DAStau_noAdj)

prepareDataset <- function(numberOfSubjects, numberOfReplicates) {
    data <-
        expand.grid(Subject = 1:numberOfSubjects,
                    rep = 1:numberOfReplicates)
    set.seed(31 * numberOfSubjects + 37 * numberOfReplicates)
    data$ranef <- rnorm(numberOfSubjects)[data$Subject]
    data$error <- rnorm(nrow(data))
    data$continuousVariable <- runif(nrow(data))
    data$binaryVariable <- as.numeric(runif(nrow(data)) < 0.5)
    data$y <-
        with(data, 1 + continuousVariable + binaryVariable + ranef + error)
    data <-
        data[, c("y", "continuousVariable", "binaryVariable", "Subject")]

    preparedDataset <-
        prepareMixedEffectDataset(
            y ~ continuousVariable + binaryVariable +
                (1 | Subject),
            data,
            overrideBeta = c(1, 1, 1),
            overrideSigma = 1,
            overrideTheta = 1
        )
    return(preparedDataset)
}

########################################################
## functions to manage parameterized simulation       ##
########################################################

createFilename <- function(numberOfSubjects,
                           numberOfReplicates,
                           randomNumberGeneratorLabel) {
    return(
        paste0(
            "datasets_convergence_",
            numberOfSubjects,
            "_",
            numberOfReplicates,
            "_",
            fs::path_sanitize(randomNumberGeneratorLabel)
        )
    )
}

generateDatasets <- function(numberOfSubjects,
                             numberOfReplicates,
                             randomNumberGenerator) {
    preparedDataset <-
        prepareDataset(numberOfSubjects, numberOfReplicates)
    datasets <-
        generateMixedEffectDatasets(
            1000,
            preparedDataset,
            errorGenerator = randomNumberGenerator,
            randomEffectGenerator = randomNumberGenerator
        )
    return(datasets)
}

createMinimalSaveFile <-
    path != system.file("simulationStudy", package = "robustlmm")

generateAndProcessDatasets <-
    Vectorize(
        function(numberOfSubjects,
                 numberOfReplicates,
                 randomNumberGenerator,
                 randomNumberGeneratorLabel) {
            datasets <- generateDatasets(numberOfSubjects,
                                         numberOfReplicates,
                                         randomNumberGenerator)
            baseFilename <-
                createFilename(numberOfSubjects,
                               numberOfReplicates,
                               randomNumberGeneratorLabel)
            results <-
                processDatasetsInParallel(
                    datasets,
                    path,
                    baseFilename,
                    fittingFunctions,
                    chunkSize = 50,
                    checkProcessed = numberOfSubjects == 5 &&
                        numberOfReplicates == 5,
                    createMinimalSaveFile = createMinimalSaveFile,
                    ncores = ncores
                )
            results$datasetIndex
            nrows <- NROW(results$datasetIndex)
            results$numberOfSubjects <-
                matrix(numberOfSubjects, nrows)
            results$numberOfReplicates <-
                matrix(numberOfReplicates, nrows)
            results$randomNumberGeneratorLabel <-
                factor(rep.int(randomNumberGeneratorLabel, nrows))
            return(results)
        },
        vectorize.args = c("numberOfSubjects", "numberOfReplicates"),
        SIMPLIFY = FALSE
    )

runAllCombinations <- function(randomNumberGenerator,
                               randomNumberGeneratorLabel) {
    return(
        outer(
            c(5, 10, 20, 50),
            c(5, 10, 20, 50),
            generateAndProcessDatasets,
            randomNumberGenerator,
            randomNumberGeneratorLabel
        )
    )
}

########################################################
## fit all datasets or load and verify results        ##
########################################################

results_N_N <- runAllCombinations(srnorm, "N/N")
results_N_N <- mergeProcessedFits(results_N_N)

results_t3_t3 <- runAllCombinations(srt3, "t3/t3")
results_t3_t3 <- mergeProcessedFits(results_t3_t3)

########################################################
## prepare data for plotting                          ##
########################################################

convertToPlotData <- function(results) {
    plotData <- cbind(
        data.frame(Method = as.factor(results$label)),
        numberOfSubjects = factor(results$numberOfSubjects),
        numberOfReplicates = factor(results$numberOfReplicates),
        results$coefficients,
        results$sigma,
        results$sigma * results$thetas,
        results$datasetIndex
    )
    levels(plotData$Method) <-
        shortenLabelsKS2022(levels(plotData$Method))
    plotData$Method <-
        factor(plotData$Method, c("lme", "RSEn", "RSEa"))
    names(plotData)[-(1:3)] <- c("Intercept",
                                 "beta.continuous",
                                 "beta.binary",
                                 "sigma",
                                 "B.sigma",
                                 "datasetIndex")
    return(plotData)
}

aggregatePlotData <- function(plotData) {
    plotDataTmp <- reshape2::melt(plotData[-ncol(plotData)], 1:3)
    plotDataAggr <-
        aggregate(plotDataTmp[["value"]], plotDataTmp[1:4], function(x)
            unlist(MASS::hubers(x, k = 1.345)))
    plotDataAggr <- cbind(plotDataAggr[1:4],
                          plotDataAggr[[5]])

    plotDataAggr <- merge(
        plotDataAggr,
        subset(
            plotDataAggr,
            Method == "lme",
            select = c("numberOfSubjects", "numberOfReplicates", "variable", "s")
        ),
        by = names(plotDataAggr)[2:4],
        suffixes = c("", ".lme")
    )
    plotDataAggr$empiricalEfficiency <-
        plotDataAggr$s.lme / plotDataAggr$s
    plotDataAggr$s.lme <- NULL
    return(plotDataAggr)
}

plotData_N_N <- convertToPlotData(results_N_N)
plotDataAggr_N_N <- aggregatePlotData(plotData_N_N)
plotData_t3_t3 <- convertToPlotData(results_t3_t3)
plotDataAggr_t3_t3 <- aggregatePlotData(plotData_t3_t3)

########################################################
## load and verify aggregated data from full results  ##
########################################################

aggregatedFile <-
    file.path(path, "datasets_convergence-aggregated.Rdata")
runningOnMinimalProcessedResults <- max(plotData_N_N$datasetIndex) == 3
if (runningOnMinimalProcessedResults) {
    if (file.exists(aggregatedFile)) {
        load(aggregatedFile)
        stopifnot(all.equal(plotData_N_N, partialPlotData_N_N, check.attributes = FALSE))
    } else {
        warning("Running on minimal processed results, ",
                "but aggregated plot data is missing.")
    }
} else if (!file.exists(aggregatedFile) && createMinimalSaveFile) {
    partialPlotData_N_N <- subset(plotData_N_N, datasetIndex <= 3)
    save(partialPlotData_N_N,
         plotDataAggr_N_N,
         plotDataAggr_t3_t3,
         file = aggregatedFile)
}

########################################################
## plot results                                       ##
########################################################

plot_convergence_N_N_bias <-
    ggplot(
        plotDataAggr_N_N,
        aes(
            numberOfSubjects,
            mu,
            linetype = Method,
            shape = Method,
            color = numberOfReplicates
        )
    ) +
    lemon::geom_pointline(aes(group = interaction(Method, numberOfReplicates))) +
    geom_hline(yintercept = 1) +
    facet_grid(numberOfReplicates ~ variable) + ylab("Location")

if (interactive()) {
    print(plot_convergence_N_N_bias)
}

plot_convergence_N_N_scale <-
    ggplot(
        plotDataAggr_N_N,
        aes(
            numberOfSubjects,
            s,
            linetype = Method,
            shape = Method,
            color = numberOfReplicates
        )
    ) +
    lemon::geom_pointline(aes(group = interaction(Method, numberOfReplicates))) +
    facet_wrap( ~ variable) + ylab("Scale")

if (interactive()) {
    print(plot_convergence_N_N_scale)
}

plot_convergence_N_N_efficiency <-
    ggplot(
        subset(plotDataAggr_N_N, Method != "lme"),
        aes(
            numberOfSubjects,
            empiricalEfficiency,
            linetype = Method,
            shape = Method,
            color = numberOfReplicates
        )
    ) +
    lemon::geom_pointline(aes(group = interaction(Method, numberOfReplicates))) +
    facet_wrap( ~ variable) + ylab("Empirical efficiency")

if (interactive()) {
    print(plot_convergence_N_N_efficiency)
}


plot_convergence_t3_t3_bias <-
    plot_convergence_N_N_bias %+% plotDataAggr_t3_t3

if (interactive()) {
    print(plot_convergence_t3_t3_bias)
}

plot_convergence_t3_t3_scale <-
    plot_convergence_N_N_scale %+% plotDataAggr_t3_t3

if (interactive()) {
    print(plot_convergence_t3_t3_scale)
}

plot_convergence_t3_t3_efficiency <-
    plot_convergence_N_N_efficiency %+% subset(plotDataAggr_t3_t3, Method != "lme")

if (interactive()) {
    print(plot_convergence_t3_t3_efficiency)
}
