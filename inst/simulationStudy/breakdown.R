########################################################
## Replicate Figure 4.5 in Koller 2013                ##
########################################################

require(robustlmm)
require(ggplot2)

## path where to find the processed results
path <- system.file("simulationStudy", package = "robustlmm")

########################################################
## generate datasets and configure fitting function   ##
########################################################

set.seed(-14561352)
datasets_base <- generateAnovaDatasets(
    numberOfDatasetsToGenerate = 1,
    numberOfLevelsInFixedFactor = 1,
    numberOfSubjects = 10,
    numberOfReplicates = 5
)
oneWay <- datasets_base$generateData(1)

baseFit <-
    processFit(fitDatasets_rlmer_DAStau_lmerNoFit(datasets_base)[[1]])

idx <-
    rep(
        datasets_base$numberOfSubjects * (1:datasets_base$numberOfReplicates - 1),
        datasets_base$numberOfSubjects
    ) +
    rep(1:datasets_base$numberOfSubjects, each = datasets_base$numberOfReplicates)

datasetList <- list(oneWay)
for (i in 1:9) {
    oneWay$y[idx[i]] <- abs(oneWay$y[idx[i]]) * 1000000
    datasetList <- c(datasetList, list(oneWay))
}

datasets <- createDatasetsFromList(
    datasetList,
    datasets_base$formula,
    trueBeta = baseFit$coefficients,
    trueSigma = baseFit$sigma,
    trueTheta = baseFit$thetas
)

fittingFunctions <- list(fitDatasets_rlmer_DAStau_lmerNoFit)

########################################################
## fit all datasets or load and verify results        ##
########################################################

filename <- "datasets_breakdown.Rdata"
results <- processFile(
    file.path(path, filename),
    fittingFunctions,
    checkProcessed = TRUE,
    datasets = datasets,
    b = 1
)

########################################################
## prepare data for plotting                          ##
########################################################

plotData <-
    data.frame(
        epsilon = (results$datasetIndex[, 1] - 1) / datasets$numberOfRows,
        beta = results$coefficients[, 1],
        sigma = results$sigma[, 1],
        theta = results$thetas[, 1],
        b1 = results$b[, 1]
    )
plotData <- reshape2::melt(plotData, 1)

plot_breakdown <- ggplot(plotData, aes(epsilon, value)) +
    geom_line() + scale_y_log10() +
    coord_cartesian(ylim = c(0.07, 100)) +
    xlab("amount of contamination epsilon") +
    facet_wrap( ~ variable)

if (interactive()) {
    print(plot_breakdown)
}
