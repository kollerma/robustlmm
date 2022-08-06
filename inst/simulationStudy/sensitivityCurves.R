########################################################
## Code to replicate Figure 4.1 in Koller 2013        ##
## and Figure 1 in Koller and Stahel 2022           ##
########################################################

require(robustlmm)
require(ggplot2)

## path where to find the processed results
path <- system.file("simulationStudy", package = "robustlmm")

########################################################
## generate base dataset and some configuration       ##
########################################################
set.seed(4)
datasets_base <- generateAnovaDatasets(
    numberOfDatasetsToGenerate = 1,
    numberOfLevelsInFixedFactor = 1,
    numberOfSubjects = 10,
    numberOfReplicates = 5
)
oneWay <- datasets_base$generateData(1)

ylim <- c(-1, 1) * 1.5 * max(abs(oneWay$y))
slim <- c(0.1, 30)
npoints <- 31

########################################################
## helper functions                                   ##
########################################################

fitAndProcess <- function(datasets) {
    basename <- deparse(substitute(datasets))
    file <- paste0(file.path(path, basename), ".Rdata")
    results <-
        processFile(
            file = file,
            fittingFunctions = list(
                fitDatasets_lmer,
                fitDatasets_rlmer_DAStau,
                fitDatasets_rlmer_DAStau_noAdj
            ),
            checkProcessed = TRUE,
            datasets = datasets,
            b = 1
        )
    return(results)
}

createPlotData <- function(results, changes) {
    method <- as.factor(results$label)
    levels(method) <- shortenLabelsKS2022(levels(method))
    method <- factor(method, c("lme", "RSEn", "RSEa"))
    plotData <- data.frame(Changes = changes[results$datasetIndex],
                           Method = method)
    plotData <- cbind(
        plotData,
        results$coefficients,
        results$sigma,
        results$thetas * results$sigma,
        results$b
    )
    names(plotData)[3:6] <- c("beta0", "sigma", "B.sigma", "b1")
    plotData <- reshape2::melt(plotData, 1:2)
    return(plotData)
}

colors <- RColorBrewer::brewer.pal(8, "Dark2")[c(6:3, 7, 1, 2, 8)]

########################################################
## Create curves for shifting the first observation   ##
########################################################
shifts <- seq(ylim[1], ylim[2], length.out = npoints)

datasets_shiftFirstObservation <- generateSensitivityCurveDatasets(
    data = oneWay,
    observationsToChange = 1,
    shifts = shifts,
    formula = datasets_base$formula
)

results_shiftFirstObservation <-
    fitAndProcess(datasets_shiftFirstObservation)

plotData_shiftFirstObservation <-
    createPlotData(results_shiftFirstObservation, shifts)

basePlot <- ggplot(plotData_shiftFirstObservation,
                   aes(Changes, value)) +
    geom_hline(
        data = subset(plotData_shiftFirstObservation, Changes == 0 &
                          Method == "lmer"),
        aes(yintercept = value),
        color = "gray"
    ) +
    geom_line(aes(color = Method, linetype = Method)) +
    facet_wrap( ~ variable, scales = "free_y", nrow = 1) +
    scale_colour_manual(
        "Method",
        values = colors[c(1, 3, 2)],
        breaks = levels(plotData_shiftFirstObservation$Method),
        labels = levels(plotData_shiftFirstObservation$Method)
    ) +
    theme(legend.position = "top")

plot_shiftFirstObservation <-
    basePlot + xlab("shift") + ggtitle("(A) shift of first observation")

if (interactive()) {
    print(plot_shiftFirstObservation)
}

########################################################
## Create curves for shifting the first group         ##
########################################################
shifts <- seq(ylim[1], ylim[2], length.out = npoints)

datasets_shiftFirstGroup <- generateSensitivityCurveDatasets(
    data = oneWay,
    observationsToChange = oneWay$Var2 == levels(oneWay$Var2)[1],
    shifts = shifts,
    formula = datasets_base$formula
)

results_shiftFirstGroup <-
    fitAndProcess(datasets_shiftFirstGroup)

plotData_shiftFirstGroup <-
    createPlotData(results_shiftFirstGroup, shifts)

plot_shiftFirstGroup <-
    basePlot %+% plotData_shiftFirstGroup +
    xlab("shift") + ggtitle("(B) shift of first group")

if (interactive()) {
    print(plot_shiftFirstGroup)
}

########################################################
## Create curves for scaling the first group          ##
########################################################
scales <- seq(slim[1], slim[2], length.out = npoints)

datasets_scaleFirstGroup <- generateSensitivityCurveDatasets(
    data = oneWay,
    observationsToChange = oneWay$Var2 == levels(oneWay$Var2)[1],
    scales = scales,
    center = datasets_base$randomEffects(1)[1],
    formula = datasets_base$formula
)

results_scaleFirstGroup <-
    fitAndProcess(datasets_scaleFirstGroup)

plotData_scaleFirstGroup <-
    createPlotData(results_scaleFirstGroup, scales)

plotData_scaleFirstGroup$value[plotData_scaleFirstGroup$value > 10] <-
    NA_real_

plot_scaleFirstGroup <-
    basePlot %+% plotData_scaleFirstGroup + scale_x_log10() +
    xlab("scale") + ggtitle("(C) scaling of first group")

if (interactive()) {
    print(plot_scaleFirstGroup)
}
