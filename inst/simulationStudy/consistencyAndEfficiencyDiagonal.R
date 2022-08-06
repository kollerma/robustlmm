########################################################
## Replicate figures 4.2, 4.3 and 4.4 in Koller 2013  ##
########################################################

require(robustlmm)
require(ggplot2)

## path where to find the processed results
path <- system.file("simulationStudy", package = "robustlmm")
## number of cores to use for fitting datasets
ncores <- parallel::detectCores()

########################################################
## generate datasets and configure fitting functions  ##
########################################################

set.seed(123)
datasets <- generateAnovaDatasets(
    numberOfDatasetsToGenerate = 1000,
    numberOfLevelsInFixedFactor = 3,
    numberOfSubjects = 10,
    numberOfReplicates = 5
)

fittingFunctions <- list(
    fitDatasets_lmer,
    fitDatasets_rlmer_DAStau_k_0_5,
    fitDatasets_rlmer_DAStau_k_0_5_noAdj,
    fitDatasets_rlmer_DAStau,
    fitDatasets_rlmer_DAStau_noAdj,
    fitDatasets_rlmer_DAStau_k_2,
    fitDatasets_rlmer_DAStau_k_2_noAdj,
    fitDatasets_rlmer_DAStau_k_5,
    fitDatasets_rlmer_DAStau_k_5_noAdj
)

########################################################
## fit all datasets or load and verify results        ##
########################################################

baseFilename <- "datasets_consistencyAndEfficiencyDiagonal"
createMinimalSaveFile <-
    path != system.file("simulationStudy", package = "robustlmm")
results <- processDatasetsInParallel(
    datasets,
    path,
    baseFilename,
    fittingFunctions,
    chunkSize = 10,
    checkProcessed = TRUE,
    createMinimalSaveFile = createMinimalSaveFile,
    ncores = ncores,
    meanB = TRUE,
    meanAbsB = TRUE
)

########################################################
## prepare data for plotting                          ##
########################################################

plotData <- cbind(
    data.frame(Method = as.factor(results$label)),
    results$coefficients,
    log(results$sigma),
    results$thetas,
    results$meanB,
    results$meanAbsB,
    results$datasetIndex
)
names(plotData)[-1] <-
    c("beta0",
      "beta1",
      "beta2",
      "log(sigma)",
      "theta",
      "meanB",
      "meanAbsB",
      "datasetIndex")
plotDataTmp <- reshape2::melt(plotData[-ncol(plotData)], 1)
plotDataAggr <-
    aggregate(plotDataTmp[["value"]], plotDataTmp[1:2], function(x)
        c(
            mean = mean(x),
            quantile(x, c(0.25, 0.75)),
            nExactZero = mean(x == 0),
            unlist(MASS::hubers(x, k = 1.345))
        ))
isNoAdj <- function(method)
    grepl("noAdj", method)
isLmer <- function(method)
    method == "fitDatasets_lmer"
plotDataAggr <- cbind(
    plotDataAggr[1:2],
    k = factor(extractPredefinedTuningParameter(
        as.character(plotDataAggr$Method), "rho.e"
    )),
    adjustment = factor(
        ifelse(
            isNoAdj(plotDataAggr$Method),
            "RSEn",
            ifelse(isLmer(plotDataAggr$Method),
                   "lme",
                   "RSEa")
        ),
        levels = c("lme",
                   "RSEn",
                   "RSEa")
    ),
    plotDataAggr[[3]]
)

plotDataLme <-
    subset(plotDataAggr, Method == "fitDatasets_lmer")
plotDataRse <- subset(plotDataAggr, Method != "fitDatasets_lmer")
plotDataTruth <-
    data.frame(
        variable = factor(
            levels(plotDataAggr$variable),
            levels(plotDataAggr$variable)
        ),
        value = c(
            datasets$trueBeta,
            log(datasets$trueSigma),
            datasets$trueTheta,
            NA_real_,
            NA_real_
        )
    )

plotDataMerged <-
    merge(
        subset(plotDataRse, variable %in% levels(variable)[1:5]),
        plotDataLme[-c(1, 3, 4)],
        by = "variable",
        suffixes = c(".RSE", ".lme")
    )
plotDataMerged <- droplevels(plotDataMerged)

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
         plotDataRse,
         plotDataTruth,
         plotDataLme,
         plotDataMerged,
         file = aggregatedFile)
}

########################################################
## Compute expected asymptotic efficiencies           ##
########################################################

## compute efficiencies
suppressWarnings(rm("asymptoticEfficiencies"))
for (method in levels(plotDataMerged$Method)) {
    if (isLmer(method)) {
        continue
    }
    k <- unname(extractPredefinedTuningParameter(method, "rho.e"))
    k.sigma <-
        unname(extractPredefinedTuningParameter(method, "rho.sigma.e"))
    rho <- chgDefaults(smoothPsi, k = k)
    efficiencyBetas <- asymptoticEfficiency(rho, "location")
    efficiencySigma <-
        asymptoticEfficiency(psi2propII(rho, k = k.sigma), "scale")
    tmp <- data.frame(
        k = k,
        adjustment = if (isNoAdj(method))
            "RSEn"
        else
            "RSEa",
        beta0 = efficiencyBetas,
        beta1 = efficiencyBetas,
        beta2 = efficiencyBetas,
        `log(sigma)` = efficiencySigma,
        theta = efficiencySigma
    )
    if (exists("asymptoticEfficiencies")) {
        asymptoticEfficiencies <- rbind(asymptoticEfficiencies, tmp)
    } else {
        asymptoticEfficiencies <- tmp
    }
}
asymptoticEfficiencies <-
    reshape2::melt(asymptoticEfficiencies, 1:2)
levels(asymptoticEfficiencies$variable)[4] <- "log(sigma)"
asymptoticEfficiencies$k <- factor(asymptoticEfficiencies$k)
asymptoticEfficiencies$adjustment <-
    factor(asymptoticEfficiencies$adjustment,
           levels(plotDataMerged$adjustment))

########################################################
## plot results                                       ##
########################################################

colors <- RColorBrewer::brewer.pal(8, "Dark2")[c(6:3, 7, 1, 2, 8)]

plot_consistencyDiagonal <-
    ggplot(plotDataRse, aes(color = adjustment)) +
    geom_hline(data = plotDataTruth, aes(yintercept = value)) +
    geom_hline(data = plotDataLme,
               aes(yintercept = `25%`, color = "lme"),
               linetype = 2) +
    geom_hline(data = plotDataLme,
               aes(yintercept = mean, color = "lme"),
               linetype = 2) +
    geom_hline(data = plotDataLme,
               aes(yintercept = `75%`, color = "lme"),
               linetype = 2) +
    geom_errorbar(
        aes(x = k, ymin = `25%`, ymax = `75%`),
        alpha = 0.6,
        width = 0.25,
        position = position_dodge(width = 0.4)
    ) +
    geom_point(aes(k, mean), position = position_dodge(width = 0.4)) +
    facet_wrap( ~ variable, scales = "free_y", nrow = 2) +
    ylab("mean and quartiles") +
    scale_colour_manual(
        "Method",
        values = colors[c(1, 3, 2)],
        breaks = levels(plotDataRse$adjustment),
        labels = levels(plotDataRse$adjustment)
    ) +
    theme(
        legend.position = c(0.95, 0.4),
        legend.justification = c("right", "top")
    )

if (interactive()) {
    print(plot_consistencyDiagonal)
}

plot_efficiencyDiagonal <-
    ggplot(plotDataMerged,
           aes(k, s.lme / s.RSE, color = adjustment, shape = adjustment)) +
    geom_hline(aes(yintercept = 1)) +
    geom_point(data = asymptoticEfficiencies, aes(y = value), color = "black") +
    lemon::geom_pointline(aes(group = adjustment)) + facet_wrap( ~ variable) +
    ylab("empirical efficiency") +
    scale_colour_manual("Method", values = colors[3:2]) +
    scale_shape_manual("Method", values = c(2, 16)) +
    theme(
        legend.position = c(0.95, 0.4),
        legend.justification = c("right", "top")
    )

if (interactive()) {
    print(plot_efficiencyDiagonal)
}
