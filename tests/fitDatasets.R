require(robustlmm)

isPackageInstalled <- function(package) {
    return(length(find.package(package, quiet = TRUE)) > 0)
}


set.seed(1)
oneWay <- generateAnovaDatasets(1, 1, 10, 4,
                                lmeFormula = y ~ 1,
                                heavyLmeRandom = ~ 1,
                                heavyLmeGroups = ~ Var2,
                                lqmmRandom = ~ 1,
                                lqmmGroup = "Var2",
                                groups = cbind(rep(1:4, each = 10), rep(1:10, 4)),
                                varcov = matrix(1, 4, 4),
                                lower = 0)
test <- function(result, expectedLength = 1) {
    stopifnot(length(result) == expectedLength,
              !is(result[[1]], "try-error"))
    print(result[[1]])
}
test(expected <- fitDatasets_lmer(oneWay))
test(fitDatasets_lmer_bobyqa(oneWay))
test(fitDatasets_lmer_Nelder_Mead(oneWay))
fitDatasets_rlmer_custom <- function(datasets) {
    return(fitDatasets_rlmer(datasets,
                             method = "DASvar",
                             tuningParameter = c(1.345, 2.28, 1.345, 2.28, 5.14, 5.14),
                             label = "fitDatasets_rlmer_custom"))
}
test(fitDatasets_rlmer_custom(oneWay))
test(fitDatasets_rlmer_DAStau(oneWay))
test(fitDatasets_rlmer_DAStau_lmerNoFit(oneWay))
test(fitDatasets_rlmer_DASvar(oneWay))
test(fitDatasets_rlmer_DAStau_noAdj(oneWay))
test(fitDatasets_rlmer_DAStau_k_0_5(oneWay))
test(fitDatasets_rlmer_DAStau_k_0_5_noAdj(oneWay))
test(fitDatasets_rlmer_DAStau_k_2(oneWay))
test(fitDatasets_rlmer_DAStau_k_2_noAdj(oneWay))
test(fitDatasets_rlmer_DAStau_k_5(oneWay))
test(fitDatasets_rlmer_DAStau_k_5_noAdj(oneWay))
if (isPackageInstalled("heavy")) {
    test(fitDatasets_heavyLme(oneWay))
}
if (isPackageInstalled("lqmm")) {
    test(fitDatasets_lqmm(oneWay))
}
## if (isPackageInstalled("rlme")) {
##     test(fitDatasets_rlme(oneWay)) # won't work as dataset is balanced
## }
if (isPackageInstalled("robustvarComp")) {
    fitDatasets_varComprob_custom <- function(datasets, postFit) {
        lcontrol <- robustvarComp::varComprob.control(lower = datasets[["lower"]])
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_custom",
                postFit = postFit
            )
        )
    }
    test(fitDatasets_varComprob_custom(oneWay))
    test(fitDatasets_varComprob_compositeTau(oneWay))
    test(fitDatasets_varComprob_compositeTau_OGK(oneWay))
    test(fitDatasets_varComprob_compositeTau_2SGS(oneWay))
    test(fitDatasets_varComprob_compositeS(oneWay))
    test(fitDatasets_varComprob_compositeS_OGK(oneWay))
    test(fitDatasets_varComprob_compositeS_2SGS(oneWay))
    test(fitDatasets_varComprob_S(oneWay))
    test(fitDatasets_varComprob_S_OGK(oneWay))
    test(fitDatasets_varComprob_S_2SGS(oneWay))
}

## test datasetIndices argument

set.seed(1)
oneWay3 <- generateAnovaDatasets(3, 1, 10, 4,
                                 lmeFormula = y ~ 1,
                                 heavyLmeRandom = ~ 1,
                                 heavyLmeGroups = ~ Var2,
                                 lqmmRandom = ~ 1,
                                 lqmmGroup = "Var2",
                                 groups = cbind(rep(1:4, each = 10), rep(1:10, 4)),
                                 varcov = matrix(1, 4, 4),
                                 lower = 0)

test(actual <- fitDatasets_lmer(oneWay3, datasetIndices = 1))
stopifnot(all.equal(expected, actual, check.attributes = FALSE))

test(fitDatasets_lmer(oneWay3, datasetIndices = 2:3), expectedLength = 2)

## test rlme implementation
if (require(rlme)) {

    data(schools)
    datasets <- createDatasetsFromList(list(schools),
                                       formula = y ~ 1 + sex + age + (1 | region) + (1 | region:school),
                                       trueBeta = rep(NA_real_, 3),
                                       trueSigma = NA_real_,
                                       trueTheta = rep(NA_real_, 2))
    fits <- fitDatasets_rlme(datasets)
    result <- processFit(fits[[1]], all = TRUE)
    stopifnot(length(result$coefficients) == 3,
              length(result$standardErrors) == 3,
              length(result$tValues) == 3,
              length(result$sigma) == 1,
              length(result$thetas) == 2,
              length(result$residuals) == NROW(schools))
}

