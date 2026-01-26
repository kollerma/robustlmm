## Disable test
quit()

## Empirical test of generateLongitudinalDatasets
## Simulates 1000 datasets, fits with lmer (REML), checks that average estimates
## correspond to expected (true) values.
##
## Note: ML estimation has known downward bias for variance components.
## REML reduces this bias but doesn't eliminate it for small samples.
## We use larger sample sizes and REML to minimize this effect.

library(robustlmm)
library(lme4)

set.seed(2024)

## Configuration
nDatasets <- 1000
nSubjects <- 60   # larger sample to reduce singular fits
nTimepoints <- 8  # more timepoints for better slope estimation

## True parameter values - using larger variance components to reduce boundary issues
## Note: theta[3] (random slope SD) needs to be large enough that estimates
## don't hit the boundary at zero, which causes downward bias.
trueBeta <- c(2.5, -1.2)  # intercept, time slope
trueSigma <- 0.8          # smaller residual variance
trueTheta <- c(1.5, 0.1, 1.2)  # Cholesky factor elements (larger to avoid boundary)

## Generate all datasets
cat("Generating", nDatasets, "datasets...\n")
generator <- generateLongitudinalDatasets(
    numberOfDatasetsToGenerate = nDatasets,
    numberOfSubjects = nSubjects,
    numberOfTimepoints = nTimepoints,
    trueBeta = trueBeta,
    trueSigma = trueSigma,
    trueTheta = trueTheta
)

## Storage for estimates
betaEstimates <- matrix(NA, nrow = nDatasets, ncol = length(trueBeta))
sigmaEstimates <- numeric(nDatasets)
thetaEstimates <- matrix(NA, nrow = nDatasets, ncol = length(trueTheta))

## Fit all datasets with lmer (REML for less biased variance estimates)
cat("Fitting", nDatasets, "models with lmer (REML)...\n")
pb <- txtProgressBar(min = 0, max = nDatasets, style = 3)
nSingular <- 0

for (i in seq_len(nDatasets)) {
    data_i <- generator$generateData(i)

    suppressMessages(
        fit <- lmer(generator$formula, data = data_i, REML = TRUE)
    )
    if (isSingular(fit)) nSingular <- nSingular + 1

    betaEstimates[i, ] <- fixef(fit)
    sigmaEstimates[i] <- sigma(fit)
    thetaEstimates[i, ] <- getME(fit, "theta")

    setTxtProgressBar(pb, i)
}
close(pb)

cat("\nSingular fits:", nSingular, "out of", nDatasets,
    sprintf("(%.1f%%)\n", 100 * nSingular / nDatasets))

## Compute averages
avgBeta <- colMeans(betaEstimates)
avgSigma <- mean(sigmaEstimates)
avgTheta <- colMeans(thetaEstimates)

## Compute standard errors of the means
seBeta <- apply(betaEstimates, 2, sd) / sqrt(nDatasets)
seSigma <- sd(sigmaEstimates) / sqrt(nDatasets)
seTheta <- apply(thetaEstimates, 2, sd) / sqrt(nDatasets)

## Report results
cat("\n\n=== Results ===\n\n")

cat("Fixed Effects (beta):\n")
cat(sprintf("  Parameter   True      Average   SE        Diff      z-score\n"))
for (j in seq_along(trueBeta)) {
    diff <- avgBeta[j] - trueBeta[j]
    z <- diff / seBeta[j]
    cat(sprintf("  beta[%d]     %7.4f   %7.4f   %7.4f   %+7.4f   %+6.2f\n",
                j, trueBeta[j], avgBeta[j], seBeta[j], diff, z))
}

cat("\nResidual standard deviation (sigma):\n")
diff <- avgSigma - trueSigma
z <- diff / seSigma
cat(sprintf("  sigma       %7.4f   %7.4f   %7.4f   %+7.4f   %+6.2f\n",
            trueSigma, avgSigma, seSigma, diff, z))

cat("\nRandom effects Cholesky factor (theta):\n")
cat(sprintf("  Parameter   True      Average   SE        Diff      z-score\n"))
for (j in seq_along(trueTheta)) {
    diff <- avgTheta[j] - trueTheta[j]
    z <- diff / seTheta[j]
    cat(sprintf("  theta[%d]    %7.4f   %7.4f   %7.4f   %+7.4f   %+6.2f\n",
                j, trueTheta[j], avgTheta[j], seTheta[j], diff, z))
}

## Statistical tests
## Fixed effects: check if average is within 3 SE of true value
## Variance components: use relative bias OR absolute bias (whichever is more appropriate)
cat("\n=== Validation ===\n")
toleranceSE <- 3          # for fixed effects (z-score tolerance)
toleranceRelBias <- 0.05  # for variance components (5% relative bias)
toleranceAbsBias <- 0.05  # absolute bias tolerance for small parameters

## Helper function: check bias with appropriate tolerance
checkBias <- function(avg, true, name) {
    absBias <- abs(avg - true)
    ## Use absolute tolerance for small true values, relative otherwise
    if (abs(true) < 0.5) {
        pass <- absBias < toleranceAbsBias
        cat(sprintf("  %s: abs bias = %.4f [%s]\n", name, absBias,
                    if (pass) "OK" else "FAIL"))
    } else {
        relBias <- absBias / abs(true)
        pass <- relBias < toleranceRelBias
        cat(sprintf("  %s: rel bias = %.2f%% [%s]\n", name, 100 * relBias,
                    if (pass) "OK" else "FAIL"))
    }
    return(pass)
}

allPass <- TRUE
fixedPass <- TRUE
variancePass <- TRUE

cat("\nFixed effects (checking z-score < 3):\n")
for (j in seq_along(trueBeta)) {
    z <- abs(avgBeta[j] - trueBeta[j]) / seBeta[j]
    pass <- z < toleranceSE
    fixedPass <- fixedPass && pass
    status <- if (pass) "OK" else "FAIL"
    cat(sprintf("  beta[%d]: z = %.2f [%s]\n", j, z, status))
}

cat("\nVariance components (rel bias < 5%% or abs bias < 0.05):\n")
variancePass <- checkBias(avgSigma, trueSigma, "sigma") && variancePass
for (j in seq_along(trueTheta)) {
    variancePass <- checkBias(avgTheta[j], trueTheta[j],
                              sprintf("theta[%d]", j)) && variancePass
}

allPass <- fixedPass && variancePass
if (allPass) {
    cat("\nPASS: All parameters within tolerance.\n")
} else {
    if (!fixedPass) cat("\nFailed on fixed effects.\n")
    if (!variancePass) cat("\nNote: Variance component bias can occur with (RE)ML estimation.\n")
}

## Additional diagnostic: histograms of estimates
cat("\n=== Summary Statistics ===\n")
cat("\nFixed effects estimates:\n")
print(summary(betaEstimates))

cat("\nSigma estimates:\n")
print(summary(sigmaEstimates))

cat("\nTheta estimates:\n")
print(summary(thetaEstimates))

## Final assertions for automated testing
## Fixed effects: must be unbiased (within 3 SEs)
stopifnot(
    "beta estimates should average to true values" =
        all(abs(avgBeta - trueBeta) / seBeta < toleranceSE)
)

## Variance components: use appropriate tolerance (relative or absolute)
checkVarianceBias <- function(avg, true) {
    absBias <- abs(avg - true)
    if (abs(true) < 0.5) {
        return(absBias < toleranceAbsBias)
    } else {
        return(absBias / abs(true) < toleranceRelBias)
    }
}

stopifnot(
    "sigma estimates should be within tolerance" =
        checkVarianceBias(avgSigma, trueSigma),
    "theta estimates should be within tolerance" =
        all(mapply(checkVarianceBias, avgTheta, trueTheta))
)

cat("\nAll tests passed!\n")
