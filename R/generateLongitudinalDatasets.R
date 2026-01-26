##' Generate balanced longitudinal datasets with random intercepts and slopes.
##' Subjects are observed at multiple time points with optional treatment groups.
##' Treatment and its interaction with time are coded as contrasts relative to
##' the first level.
##'
##' The generated data follows the model:
##' \deqn{y_{ij} = \beta_0 + \beta_1 \cdot \text{time}_{ij} +
##' \sum_{k=1}^{K-1} \beta_{1+k} \cdot \text{treatment}_{k,i} +
##' \sum_{k=1}^{K-1} \beta_{K+k} \cdot \text{treatment}_{k,i} \cdot
##' \text{time}_{ij} + b_{0i} + b_{1i} \cdot \text{time}_{ij} +
##' \epsilon_{ij}}{y_ij = beta_0 + beta_1 * time_ij +
##'                sum_k beta_{1+k} * treat_{k,i} +
##'                sum_k beta_{K+k} * treat_{k,i} * time_ij +
##'                b_0i + b_1i * time_ij + eps_ij}
##'
##' where \eqn{K} is the number of treatment groups,
##'  \eqn{b_i = (b_{0i}, b_{1i})^T \sim N(0, \sigma^2 \Lambda \Lambda^T)}{b_i ~ N(0, sigma^2 * Lambda * Lambda^T)}
##' with \eqn{\Lambda}{\Lambda} being the lower-triangular Cholesky factor
##' reconstructed from the theta vector.
##'
##' The theta parameterization follows lme4 conventions:
##' \itemize{
##'   \item For a 2x2 random-effects covariance structure (intercept and slope),
##'         theta has 3 elements:
##'         \eqn{\theta = (\lambda_{11}, \lambda_{21}, \lambda_{22})}{theta = (lambda_11, lambda_21, lambda_22)}
##'   \item The Cholesky factor is:
##'         \eqn{\Lambda = \begin{pmatrix} \lambda_{11} & 0 \\ \lambda_{21} & \lambda_{22} \end{pmatrix}}{Lambda = [lambda_11, 0; lambda_21, lambda_22]}
##' }
##'
##' In order to save memory, only the generated random effects and the errors
##' are stored. The dataset is only created on demand when the method
##' \code{generateData} in the returned list is evaluated.
##'
##' The random variables are generated in a way that one can simulate more
##' datasets easily. When starting from the same seed, the first generated
##' datasets will be the same as for a previous call of
##' \code{generateLongitudinalDatasets} with a smaller number of datasets to
##' generate, see examples.
##'
##' @title Generate Longitudinal Datasets
##' @param numberOfDatasetsToGenerate number of datasets to generate.
##' @param numberOfSubjects number of subjects per dataset.
##' @param numberOfTimepoints number of observation time points per subject.
##' @param numberOfTreatmentLevels number of treatment levels. Default: 1
##'   (no treatment effect, intercept and time only).
##' @param timeRange numeric vector of length 2, range of time values (min, max).
##'   Default: \code{c(0, 1)}.
##' @param errorGenerator random number generator used for the errors. Called as
##'   \code{errorGenerator(n) * trueSigma}.
##' @param randomEffectGenerator random number generator used for the spherical
##'   random effects. Called as \code{randomEffectGenerator(n) * trueSigma}.
##' @param trueBeta scalar or vector with the true values of the fixed effects
##'   coefficients. Can be of length one in which case it will be replicated to
##'   the required length. The order is: intercept, time, treatment contrasts
##'   (if any), treatment-by-time interactions (if any).
##' @param trueSigma scalar with the true value of the error scale.
##' @param trueTheta numeric vector of length 3 with the true values for the
##'   Cholesky factor of the random effects covariance matrix (lme4 convention).
##'   Default: \code{c(1, 0, 1)} (independent random intercepts and slopes).
##' @param contamFun optional contamination function. If provided, it receives
##'   the full dataset (a data frame with columns id, time, treatment, y) and
##'   an info list, and must return the (possibly modified) data frame. This
##'   allows arbitrary contamination including changing group assignments. See
##'   Details for the contents of the info list.
##' @param ... all additional arguments are added to the returned list.
##' @author Manuel Koller
##' @return list with generators and the original arguments
##' \item{\code{generateData}: }{function to generate data taking one argument,
##'   the dataset index.}
##' \item{\code{createXMatrix}: }{function to generate X matrix taking one
##'   argument, the result of \code{generateData}.}
##' \item{\code{createZMatrix}: }{function to generate Z matrix taking one
##'   argument, the result of \code{generateData}.}
##' \item{\code{createLambdaMatrix}: }{function to generate Lambda matrix taking
##'   one argument, the result of \code{generateData}.}
##' \item{\code{randomEffects}: }{function to return the generated random effects
##'   taking one argument, the dataset index.}
##' \item{\code{sphericalRandomEffects}: }{function to return the generated
##'   spherical random effects taking one argument, the dataset index.}
##' \item{\code{errors}: }{function to return the generated errors taking one
##'   argument, the dataset index.}
##' \item{\code{allRandomEffects}: }{function without arguments that returns the
##'   matrix of all generated random effects.}
##' \item{\code{allErrors}: }{function without arguments that returns the matrix
##'   of all generated errors.}
##' \item{\code{numberOfDatasets}: }{\code{numberOfDatasetsToGenerate} as supplied}
##' \item{\code{numberOfSubjects}: }{\code{numberOfSubjects} as supplied}
##' \item{\code{numberOfTimepoints}: }{\code{numberOfTimepoints} as supplied}
##' \item{\code{numberOfTreatmentLevels}: }{\code{numberOfTreatmentLevels} as supplied}
##' \item{\code{numberOfRows}: }{number of rows in the generated dataset}
##' \item{\code{trueBeta}: }{true values used for beta}
##' \item{\code{trueSigma}: }{true value used for sigma}
##' \item{\code{trueTheta}: }{true values used for theta}
##' \item{\code{formula}: }{formula to fit the model using \code{lmer}}
##' \item{\code{...}: }{additional arguments passed via \code{...}}
##'
##' @details
##' \strong{Treatment Assignment:}
##' Subjects are assigned to treatment groups in a balanced, deterministic manner.
##' Subject i is assigned to treatment \code{(i - 1) mod numberOfTreatmentLevels + 1}.
##'
##' \strong{Contamination Function:}
##' If \code{contamFun} is provided, it is called as \code{contamFun(data, info)}
##' after the response \code{y} is computed. The \code{info} list contains:
##' \itemize{
##'   \item \code{datasetIndex}: the dataset index
##'   \item \code{randomEffects}: the random effects vector for this dataset
##'   \item \code{errors}: the error vector for this dataset
##'   \item \code{trueBeta}: as passed to \code{generateLongitudinalDatasets}
##'   \item \code{trueSigma}: as passed to \code{generateLongitudinalDatasets}
##'   \item \code{trueTheta}: as passed to \code{generateLongitudinalDatasets}
##' }
##' The function must return a data frame with the same structure (columns id,
##' time, treatment, y). This allows arbitrary modifications including:
##' \itemize{
##'   \item Modifying the response \code{y} (e.g., adding outliers)
##'   \item Changing group assignments (e.g., moving subjects between treatments)
##'   \item Modifying time values
##'   \item Any combination of the above
##' }
##'
##' @seealso \code{\link{generateAnovaDatasets}}, \code{\link{generateMixedEffectDatasets}}
##' @examples
##'   oneGroup <- generateLongitudinalDatasets(2, 10, 5)
##'   head(oneGroup$generateData(1))
##'   head(oneGroup$generateData(2))
##'   oneGroup$formula
##'
##'   twoGroups <- generateLongitudinalDatasets(2, 20, 5, numberOfTreatmentLevels = 2)
##'   head(twoGroups$generateData(1))
##'   twoGroups$formula
##'
##'   ## illustration how to generate more datasets
##'   set.seed(1)
##'   datasets1 <- generateLongitudinalDatasets(2, 10, 5)
##'   set.seed(1)
##'   datasets2 <- generateLongitudinalDatasets(3, 10, 5)
##'   stopifnot(all.equal(datasets1$generateData(1), datasets2$generateData(1)),
##'             all.equal(datasets1$generateData(2), datasets2$generateData(2)))
##'
##'   ## contamination example: add outliers to 10% of observations
##'   set.seed(42)
##'   contam <- generateLongitudinalDatasets(
##'     numberOfDatasetsToGenerate = 5,
##'     numberOfSubjects = 20,
##'     numberOfTimepoints = 5,
##'     contamFun = function(data, info) {
##'       n <- nrow(data)
##'       idx <- sample(n, size = ceiling(0.1 * n))
##'       data$y[idx] <- data$y[idx] + 10
##'       data
##'     }
##'   )
##'   head(contam$generateData(1))
##'
##'   ## contamination example: reassign some subjects to different treatment
##'   set.seed(42)
##'   contamGroup <- generateLongitudinalDatasets(
##'     numberOfDatasetsToGenerate = 5,
##'     numberOfSubjects = 20,
##'     numberOfTimepoints = 5,
##'     numberOfTreatmentLevels = 2,
##'     contamFun = function(data, info) {
##'       ## move first subject from T1 to T2
##'       data$treatment[data$id == 1] <- "T2"
##'       data
##'     }
##'   )
##'   head(contamGroup$generateData(1), 10)
##'
##'   ## medsim: simulation inspired by the medication dataset from confintROB
##'   ## Two subjects from treatment are mislabeled as control, and responses
##'   ## are truncated at a measurement floor of 100.
##'   contaminateMedsim <- function(data, info) {
##'     data$y <- pmax(data$y, 100)  # measurement floor
##'     data$treatment[data$id %in% c("2", "4")] <- "T1"
##'     data
##'   }
##'   set.seed(2000)
##'   medsim <- generateLongitudinalDatasets(
##'     numberOfDatasetsToGenerate = 100,
##'     numberOfSubjects = 60,
##'     numberOfTimepoints = 7,
##'     numberOfTreatmentLevels = 2,
##'     timeRange = c(0, 18),
##'     trueBeta = c(240, -3.11, -2.42, 4.00),
##'     trueSigma = sqrt(1229.93),
##'     trueTheta = c(1.310266, -0.07547461, 0.2147735),
##'     contamFun = contaminateMedsim
##'   )
##'   head(medsim$generateData(1))
##'
##' @export
##' @importFrom stats rnorm as.formula
generateLongitudinalDatasets <- function(numberOfDatasetsToGenerate,
                                         numberOfSubjects,
                                         numberOfTimepoints,
                                         numberOfTreatmentLevels = 1L,
                                         timeRange = c(0, 1),
                                         errorGenerator = rnorm,
                                         randomEffectGenerator = rnorm,
                                         trueBeta = 0,
                                         trueSigma = 1,
                                         trueTheta = c(1, 0, 1),
                                         contamFun = NULL,
                                         ...) {
    checkIsScalar(numberOfDatasetsToGenerate, "numberOfDatasetsToGenerate")
    checkIsScalar(numberOfSubjects, "numberOfSubjects")
    checkIsScalar(numberOfTimepoints, "numberOfTimepoints")
    checkIsScalar(numberOfTreatmentLevels, "numberOfTreatmentLevels")
    checkIsScalar(trueSigma, "trueSigma")
    checkLongitudinalTreatmentLevels(numberOfTreatmentLevels)
    checkLongitudinalSubjects(numberOfSubjects)
    checkLongitudinalTimepoints(numberOfTimepoints)
    checkLongitudinalTimeRange(timeRange)
    checkLongitudinalTheta(trueTheta)

    numberOfRows <- numberOfSubjects * numberOfTimepoints
    numberOfFixedEffectCoefficients <- computeNumberOfFixedEffects(numberOfTreatmentLevels)

    trueBeta <- checkLengthOrReplicate(trueBeta, "trueBeta",
                                       numberOfFixedEffectCoefficients)
    trueTheta <- checkLengthOrReplicate(trueTheta, "trueTheta", 3L)

    randomVariables <-
        generateAndAssembleLongitudinalRandomVariables(
            numberOfDatasetsToGenerate,
            numberOfSubjects,
            numberOfRows,
            trueSigma,
            trueTheta,
            randomEffectGenerator,
            errorGenerator
        )

    return(list(
        generateData = function(datasetIndex) {
            data <- generateLongitudinalData(numberOfSubjects,
                                             numberOfTimepoints,
                                             numberOfTreatmentLevels,
                                             timeRange)
            X <- createLongitudinalXMatrix(data, numberOfTreatmentLevels)
            Z <- createLongitudinalZMatrix(data)
            data$y <- as(
                X %*% trueBeta +
                    Z %*% randomVariables[["randomEffects"]][, datasetIndex] +
                    randomVariables[["errors"]][, datasetIndex],
                "vector"
            )
            if (!is.null(contamFun)) {
                contamInfo <- list(
                    datasetIndex = datasetIndex,
                    randomEffects = randomVariables[["randomEffects"]][, datasetIndex],
                    errors = randomVariables[["errors"]][, datasetIndex],
                    trueBeta = trueBeta,
                    trueSigma = trueSigma,
                    trueTheta = trueTheta
                )
                data <- contamFun(data, contamInfo)
            }
            attr(data, "datasetIndex") <-
                datasetIndex + randomVariables[["datasetIndexOffset"]]
            return(data)
        },
        createXMatrix = function(data) {
            return(createLongitudinalXMatrix(data, numberOfTreatmentLevels))
        },
        createZMatrix = function(data) {
            return(createLongitudinalZMatrix(data))
        },
        createLambdaMatrix = function(data) {
            return(createLongitudinalLambda(trueTheta, numberOfSubjects))
        },
        randomEffects = function(datasetIndex) {
            return(randomVariables[["randomEffects"]][, datasetIndex])
        },
        sphericalRandomEffects = function(datasetIndex) {
            Lambda <- createLongitudinalLambda(trueTheta, numberOfSubjects)
            sphericalRandomEffects <-
                solve(Lambda, randomVariables[["randomEffects"]][, datasetIndex])
            return(sphericalRandomEffects)
        },
        errors = function(datasetIndex) {
            return(randomVariables[["errors"]][, datasetIndex])
        },
        allRandomEffects = function() {
            return(randomVariables[["randomEffects"]])
        },
        allErrors = function() {
            return(randomVariables[["errors"]])
        },
        numberOfDatasets = numberOfDatasetsToGenerate,
        numberOfSubjects = numberOfSubjects,
        numberOfTimepoints = numberOfTimepoints,
        numberOfTreatmentLevels = numberOfTreatmentLevels,
        numberOfRows = numberOfRows,
        timeRange = timeRange,
        trueBeta = trueBeta,
        trueSigma = trueSigma,
        trueTheta = trueTheta,
        formula = generateLongitudinalFormula(numberOfTreatmentLevels),
        ...
    ))
}

checkLongitudinalTreatmentLevels <- function(numberOfTreatmentLevels) {
    if (numberOfTreatmentLevels < 1L) {
        stop("Argument 'numberOfTreatmentLevels' must be at least 1.")
    }
}

checkLongitudinalSubjects <- function(numberOfSubjects) {
    if (numberOfSubjects < 2L) {
        stop("Argument 'numberOfSubjects' must be at least 2.")
    }
}

checkLongitudinalTimepoints <- function(numberOfTimepoints) {
    if (numberOfTimepoints < 2L) {
        stop("Argument 'numberOfTimepoints' must be at least 2.")
    }
}

checkLongitudinalTimeRange <- function(timeRange) {
    if (length(timeRange) != 2L) {
        stop("Argument 'timeRange' must be a vector of length 2.")
    }
    if (timeRange[1] >= timeRange[2]) {
        stop("Argument 'timeRange' must have timeRange[1] < timeRange[2].")
    }
}

checkLongitudinalTheta <- function(trueTheta) {
    if (length(trueTheta) != 3L && length(trueTheta) != 1L) {
        stop("Argument 'trueTheta' must have length 1 or 3 for random intercept + slope model.")
    }
}

computeNumberOfFixedEffects <- function(numberOfTreatmentLevels) {
    ## intercept + time + (K-1) treatment contrasts + (K-1) interactions
    1L + 1L + max(0L, numberOfTreatmentLevels - 1L) +
        max(0L, numberOfTreatmentLevels - 1L)
}

generateLongitudinalData <- function(numberOfSubjects,
                                     numberOfTimepoints,
                                     numberOfTreatmentLevels,
                                     timeRange) {
    subjectIds <- rep(seq_len(numberOfSubjects), each = numberOfTimepoints)
    times <- rep(seq(timeRange[1], timeRange[2], length.out = numberOfTimepoints),
                 times = numberOfSubjects)
    treatmentAssignment <- rep(
        ((seq_len(numberOfSubjects) - 1L) %% numberOfTreatmentLevels) + 1L,
        each = numberOfTimepoints
    )
    data.frame(
        id = factor(subjectIds),
        time = times,
        treatment = factor(treatmentAssignment,
                           levels = seq_len(numberOfTreatmentLevels),
                           labels = paste0("T", seq_len(numberOfTreatmentLevels)))
    )
}

createLongitudinalXMatrix <- function(data, numberOfTreatmentLevels) {
    if (numberOfTreatmentLevels == 1L) {
        X <- cbind(
            Intercept = 1,
            time = data$time
        )
    } else {
        treatMatrix <- model.matrix(~ treatment - 1, data = data)
        treatContrasts <- treatMatrix[, -1, drop = FALSE]
        colnames(treatContrasts) <- paste0("treatment",
                                           seq(2L, numberOfTreatmentLevels))
        treatByTimeInteractions <- treatContrasts * data$time
        colnames(treatByTimeInteractions) <- paste0("time:treatment",
                                                    seq(2L, numberOfTreatmentLevels))
        X <- cbind(
            Intercept = 1,
            time = data$time,
            treatContrasts,
            treatByTimeInteractions
        )
    }
    return(X)
}

createLongitudinalZMatrix <- function(data) {
    subjects <- unique(data$id)
    numberOfSubjects <- length(subjects)
    numberOfRows <- nrow(data)
    Z <- matrix(0, nrow = numberOfRows, ncol = 2L * numberOfSubjects)
    for (i in seq_len(numberOfSubjects)) {
        rowsI <- which(data$id == subjects[i])
        Z[rowsI, 2L * i - 1L] <- 1
        Z[rowsI, 2L * i] <- data$time[rowsI]
    }
    return(Z)
}

createLongitudinalLambda <- function(trueTheta, numberOfSubjects) {
    dim <- 2L * numberOfSubjects
    Lambda <- matrix(0, nrow = dim, ncol = dim)
    cholBlock <- matrix(
        c(trueTheta[1], 0,
          trueTheta[2], trueTheta[3]),
        nrow = 2L,
        ncol = 2L,
        byrow = TRUE
    )
    for (i in seq_len(numberOfSubjects)) {
        rows <- (2L * i - 1L):(2L * i)
        Lambda[rows, rows] <- cholBlock
    }
    return(Lambda)
}

generateAndAssembleLongitudinalRandomVariables <-
    function(numberOfDatasetsToGenerate,
             numberOfSubjects,
             numberOfRows,
             trueSigma,
             trueTheta,
             randomEffectGenerator,
             errorGenerator) {
        Lambda <- createLongitudinalLambda(trueTheta, numberOfSubjects)
        numberOfRandomEffects <- 2L * numberOfSubjects
        generatedValues <-
            replicate(
                numberOfDatasetsToGenerate,
                generateOneLongitudinalSetOfRandomVariables(
                    numberOfRandomEffects,
                    numberOfRows,
                    trueSigma,
                    Lambda,
                    randomEffectGenerator,
                    errorGenerator
                ),
                simplify = FALSE
            )
        return(list(
            randomEffects = sapply(generatedValues, `[[`, "randomEffects"),
            errors = sapply(generatedValues, `[[`, "errors"),
            datasetIndexOffset = 0
        ))
    }

generateOneLongitudinalSetOfRandomVariables <-
    function(numberOfRandomEffects,
             numberOfRows,
             trueSigma,
             Lambda,
             randomEffectGenerator,
             errorGenerator) {
        sphericalRandomEffects <- randomEffectGenerator(numberOfRandomEffects) * trueSigma
        randomEffects <- drop(Lambda %*% sphericalRandomEffects)
        errors <- errorGenerator(numberOfRows) * trueSigma
        return(list(
            randomEffects = randomEffects,
            errors = errors
        ))
    }

generateLongitudinalFormula <- function(numberOfTreatmentLevels) {
    if (numberOfTreatmentLevels == 1L) {
        as.formula("y ~ 1 + time + (1 + time | id)")
    } else {
        as.formula("y ~ 1 + time + treatment + time:treatment + (1 + time | id)")
    }
}
