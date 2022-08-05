##' Generate balanced datasets with multiple factors. All combinations of all
##' factor variables are generated, i.e., a fully crossed dataset will be
##' generated. \code{numberOfReplicates} specifies the number of replications
##' per unique combination.
##'
##' \code{numberOfLevelsInFixedFactor} can either be a scalar or a vector with
##' the number of levels for each fixed effects group. If
##' \code{numberOfLevelsInFixedFactor} is a scalar, the value of \code{1} is
##' allowed. This can be used to generate a dataset with an intercept only. If
##' \code{numberOfLevelsInFixedFactor} is a vector with more than one entry,
##' then all the values need to be larger than one.
##'
##' \code{numberOfSubjects} can also be a scalar of a vector with the number of
##' levels for each variance component. Each group needs to have more than one
##' level. The vector is sorted descending before the names are assigned. This
##' ensures that, when running \code{lmer}, the order of the random effects does
##' not change. \code{lmer} also sorts the random effects by decending number of
##' levels.
##'
##' In order to save memory, only the generated random effects and the errors
##' are stored. The dataset is only created on demand when the method
##' \code{generateData} in the returned list is evaluated.
##'
##' The random variables are generated in a way that one can simulate more
##' datasets easily. When starting from the same seed, the first generated
##' datasets will be the same as for the a previous call of
##' \code{generateAnovaDatasets} with a smaller number of datasets to generate,
##' see examples.
##'
##' The package \code{MatrixModels} needs to be installed for
##' \code{generateAnovaDatasets} to work.
##' @title Generate ANOVA type datasets
##' @param numberOfDatasetsToGenerate number of datasets to generate.
##' @param numberOfLevelsInFixedFactor scalar or vector with the number of
##'   levels per fixed factor or grouping variable.
##' @param numberOfSubjects scalar or vector with the number of levels per
##'   variance component.
##' @param numberOfReplicates number of replicates per unique combination of
##'   fixed factor and variance component.
##' @param errorGenerator random number generator used for the errors.
##' @param randomEffectGenerator random number generator used for the spherical
##'   random effects.
##' @param trueBeta scalar or vector with the true values of the fixed effects
##'   coefficients. Can be of length one in which case it will be replicated to
##'   the required length if needed.
##' @param trueSigma scalar with the true value of the error scale.
##' @param trueTheta scalar of vector with the true values for the variance
##'   component coefficients, not including sigma. Can be of length one in which
##'   case it will be replicated to the required length if needed.
##' @param ... all additional arguments are added to the returned list.
##' @param arrange If \code{TRUE}, the observations in the dataset are arranged
##'   such that the call to \code{\link[dplyr]{arrange}} in
##'   \code{\link[robustvarComp]{varComprob}} does not break the observation-
##'   group relationship. This requires package dplyr to be installed.
##' @author Manuel Koller
##' @return list with generators and the original arguments
##' \describe{
##' \item{\code{generateData}: }{function to generate data taking one argument, the dataset index.}
##' \item{\code{createXMatrix}: }{function to generate X matrix taking one argument,
##'   the result of \code{generateData}.}
##' \item{\code{createZMatrix}: }{function to generate Z matrix taking one argument,
##'   the result of \code{generateData}.}
##' \item{\code{createLambdaMatrix}: }{function to generate Lambda matrix taking one
##'   argument, the result of \code{generateData}.}
##' \item{\code{randomEffects}: }{function to return the generated random effects
##'   taking one argument, the dataset index.}
##' \item{\code{sphericalRandomeffects}: }{function to return the generated spherical random effects
##'   taking one argument, the dataset index.}
##' \item{\code{errors}: }{function to return the generated errors taking one argument,
##'   the dataset index.}
##' \item{\code{allRandomEffects}: }{function without arguments that returns the matrix of
##'   all generated random effects.}
##' \item{\code{allErrors}: }{function without arguments that returns the matrix of all
##'   generated errors.}
##' \item{\code{numberOfDatasets}: }{\code{numberOfDatasetsToGenerate} as supplied}
##' \item{\code{numberOfLevelsInFixedFactor}: }{\code{numberOfLevelsInFixedFactor} as supplied}
##' \item{\code{numberOfSubjects}: }{\code{numberOfSubjects} sorted.}
##' \item{\code{numberOfReplicates}: }{\code{numberOfReplicates} as supplied}
##' \item{\code{numberOfRows}: }{number of rows in the generated dataset}
##' \item{\code{trueBeta}: }{true values used for beta}
##' \item{\code{trueSigma}: }{true value used for sigma}
##' \item{\code{trueTheta}: }{true values used for theta}
##' \item{\code{formula}: }{formula to fit the model using \code{lmer}}
##' \item{\code{...}: }{additional arguments passed via \code{...}}
##' }
##' @seealso \code{\link{generateMixedEffectDatasets}} and
##'   \code{\link{createDatasetsFromList}}
##' @examples
##'   oneWay <- generateAnovaDatasets(2, 1, 5, 4)
##'   head(oneWay$generateData(1))
##'   head(oneWay$generateData(2))
##'   oneWay$formula
##'   head(oneWay$randomEffects(1))
##'   head(oneWay$sphericalRandomEffects(1))
##'   head(oneWay$errors(1))
##'
##'   twoWayFixedRandom <- generateAnovaDatasets(2, 3, 5, 4)
##'   head(twoWayFixedRandom$generateData(1))
##'   twoWayFixedRandom$formula
##'
##'   twoWayRandom <- generateAnovaDatasets(2, 1, c(3, 5), 4)
##'   head(twoWayRandom$generateData(1))
##'   twoWayRandom$formula
##'
##'   large <- generateAnovaDatasets(2, c(10, 15), c(20, 30), 5)
##'   head(large$generateData(1))
##'   large$formula
##'
##'   ## illustration how to generate more datasets
##'   set.seed(1)
##'   datasets1 <- generateAnovaDatasets(2, 1, 5, 4)
##'   set.seed(1)
##'   datasets2 <- generateAnovaDatasets(3, 1, 5, 4)
##'   stopifnot(all.equal(datasets1$generateData(1), datasets2$generateData(1)),
##'             all.equal(datasets1$generateData(2), datasets2$generateData(2)))
##' @export
##' @importFrom stats rnorm
generateAnovaDatasets <- function(numberOfDatasetsToGenerate,
                                  numberOfLevelsInFixedFactor,
                                  numberOfSubjects,
                                  numberOfReplicates,
                                  errorGenerator = rnorm,
                                  randomEffectGenerator = rnorm,
                                  trueBeta = 1,
                                  trueSigma = 4,
                                  trueTheta = 1,
                                  ...,
                                  arrange = FALSE) {
    checkIsScalar(numberOfDatasetsToGenerate,
                  "numberOfDatasetsToGenerate")
    checkLevelsFixedEffects(numberOfLevelsInFixedFactor)
    checkLevelsVarianceComponents(numberOfSubjects)
    checkReplicates(numberOfReplicates)
    checkArrange(arrange, numberOfLevelsInFixedFactor, numberOfSubjects)

    numberOfRows <-
        prod(c(numberOfLevelsInFixedFactor, numberOfSubjects)) * numberOfReplicates
    numberOfSubjects <-
        sortRandomEffectsLikeLmer(numberOfSubjects)
    numberOfFixedFactors <- length(numberOfLevelsInFixedFactor)
    numberOfFixedEffectCoefficients <-
        sum(numberOfLevelsInFixedFactor) - numberOfFixedFactors + 1
    numberOfRandomEffectCoefficients <- length(numberOfSubjects)

    trueBeta <-
        checkLengthOrReplicate(trueBeta, "trueBeta", numberOfFixedEffectCoefficients)
    checkIsScalar(trueSigma, "trueSigma")
    trueTheta <-
        checkLengthOrReplicate(trueTheta, "trueTheta", numberOfRandomEffectCoefficients)

    randomVariables <-
        generateAndAssembleRandomVariables(
            numberOfDatasetsToGenerate,
            numberOfSubjects,
            numberOfRows,
            trueSigma,
            trueTheta,
            randomEffectGenerator,
            errorGenerator
        )

    return(
        list(
            generateData = function(datasetIndex) {
                data <-
                    generateData(numberOfLevelsInFixedFactor,
                                 numberOfSubjects,
                                 numberOfRows,
                                 arrange)
                X <-
                    createXMatrix(data, numberOfLevelsInFixedFactor, numberOfRows)
                Z <-
                    createZMatrix(data, numberOfLevelsInFixedFactor, numberOfSubjects)
                data$y <-
                    as(X %*% trueBeta +
                           Z %*% randomVariables[["randomEffects"]][, datasetIndex] +
                           randomVariables[["errors"]][, datasetIndex],
                       'vector')
                attr(data, "datasetIndex") <-
                    datasetIndex + randomVariables[["datasetIndexOffset"]]
                return(data)
            },
            createXMatrix = function(data) {
                return(createXMatrix(data, numberOfLevelsInFixedFactor, numberOfRows))
            },
            createZMatrix = function(data) {
                return(createZMatrix(data, numberOfLevelsInFixedFactor, numberOfSubjects))
            },
            createLambdaMatrix = function(data) {
                return(createLambda(trueTheta, numberOfSubjects))
            },
            randomEffects = function(datsetIndex) {
                return(randomVariables[["randomEffects"]][, datsetIndex])
            },
            sphericalRandomEffects = function(datsetIndex) {
                Lambda <- createLambda(trueTheta, numberOfSubjects)
                sphericalRandomEffects <-
                    solve(Lambda, randomVariables[["randomEffects"]][, datsetIndex])
                return(sphericalRandomEffects)
            },
            errors = function(datsetIndex) {
                return(randomVariables[["errors"]][, datsetIndex])
            },
            allRandomEffects = function() {
                return(randomVariables[["randomEffects"]])
            },
            allErrors = function() {
                return(randomVariables[["errors"]])
            },
            numberOfDatasets = numberOfDatasetsToGenerate,
            numberOfLevelsInFixedFactor = numberOfLevelsInFixedFactor,
            numberOfSubjects = numberOfSubjects,
            numberOfReplicates = numberOfReplicates,
            numberOfRows = numberOfRows,
            trueBeta = trueBeta,
            trueSigma = trueSigma,
            trueTheta = trueTheta,
            formula = generateFormula(numberOfLevelsInFixedFactor, numberOfSubjects),
            ...
        )
    )
}

checkLevelsFixedEffects <- function(numberOfLevelsInFixedFactor) {
    if (length(numberOfLevelsInFixedFactor) > 1 &&
        any(numberOfLevelsInFixedFactor == 1)) {
        stop(
            "Argument 'numberOfLevelsInFixedFactor' contains factors with just one level. ",
            "This is only allowed if the argument is of length one."
        )
    }
    checkLargerThanZero(numberOfLevelsInFixedFactor,
                        "numberOfLevelsInFixedFactor")
}

checkLargerThanZero <- function(arg, argName) {
    if (any(arg < 1)) {
        stop("Argument '",
             argName,
             "' contains values less than one. This is not allowed.")
    }
}

checkLevelsVarianceComponents <- function(numberOfSubjects) {
    if (any(numberOfSubjects == 1)) {
        stop(
            "Argument 'numberOfSubjects' contains factors with just one level. ",
            "This is not allowed, please remove them."
        )
    }
    checkLargerThanZero(numberOfSubjects, "numberOfSubjects")
}

checkReplicates <- function(numberOfReplicates) {
    checkIsScalar(numberOfReplicates, "numberOfReplicates")
    if (numberOfReplicates < 1) {
        stop("Argument 'numberOfReplicates' needs to be larger than 0.")
    }
}

checkArrange <-
    function(arrange,
             numberOfLevelsInFixedFactor,
             numberOfSubjects) {
        if (arrange) {
            postfix <- " if arrange is set to TRUE"
            checkIsScalar(numberOfLevelsInFixedFactor,
                          "numberOfLevelsInFixedFactor",
                          postfix)
            checkIsScalar(numberOfSubjects, "numberOfSubjects", postfix)
        }
    }

checkIsScalar <- function(arg, argName, postfix = "") {
    if (length(arg) != 1) {
        stop(
            "Argument '",
            argName,
            "' is of length ",
            length(arg),
            " but it should be a scalar (i.e., a vector of length one)",
            postfix,
            "."
        )
    }
}

sortRandomEffectsLikeLmer <- function(numberOfSubjects) {
    return(sort(numberOfSubjects, decreasing = TRUE))
}

checkLengthOrReplicate <- function(arg, argName, expectedLength) {
    if (length(arg) == 1) {
        arg <- rep(arg, expectedLength)
    } else if (length(arg) != expectedLength) {
        stop(
            "Argument '",
            argName,
            "' is of length ",
            length(arg),
            " but it needs to be ",
            "either be of length 1 or ",
            expectedLength,
            "."
        )
    }
    return(arg)
}

generateAndAssembleRandomVariables <-
    function(numberOfDatasetsToGenerate,
             numberOfSubjects,
             numberOfRows,
             trueSigma,
             trueTheta,
             randomEffectGenerator,
             errorGenerator) {
        Lambda <- createLambda(trueTheta, numberOfSubjects)
        numberOfRandomEffects <- sum(numberOfSubjects)
        return(
            generateAndAssembleRandomVariablesUsingLambda(
                numberOfDatasetsToGenerate,
                numberOfRandomEffects,
                numberOfRows,
                trueSigma,
                Lambda,
                randomEffectGenerator,
                errorGenerator
            )
        )

    }

generateAndAssembleRandomVariablesUsingLambda <-
    function(numberOfDatasetsToGenerate,
             numberOfRandomEffects,
             numberOfRows,
             trueSigma,
             Lambda,
             randomEffectGenerator,
             errorGenerator) {
        generatedValues <-
            replicate(
                numberOfDatasetsToGenerate,
                generateOneSetOfRandomVariables(
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

generateOneSetOfRandomVariables <- function(numberOfRandomEffects,
                                            numberOfRows,
                                            trueSigma,
                                            Lambda,
                                            randomEffectGenerator,
                                            errorGenerator) {
    return(list(
        randomEffects = generateRandomEffects(
            numberOfRandomEffects,
            trueSigma,
            Lambda,
            randomEffectGenerator
        ),
        errors = generateErrors(trueSigma, numberOfRows, errorGenerator)
    ))
}

generateRandomEffects <-
    function(numberOfRandomEffects,
             trueSigma,
             Lambda,
             randomEffectGenerator) {
        sphericalRandomEffects <-
            randomEffectGenerator(numberOfRandomEffects) * trueSigma
        randomEffects <- drop(Lambda %*% sphericalRandomEffects)
        return(randomEffects)
    }

createLambda <- function(trueTheta, numberOfSubjects) {
    replicationUnits <-
        as.vector(unlist(mapply(
            rep, times = numberOfSubjects,
            x = seq_along(numberOfSubjects)
        )))
    Lambda <- diag(trueTheta[replicationUnits])
    return(Lambda)
}

generateErrors <-
    function(trueSigma, numberOfRows, errorGenerator) {
        errors <- errorGenerator(numberOfRows) * trueSigma
        return(errors)
    }

createXMatrix <-
    function(data,
             numberOfLevelsInFixedFactor,
             numberOfRows) {
        if (length(numberOfLevelsInFixedFactor) == 1 &&
            numberOfLevelsInFixedFactor == 1) {
            X <- matrix(1, numberOfRows, 1)
        } else {
            X <-
                MatrixModels::model.Matrix(~ ., data[, seq_along(numberOfLevelsInFixedFactor), drop = FALSE])
        }
        return(X)
    }

createZMatrix <-
    function(data,
             numberOfLevelsInFixedFactor,
             numberOfSubjects) {
        Ztmp <-
            apply(data[, length(numberOfLevelsInFixedFactor) + seq_along(numberOfSubjects), drop = FALSE],
                  2, function(componentData)
                      MatrixModels::model.Matrix(~ . - 1, data.frame(componentData)),
                  simplify = FALSE)
        Z <- Ztmp[[1]]
        if (length(numberOfSubjects) > 1) {
            for (i in 2:length(numberOfSubjects)) {
                Z <- cbind2(Z, Ztmp[[i]])
            }
        }
        return(Z)
    }

##' @importFrom utils as.roman
generateData <-
    function(numberOfLevelsInFixedFactor,
             numberOfSubjects,
             numberOfRows,
             arrange) {
        data <-
            expand.grid(c(
                lapply(numberOfLevelsInFixedFactor, function(x)
                    tolower(as.roman(1:x))),
                lapply(numberOfSubjects, function(x)
                    as.roman(1:x))
            ))
        data <-
            data.frame(apply(data, 2, rep, length.out = numberOfRows))
        data <- as.data.frame(lapply(data, factor))
        if (arrange) {
            data <- arrangeData(data)
        }
        return(data)
    }

globalVariables(c("Var2", "Var1"))

arrangeData <- function(data) {
    if (isPackageInstalled("dplyr")) {
        data <- dplyr::arrange(data, Var2, Var1)
    } else {
        warn("Cannot arrange data as package 'dplyr' is not installed.")
    }
    return(data)
}

generateFormula <-
    function(numberOfLevelsInFixedFactor,
             numberOfSubjects) {
        formulaString <- "y ~ "
        variableCounter <- 1
        for (numberOfLevels in numberOfLevelsInFixedFactor) {
            formulaString <-
                generateFormulaPart(formulaString,
                                    numberOfLevels,
                                    variableCounter,
                                    "Var",
                                    " ")
            variableCounter <- variableCounter + 1
        }
        for (numberOfLevels in numberOfSubjects) {
            formulaString <-
                generateFormulaPart(formulaString,
                                    numberOfLevels,
                                    variableCounter,
                                    "(1 | Var",
                                    ") ")
            variableCounter <- variableCounter + 1
        }
        return(as.formula(formulaString))
    }

generateFormulaPart <-
    function(formulaString,
             numberOfLevels,
             variableCounter,
             prefix,
             postfix) {
        if (numberOfLevels > 1) {
            if (substring(formulaString, nchar(formulaString) - 1) != "~ ") {
                formulaString <- paste0(formulaString, "+ ")
            }
            formulaString <-
                paste0(formulaString, prefix, variableCounter, postfix)
        }
        return(formulaString)
    }
