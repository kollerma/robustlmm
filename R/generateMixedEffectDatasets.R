##' This function runs \code{\link[lme4]{lmer}} and extracts all information needed to
##' generate new datasets using parametric bootstrap later.
##' @title Prepare Dataset for Parametric Bootstrap
##' @param formula passed on to \code{\link[lme4]{lmer}}
##' @param data passed on to \code{\link[lme4]{lmer}}
##' @param REML passed on to \code{\link[lme4]{lmer}}
##' @param overrideBeta use to override beta used to simulate new datasets, by
##'   default \code{\link{getME}}\code{(fm, "beta")} where \code{fm} is the
##'   fitted model returned by \code{\link[lme4]{lmer}}.
##' @param overrideSigma use to override sigma used to simulate new datasets, by
##'   default \code{\link{getME}}\code{(fm, "sigma")} where \code{fm} is the
##'   fitted model returned by \code{\link[lme4]{lmer}}.
##' @param overrideTheta use to override theta used to simulate new datasets, by
##'   default \code{\link{getME}}\code{(fm, "theta")} where \code{fm} is the
##'   fitted model returned by \code{\link[lme4]{lmer}}.
##' @param ... all additional arguments are added to the returned list.
##' @author Manuel Koller
##' @return List that can be passed to
##'   \code{\link{generateMixedEffectDatasets}}. \item{\code{data}:
##'   }{the original dataset} \item{\code{X}: }{the X matrix as returned by
##'   \code{\link{getME}}} \item{\code{Z}: }{the Z matrix as returned by
##'   \code{\link{getME}}} \item{\code{Lambda}: }{the Lambda matrix as returned
##'   by \code{\link{getME}}} \item{\code{numberOfFixedEffects}: }{the number of
##'   fixed effects coefficients} \item{\code{numberOfRandomEffects}: }{the
##'   number of random effects} \item{\code{numberOfRows}: }{number of rows in
##'   the generated dataset} \item{\code{trueBeta}: }{true values used for beta}
##'   \item{\code{trueSigma}: }{true value used for sigma}
##'   \item{\code{trueTheta}: }{true values used for theta}
##'   \item{\code{formula}: }{formula to fit the model using \code{lmer}}
##'   \item{\code{...}: }{additional arguments passed via \code{...}}
##' @examples
##'   preparedDataset <- prepareMixedEffectDataset(Reaction ~ Days + (Days|Subject), sleepstudy)
##'   str(preparedDataset)
##' @export
prepareMixedEffectDataset <-
    function(formula,
             data,
             REML = TRUE,
             overrideBeta,
             overrideSigma,
             overrideTheta,
             ...) {
        fm <- lmer(formula, data, REML)

        beta <- getME(fm, "beta")
        sigma <- getME(fm, "sigma")
        theta <- getME(fm, "theta")
        Lambda <- getME(fm, "Lambda")

        numberOfFixedEffects <- length(beta)

        if (!missing(overrideBeta)) {
            checkLength(overrideBeta, numberOfFixedEffects, "overrideBeta")
            beta <- overrideBeta
        }
        if (!missing(overrideSigma)) {
            checkIsScalar(overrideSigma, "overrideSigma")
            sigma <- overrideSigma
        }
        if (!missing(overrideTheta)) {
            checkLength(overrideTheta, length(theta), "overrideTheta")
            theta <- overrideTheta
            Lind <- getME(fm, "Lind")
            Lambda@x <- theta[Lind]
        }

        return(
            list(
                data = data,
                X = getME(fm, "X"),
                Z = getME(fm, "Z"),
                Lambda = Lambda,
                numberOfFixedEffects = numberOfFixedEffects,
                numberOfRandomEffects = length(getME(fm, "b")),
                numberOfRows = nobs(fm),
                trueBeta = beta,
                trueSigma = sigma,
                trueTheta = theta,
                formula = formula,
                ...
            )
        )
    }

checkLength <- function(arg, expectedLength, argName) {
    if (length(arg) != expectedLength) {
        stop(
            "Argument '",
            argName,
            "' is of length ",
            length(arg),
            " but it should be of length ",
            expectedLength,
            "."
        )
    }
}

##' Generates mixed effects datasets using parametric bootstrap.
##' @title Generate Mixed Effects Datasets
##' @param numberOfDatasetsToGenerate number of datasets to generate.
##' @param preparedDataset dataset as prepared by
##'   \code{\link{prepareMixedEffectDataset}}.
##' @param errorGenerator random number generator used for the errors.
##' @param randomEffectGenerator random number generator used for the spherical
##'   random effects.
##' @author Manuel Koller
##' @return list with generators and the contents of the prepared dataset. See
##'   \code{\link{prepareMixedEffectDataset}} and
##'   \code{\link{generateAnovaDatasets}} for a description of the contents.
##' @seealso \code{\link{generateAnovaDatasets}},
##'   \code{\link{prepareMixedEffectDataset}} and
##'   \code{\link{createDatasetsFromList}}
##' @examples
##'   preparedDataset <- prepareMixedEffectDataset(Reaction ~ Days + (Days|Subject), sleepstudy)
##'   datasets <- generateMixedEffectDatasets(2, preparedDataset)
##'   head(datasets$generateData(1))
##'   head(datasets$generateData(2))
##'   datasets$formula
##'   head(datasets$randomEffects(1))
##'   head(datasets$sphericalRandomEffects(1))
##'   head(datasets$errors(1))
##' @export
generateMixedEffectDatasets <- function(numberOfDatasetsToGenerate,
                                        preparedDataset,
                                        errorGenerator = rnorm,
                                        randomEffectGenerator = rnorm) {
    randomVariables <-
        generateAndAssembleRandomVariablesUsingLambda(
            numberOfDatasetsToGenerate,
            preparedDataset[["numberOfRandomEffects"]],
            preparedDataset[["numberOfRows"]],
            preparedDataset[["trueSigma"]],
            preparedDataset[["Lambda"]],
            randomEffectGenerator,
            errorGenerator
        )

    return(c(
        list(
            generateData = function(datasetIndex) {
                resp <-
                    as.character(preparedDataset[["formula"]])[2]
                data <- preparedDataset[["data"]]
                data[[resp]] <-
                    with(
                        preparedDataset,
                        as(
                            X %*% trueBeta +
                                Z %*% randomVariables[["randomEffects"]][, datasetIndex] +
                                randomVariables[["errors"]][, datasetIndex],
                            'vector'
                        )
                    )
                attr(data, "datasetIndex") <-
                    datasetIndex + randomVariables[["datasetIndexOffset"]]
                return(data)
            },
            createXMatrix = function(data) {
                return(preparedDataset[["X"]])
            },
            createZMatrix = function(data) {
                return(preparedDataset[["Z"]])
            },
            createLambdaMatrix = function(data) {
                return(preparedDataset[["Lambda"]])
            },
            randomEffects = function(datasetIndex) {
                return(randomVariables[["randomEffects"]][, datasetIndex])
            },
            sphericalRandomEffects = function(datasetIndex) {
                Lambda <- preparedDataset[["Lambda"]]
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
            numberOfDatasets = numberOfDatasetsToGenerate
        ),
        preparedDataset
    ))
}
