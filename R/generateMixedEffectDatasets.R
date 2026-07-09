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

##' Generate balanced repeated-measures datasets with a structured
##' random-effect covariance.
##'
##' Each subject is observed once per level of a within-subject factor
##' \code{visit} (with \code{numberOfVisits} levels), optionally replicated
##' \code{numberOfReplicates} times, giving a \code{numberOfVisits}-dimensional
##' random effect per subject through the term \code{(0 + visit | subject)}.
##' The random-effect covariance follows the requested \code{structure}:
##' \describe{
##'   \item{\code{"unstructured"}}{an arbitrary covariance (the
##'     \code{(0 + visit | subject)} term).}
##'   \item{\code{"cs"}}{compound symmetry: a common \code{correlation}
##'     between all visits (\code{cs(0 + visit | subject)}).}
##'   \item{\code{"ar1"}}{autoregressive: \eqn{\mathrm{Cor}(i, j) =
##'     \code{correlation}^{|i - j|}} (\code{ar1(0 + visit | subject)}).}
##'   \item{\code{"diag"}}{uncorrelated visits (\code{diag(0 + visit |
##'     subject)}).}
##' }
##' The data are simulated from the chosen covariance using the same
##' machinery as \code{\link{generateMixedEffectDatasets}}; the returned
##' object has the identical interface (including \code{generateData},
##' \code{formula} and the \code{fitDatasets_*} compatibility) and stores the
##' structured covariance in \code{trueTheta}. Structured covariances require
##' \pkg{lme4} >= 2.0-0.
##'
##' @title Generate Repeated-Measures Datasets With Structured Covariance
##' @param numberOfDatasetsToGenerate number of datasets to generate.
##' @param numberOfSubjects number of subjects (grouping levels).
##' @param numberOfVisits number of levels of the within-subject factor
##'   (the dimension of the random effect).
##' @param numberOfReplicates number of replicates per subject-by-visit cell.
##' @param structure random-effect covariance structure, one of
##'   \code{"unstructured"}, \code{"cs"}, \code{"ar1"}, \code{"diag"}.
##' @param marginalSd marginal standard deviation(s) of the random effects,
##'   recycled to length \code{numberOfVisits}.
##' @param correlation the structure's single correlation parameter (the
##'   common correlation for \code{"cs"}, the lag-1 correlation for
##'   \code{"ar1"}); ignored for \code{"diag"}.
##' @param trueBeta the true intercept (the only fixed effect).
##' @param trueSigma the true residual standard deviation.
##' @param errorGenerator,randomEffectGenerator functions used to draw the
##'   errors and spherical random effects, see
##'   \code{\link{generateMixedEffectDatasets}}.
##' @return A list with the same structure as the return value of
##'   \code{\link{generateMixedEffectDatasets}}.
##' @seealso \code{\link{generateMixedEffectDatasets}},
##'   \code{\link{generateLongitudinalDatasets}}
##' @examples
##'   if (packageVersion("lme4") >= "2.0.0") {
##'     datasets <- generateRepeatedMeasuresDatasets(
##'         1, numberOfSubjects = 30, numberOfVisits = 3, numberOfReplicates = 4,
##'         structure = "cs", marginalSd = c(2, 1.5, 1.2), correlation = 0.5)
##'     fit <- rlmer(datasets$formula, datasets$generateData(1), method = "DASvar")
##'     VarCorr(fit)
##'   }
##' @export
generateRepeatedMeasuresDatasets <-
    function(numberOfDatasetsToGenerate,
             numberOfSubjects,
             numberOfVisits,
             numberOfReplicates = 1L,
             structure = c("unstructured", "cs", "ar1", "diag"),
             marginalSd = 1,
             correlation = 0,
             trueBeta = 0,
             trueSigma = 1,
             errorGenerator = rnorm,
             randomEffectGenerator = rnorm) {
        structure <- match.arg(structure)
        checkIsScalar(numberOfDatasetsToGenerate, "numberOfDatasetsToGenerate")
        checkIsScalar(numberOfSubjects, "numberOfSubjects")
        checkIsScalar(numberOfVisits, "numberOfVisits")
        checkIsScalar(numberOfReplicates, "numberOfReplicates")
        checkIsScalar(correlation, "correlation")
        checkIsScalar(trueBeta, "trueBeta")
        checkIsScalar(trueSigma, "trueSigma")
        if (structure != "unstructured" && packageVersion("lme4") < "2.0.0") {
            stop("structured covariances ('", structure,
                 "') require lme4 >= 2.0-0.")
        }
        q <- as.integer(numberOfVisits)
        sds <- rep_len(marginalSd, q)
        R <- switch(structure,
                    cs = {
                        m <- matrix(correlation, q, q)
                        diag(m) <- 1
                        m
                    },
                    ar1 = correlation ^ abs(outer(seq_len(q), seq_len(q), "-")),
                    diag(q))   # diag and unstructured both start from identity
        Sigma <- outer(sds, sds) * R
        reTerm <- if (structure == "unstructured") {
            "(0 + visit | subject)"
        } else {
            sprintf("%s(0 + visit | subject)", structure)
        }
        formula <- stats::as.formula(paste("y ~ 1 +", reTerm))
        design <- expand.grid(
            visit = gl(q, 1L, labels = paste0("v", seq_len(q))),
            replicate = seq_len(numberOfReplicates),
            subject = gl(numberOfSubjects, 1L))
        design[["y"]] <- 0
        preparedDataset <-
            prepareMixedEffectDataset(formula, design, REML = TRUE,
                                      overrideBeta = trueBeta,
                                      overrideSigma = trueSigma)
        ## set the relative covariance factor (Sigma / sigma^2 = Lambda
        ## Lambda^T) to the requested structure exactly by overwriting the
        ## nonzero entries of the (block-diagonal, per-subject) Lambda.
        relativeChol <- t(chol(Sigma / trueSigma ^ 2))
        nonZerosPerBlock <-
            length(preparedDataset[["Lambda"]]@x) / numberOfSubjects
        block <- if (nonZerosPerBlock == q) {
            diag(relativeChol)
        } else {
            relativeChol[lower.tri(relativeChol, diag = TRUE)]
        }
        preparedDataset[["Lambda"]]@x <- rep(block, numberOfSubjects)
        preparedDataset[["trueTheta"]] <- block
        return(
            generateMixedEffectDatasets(
                numberOfDatasetsToGenerate,
                preparedDataset,
                errorGenerator = errorGenerator,
                randomEffectGenerator = randomEffectGenerator
            )
        )
    }
