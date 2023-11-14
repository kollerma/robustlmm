##' This method can be used to bind multiple datasets generated using different
##' random genrators into one large dataset. The underlying dataset needs to be
##' the same.
##' @title Bind Generated Datasets
##' @param ... multiple datasets to be bound together
##' @param datasetList list of datasets created with one of the generate dataset
##'   functions
##' @author Manuel Koller
##' @return merged list with generators and the contents of the prepared
##'   dataset. See `\code{\link{prepareMixedEffectDataset}} and
##'   \code{\link{generateAnovaDatasets}} for a description of the contents.
##' @seealso \code{\link{splitDatasets}}
##' @examples
##'   datasets1 <- generateAnovaDatasets(2, 4, 4, 4)
##'   datasets2 <- generateAnovaDatasets(2, 4, 4, 4)
##'   datasets <- bindDatasets(datasets1, datasets2)
##'   data <- datasets$generateData(1)
##'   stopifnot(data$numberOfDatasets == 4,
##'             all.equal(datasets2$generateData(1), datasets$generateData(3),
##'                       check.attributes = FALSE),
##'             all.equal(datasets2$sphericalRandomEffects(1), datasets$sphericalRandomEffects(3)),
##'             all.equal(datasets2$createXMatrix(data), datasets$createXMatrix(data)),
##'             all.equal(datasets2$createZMatrix(data), datasets$createZMatrix(data)))
##'
##'   preparedDataset <- prepareMixedEffectDataset(Reaction ~ Days + (Days|Subject), sleepstudy)
##'   datasets1 <- generateMixedEffectDatasets(2, preparedDataset)
##'   datasets2 <- generateMixedEffectDatasets(2, preparedDataset)
##'   datasets <- bindDatasets(datasets1, datasets2)
##'   data <- datasets$generateData(1)
##'   stopifnot(data$numberOfDatasets == 4,
##'             all.equal(datasets2$generateData(1), datasets$generateData(3),
##'                       check.attributes = FALSE),
##'             all.equal(datasets2$sphericalRandomEffects(1), datasets$sphericalRandomEffects(3)),
##'             all.equal(datasets2$createXMatrix(data), datasets$createXMatrix(data)),
##'             all.equal(datasets2$createZMatrix(data), datasets$createZMatrix(data)))
##' @export
bindDatasets <- function(..., datasetList = list(...)) {
    checkAdditionalInformation(datasetList)
    randomVariables <- bindRandomVariables(datasetList)
    result <- datasetList[[1]]
    environmentClone <- cloneEnvironment(result$generateData)
    environmentClone[["randomVariables"]] <- randomVariables
    result <- lapply(result, setEnvironment, environmentClone)
    result$numberOfDatasets <-
        sum(sapply(datasetList, `[[`, "numberOfDatasets"))
    return(result)
}

checkAdditionalInformation <- function(datasetList) {
    additionalInformation <- datasetList[[1]][-(1:9)]
    for (dataIndex in seq_along(datasetList)[-1]) {
        otherAdditionalInformation <- datasetList[[dataIndex]][-(1:9)]
        result <-
            all.equal(additionalInformation, additionalInformation)
        if (!isTRUE(result)) {
            stop(
                "Can only bind datasets prepared using the same arguments. Prepared dataset ",
                dataIndex,
                " appears to be different: ",
                result
            )
        }
    }
    return(invisible())
}

bindRandomVariables <- function(datasetList) {
    return(list(
        randomEffects = do.call(
            cbind,
            lapply(datasetList, function(dataset)
                dataset$allRandomEffects())
        ),
        errors = do.call(cbind, lapply(datasetList, function(dataset)
            dataset$allErrors())),
        datasetIndexOffset = 0
    ))
}

cloneEnvironment <- function(f) {
    env <- environment(f)
    clone <- as.environment(as.list(env, all.names = TRUE))
    parent.env(clone) <- parent.env(env)
    return(clone)
}

setEnvironment <- function(f, environmentClone) {
    if (is.function(f)) {
        environment(f) <- environmentClone
    }
    return(f)
}

##' Method that splits up dataset objects into smaller chunks, so that they can
##' be processed separately.
##' @title Split Datasets Into Chunks
##' @param datasets dataset object to split into chunks
##' @param chunkSize number of datasets to keep in one chunk
##' @author Manuel Koller
##' @return list of dataset lists with generators and the contents of the
##'   original dataset. See \code{\link{prepareMixedEffectDataset}} and
##'   \code{\link{generateAnovaDatasets}} for a description of the contents.
##'   There is one additional entry in the list:
##'   \item{\code{chunkIndex}: }{index of the chunk}
##' @seealso \code{\link{bindDatasets}}
##' @examples
##'   oneWay <- generateAnovaDatasets(18, 1, 5, 4)
##'   datasetList <- splitDatasets(oneWay, 5)
##'   data <- datasetList[[4]]$generateData(1)
##'   stopifnot(all.equal(oneWay$generateData(16), datasetList[[4]]$generateData(1),
##'                       check.attributes = TRUE),
##'             all.equal(oneWay$sphericalRandomEffects(16),
##'                       datasetList[[4]]$sphericalRandomEffects(1)),
##'             all.equal(oneWay$createXMatrix(data), datasetList[[4]]$createXMatrix(data)),
##'             all.equal(oneWay$createZMatrix(data), datasetList[[4]]$createZMatrix(data)))
##' @export
splitDatasets <- function(datasets, chunkSize = 50) {
    numberOfChunks <-
        ceiling(datasets[["numberOfDatasets"]] / chunkSize)
    return(lapply(1:numberOfChunks, function(chunkIndex)
        createDatasetChunk(chunkIndex, datasets, chunkSize)))
}

createDatasetChunk <- function(chunkIndex, datasets, chunkSize) {
    firstDatasetIndex <- (chunkIndex - 1) * chunkSize + 1
    lastDatasetIndex <-
        min(chunkIndex * chunkSize, datasets[["numberOfDatasets"]])
    datasetIndices <- firstDatasetIndex:lastDatasetIndex
    randomVariables <-
        splitRandomVariables(datasets, datasetIndices, firstDatasetIndex - 1)
    result <- datasets
    environmentClone <- cloneEnvironment(result$generateData)
    environmentClone[["randomVariables"]] <- randomVariables
    result <- lapply(result, setEnvironment, environmentClone)
    result$numberOfDatasets <- length(datasetIndices)
    result$chunkIndex = chunkIndex
    return(result)
}

splitRandomVariables <-
    function(datasets,
             datasetIndices,
             datasetIndexOffset) {
        return(
            list(
                randomEffects = datasets$allRandomEffects()[, datasetIndices, drop = FALSE],
                errors = datasets$allErrors()[, datasetIndices, drop = FALSE],
                datasetIndexOffset = datasetIndexOffset
            )
        )
    }

##' Convert a list of datasets to a dataset list similar to the ones created by
##' \code{\link{generateAnovaDatasets}} and
##' \code{\link{generateMixedEffectDatasets}}.
##'
##' The returned list can be passed to \code{\link{processFit}} and to any of
##' the \code{\link{fitDatasets}} functions. Splitting and binding of datasets
##' using \code{\link{splitDatasets}} and \code{\link{bindDatasets}} is not
##' supported.
##' @title Create Dataset List From List of Data Objects
##' @param datasetList list of data objects, usually of type \code{data.frame}.
##' @param formula formula to fit the model using \code{lmer}.
##' @param trueBeta scalar or vector with the true values of the fixed effects
##'   coefficients. Can be of length one in which case it will be replicated to
##'   the required length if needed.
##' @param trueSigma scalar with the true value of the error scale.
##' @param trueTheta scalar or vector with the true values for the variance
##'   component coefficients, not including sigma. Can be of length one in which
##'   case it will be replicated to the required length if needed.
##' @param ... all additional arguments are added to the returned list.
##' @return list that can be passed to \code{\link{processFit}} and to any of
##'   the \code{\link{fitDatasets}} functions. Only \code{generateData} is
##'   implemented, all the other functions return an error if called.
##' @seealso \code{\link{generateAnovaDatasets}} and
##'   \code{\link{generateMixedEffectDatasets}}
##' @examples
##'   data(sleepstudy)
##'   sleepstudy2 <- sleepstudy
##'   sleepstudy2[1, "Reaction"] <- sleepstudy2[1, "Reaction"] + 10
##'   fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##'   datasets <- createDatasetsFromList(list(sleepstudy, sleepstudy2),
##'                                      formula = Reaction ~ Days + (Days|Subject),
##'                                      trueBeta = getME(fm1, "beta"),
##'                                      trueSigma = sigma(fm1),
##'                                      trueTheta = getME(fm1, "theta"))
##'   fitDatasets_lmer(datasets)
##' @export
createDatasetsFromList <-
    function(datasetList,
             formula,
             trueBeta,
             trueSigma,
             trueTheta,
             ...) {
        notImplemented <- function(...)
            stop("Function not implemented")
        return(
            list(
                generateData = function(i) {
                    result <- datasetList[[i]]
                    attr(result, "datasetIndex") <- i
                    return(result)
                },
                createXMatrix = function(data) {
                    fm <- lmerNoFit(formula, data)
                    return(getME(fm, "X"))
                },
                createZMatrix = function(data) {
                    fm <- lmerNoFit(formula, data)
                    return(getME(fm, "Z"))
                },
                createLambdaMatrix = function(data) {
                    fm <- lmerNoFit(formula, data)
                    return(getME(fm, "Lambda"))
                },
                randomEffects = notImplemented,
                sphericalRandomEffects = notImplemented,
                errors = notImplemented,
                allRandomEffects = notImplemented,
                allErrors = notImplemented,
                numberOfDatasets = length(datasetList),
                numberOfRows = NROW(datasetList[[1]]),
                trueBeta = trueBeta,
                trueSigma = trueSigma,
                trueTheta = trueTheta,
                formula = formula,
                ...
            )
        )
    }

##' Saves dataset to one or more files.
##'
##' The file will be saved to \code{path/filename.Rdata}.
##'
##' If \code{chunkSize} is not missing, the filename is interpreted as format
##' specifier and passed onto \code{\link{sprintf}}. One argument is given, the
##' index of the chunk.
##' @title Save datasets
##' @param datasets dataset list generated by one of the generate functions.
##' @param path path to save the datasets to.
##' @param file filename to use, without extension.
##' @param chunkSize if provided, datasets are split into \code{chunkSize}
##'   chunks and then saved.
##' @author Manuel Koller
##' @return filename or vector of filenames.
##' @export
saveDatasets <-
    function(datasets, path = getwd(), file, chunkSize) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (missing(chunkSize)) {
            return(saveChunk(datasets, path, file))
        }

        datasetList <- splitDatasets(datasets, chunkSize)
        ## TODO Save base dataset as full, then save only
        ## randomVariables for subsequent datasets.
        files <-
            sapply(datasetList, saveChunk, path = path, file = file)
        return(files)
    }

saveChunk <- function(datasets, path, file) {
    if (exists("chunkIndex", datasets)) {
        file <- sprintf(file, datasets[["chunkIndex"]])
    }
    file <- paste0(file.path(path, file), ".Rdata")
    datasets <- cleanSrcref(datasets)
    save(datasets, file = file)
    return(file)
}

cleanSrcref <- function(datasets) {
    return(lapply(datasets, function(x) {
        attr(x, "srcref") <- NULL
        return(x)
    }))
}
