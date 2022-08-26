globalVariables("datasets")

fitOrLoadFits <-
    function(file,
             fittingFunctions,
             saveFitted,
             datasets) {
        fittedFile <- sub(".Rdata", "-fitted.Rdata", file)
        if (!file.exists(fittedFile)) {
            if (missing(datasets)) {
                contents <- load(file)
                checkContents(contents, file, "datasets")
            }
            fits <-
                unlist(lapply(fittingFunctions, function(fun)
                    fun(datasets)), recursive = FALSE)
            attr(fits, "sessionInfo") <- createSessionInfoString()
            if (saveFitted) {
                save(fits, file = fittedFile)
            }
        } else {
            contents <- load(fittedFile)
            checkContents(contents, fittedFile, "fits")
        }
        return(fits)
    }

checkContents <- function(contents, file, expected) {
    if (!(expected %in% contents)) {
        stop(
            "Expected file '",
            file,
            "' to contain object '",
            expected ,
            "' but got ",
            paste(contents, collapse = ", "),
            "."
        )
    }
}

##' @importFrom stats sessionInfo toLatex
createSessionInfoString <- function() {
    return(toLatex(sessionInfo(), locale = FALSE))
}

##' Combine list of processed fits into one list in matrix form.
##' @title Merge Processed Fits
##' @param processedFitList list of processed fits as produced by
##'   \code{\link{processFit}}.
##' @return similar list as returned by \code{\link{processFit}} just with
##'   matrix entries instead of vectors.
##' @examples
##'   preparedDataset <-
##'       prepareMixedEffectDataset(Reaction ~ Days + (Days|Subject),
##'                                 sleepstudy)
##'   set.seed(1)
##'   datasets <- generateMixedEffectDatasets(2, preparedDataset)
##'
##'   fits <- fitDatasets_lmer(datasets)
##'   processedFits <- lapply(fits, processFit, all = TRUE)
##'   merged <- mergeProcessedFits(processedFits)
##'   str(merged)
##' @export
mergeProcessedFits <- function(processedFitList) {
    processedFitList <- findValidEntries(processedFitList)
    allNames <- unique(unlist(lapply(processedFitList, names)))
    results <-
        lapply(allNames, extractValuesAndBind, processedFitList)
    names(results) <- allNames
    return(results)
}

findValidEntries <- function(processedFitList) {
    validEntries <- sapply(processedFitList, function(entry)
        ! is.null(entry[["label"]]) &&
            !is.null(entry[["datasetIndex"]]))
    return(processedFitList[validEntries])
}

extractValuesAndBind <- function(name, processedFitList) {
    expectedNumberOfColumns <-
        determineExpectedNumberOfColumns(processedFitList, name)
    valuesList <-
        lapply(processedFitList,
               getValuesOrNa,
               name,
               expectedNumberOfColumns)
    if (name == "label") {
        values <- factor(unlist(valuesList))
    } else {
        values <- do.call(rbind, valuesList)
    }
    return(values)
}

determineExpectedNumberOfColumns <-
    function(processedFitList, name) {
        maximumNumberOfColumns <-
            max(sapply(processedFitList, extractNumberOfColumns, name))
        return(maximumNumberOfColumns)
    }

extractNumberOfColumns <- function(processedFit, name) {
    value <- processedFit[[name]]
    if (is.null(value)) {
        return(0)
    } else if (is.matrix(value)) {
        return(NCOL(value))
    } else if (is.factor(value)) {
        return(1)
    } else {
        return(length(value))
    }
}

getValuesOrNa <-
    function(processedFit,
             name,
             expectedNumberOfColumns) {
        values <- processedFit[[name]]
        if (name == "label") {
            return(as.character(values))
        }
        if (is.null(values)) {
            values <- NA_real_
        }
        if (!is.matrix(values)) {
            colnames <- names(values)
            if (is.null(colnames)) {
                colnames <- name
                if (endsWith(colnames, "s") &&
                    name != "numberOfWarnings") {
                    colnames <- substr(colnames, 1, nchar(colnames) - 1)
                }
                if (expectedNumberOfColumns > 1) {
                    colnames <- paste(colnames, 1:expectedNumberOfColumns, sep = "")
                }
            }
            values <-
                matrix(
                    values,
                    ncol = expectedNumberOfColumns,
                    byrow = TRUE,
                    dimnames = list(NULL, colnames)
                )

        }
        if (NCOL(values) < expectedNumberOfColumns) {
            values <-
                cbind(values,
                      matrix(
                          NA_real_,
                          NROW(values),
                          expectedNumberOfColumns - NCOL(values)
                      ))
        }
        return(values)
    }


##' Call this function for each file stored using \code{\link{saveDatasets}}. If
##' a file hasn't been processed yet, then it is processed and a new file with
##' the postfix \dQuote{processed} is created containing the results.
##'
##' In case the raw fits may have to be inspected or \code{\link{processFit}}
##' may be called with another set of arguments, then set \code{saveFitted} to
##' TRUE. In that case, another file with the postfix \dQuote{fitted} is
##' created. Remove the files with postfix \dQuote{processed} and run
##' \code{processFile} again. The fits will not be re-done but instead loaded
##' from the file with postfix \dQuote{fitted}.
##'
##' @title Process File of Stored Datasets
##' @param file file saved by \code{\link{saveDatasets}}.
##' @param fittingFunctions vector of \code{\link{fitDatasets}} functions that
##'   should be applied to each dataset.
##' @param saveFitted logical, if true, the raw fits are also stored.
##' @param checkProcessed logical, if true, will check whether the contents of
##'   the processed output is reproduced for the first dataset. This is useful
##'   to ensure that everything is still working as expected without having to
##'   re-run the whole simulation study.
##' @param createMinimalSaveFile logical, if true, will create a file with the
##'   processed results of the first three datasets. This is helpful if one
##'   wants to store only the final aggregated results but still wants to make
##'   sure that the full code works as expected.
##' @param datasets optional, datasets as stored in \code{file}, to avoid doing
##'   a detour of saving and loading the file.
##' @param ... passed on to \code{\link{processFit}}. Use this to control what
##'   to save.
##' @return The list of all processed results merged together.
##'
##'   To help reproduciblility, the output of \code{toLatex(sessionInfo(),
##'   locale = FALSE)} is stored in the \code{sessionInfo} attribute.
##' @author Manuel Koller
##' @export
processFile <-
    function(file,
             fittingFunctions,
             saveFitted = FALSE,
             checkProcessed = FALSE,
             createMinimalSaveFile = FALSE,
             datasets,
             ...) {
        resultsFile <-
            processedFile <- sub(".Rdata", "-processed.Rdata", file)
        minimalFile <- sub(".Rdata", "-minimal.Rdata", file)
        processedResultsExist <- file.exists(processedFile)
        if (!processedResultsExist && file.exists(minimalFile)) {
            resultsFile <- minimalFile
            processedResultsExist <- TRUE
        }
        if (processedResultsExist) {
            if (checkProcessed) {
                results <-
                    checkProcessed(file,
                                   resultsFile,
                                   fittingFunctions,
                                   datasets,
                                   ...)
            } else {
                contents <- load(resultsFile)
                checkContents(contents, resultsFile, "results")
            }
            if (createMinimalSaveFile &&
                !file.exists(minimalFile)) {
                doCreateMinimalSaveFile(results, minimalFile)
            }
            return(results)
        }
        fits <-
            fitOrLoadFits(file, fittingFunctions, saveFitted, datasets)
        processedFitList <- lapply(fits, processFit, ...)
        results <- mergeProcessedFits(processedFitList)
        attr(results, "sessionInfo") <- attr(fits, "sessionInfo")
        save(results, file = processedFile)
        if (createMinimalSaveFile) {
            doCreateMinimalSaveFile(results, minimalFile)
        }
        return(results)
    }

doCreateMinimalSaveFile <- function(fullResults, minimalFile) {
    index <- fullResults$datasetIndex <= 3
    results <- list()
    for (name in names(fullResults)) {
        storedItem <- fullResults[[name]]
        if (length(dim(storedItem)) == 2) {
            extractedItem <- storedItem[index, , drop = FALSE]
        } else {
            extractedItem <- storedItem[index, drop = FALSE]
        }
        extractedItemList <- list(extractedItem)
        names(extractedItemList) <- name
        results <- c(results, extractedItemList)
    }
    attributes(results) <- attributes(fullResults)
    save(results, file = minimalFile)
}

globalVariables(c("datasets", "results"))

checkProcessed <-
    function(file,
             processedFile,
             fittingFunctions,
             datasets,
             ...) {
        if (missing(datasets)) {
            contents <- load(file)
            checkContents(contents, file, "datasets")
        }
        contents <- load(processedFile)
        checkContents(contents, processedFile, "results")
        datasets[["numberOfDatasets"]] <- 1
        firstFits <-
            unlist(lapply(fittingFunctions, function(fun)
                fun(datasets)), recursive = FALSE)
        for (fit in firstFits) {
            checkEqualsStoredResult(fit, results, ...)
        }
        return(results)
    }

checkEqualsStoredResult <- function(fit, storedResults, ...) {
    if (is(fit, "package-not-installed")) {
        return("Skipping check as required package is not installed.")
    }
    if (is(fit, "try-error")) {
        stopOnFailingFit(fit,
                         "Running the function produced the following error:\n",
                         fit)
    }
    expected <- extractFromStoredResults(storedResults, fit)
    if (!is.list(expected)) {
        stopOnFailingFit(fit,
                         "Could find corresponding results in stored results:\n",
                         expected)
    }
    actual <- processFit(fit, ...)
    actual <- actual[names(expected)]
    tolerance <- sqrt(.Machine$double.eps)
    if (is(fit, "lqmm")) {
        tolerance <- 0.01 ## summary.lqmm uses bootstrap.
        ## The result are sometimes different even though we specify the rng seed.
    }
    check <-
        all.equal(expected,
                  actual,
                  check.attributes = FALSE,
                  tolerance = tolerance)
    if (isTRUE(check)) {
        return(check)
    } else {
        stopOnFailingFit(fit,
                         "The processed results differed from the stored results:\n",
                         check)
    }
}

stopOnFailingFit <- function(fit, description, errorMessage = "") {
    label <- attr(fit, "label")
    if (is.null(label)) {
        label <- "<unknown function>"
    }
    warnings <- attr(fit, "warnings")
    if (is.null(warnings)) {
        warnings <- ""
    } else {
        warnings <- paste0("\nAnd some warnings were issued:\n", warnings)
    }
    stop("Error calling function '",
         label,
         "'. ",
         description,
         errorMessage,
         warnings)
}

extractFromStoredResults <- function(storedResults, fit) {
    label <- attr(fit, "label")
    if (is.null(label)) {
        return("Attribute 'label' is null for the fit, but a non-null value is required.")
    }
    datasetIndex <- attr(fit, "datasetIndex")
    if (is.null(datasetIndex)) {
        return("Attribute 'datasetIndex' is null for the fit, but a non-null value is required.")
    }
    index <-
        storedResults[["label"]] == label &
        storedResults[["datasetIndex"]] == datasetIndex
    if (sum(index) == 0) {
        return(
            paste0(
                "Could not find processed results matching label '",
                label,
                "' and datasetIndex '",
                datasetIndex,
                "'."
            )
        )
    } else if (sum(index) > 1) {
        return(
            paste0(
                "Expected just one match for '",
                label,
                "' and datasetIndex '",
                datasetIndex,
                "' but got ",
                sum(index),
                "."
            )
        )
    }
    extracted <- list(label = label, datasetIndex = datasetIndex)
    for (name in names(storedResults)) {
        if (!name %in% c("label", "datasetIndex")) {
            extractedItem <- list(storedResults[[name]][index, drop = FALSE])
            names(extractedItem) <- name
            extracted <- c(extracted, extractedItem)
        }
    }
    return(extracted)
}

createSessionInfoString <- function() {
    return(toLatex(sessionInfo(), locale = FALSE))
}

##' Convenience function to run simulation study in parallel on a single
##' machine.
##'
##' The merged results are saved in a file taking the name
##' \code{<path>/<baseFilename>-processed.Rdata}. You can delete the
##' intermediate result files with the numbers (the chunk index) in the name.
##'
##' To run on multiple machines, use \code{\link{saveDatasets}} to save datasets
##' into multiple files. Then call \code{\link{processFile}} on each of them on
##' the designated machine. Finally, load and merge the results together using
##' \code{\link{loadAndMergePartialResults}}.
##'
##' @title Process Datasets in Parallel
##' @param datasets dataset list generated by one of the generate functions.
##' @param path path to save the datasets to.
##' @param baseFilename filename to use, without extension.
##' @param fittingFunctions vector of \code{\link{fitDatasets}} functions that
##'   should be applied to each dataset.
##' @param chunkSize number of datasets to process together in a single job.
##' @param saveFitted logical, if true, the raw fits are also stored.
##' @param checkProcessed logical, if true, will check whether the contents of
##'   the processed output is reproduced for the first dataset. This is useful
##'   to ensure that everything is still working as expected without having to
##'   re-run the whole simulation study.
##' @param createMinimalSaveFile logical, if true, will create a file with the
##'   processed results of the first three datasets. This is helpful if one
##'   wants to store only the final aggregated results but still wants to make
##'   sure that the full code works as expected.
##' @param ncores number of cores to use in processing, if set to 1, datasets
##'   are processed in the current R session. Use
##'   \code{\link[parallel]{detectCores}} to find out how many cores are
##'   available on your machine.
##' @param clusterType type of cluster to be created, passed to
##'   \code{\link[parallel]{makeCluster}}.
##' @param ... passed on to \code{\link{processFit}}. Use this to control what
##'   to save.
##' @return The list of all processed results merged together.
##'
##'   To help reproduciblility, the output of \code{toLatex(sessionInfo(),
##'   locale = FALSE)} is stored in the \code{sessionInfo} attribute.
##' @seealso \code{\link{saveDatasets}}, \code{\link{processFile}}
##' @author Manuel Koller
##' @export
processDatasetsInParallel <- function(datasets,
                                      path,
                                      baseFilename,
                                      fittingFunctions,
                                      chunkSize,
                                      saveFitted = FALSE,
                                      checkProcessed = FALSE,
                                      createMinimalSaveFile = FALSE,
                                      ncores = 1,
                                      clusterType = "PSOCK",
                                      ...) {
    file <- paste0(file.path(path, baseFilename), ".Rdata")
    processedFile <- sub(".Rdata", "-processed.Rdata", file)
    minimalFile <- sub(".Rdata", "-minimal.Rdata", file)

    if (file.exists(processedFile) || file.exists(minimalFile)) {
        results <-
            processFile(
                file,
                fittingFunctions,
                checkProcessed = checkProcessed,
                createMinimalSaveFile = createMinimalSaveFile,
                datasets = datasets,
                ...
            )
        return(results)
    }
    files <-
        saveDatasets(datasets,
                     path,
                     file = paste0(baseFilename, "-%05d"),
                     chunkSize = chunkSize)
    if (ncores > 1) {
        ncores <- min(ncores, length(files))
    }
    if (ncores > 1) {
        cl <- parallel::makeCluster(ncores, type = clusterType)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterEvalQ(cl, require(robustlmm))
    }
    done <- listIntermediateResults(path, baseFilename)
    done <- sub("-processed", "", done)
    files <- setdiff(files, done)
    if (ncores > 1) {
        parallel::clusterApplyLB(
            cl,
            files,
            processFile,
            fittingFunctions,
            saveFitted = saveFitted,
            checkProcessed = checkProcessed,
            createMinimalSaveFile = FALSE,
            ...
        )
    } else {
        sapply(files, function(x) {
            processFile(
                x,
                fittingFunctions,
                saveFitted = saveFitted,
                checkProcessed = checkProcessed,
                createMinimalSaveFile = FALSE,
                ...
            )
        })
    }
    done <- listIntermediateResults(path, baseFilename)
    results <- loadAndMergePartialResults(done)
    save(results, file = processedFile)
    if (createMinimalSaveFile) {
        doCreateMinimalSaveFile(results, minimalFile)
    }
    return(results)
}

listIntermediateResults <- function(path, baseFilename) {
    return(list.files(
        path,
        paste0(baseFilename, "-\\d{1,5}-processed.Rdata"),
        full.names = TRUE
    ))
}

##' Convenience function that loads the results stored in each of the files and
##' then calls \code{\link{mergeProcessedFits}} to merge them.
##' @title Load And Merge Partial Results
##' @param files vector of filenames (including paths) of files containing the
##'   processed results
##' @seealso \code{\link{processDatasetsInParallel}}
##' @author Manuel Koller
##' @export
loadAndMergePartialResults <- function(files) {
    partialResults <- lapply(files, function(file) {
        contents <- load(file)
        checkContents(contents, file, "results")
        return(results)
    })
    results <- mergeProcessedFits(partialResults)
    attr(results, "sessionInfo") <-
        attr(partialResults[[1]], "sessionInfo")
    return(results)
}
