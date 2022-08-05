##' This method creates a list of datasets that can be used to create
##' sensitivity curves. The response of the dataset is modified according to the
##' supplied arguments.
##'
##' Either \code{shifts} or \code{scales} need to be provided. Both are also
##' possible.
##'
##' The argument \code{shifts} contains all the values that shall be added to
##' each of the observations that should be changed. One value per generated
##' dataset.
##'
##' The argument \code{scales} contains all the values that shall be used to
##' move observations away from their center. If \code{scales} is provided, then
##' \code{observationsToChange} needs to select more than one observation.
##'
##' The returned list can be passed to \code{\link{processFit}} and to any of
##' the \code{\link{fitDatasets}} functions. Splitting and binding of datasets
##' using \code{\link{splitDatasets}} and \code{\link{bindDatasets}} is not
##' supported.
##' @title Generate Datasets To Create Sensitivity Curves
##' @param data dataset to be modified.
##' @param observationsToChange index or logical vector indicating which
##'   observations should be modified.
##' @param shifts vector of shifts that should be applied one by one to each of
##'   the modified observations.
##' @param scales vector scales that should be used to scale the observations
##'   around their original center.
##' @param center optional scalar used to define the center from which the
##'   observations are scaled from. If missing, the mean of all the changed
##'   observations is used.
##' @param formula formula to fit the model using \code{lmer}.
##' @param ... all additional arguments are added to the returned list.
##' @return list that can be passed to \code{\link{processFit}} and to any of
##'   the \code{\link{fitDatasets}} functions. Only \code{generateData} is
##'   implemented, all the other functions return an error if called.
##' @seealso \code{\link{generateAnovaDatasets}}
##' @examples
##'   oneWay <- generateAnovaDatasets(1, 1, 10, 5)
##'   datasets <-
##'       generateSensitivityCurveDatasets(oneWay$generateData(1),
##'                                        observationsToChange = 1:5,
##'                                        shifts = -10:10,
##'                                        formula = oneWay$formula)
##'   datasets$generateData(1)
##' @export
generateSensitivityCurveDatasets <-
    function(data,
             observationsToChange,
             shifts,
             scales,
             center,
             formula,
             ...) {
        if (length(observationsToChange) < 1) {
            stop("Argument \"observationsToChange\" needs to be of length >= 1.")
        }
        if (is.logical(observationsToChange)) {
            if (length(observationsToChange) != NROW(data)) {
                stop(
                    "If argument \"observationsToChange\" is logical, then it needs to",
                    " be of the same length as the number of rows in data."
                )
            }
            observationsToChange <- which(observationsToChange)
        }
        if (missing(shifts) && missing(scales)) {
            stop("Both arguments \"shifts\" and \"scales\" missing, need at least one.")
        }
        if (!missing(scales) && length(observationsToChange) == 1) {
            stop(
                "Argument \"scales\" only makes sense if more than one observation",
                " is changed."
            )
        }
        if (!missing(shifts) &&
            !missing(scales) && length(shifts) != length(scales)) {
            stop(
                "If both arguments \"shifts\" and \"scales\" are given, then they",
                " need to be of the same length."
            )
        }
        if (missing(shifts)) {
            shifts <- rep(0, length(scales))
        }
        fm <- lmer(formula, data)
        response <- which(colnames(data) == as.character(formula)[2])
        scaleFromCenter <- length(observationsToChange) > 1 && !missing(scales)
        if (scaleFromCenter && missing(center)) {
            center <- mean(data[observationsToChange, response])
        }
        trueBeta <- fixef(fm)
        trueSigma <- sigma(fm)
        trueTheta <- getME(fm, "theta")
        rm(fm)
        notImplemented <- function(...)
            stop("Function not implemented")
        return(
            list(
                generateData = function(datasetIndex) {
                    ldata <- data
                    if (scaleFromCenter) {
                        ldata[observationsToChange, response] <-
                            center +
                            (ldata[observationsToChange, response] - center) * scales[datasetIndex] +
                            shifts[datasetIndex]
                    } else {
                        ldata[observationsToChange, response] <-
                            ldata[observationsToChange, response] + shifts[datasetIndex]
                    }
                    attr(ldata, "datasetIndex") <- datasetIndex
                    return(ldata)
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
                numberOfDatasets = length(shifts),
                numberOfRows = NROW(data),
                trueBeta = trueBeta,
                trueSigma = trueSigma,
                trueTheta = trueTheta,
                formula = formula,
                ...
            )
        )
    }
