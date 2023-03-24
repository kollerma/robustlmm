##' Apply function for all generated datasets.
##'
##' @title Lapply for generated datasets
##' @param datasets Datasets list to be used to generate datasets.
##' @param FUN the function to be applied to each generated dataset. The
##'   function will be called like \code{FUN(data, ...)}.
##' @param ... optional arguments to \code{FUN}.
##' @param label optional parameter, if present, each result is added an
##'   attribute named \emph{label} with the value of \code{label}.
##' @param POST_FUN function to be applied to the result of \code{FUN}. While
##'   one could just modify \code{FUN} instead, this additional argument makes
##'   it a bit easier to combine different kinds of methods together.
##' @param datasetIndices optional vector of dataset indices to fit, useful to
##'   try only a few datasets instead of all of them. Use \code{"all"} to
##'   process all datasets (default).
##' @author Manuel Koller
##' @return list of results. The items in the resulting list will have two
##'   additional attributes: \code{datasetIndex} and \code{proc.time}. If
##'   \code{FUN} failed for an item, then the item  will be the error as
##'   returned by try, i.e., it ill be of class \code{try-error}.
##' @examples
##'   oneWay <- generateAnovaDatasets(2, 1, 5, 4)
##'   lapplyDatasets(oneWay, function(data) sum(data$y))
##'   lapplyDatasets(oneWay, function(data) sum(data$y), POST_FUN = function(x) x^2)
##' @export
lapplyDatasets <- function(datasets,
                           FUN,
                           ...,
                           label,
                           POST_FUN,
                           datasetIndices = "all") {
    haveLabel <- !missing(label)
    if (missing(POST_FUN)) {
        POST_FUN <- identity
    }
    datasetIndices <- checkDatasetIndices(datasetIndices, datasets)
    return(lapply(datasetIndices, function(datasetIndex) {
        data <- datasets$generateData(datasetIndex)
        warnings <- NULL
        time <-
            system.time(result <- withCallingHandlers(
                try(POST_FUN(FUN(data, ...)))
                ,
                warning = function(w) {
                    warnings <<- append(warnings, list(w))
                    invokeRestart("muffleWarning")
                }
            ))
        attr(result, "datasetIndex") <- attr(data, "datasetIndex")
        attr(result, "proc.time") <- time
        if (!is.null(warnings)) {
            attr(result, "warnings") <- warnings
        }
        if (haveLabel) {
            attr(result, "label") <- label
        }
        return(result)
    }))
}

checkDatasetIndices <- function(datasetIndices, datasets) {
    if (length(datasetIndices) == 1L) {
        if (datasetIndices == "all") {
            datasetIndices <- 1:datasets[["numberOfDatasets"]]
        } else {
            datasetIndices <- as.integer(datasetIndices)
        }
    } else if (!is.integer(datasetIndices)) {
        stop("Expected a vector of integers for 'datasetIndices' but got ",
             datasetIndices)
    }
    return(datasetIndices)
}

##' Methods to fit various mixed effects estimators to all generated datasets.
##'
##' @details Existing fitting functions are:
##' @title Fitting Functions
##' @param datasets Datasets list to be used to generate datasets.
##' @param control a list (of correct class for the respective fitting function)
##'   containing control parameters to be passed through.
##' @param label a string used to identify which fits have been created by which
##'   function.
##' @param postFit a function, taking one argument, the resulting fit. This
##'   makes it easy to add an additional step after fitting.
##' @param datasetIndices optional vector of dataset indices to fit, useful to
##'   try only a few datasets instead of all of them.
##' @author Manuel Koller
##' @return list of fitted models. See also \code{\link{lapplyDatasets}} which
##'   is called internally.
##' @examples
##'   set.seed(1)
##'   oneWay <- generateAnovaDatasets(1, 1, 10, 4,
##'                                   lmeFormula = y ~ 1,
##'                                   heavyLmeRandom = ~ 1,
##'                                   heavyLmeGroups = ~ Var2,
##'                                   lqmmRandom = ~ 1,
##'                                   lqmmGroup = "Var2",
##'                                   groups = cbind(rep(1:4, each = 10), rep(1:10, 4)),
##'                                   varcov = matrix(1, 4, 4),
##'                                   lower = 0)
##'   fitDatasets_lmer(oneWay)
##' @aliases fitDatasets
##' @details \code{fitDatasets_lmer}: Fits datasets using \code{\link{lmer}}
##'   using its default options.
##' @rdname fitDatasets
##' @export
##' @importFrom lmer lmerControl
fitDatasets_lmer <-
    function(datasets,
             control,
             label,
             postFit,
             datasetIndices = "all") {
        if (missing(control)) {
            control <- lmerControl()
        }
        if (missing(label)) {
            label <- "fitDatasets_lmer"
        }
        fun <- function(data) {
            return(eval(substitute(
                lmer(formula, data = data, control = control),
                list(formula = datasets[["formula"]])
            )))
        }
        return(
            lapplyDatasets(
                datasets,
                fun,
                label = label,
                POST_FUN = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_lmer_bobyqa}: Fits datasets using \code{\link{lmer}} using
##'   the \code{bobyqa} optimizer.
##' @rdname fitDatasets
##' @export
fitDatasets_lmer_bobyqa <-
    function(datasets, postFit, datasetIndices = "all") {
        return(
            fitDatasets_lmer(
                datasets,
                lmerControl(optimizer = "bobyqa"),
                "fitDatasets_lmer_bobyqa",
                postFit,
                datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_lmer_Nelder_Mead}: Fits datasets using
##'   \code{\link{lmer}} using the \code{Nelder Mead} optimizer.
##' @rdname fitDatasets
##' @export
fitDatasets_lmer_Nelder_Mead <-
    function(datasets, postFit, datasetIndices = "all") {
        return(
            fitDatasets_lmer(
                datasets,
                lmerControl(optimizer = "Nelder_Mead"),
                "fitDatasets_lmer_Nelder_Mead",
                postFit,
                datasetIndices
            )
        )
    }

##' @importFrom Matrix isDiagonal
isDiagonalDataset <- function(datasets) {
    data <- datasets$generateData(1)
    Lambda <- datasets$createLambdaMatrix(data)
    return(isDiagonal(Lambda))
}

##' @param method argument passed on to \code{\link{rlmer}}.
##' @param tuningParameter argument passed on to
##'   \code{\link{extractTuningParameter}}.
##' @param init optional argument passed on to \code{\link{rlmer}}.
##' @param ... argument passed on to \code{\link{createRhoFunction}}.
##' @details \code{fitDatasets_rlmer}: Fits datasets using \code{\link{rlmer}}
##'   using a custom configuration. The argument 'tuningParameter' is passed to
##'   \code{\link{extractTuningParameter}}, details are documented there.
##' @examples
##'   ## call rlmer with custom arguments
##'   fitDatasets_rlmer_custom <- function(datasets) {
##'     return(fitDatasets_rlmer(datasets,
##'                              method = "DASvar",
##'                              tuningParameter = c(1.345, 2.28, 1.345, 2.28, 5.14, 5.14),
##'                              label = "fitDatasets_rlmer_custom"))
##'   }
##'   fitDatasets_rlmer_custom(oneWay)
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer <-
    function(datasets,
             method,
             tuningParameter,
             label,
             postFit,
             datasetIndices = "all",
             ...,
             init) {
        if (missing(method)) {
            method <- "DAStau"
        }
        if (missing(tuningParameter)) {
            tuningParameter <- c(1.345, 2.28, 1.345, 2.28, 5.14, 5.14)
        }
        if (missing(label)) {
            label <- "fitDatasets_rlmer_DAStau"
        }
        rho.e <- createRhoFunction(tuningParameter, "rho.e", ...)
        rho.sigma.e <-
            createRhoFunction(tuningParameter, "rho.sigma.e", ...)
        if (isDiagonalDataset(datasets)) {
            rho.b <- createRhoFunction(tuningParameter, "rho.b.diagonal", ...)
            rho.sigma.b <-
                createRhoFunction(tuningParameter, "rho.sigma.b.diagonal", ...)
        } else {
            rho.b <-
                createRhoFunction(tuningParameter, "rho.b.blockDiagonal", ...)
            rho.sigma.b <-
                createRhoFunction(tuningParameter, "rho.sigma.b.blockDiagonal", ...)
        }
        if (missing(init)) {
            fun <- function(data) {
                eval(substitute(
                    rlmer(
                        formula,
                        data = data,
                        method = method,
                        rho.e = rho.e,
                        rho.b = rho.b,
                        rho.sigma.e = rho.sigma.e,
                        rho.sigma.b = rho.sigma.b
                    ),
                    list(formula = datasets[["formula"]])
                ))
            }
        } else {
            fun <- function(data) {
                eval(substitute(
                    rlmer(
                        formula,
                        data = data,
                        method = method,
                        rho.e = rho.e,
                        rho.b = rho.b,
                        rho.sigma.e = rho.sigma.e,
                        rho.sigma.b = rho.sigma.b,
                        init = init
                    ),
                    list(formula = datasets[["formula"]])
                ))
            }
        }
        return(
            lapplyDatasets(
                datasets,
                fun,
                label = label,
                POST_FUN = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' Convenience function to create rho-functions with custom tuning parameter.
##'
##' 'rho.b.diagonal' denotes the tuning parameter to be used for 'rho.b' for
##' models with diagonal random effects covariance matrix. 'rho.b.blockDiagonal'
##' is the tuning parameter to be used in the block diagonal case, respectively.
##'
##' For arguments \code{rho.sigma.e} (and \code{rho.sigma.b.diagonal}), the
##' Proposal 2 variant of the function specified for \code{rho.e} (and
##' \code{rho.b}) is used.
##' @title Create Rho-Functions With Custom Tuning Parameter
##' @param tuningParameter argument passed on to
##'   \code{\link{extractTuningParameter}}. See its documentation for details.
##' @param which string specifiying which tuning parameter should be extracted.
##' @param rho.e \code{\link{PsiFunction}} to be used for \code{rho.e}.
##' @param rho.sigma.e \code{\link{PsiFunction}} to be used for
##'   \code{rho.sigma.e}.
##' @param rho.b.diagonal \code{\link{PsiFunction}} to be used for \code{rho.b}
##'   for models with diagonal random effects covariance matrix.
##' @param rho.sigma.b.diagonal \code{\link{PsiFunction}} to be used for
##'   \code{rho.sigma.b} for models with diagonal random effects covariance
##'   matrix.
##' @param rho.b.blockDiagonal \code{\link{PsiFunction}} to be used for
##'   \code{rho.b} for models with block-diagonal random effects covariance
##'   matrix.
##' @param rho.sigma.b.blockDiagonal \code{\link{PsiFunction}} to be used for
##'   \code{rho.sigma.b} for models with block-diagonal random effects
##'   covariance matrix.
##' @param ... passed on to \code{\link{chgDefaults}}.
##' @author Manuel Koller
##' @examples
##'   createRhoFunction(c(1.345, 2.28, 1.345, 2.28, 5.14, 5.14), "rho.sigma.e")
##' @export
createRhoFunction <-
    function(tuningParameter,
             which = c(
                 "rho.e",
                 "rho.sigma.e",
                 "rho.b.diagonal",
                 "rho.sigma.b.diagonal",
                 "rho.b.blockDiagonal",
                 "rho.sigma.b.blockDiagonal"
             ),
             rho.e = smoothPsi,
             rho.sigma.e = psi2propII(rho.e),
             rho.b.diagonal = rho.e,
             rho.sigma.b.diagonal = psi2propII(rho.b.diagonal),
             rho.b.blockDiagonal = rho.e,
             rho.sigma.b.blockDiagonal = rho.b.blockDiagonal,
             ...) {
        which <- match.arg(which)
        k <- extractTuningParameter(tuningParameter, which)
        psiFunction <- switch (
            which,
            rho.e = rho.e,
            rho.sigma.e = rho.sigma.e,
            rho.b.diagonal = rho.b.diagonal,
            rho.sigma.b.diagonal = rho.sigma.b.diagonal,
            rho.b.blockDiagonal = rho.b.blockDiagonal,
            rho.sigma.b.blockDiagonal = rho.sigma.b.blockDiagonal
        )
        return(chgDefaults(psiFunction, k = k, ...))
    }


##' Methods to extract which tuning parameters have been used for fitting
##' models. Use \code{extractTuningParameter} for custom configurations and
##' \code{extractPredefinedTuningParameter} for predefined configurations
##' provided in this package.
##' @title Extract Tuning Parameters Used In Fitting
##' @param tuningParameter vector of tuning parameters. The vector is expected
##'   to be of length 6, containing the tuning parameters for rho.e,
##'   rho.sigma.e, rho.b.diagonal, rho.sigma.b.diagonal, rho.b.blockDiagonal and
##'   rho.sigma.b.blockDiagonal. 'rho.b.diagonal' denotes the tuning parameter
##'   to be used for 'rho.b' for models with diagonal random effects covariance
##'   matrix. Names are optional.
##' @param which string specifiying which tuning parameter should be extracted.
##' @author Manuel Koller
##' @return scalar tuning parameter
##' @export
extractTuningParameter <-
    function(tuningParameter,
             which = c(
                 "rho.e",
                 "rho.sigma.e",
                 "rho.b.diagonal",
                 "rho.sigma.b.diagonal",
                 "rho.b.blockDiagonal",
                 "rho.sigma.b.blockDiagonal"
             )) {
        which <- match.arg(which)
        expectedNames <- c(
            "rho.e",
            "rho.sigma.e",
            "rho.b.diagonal",
            "rho.sigma.b.diagonal",
            "rho.b.blockDiagonal",
            "rho.sigma.b.blockDiagonal"
        )
        if (length(tuningParameter) != 6) {
            stop(
                "Argument 'tuningParameter' needs to be of length 6. ",
                "The contents are assumed to be the tuning parameters for ",
                "the rho functions in the order of [",
                paste(expectedNames, collapse = ", "),
                "]."
            )
        }
        if (is.null(names(tuningParameter))) {
            names(tuningParameter) <- expectedNames
        } else {
            diff <- setdiff(expectedNames, names(tuningParameter))
            if (length(diff) > 0) {
                stop(
                    "Expected names of argument 'tuningParameter' to be [",
                    paste(expectedNames, collapse = ", "),
                    "] but got [",
                    paste(names(tuningParameter), collapse = ", "),
                    "] instead."
                )
            }
        }
        return(tuningParameter[[which]])
    }

predefinedTuningParameter <- list()

##' @param label label or vector of labels in results. Only predefined labels of
##'   the form 'fitDatasets_rlmer_...' are supported (for others NA is
##'   returned).
##' @examples
##'   extractPredefinedTuningParameter("fitDatasets_rlmer_DAStau", "rho.e")
##' @rdname extractTuningParameter
##' @export
extractPredefinedTuningParameter <- function(label, which) {
    extractForOneLabel <- function(oneLabel) {
        tuningParameter <- predefinedTuningParameter[[oneLabel]]
        if (is.null(tuningParameter)) {
            return(NA_real_)
        } else {
            return(extractTuningParameter(tuningParameter, which))
        }
    }
    return(sapply(label, extractForOneLabel))
}

predefinedTuningParameter[["fitDatasets_rlmer_DAStau"]] <-
    c(1.345, 2.28, 1.345, 2.28, 5.14, 5.14)

##' @details \code{fitDatasets_rlmer_DAStau}: Fits datasets using
##'   \code{\link{rlmer}} using method DAStau and \code{\link{smoothPsi}} for
##'   the rho functions. The tuning parameters are k = 1.345 for \code{rho.e}.
##'   For \code{rho.sigma.e}, the Proposal 2 variant is used using k = 2.28.
##'   The choices for \code{rho.b} and \code{rho.sigma.b} depend on whether the
##'   model uses a diagonal or a block diagonal matrix for Lambda. In the former
##'   case, the same psi functions and tuning parameters are use as for
##'   \code{rho.e} and \code{rho.sigma.b}. In the block diagonal case,
##'   \code{rho.b} and \code{rho.sigma.b} both use \code{\link{smoothPsi}} using
##'   a tuning parameter k = 5.14 (assuming blocks of dimension 2).
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DAStau <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DAStau"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DAStau",
                tuningParameter = predefinedTuningParameter[[label]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_rlmer_DAStau_lmerNoFit}: Fits datasets using
##'   \code{\link{rlmer}} using the same configuration as
##'   \code{fitDatasets_rlmer_DAStau} except for that it is using
##'   \code{\link{lmerNoFit}} as initial estimator.
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DAStau_lmerNoFit <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DAStau_lmerNoFit"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DAStau",
                tuningParameter = predefinedTuningParameter[["fitDatasets_rlmer_DAStau"]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices,
                init = lmerNoFit
            )
        )
    }

predefinedTuningParameter[["fitDatasets_rlmer_DASvar"]] <-
    c(1.345, 2.28, 1.345, 2.28, 5.14, 5.14)

##' @details \code{fitDatasets_rlmer_DASvar}: Fits datasets using
##'   \code{\link{rlmer}} using method DASvar. The same rho functions and tuning
##'   parameters are used as for \code{fitDatasets_rlmer_DAStau}.
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DASvar <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DASvar"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DASvar",
                tuningParameter = predefinedTuningParameter[[label]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

predefinedTuningParameter[["fitDatasets_rlmer_DAStau_noAdj"]] <-
    c(1.345, 1.345, 1.345, 1.345, 5.14, 5.14)

##' @details \code{fitDatasets_rlmer_DAStau_noAdj}: Fits datasets using
##'   \code{\link{rlmer}} using method DAStau. The same rho functions and tuning
##'   parameters are used as for \code{fitDatasets_rlmer_DAStau}, except for
##'   \code{rho.sigma.e} (and \code{rho.sigma.b} in the diagonal case) for which
##'   the Proposal 2 variant of \code{smoothPsi} using k = 1.345 is used.
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DAStau_noAdj <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DAStau_noAdj"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DAStau",
                tuningParameter = predefinedTuningParameter[[label]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices,
            )
        )
    }

predefinedTuningParameter[["fitDatasets_rlmer_DAStau_k_0_5"]] <-
    c(0.5, 1.47, 0.5, 1.47, 2.17, 2.17)

##' @details \code{fitDatasets_rlmer_DAStau_k_0_5}: Fits datasets using
##'   \code{\link{rlmer}} using method DAStau. Use \code{\link{smoothPsi}}
##'   psi-function with tuning parameter \code{k = 0.5} for \code{rho.e} and
##'   \code{k = 1.47} for \code{rho.sigma.e}, the latter adjusted to reach the
##'   same asymptotic efficiency. In the diagonal case, the same are used for
##'   \code{rho.b} and \code{rho.sigma.b} as well. In the block-diagonal case,
##'   the tuning parameter \code{k = 2.17} is used for \code{rho.b} and
##'   \code{rho.sigma.b}. The tuning parameter is chosen to reach about the same
##'   asymptotic efficiency for theta as for the fixed effects.
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DAStau_k_0_5 <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DAStau_k_0_5"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DAStau",
                tuningParameter = predefinedTuningParameter[[label]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

predefinedTuningParameter[["fitDatasets_rlmer_DAStau_k_0_5_noAdj"]] <-
    c(0.5, 0.5, 0.5, 0.5, 2.17, 2.17)

##' @details \code{fitDatasets_rlmer_DAStau_k_0_5_noAdj}: Fits datasets using
##'   \code{\link{rlmer}} using method DAStau. Use \code{\link{smoothPsi}}
##'   psi-function with tuning parameter \code{k = 0.5} for \code{rho.e} and
##'   \code{rho.sigma.e}. In the diagonal case, the same are used for
##'   \code{rho.b} and \code{rho.sigma.b} as well. In the block-diagonal case,
##'   the tuning parameter \code{k = 2.17} is used for \code{rho.b} and
##'   \code{rho.sigma.b}. The tuning parameter is chosen to reach about the same
##'   asymptotic efficiency for theta as for the fixed effects.
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DAStau_k_0_5_noAdj <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DAStau_k_0_5_noAdj"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DAStau",
                tuningParameter = predefinedTuningParameter[[label]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

predefinedTuningParameter[["fitDatasets_rlmer_DAStau_k_2"]] <-
    c(2, 2.9, 2, 2.9, 8.44, 8.44)

##' @details \code{fitDatasets_rlmer_DAStau_k_2}: Fits datasets using
##'   \code{\link{rlmer}} using method DAStau. Use \code{\link{smoothPsi}}
##'   psi-function with tuning parameter \code{k = 2} for \code{rho.e} and
##'   \code{k = 2.9} \code{rho.sigma.e}, the latter adjusted to reach the same
##'   asymptotic efficiency. In the diagonal case, the same are used for
##'   \code{rho.b} and \code{rho.sigma.b} as well. In the block-diagonal case,
##'   the tuning parameter \code{k = 8.44} is used for \code{rho.b} and
##'   \code{rho.sigma.b}. The tuning parameter is chosen to reach about the same
##'   asymptotic efficiency for theta as for the fixed effects.
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DAStau_k_2 <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DAStau_k_2"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DAStau",
                tuningParameter = predefinedTuningParameter[[label]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

predefinedTuningParameter[["fitDatasets_rlmer_DAStau_k_2_noAdj"]] <-
    c(2, 2, 2, 2, 8.44, 8.44)

##' @details \code{fitDatasets_rlmer_DAStau_k_2_noAdj}: Fits datasets using
##'   \code{\link{rlmer}} using method DAStau. Use \code{\link{smoothPsi}}
##'   psi-function with tuning parameter \code{k = 2} for \code{rho.e} and
##'   \code{rho.sigma.e}. In the diagonal case, the same are used for
##'   \code{rho.b} and \code{rho.sigma.b} as well. In the block-diagonal case,
##'   the tuning parameter \code{k = 8.44} is used for \code{rho.b} and
##'   \code{rho.sigma.b}. The tuning parameter is chosen to reach about the same
##'   asymptotic efficiency for theta as for the fixed effects.
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DAStau_k_2_noAdj <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DAStau_k_2_noAdj"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DAStau",
                tuningParameter = predefinedTuningParameter[[label]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

predefinedTuningParameter[["fitDatasets_rlmer_DAStau_k_5"]] <-
    c(5, 5.03, 5, 5.03, 8.44, 8.44)

##' @details \code{fitDatasets_rlmer_DAStau_k_5}: Fits datasets using
##'   \code{\link{rlmer}} using method DAStau. Use \code{\link{smoothPsi}}
##'   psi-function with tuning parameter \code{k = 5} for \code{rho.e} and
##'   \code{k = 5.03} \code{rho.sigma.e}, the latter adjusted to reach the same
##'   asymptotic efficiency. In the diagonal case, the same are used for
##'   \code{rho.b} and \code{rho.sigma.b} as well. In the block-diagonal case,
##'   the tuning parameter \code{k = 34.21} is used for \code{rho.b} and
##'   \code{rho.sigma.b}. The tuning parameter is chosen to reach about the same
##'   asymptotic efficiency for theta as for the fixed effects.
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DAStau_k_5 <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DAStau_k_5"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DAStau",
                tuningParameter = predefinedTuningParameter[[label]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

predefinedTuningParameter[["fitDatasets_rlmer_DAStau_k_5_noAdj"]] <-
    c(5, 5, 5, 5, 8.44, 8.44)

##' @details \code{fitDatasets_rlmer_DAStau_k_5_noAdj}: Fits datasets using
##'   \code{\link{rlmer}} using method DAStau. Use \code{\link{smoothPsi}}
##'   psi-function with tuning parameter \code{k = 5} for \code{rho.e} and
##'   \code{rho.sigma.e}. In the diagonal case, the same are used for
##'   \code{rho.b} and \code{rho.sigma.b} as well. In the block-diagonal case,
##'   the tuning parameter \code{k = 34.21} is used for \code{rho.b} and
##'   \code{rho.sigma.b}. The tuning parameter is chosen to reach about the same
##'   asymptotic efficiency for theta as for the fixed effects.
##' @rdname fitDatasets
##' @export
fitDatasets_rlmer_DAStau_k_5_noAdj <-
    function(datasets, postFit, datasetIndices = "all") {
        label <- "fitDatasets_rlmer_DAStau_k_5_noAdj"
        return(
            fitDatasets_rlmer(
                datasets,
                method = "DAStau",
                tuningParameter = predefinedTuningParameter[[label]],
                label = label,
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_heavyLme}: Fits datasets using
##'   \code{heavyLme} from package \code{heavy}. Additional
##'   required arguments are: \code{lmeFormula}, \code{heavyLmeRandom} and
##'   \code{heavyLmeGroups}. They are passed to the \code{formula},
##'   \code{random} and \code{groups} arguments of \code{heavyLme}.
##' @rdname fitDatasets
##' @export
fitDatasets_heavyLme <-
    function(datasets, postFit, datasetIndices = "all") {
        packageInstalled <- isPackageInstalled("heavy")
        fun <- function(data) {
            if (!packageInstalled) {
                return(createPackageMissingReturnValue("heavy", "fitDatasets_heavyLme"))
            }
            fun <- eval(parse(text="heavy::heavyLme"))
            eval(substitute(
                fun(
                    formula,
                    random = random,
                    groups = groups,
                    data = data
                ),
                list(
                    formula = datasets[["lmeFormula"]],
                    random = datasets[["heavyLmeRandom"]],
                    groups = datasets[["heavyLmeGroups"]]
                )
            ))
        }
        return(
            lapplyDatasets(
                datasets,
                fun,
                label = "fitDatasets_heavyLme",
                POST_FUN = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

isPackageInstalled <- function(package) {
    return(length(find.package(package, quiet = TRUE)) > 0)
}

createPackageMissingReturnValue <- function(package, functionName) {
    value <-
        paste0(
            "Package '",
            package,
            "' is not installed. Please install it to run '",
            functionName,
            "'."
        )
    attr(value, "class") <- "package-not-installed"
    return(value)
}

##' @details \code{fitDatasets_lqmm}: Fits datasets using
##'   \code{\link[lqmm]{lqmm}} from package \code{lqmm}. Additional required
##'   arguments are: \code{lmeFormula}, \code{lqmmRandom}, \code{lqmmGroup} and
##'   \code{lqmmCovariance}. They are passed to the \code{formula},
##'   \code{random}, \code{groups} and \code{covariance} arguments of
##'   \code{lqmm}. \code{lqmmCovariance} is optional, if omitted \code{pdDiag}
##'   is used.
##' @rdname fitDatasets
##' @export
fitDatasets_lqmm <-
    function(datasets, postFit, datasetIndices = "all") {
        packageInstalled <- isPackageInstalled("lqmm")
        covariance <- datasets[["lqmmCovariance"]]
        if (is.null(covariance))
            covariance <- "pdDiag"
        fun <- function(data) {
            if (!packageInstalled) {
                return(createPackageMissingReturnValue("lqmm", "fitDatasets_lqmm"))
            }
            eval(substitute(
                lqmm::lqmm(
                    formula,
                    random = random,
                    group = group,
                    data = data,
                    covariance = covariance
                ),
                list(
                    formula = datasets[["lmeFormula"]],
                    random = datasets[["lqmmRandom"]],
                    group = datasets[["lqmmGroup"]],
                    covariance = covariance
                )
            ))
        }
        return(
            lapplyDatasets(
                datasets,
                fun,
                label = "fitDatasets_lqmm",
                POST_FUN = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_rlme}: Fits datasets using
##'   \code{\link[rlme]{rlme}} from package \code{rlme}.
##' @rdname fitDatasets
##' @export
fitDatasets_rlme <-
    function(datasets, postFit, datasetIndices = "all") {
        packageInstalled <- isPackageInstalled("rlme")
        fun <- function(data) {
            if (!packageInstalled) {
                return(createPackageMissingReturnValue("rlme", "fitDatasets_rlme"))
            }
            return(rlme::rlme(datasets[["formula"]], data))
        }
        return(
            lapplyDatasets(
                datasets,
                fun,
                label = "fitDatasets_rlme",
                POST_FUN = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_varComprob}: Prototype method to fit datasets
##'   using \code{varComprob} from package
##'   \code{robustvarComp}. Additional required items in \code{datasets} are:
##'   \code{lmeFormula}, \code{groups}, \code{varcov} and \code{lower}. They are
##'   passed to the \code{fixed}, \code{groups}, \code{varcov} and \code{lower}
##'   arguments of \code{varComprob}. The running of this method produces many
##'   warnings of the form "passing a char vector to .Fortran is not portable"
##'   which are suppressed.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob <-
    function(datasets,
             control,
             label,
             postFit,
             datasetIndices = "all") {
        packageInstalled <- isPackageInstalled("robustvarComp")
        fun <- function(data) {
            if (!packageInstalled) {
                return(createPackageMissingReturnValue("robustvarComp", label))
            }
            warnings <- NULL
            set.seed(attr(data, "datasetIndex"))
            result <- withCallingHandlers(
                eval(substitute(
                    robustvarComp::varComprob(
                        fixed = formula,
                        groups = groups,
                        data = data,
                        varcov = varcov,
                        control = control,
                        model = FALSE,
                        X = FALSE,
                        Y = FALSE,
                        K = FALSE
                    ),
                    list(
                        formula = datasets[["lmeFormula"]],
                        groups = datasets[["groups"]],
                        varcov = datasets[["varcov"]]
                    )
                )),
                warning = function(w) {
                    warnings <<- append(warnings, list(w))
                    invokeRestart("muffleWarning")
                }
            )
            for (condition in warnings) {
                if (!grepl(
                    "passing a char vector to .Fortran is not portable",
                    conditionMessage(condition),
                    fixed = TRUE
                )) {
                    warning(condition)
                }
            }
            return(result)
        }
        return(
            lapplyDatasets(
                datasets,
                fun,
                label = label,
                POST_FUN = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_varComprob_compositeTau}: Fits datasets with the
##'   composite Tau method using \code{varComprob} from
##'   package \code{robustvarComp}. See \code{fitDatasets_varComprob} for
##'   additional details.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob_compositeTau <-
    function(datasets, postFit, datasetIndices = "all") {
        lcontrol <-
            varComrob.control.or.null(lower = datasets[["lower"]])
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_compositeTau",
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

varComrob.control.or.null <- function(...) {
    if (isPackageInstalled("robustvarComp")) {
        return(robustvarComp::varComprob.control(...))
    } else {
        return(NULL)
    }
}

##' @details \code{fitDatasets_varComprob_compositeTau_OGK}: Similar to
##'   \code{fitDatasets_varComprob_compositeTau} but using \code{covOGK} as initial
##'   covariance matrix estimator.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob_compositeTau_OGK <-
    function(datasets, postFit, datasetIndices = "all") {
        lcontrol <-
            varComrob.control.or.null(lower = datasets[["lower"]],
                                      cov.init = "covOGK")
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_compositeTau_OGK",
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_varComprob_compositeTau_2SGS}: Similar to
##'   \code{fitDatasets_varComprob_compositeTau} but using \code{2SGS} as initial covariance
##'   matrix estimator.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob_compositeTau_2SGS <-
    function(datasets, postFit, datasetIndices = "all") {
        lcontrol <-
            varComrob.control.or.null(lower = datasets[["lower"]],
                                      cov.init = "2SGS")
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_compositeTau_2SGS",
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_varComprob_compositeS}: Similar to
##'   \code{fitDatasets_varComprob_compositeTau} but using method composite S.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob_compositeS <-
    function(datasets, postFit, datasetIndices = "all") {
        lcontrol <-
            varComrob.control.or.null(lower = datasets[["lower"]],
                                      method = "compositeS")
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_compositeS",
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_varComprob_compositeS_OGK}: Similar to
##'   \code{fitDatasets_varComprob_compositeS} but using \code{covOGK} as
##'   initial covariance matrix estimator.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob_compositeS_OGK <-
    function(datasets, postFit, datasetIndices = "all") {
        lcontrol <-
            varComrob.control.or.null(lower = datasets[["lower"]],
                                      method = "compositeS",
                                      cov.init = "covOGK")
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_compositeS_OGK",
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_varComprob_compositeS_2SGS}: Similar to
##'   \code{fitDatasets_varComprob_compositeS} but using \code{2SGS} as initial
##'   covariance matrix estimator.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob_compositeS_2SGS <-
    function(datasets, postFit, datasetIndices = "all") {
        lcontrol <-
            varComrob.control.or.null(lower = datasets[["lower"]],
                                      method = "compositeS",
                                      cov.init = "2SGS")
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_compositeS_2SGS",
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_varComprob_S}: Similar to
##'   \code{fitDatasets_varComprob_compositeTau} but using method S and the
##'   Rocke psi-function.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob_S <-
    function(datasets, postFit, datasetIndices = "all") {
        lcontrol <-
            varComrob.control.or.null(lower = datasets[["lower"]],
                                      method = "S",
                                      psi = "rocke")
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_S",
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_varComprob_S_OGK}: Similar to
##'   \code{fitDatasets_varComprob_S} but using \code{covOGK} as initial
##'   covariance matrix estimator.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob_S_OGK <-
    function(datasets, postFit, datasetIndices = "all") {
        lcontrol <-
            varComrob.control.or.null(
                lower = datasets[["lower"]],
                method = "S",
                psi = "rocke",
                cov.init = "covOGK"
            )
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_S_OGK",
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' @details \code{fitDatasets_varComprob_S_2SGS}: Similar to
##'   \code{fitDatasets_varComprob_S} but using \code{2SGS} as initial
##'   covariance matrix estimator.
##' @rdname fitDatasets
##' @export
fitDatasets_varComprob_S_2SGS <-
    function(datasets, postFit, datasetIndices = "all") {
        lcontrol <-
            varComrob.control.or.null(
                lower = datasets[["lower"]],
                method = "S",
                psi = "rocke",
                cov.init = "2SGS"
            )
        return(
            fitDatasets_varComprob(
                datasets,
                lcontrol,
                "fitDatasets_varComprob_S_2SGS",
                postFit = postFit,
                datasetIndices = datasetIndices
            )
        )
    }

##' Shorten labels created by the various \code{\link{fitDatasets}} functions,
##' for use in plotting, etc.
##'
##' The labels are shortened as they are in the simulation study published in
##' Koller and Stahel (2022).
##' @title Shorten Labels
##' @param labels vector of labels as assigned by \code{\link{fitDatasets}}
##' @author Manuel Koller
##' @return Vector of shortened labels
##' @references Koller M, Stahel WA (2022). "Robust Estimation of General Linear
##'  Mixed Effects Models.” In PM Yi, PK Nordhausen (eds.), Robust and Multivariate
##'  Statistical Methods, Springer Nature Switzerland AG.
##' @examples
##'   labels <- c("fitDatasets_lmer", "fitDatasets_rlmer_DAStau",
##'               "fitDatasets_rlmer_DAStau_noAdj",
##'               "fitDatasets_varComprob_compositeTau_OGK",
##'               "fitDatasets_varComprob_S_OGK",
##'               "fitDatasets_heavyLme",
##'               "fitDatasets_lqmm")
##'   shortenLabelsKS2022(labels)
##' @export
shortenLabelsKS2022 <- function(labels) {
    idx <- labels == "fitDatasets_lmer"
    labels[idx] <- "lme"
    idx <- labels == "fitDatasets_rlmer_DAStau"
    labels[idx] <- "RSEa"
    idx <- labels == "fitDatasets_rlmer_DAStau_noAdj"
    labels[idx] <- "RSEn"
    idx <- grepl("fitDatasets_varComprob_compositeTau", labels)
    labels[idx] <- "cTau"
    idx <- grepl("fitDatasets_varComprob_S", labels)
    labels[idx] <- "S"
    idx <- grepl("fitDatasets_heavyLme", labels)
    labels[idx] <- "t4"
    idx <- grepl("fitDatasets_lqmm", labels)
    labels[idx] <- "med"
    return(labels)
}
