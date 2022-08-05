##' Methods to process fitted objects and convert into a data structure that is
##' useful in post-processing.
##' @title Process Fitted Objects
##' @param obj object returned by the fitting method.
##' @param all logical, shorthand to enable all exports.
##' @param coefs logical, if true coefficients are added to export.
##' @param stdErrors logical, if true, standard errors are added to export.
##' @param tValues logical, if true, t-values are added to export.
##' @param sigma logical, if true, sigma is added to export.
##' @param thetas logical, if true, thetas are added to export.
##' @param b scalar logical or index vector, if true, all random effects are
##'   added to export. If an index vector is given, then only the corresponding
##'   random effects are added to the export. The same order as in \code{lmer}
##'   is used for all methods.
##' @param meanB logical, if true, the mean of the random effects is added to
##'   the export.
##' @param meanAbsB logical, if true, the mean of the absolute value of the
##'   random effects is added to the export.
##' @param residuals scalar logical or index vector, similar to argument
##'   \code{b}, just returning the residuals.
##' @param converged logical, if true, convergence code is added to export.
##' @param numWarnings logical, if true, the number of warnings generated during
##'   the fitting process is added to export.
##' @param procTime logical, if true, time needed to fit object is added to
##'   export.
##' @param ... optional parameters used for some implementations.
##' @return List with extracted values, most items can be suppressed
##'   to save disk space.
##' \describe{
##' \item{label}{Name of fitting method used to create the fit}
##' \item{datasetIndex}{Index of the dataset in the dataset list}
##' \item{coefficients}{Vector of estimated fixed-effects coefficients of the fitted model}
##' \item{standardErrors}{Vector of estimated standard errors of the fixed-effects coefficients}
##' \item{tValues}{Vector of t-Values (or z-Values depending on fitting method)
##'   of the fixed-effects coefficients}
##' \item{sigma}{Estimated residual standard error}
##' \item{thetas}{Vector of random-effects parameter estimates. As parameterized as by
##'   \code{\link{lmer}} and \code{\link{rlmer}}.}
##' \item{b}{Vector of requested predicted random-effects.}
##' \item{meanB}{Vector of means of the predicted random-effects.}
##' \item{meanAbsB}{Vector of means of the absolute values of the  predicted random-effects.}
##' \item{residuals}{Vector of requested residuals.}
##' \item{converged}{Convergence status as reported by the fitting method. \code{0} means converged.
##'   If not available, \code{NA} is used. Other values are to be interpreted carefully as codes
##'   vary from method to method.}
##' \item{numberOfWarnings}{the number of warnings generated during the fitting process.}
##' \item{proc.time}{Vector of times (user, system, elapsed) as reported by \code{\link{proc.time}}
##'   required to fit the model.}
##' }
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
##' @export
processFit <- function(obj,
                       all = FALSE,
                       coefs = TRUE,
                       stdErrors = all,
                       tValues = all,
                       sigma = TRUE,
                       thetas = TRUE,
                       b = all,
                       meanB = all,
                       meanAbsB = all,
                       residuals = all,
                       converged = TRUE,
                       numWarnings = all,
                       procTime = all,
                       ...)
    UseMethod("processFit")

##' @export
`processFit.try-error` <-
    function(obj,
             all = FALSE,
             coefs = TRUE,
             stdErrors = all,
             tValues = all,
             sigma = TRUE,
             thetas = TRUE,
             b = all,
             meanB = all,
             meanAbsB = all,
             residuals = all,
             converged = TRUE,
             numWarnings = all,
             procTime = all,
             ...)
        list()

##' @export
`processFit.package-not-installed` <- `processFit.try-error`

##' @examples
##'   processFit(fitDatasets_lmer(oneWay)[[1]], all = TRUE)
##' @rdname processFit
##' @export
processFit.lmerMod <- function(obj,
                               all = FALSE,
                               coefs = TRUE,
                               stdErrors = all,
                               tValues = all,
                               sigma = TRUE,
                               thetas = TRUE,
                               b = all,
                               meanB = all,
                               meanAbsB = all,
                               residuals = all,
                               converged = TRUE,
                               numWarnings = all,
                               procTime = all,
                               ...) {
    smry <- summary(obj)
    coef <- coef(smry)

    ret <- list(label = attr(obj, "label"),
                datasetIndex = attr(obj, "datasetIndex"))

    if (coefs) {
        ret[["coefficients"]] <- coef[, 1]
    }
    if (stdErrors) {
        ret[["standardErrors"]] <- coef[, 2]
    }
    if (tValues) {
        ret[["tValues"]] <- coef[, 3]
    }
    if (sigma) {
        ret[["sigma"]] <- sigma(obj)
    }
    if (thetas) {
        theta <- getME(obj, "theta")
        if (length(theta) == 1) {
            theta <- unname(theta)
        }
        ret[["thetas"]] <- theta
    }
    if (!isFALSE(b)) {
        randomEffects <- as.vector(getME(obj, "b"))
        if (length(b) > 1 || !is.logical(b)) {
            randomEffects <- randomEffects[b]
        }
        ret[["b"]] <- randomEffects
    }
    if (meanB) {
        ret[["meanB"]] <- unlist(lapply(ranef(obj), colMeans))
    }
    if (meanAbsB) {
        ret[["meanAbsB"]] <- unlist(lapply(ranef(obj), function(x) colMeans(abs(x))))
    }
    if (!isFALSE(residuals)) {
        resid <- resid(obj)
        if (length(residuals) > 1 || !is.logical(residuals)) {
            resid <- resid[residuals]
        }
        ret[["residuals"]] <- resid
    }
    if (converged) {
        opt <- length(obj@optinfo$conv$lme4)
        if (is.null(opt) || opt == 0) {
            opt <- obj@optinfo$conv$opt
        }
        ret[["converged"]] <- if (is.null(opt))
            NA_integer_
        else
            opt
    }
    if (numWarnings) {
        warnings <- attr(obj, "warnings")
        if (is.null(warnings)) {
            ret[["numberOfWarnings"]] <- 0
        } else {
            ret[["numberOfWarnings"]] <- length(warnings)
        }
    }
    if (procTime) {
        ret[["proc.time"]] <- attr(obj, "proc.time")
    }
    return(ret)
}

##' @examples
##'   processFit(fitDatasets_rlmer_DASvar(oneWay)[[1]], all = TRUE)
##' @rdname processFit
##' @export
processFit.rlmerMod <- processFit.lmerMod

##' @examples
##'   if (require(heavy)) {
##'     processFit(fitDatasets_heavyLme(oneWay)[[1]], all = TRUE)
##'   }
##' @rdname processFit
##' @export
processFit.heavyLme <- function(obj,
                                all = FALSE,
                                coefs = TRUE,
                                stdErrors = all,
                                tValues = all,
                                sigma = TRUE,
                                thetas = TRUE,
                                b = all,
                                meanB = all,
                                meanAbsB = all,
                                residuals = all,
                                converged = TRUE,
                                numWarnings = all,
                                procTime = all,
                                ...) {
    smry <- summary(obj)
    coef <- coef(smry)

    ret <- list(label = attr(obj, "label"),
                datasetIndex = attr(obj, "datasetIndex"))

    if (coefs) {
        ret[["coefficients"]] <- coef[, 1]
    }
    if (stdErrors) {
        ret[["standardErrors"]] <- coef[, 2]
    }
    if (tValues) {
        ret[["tValues"]] <- coef[, 3]
    }
    if (sigma) {
        ret[["sigma"]] <- obj$scale
    }
    if (thetas) {
        theta <- obj$theta / obj$scale
        if (length(theta) == 1) {
            theta <- unname(theta)
        }
        ret[["thetas"]] <- theta[upper.tri(theta, diag = TRUE)]
    }
    if (!isFALSE(b)) {
        randomEffects <- c(t(obj$ranef))
        if (length(b) > 1 || !is.logical(b)) {
            randomEffects <- randomEffects[b]
        }
        ret[["b"]] <- randomEffects
    }
    if (meanB) {
        ret[["meanB"]] <- colMeans(obj$ranef)
    }
    if (meanAbsB) {
        ret[["meanAbsB"]] <- colMeans(abs(obj$ranef))
    }
    if (!isFALSE(residuals)) {
        resid <- obj$Resid$marginal
        if (length(residuals) > 1 || !is.logical(residuals)) {
            resid <- resid[residuals]
        }
        ret[["residuals"]] <- resid
    }
    if (converged) {
        ret[["converged"]] <- as.numeric(!obj$converged)
    }
    if (numWarnings) {
        warnings <- attr(obj, "warnings")
        if (is.null(warnings)) {
            ret[["numberOfWarnings"]] <- 0
        } else {
            ret[["numberOfWarnings"]] <- length(warnings)
        }
    }
    if (procTime) {
        ret[["proc.time"]] <- attr(obj, "proc.time")
    }
    return(ret)
}

##' @examples
##'   if (require(lqmm)) {
##'     processFit(fitDatasets_lqmm(oneWay)[[1]], all = TRUE)
##'   }
##' @rdname processFit
##' @export
processFit.lqmm <- function(obj,
                            all = FALSE,
                            coefs = TRUE,
                            stdErrors = all,
                            tValues = all,
                            sigma = TRUE,
                            thetas = TRUE,
                            b = all,
                            meanB = all,
                            meanAbsB = all,
                            residuals = all,
                            converged = TRUE,
                            numWarnings = all,
                            procTime = all,
                            ...) {
    smry <- summary(obj, seed = attr(obj$mfArgs$data, "datasetIndex"))
    coef <- smry$tTable

    ret <- list(label = attr(obj, "label"),
                datasetIndex = attr(obj, "datasetIndex"))

    if (coefs) {
        ret[["coefficients"]] <- coef[, 1]
    }
    if (stdErrors) {
        ret[["standardErrors"]] <- coef[, 2]
    }
    if (tValues) {
        ret[["tValues"]] <- coef[, 1] / coef[, 2]
    }
    if (sigma) {
        ret[["sigma"]] <- obj$scale
    }
    if (thetas) {
        theta <- lqmm::covHandling(
            obj$theta_z,
            n = obj$dim_theta[2],
            cov_name = obj$cov_name,
            quad_type = obj$type
        )
        if (is.matrix(theta)) {
            theta <- chol(theta)
            theta <- theta[upper.tri(theta, diag = TRUE)]
        } else {
            theta <- sqrt(theta)
        }
        ret[["thetas"]] <- theta / obj$scale
    }
    if (!isFALSE(b)) {
        randomEffects <- c(t(ranef(obj)))
        if (length(b) > 1 || !is.logical(b)) {
            randomEffects <- randomEffects[b]
        }
        ret[["b"]] <- randomEffects
    }
    if (!isFALSE(residuals)) {
        resid <- resid(obj)
        if (length(residuals) > 1 || !is.logical(residuals)) {
            resid <- resid[residuals]
        }
        ret[["residuals"]] <- resid
    }
    if (meanB) {
        ret[["meanB"]] <- colMeans(ranef(obj))
    }
    if (meanAbsB) {
        ret[["meanAbsB"]] <- colMeans(abs(ranef(obj)))
    }
    if (converged) {
        ret[["converged"]] <- NA_real_
    }
    if (numWarnings) {
        warnings <- attr(obj, "warnings")
        if (is.null(warnings)) {
            ret[["numberOfWarnings"]] <- 0
        } else {
            ret[["numberOfWarnings"]] <- length(warnings)
        }
    }
    if (procTime) {
        ret[["proc.time"]] <- attr(obj, "proc.time")
    }
    return(ret)
}

##' @rdname processFit
##' @export
processFit.rlme <- function(obj,
                            all = FALSE,
                            coefs = TRUE,
                            stdErrors = all,
                            tValues = all,
                            sigma = TRUE,
                            thetas = TRUE,
                            b = all,
                            meanB = all,
                            meanAbsB = all,
                            residuals = all,
                            converged = TRUE,
                            numWarnings = all,
                            procTime = all,
                            ...) {
    coef <- obj$fixed.effects

    ret <- list(label = attr(obj, "label"),
                datasetIndex = attr(obj, "datasetIndex"))

    if (coefs) {
        ret[["coefficients"]] <- coef[, 2]
    }
    if (stdErrors) {
        ret[["standardErrors"]] <- coef[, 3]
    }
    if (tValues) {
        ret[["tValues"]] <- unname(obj$t.value)
    }
    haveRandomEffects <- "random.effects" %in% attributes(obj)$names
    if (sigma && haveRandomEffects) {
        idx <- NROW(obj$random.effects)
        ret[["sigma"]] <- sqrt(obj$random.effects[idx, 3])
    }
    if (thetas && haveRandomEffects) {
        idx <- 1:(NROW(obj$random.effects) - 1)
        ret[["thetas"]] <- sqrt(obj$random.effects[idx, 3]) / ret[["sigma"]]
    }
    if (!isFALSE(b)) {
        ret[["b"]] <- rep.int(NA_real_, length(b))
    }
    if (!isFALSE(residuals)) {
        resid <- obj$standard.residual
        if (length(residuals) > 1 || !is.logical(residuals)) {
            resid <- resid[residuals]
        }
        ret[["residuals"]] <- resid
    }
    if (meanB) {
        ret[["meanB"]] <- NA_real_
    }
    if (meanAbsB) {
        ret[["meanAbsB"]] <- NA_real_
    }
    if (converged) {
        ret[["converged"]] <- NA_real_
    }
    if (numWarnings) {
        warnings <- attr(obj, "warnings")
        if (is.null(warnings)) {
            ret[["numberOfWarnings"]] <- 0
        } else {
            ret[["numberOfWarnings"]] <- length(warnings)
        }
    }
    if (procTime) {
        ret[["proc.time"]] <- attr(obj, "proc.time")
    }
    return(ret)
}

##' @param isInterceptCorrelationSlopeModel optional logical, can be used to
##'   override the assumption that a model with three variance components can be
##'   interpreted as having intercept, correlation and slope.
##' @details Warning. \code{processFit.varComprob} uses simplistic logic to
##'   convert from the parameterisation used in the robustvarComp package to
##'   \code{theta} as used in \code{\link{lmer}} and \code{\link{rlmer}}. If
##'   there are three variance components, the code assumes that they are
##'   intercept, correlation and slope. Otherwise the code assumes that the
##'   variance components are independent. Exports \code{b} and \code{residuals}
##'   are not supported.
##' @examples
##'   if (require(robustvarComp)) {
##'     processFit(fitDatasets_varComprob_compositeTau(oneWay)[[1]], all = TRUE)
##'   }
##' @rdname processFit
##' @export
processFit.varComprob <- function(obj,
                                  all = FALSE,
                                  coefs = TRUE,
                                  stdErrors = all,
                                  tValues = all,
                                  sigma = TRUE,
                                  thetas = TRUE,
                                  b = all,
                                  meanB = all,
                                  meanAbsB = all,
                                  residuals = all,
                                  converged = TRUE,
                                  numWarnings = all,
                                  procTime = all,
                                  isInterceptCorrelationSlopeModel,
                                  ...) {
    ret <- list(label = attr(obj, "label"),
                datasetIndex = attr(obj, "datasetIndex"))

    if (coefs) {
        ret[["coefficients"]] <- obj$beta
    }
    if (stdErrors || tValues) {
        if (length(obj$beta) == 1) {
            standardErrors <- sqrt(obj$vcov.beta)
        } else {
            standardErrors <- sqrt(diag(obj$vcov.beta))
        }
    }
    if (stdErrors) {
        ret[["standardErrors"]] <- standardErrors
    }
    if (tValues) {
        ret[["tValues"]] <- obj$beta / standardErrors
    }
    if (sigma) {
        ret[["sigma"]] <- sqrt(obj$eta0)
    }
    if (thetas) {
        v <- obj$eta / obj$sigma2
        if (missing(isInterceptCorrelationSlopeModel)) {
            isInterceptCorrelationSlopeModel <- length(v) == 3
        }
        if (isInterceptCorrelationSlopeModel) {
            theta <- c(sqrt(v[1]), v[2] / sqrt(v[1]),
                       sqrt(v[3] - v[2] ^ 2 / v[1]))
        } else {
            theta <- sqrt(v)
        }
        if (length(theta) > 1) {
            names(theta) <- obj$random.labels
        }
        ret[["thetas"]] <- theta
    }
    if (!isFALSE(b)) {
        ret[["b"]] <- rep.int(NA_real_, length(b))
    }
    if (!isFALSE(residuals)) {
        nobs <- length(residuals)
        if (nobs == 1) {
            nobs <- obj$nobs
        }
        ret[["residuals"]] <- rep.int(NA_real_, nobs)
    }
    if (meanB) {
        ret[["meanB"]] <- NA_real_
    }
    if (meanAbsB) {
        ret[["meanAbsB"]] <- NA_real_
    }
    if (converged) {
        ret[["converged"]] <-
            as.numeric(obj$iterations == obj$control$max.it)
    }
    if (numWarnings) {
        warnings <- attr(obj, "warnings")
        if (is.null(warnings)) {
            ret[["numberOfWarnings"]] <- 0
        } else {
            ret[["numberOfWarnings"]] <- length(warnings)
        }
    }
    if (procTime) {
        ret[["proc.time"]] <- attr(obj, "proc.time")
    }
    return(ret)
}
