globalVariables(".data")

##' Plot longitudinal data with robustness-weight colored lines
##'
##' Creates a visualization of longitudinal data with one facet per treatment
##' group. Subject trajectories are colored by their robustness weight from a
##' robust mixed-effects model fit, with darker lines indicating lower weights
##' (potential outliers). Fixed-effect predictions are overlaid as reference
##' lines.
##'
##' @param data A data frame containing longitudinal data. Must have columns
##'   for subject ID, time, treatment group, and response variable.
##' @param formula A formula for the mixed-effects model. Default is
##'   \code{y ~ treatment * time + (1 + time | id)} where \code{y},
##'   \code{treatment}, \code{time}, and \code{id} refer to the standardized
##'   internal column names (mapped from \code{responseVar}, \code{treatmentVar},
##'   \code{timeVar}, and \code{idVar}).
##' @param idVar Character string naming the subject ID column in \code{data}.
##'   Default: \code{"id"}.
##' @param timeVar Character string naming the time column in \code{data}.
##'   Default: \code{"time"}.
##' @param treatmentVar Character string naming the treatment column in
##'   \code{data}. Default: \code{"treatment"}.
##' @param responseVar Character string naming the response column in
##'   \code{data}. Default: \code{"y"}.
##' @param rlmerArgs A list of additional arguments passed to
##'   \code{\link{rlmer}}.
##' @param lineAlpha Numeric in [0, 1]. Transparency of subject lines.
##'   Default: 0.6.
##' @param fixedLineWidth Numeric. Width of fixed-effect overlay lines.
##'   Default: 1.2.
##' @param lowColor Color for low robustness weights (potential outliers).
##'   Default: \code{"black"}.
##' @param highColor Color for high robustness weights (typical observations).
##'   Default: \code{"lightgray"}.
##' @param fixedLineColor Color for the fixed-effect prediction lines.
##'   Default: \code{"firebrick"}.
##' @param fixedLinetype Linetype for the fixed-effect prediction lines. Can be
##'   a single value (e.g., \code{"solid"}, \code{"dashed"}) applied to all
##'   lines, or \code{"byTreatment"} to use different linetypes for each
##'   treatment group. Default: \code{"solid"}.
##' @param title Optional plot title.
##' @param xlab Label for x-axis. If \code{NULL} (default), uses the value of
##'   \code{timeVar}.
##' @param ylab Label for y-axis. If \code{NULL} (default), uses the value of
##'   \code{responseVar}.
##'
##' @return A \code{ggplot} object.
##'
##' @details
##' The function fits a robust linear mixed-effects model using
##' \code{\link{rlmer}} and extracts the robustness weights for the random
##' effects. Subjects with low weights (shown in darker colors) are those whose
##' random effects deviate substantially from the assumed distribution.
##'
##' The fixed-effect prediction lines show the population-average trajectory
##' for each treatment group, ignoring random effects.
##'
##' @seealso \code{\link{rlmer}}, \code{\link{generateLongitudinalDatasets}}
##'
##' @examples
##' \dontrun{
##'   ## Using the medication dataset from confintROB
##'   library(confintROB)
##'   plotLongitudinalBySubject(
##'     medication,
##'     idVar = "id",
##'     treatmentVar = "treat",
##'     responseVar = "pos"
##'   )
##'
##'   ## Using simulated data
##'   set.seed(123)
##'   simdat <- generateLongitudinalDatasets(
##'     numberOfDatasetsToGenerate = 1,
##'     numberOfSubjects = 40,
##'     numberOfTimepoints = 7,
##'     numberOfTreatmentLevels = 2,
##'     timeRange = c(0, 18),
##'     trueBeta = c(200, -2, -5, 3),
##'     trueSigma = 30
##'   )
##'   plotLongitudinalBySubject(simdat$generateData(1))
##' }
##'
##' @export
##' @importFrom ggplot2 ggplot aes geom_line scale_color_gradient facet_wrap
##' @importFrom ggplot2 labs theme_bw theme .data scale_linetype_manual
##' @importFrom stats predict as.formula
plotLongitudinalBySubject <- function(data,
                                      formula = NULL,
                                      idVar = "id",
                                      timeVar = "time",
                                      treatmentVar = "treatment",
                                      responseVar = "y",
                                      rlmerArgs = list(),
                                      lineAlpha = 0.6,
                                      fixedLineWidth = 1.2,
                                      lowColor = "black",
                                      highColor = "lightgray",
                                      fixedLineColor = "firebrick",
                                      fixedLinetype = "solid",
                                      title = NULL,
                                      xlab = NULL,
                                      ylab = NULL) {
    ## Input validation
    if (!is.data.frame(data)) {
        stop("Argument 'data' must be a data frame.")
    }
    requiredVars <- c(idVar, timeVar, treatmentVar, responseVar)
    missingVars <- setdiff(requiredVars, names(data))
    if (length(missingVars) > 0L) {
        stop("Missing required columns in 'data': ",
             paste(missingVars, collapse = ", "))
    }
    if (!is.numeric(lineAlpha) || lineAlpha < 0 || lineAlpha > 1) {
        stop("Argument 'lineAlpha' must be a number between 0 and 1.")
    }
    if (!is.numeric(fixedLineWidth) || fixedLineWidth <= 0) {
        stop("Argument 'fixedLineWidth' must be a positive number.")
    }

    ## Create standardized data frame for internal use
    df <- data.frame(
        id = as.character(data[[idVar]]),
        time = as.numeric(data[[timeVar]]),
        treatment = factor(data[[treatmentVar]]),
        y = as.numeric(data[[responseVar]]),
        stringsAsFactors = FALSE
    )

    ## Check for complete cases
    complete_rows <- stats::complete.cases(df)
    if (!all(complete_rows)) {
        n_removed <- sum(!complete_rows)
        warning("Removed ", n_removed, " rows with missing values.")
        df <- df[complete_rows, , drop = FALSE]
    }

    ## Build formula if not provided
    if (is.null(formula)) {
        formula <- as.formula("y ~ treatment * time + (1 + time | id)")
    }

    ## Fit robust linear mixed-effects model
    fit <- suppressWarnings(do.call(rlmer, c(
        list(formula = formula, data = df), rlmerArgs
    )))

    ## Extract random-effect robustness weights
    w_b <- getME(fit, "w_b")
    wb_df <- w_b[[1L]]

    ## Get intercept weights (named "(Intercept)" or first column)
    if ("(Intercept)" %in% names(wb_df)) {
        wb_col <- wb_df[["(Intercept)"]]
    } else {
        wb_col <- wb_df[[1L]]
    }

    ## Create weight lookup and merge with data
    w_b_lookup <- stats::setNames(wb_col, row.names(ranef(fit)[[1L]]))
    df$robust_w <- w_b_lookup[df$id]

    ## Create fixed-effect prediction data
    time_range <- range(df$time, na.rm = TRUE)
    newdat <- expand.grid(
        treatment = levels(df$treatment),
        time = seq(time_range[1L], time_range[2L], length.out = 100L),
        KEEP.OUT.ATTRS = FALSE
    )
    newdat$y_fixed <- predict(fit, newdat, re.form = NA)

    ## Set axis labels
    if (is.null(xlab)) xlab <- timeVar
    if (is.null(ylab)) ylab <- responseVar

    ## Build plot using .data pronoun for robust column references
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$y,
                         group = .data$id, color = .data$robust_w)) +
        ggplot2::geom_line(alpha = lineAlpha) +
        ggplot2::scale_color_gradient(
            name = "Robustness weight",
            low = lowColor,
            high = highColor,
            limits = c(0, 1)
        ) +
        (if (identical(fixedLinetype, "byTreatment")) {
            ggplot2::geom_line(
                data = newdat,
                mapping = ggplot2::aes(x = .data$time, y = .data$y_fixed,
                                       group = .data$treatment,
                                       linetype = .data$treatment),
                linewidth = fixedLineWidth,
                color = fixedLineColor,
                inherit.aes = FALSE
            )
        } else {
            ggplot2::geom_line(
                data = newdat,
                mapping = ggplot2::aes(x = .data$time, y = .data$y_fixed,
                                       group = .data$treatment),
                linewidth = fixedLineWidth,
                color = fixedLineColor,
                linetype = fixedLinetype,
                inherit.aes = FALSE
            )
        }) +
        (if (identical(fixedLinetype, "byTreatment")) {
            ggplot2::scale_linetype_manual(
                values = c("solid", "dotted"),
                guide = "none"
            )
        } else {
            NULL
        }) +
        ggplot2::facet_wrap(~treatment) +
        ggplot2::labs(x = xlab, y = ylab, title = title) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
    p
}
