## rlmerControl(): control parameters of rlmer() in one object.

##' Control of the robust linear mixed model fit
##'
##' Collects the tuning parameters of \code{\link{rlmer}}'s fitting
##' algorithm in one object for its argument \code{control}.
##'
##' The arguments \code{rel.tol} and \code{max.iter} of \code{rlmer} keep
##' working. A setting may be given either there or in \code{control}, not
##' in both. To change a setting of \code{control} with \code{update},
##' update \code{control} itself.
##'
##' @section Retry:
##'   With \code{retry = TRUE}, a fit that does not converge (see the
##'   section \dQuote{Convergence} of \code{\link{rlmer}}) is fitted again
##'   along a chain of fallbacks, and the first converged fit is returned:
##'   \enumerate{
##'     \item \code{"iterations"}: continue from the unconverged estimate;
##'     \item \code{"DASvar"}: fit with \code{method = "DASvar"}, and start
##'       the requested method from that fit (skipped if the method already
##'       is \code{"DASvar"});
##'     \item \code{"robust"}: start from the high-breakdown start,
##'       \code{init = "robust"}.
##'   }
##'   Every attempt may use \code{retry.max.iter.factor} times as many
##'   iterations as the first fit.
##'   Which attempt succeeded is recorded in \code{object@optinfo$retry}
##'   and shown in the convergence block of \code{print} and
##'   \code{summary}. If no attempt converges, the first fit is returned
##'   with its warnings.
##'
##' @param rel.tol relative tolerance used as criterion in the fitting
##'   algorithms. Default \code{1e-8}.
##' @param max.iter maximum number of iterations. Default \code{NULL}:
##'   \eqn{40 (r + 1)^2} for \eqn{r} variance parameters.
##' @param dastau.mc,dasmc.nsim,dasmc.seed settings of the experimental
##'   Monte-Carlo calibration of the DAS-tau method for this fit. Default
##'   \code{NULL}: use the options \code{robustlmm.dastau.mc},
##'   \code{robustlmm.dasmc.nsim} and \code{robustlmm.dasmc.seed}, see
##'   \code{\link{robustlmm-options}}.
##' @param retry fit again if the fit does not converge? Default
##'   \code{FALSE}.
##' @param retry.steps the fallbacks to try, in this order; a subset of
##'   \code{c("iterations", "DASvar", "robust")}.
##' @param retry.max.iter.factor factor by which \code{max.iter} is
##'   multiplied for the retry attempts.
##' @param lmerControl a control object of \code{\link[lme4]{lmerControl}}
##'   for the classical fit that provides the starting values, or
##'   \code{NULL}.
##' @return A list of class \code{"rlmerControl"}.
##' @seealso \code{\link{rlmer}}, \code{\link{robustlmm-options}}
##' @examples
##' rlmerControl(max.iter = 200, retry = TRUE)
##' \dontrun{
##'   rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
##'         control = rlmerControl(retry = TRUE))
##' }
##' @section Lifecycle:
##'   Stable. Experimental, and subject to change: the retry
##'   (\code{retry}, \code{retry.steps}, \code{retry.max.iter.factor}) and the
##'   Monte-Carlo DAS-tau settings.
##' @export
rlmerControl <- function(rel.tol = 1e-8, max.iter = NULL,
                         dastau.mc = NULL, dasmc.nsim = NULL,
                         dasmc.seed = NULL, retry = FALSE,
                         retry.steps = c("iterations", "DASvar", "robust"),
                         retry.max.iter.factor = 5,
                         lmerControl = NULL) {
    .check_number(rel.tol, "rel.tol", lower = 0, open = TRUE)
    if (!is.null(max.iter))
        .check_number(max.iter, "max.iter", lower = 0, integer = TRUE)
    if (!is.null(dastau.mc)) .check_flag(dastau.mc, "dastau.mc")
    if (!is.null(dasmc.nsim))
        .check_number(dasmc.nsim, "dasmc.nsim", lower = 100, integer = TRUE)
    if (!is.null(dasmc.seed))
        .check_number(dasmc.seed, "dasmc.seed", integer = TRUE)
    .check_flag(retry, "retry")
    steps <- c("iterations", "DASvar", "robust")
    if (!is.character(retry.steps) || !length(retry.steps) ||
        anyDuplicated(retry.steps) || !all(retry.steps %in% steps))
        stop("'retry.steps' must be a subset of ",
             paste0("\"", steps, "\"", collapse = ", "), call. = FALSE)
    .check_number(retry.max.iter.factor, "retry.max.iter.factor", lower = 1)
    if (!is.null(lmerControl) && !inherits(lmerControl, "lmerControl"))
        stop("'lmerControl' must be NULL or the result of ",
             "lme4::lmerControl()", call. = FALSE)
    structure(list(rel.tol = rel.tol, max.iter = max.iter,
                   dastau.mc = dastau.mc, dasmc.nsim = dasmc.nsim,
                   dasmc.seed = dasmc.seed, retry = retry,
                   retry.steps = retry.steps,
                   retry.max.iter.factor = retry.max.iter.factor,
                   lmerControl = lmerControl),
              given = setdiff(names(match.call())[-1L], ""),
              class = "rlmerControl")
}

##' @export
print.rlmerControl <- function(x, ...) {
    cat("rlmer control parameters\n")
    v <- unclass(x); attr(v, "given") <- NULL
    v$lmerControl <- if (is.null(v$lmerControl)) NULL else "<lmerControl>"
    utils::str(v, give.head = FALSE, no.list = TRUE)
    invisible(x)
}

## Resolve the `control` argument of rlmer(): an rlmerControl object, or
## (as before this argument existed) an lmerControl object meant for the
## classical start, which rlmer() used to pass on to lmer() through `...`.
.resolveControl <- function(control) {
    if (inherits(control, "rlmerControl")) return(control)
    if (inherits(control, "lmerControl")) {
        res <- rlmerControl(lmerControl = control)
        attr(res, "given") <- character(0)
        return(res)
    }
    stop("'control' must be the result of rlmerControl() (or of ",
         "lme4::lmerControl(), which is passed to the classical start)",
         call. = FALSE)
}

## A setting given both as an argument of rlmer() and in `control`
.checkControlClash <- function(control, given) {
    both <- intersect(given, attr(control, "given"))
    if (length(both))
        stop("'", both[1L], "' is given both as an argument and in ",
             "'control'; give it in one place only", call. = FALSE)
    invisible(NULL)
}

## options() in force during a fit with these control settings
.controlOptions <- function(control) {
    op <- list(robustlmm.dastau.mc = control$dastau.mc,
               robustlmm.dasmc.nsim = control$dasmc.nsim,
               robustlmm.dasmc.seed = control$dasmc.seed)
    op[!vapply(op, is.null, logical(1))]
}
