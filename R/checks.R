## Small argument checks for the exported entry points: clear, early
## errors that name the argument.

.check_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x))
        stop("'", name, "' must be TRUE or FALSE", call. = FALSE)
    invisible(x)
}

.check_number <- function(x, name, lower = -Inf, upper = Inf,
                          integer = FALSE, open = FALSE) {
    ok <- is.numeric(x) && length(x) == 1L && !is.na(x) &&
        (if (open) x > lower && x < upper else x >= lower && x <= upper) &&
        (!integer || x == round(x))
    if (!ok) {
        rng <- if (is.finite(lower) && is.finite(upper))
                   paste0(" between ", lower, " and ", upper,
                          if (open) " (exclusive)")
               else if (is.finite(lower))
                   paste0(if (open) " greater than " else " at least ", lower)
               else if (is.finite(upper))
                   paste0(if (open) " less than " else " at most ", upper)
               else ""
        stop("'", name, "' must be a single ",
             if (integer) "whole " else "", "number", rng, call. = FALSE)
    }
    invisible(x)
}

.check_level <- function(x, name = "level")
    .check_number(x, name, lower = 0, upper = 1, open = TRUE)

## Guard against long-running computations: a message (never a prompt;
## the package runs in batch simulations) when the dimension-only workload
## of a computation exceeds a threshold option. `work` is in the units of
## .satterthwaite_df_workload(); `what` names the computation and `avoid`
## how to avoid it.
.slowMessage <- function(work, option, default, what, avoid) {
    max <- getOption(option, default)
    if (!is.numeric(max) || length(max) != 1L || is.na(max) || work <= max)
        return(invisible(FALSE))
    message(what, ": workload ", format(signif(work, 2), big.mark = ","),
            " (above the message threshold ",
            format(max, big.mark = ",", scientific = FALSE),
            "); this can take a long time. ", avoid,
            " Silence this message with options(", option, " = Inf).")
    invisible(TRUE)
}
