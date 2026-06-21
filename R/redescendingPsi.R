##' Construct a redescending Tukey bisquare psi function.
##'
##' Returns a \code{psi_func_rcpp} object whose psi, rho, weight and
##' derivative slots implement Tukey's bisquare (biweight) function
##' \deqn{\psi(x) = x (1 - (x/c)^2)^2 \quad \text{for } |x| \le c,\quad 0 \text{ otherwise}}
##' suitable for use as \code{rho.e} (or \code{rho.b}, \code{rho.sigma.e},
##' \code{rho.sigma.b}) in \code{\link{rlmer}}.
##'
##' The Fisher consistency expectations (\code{Erho}, \code{Epsi2},
##' \code{EDpsi}) are computed by numerical integration against the
##' standard normal.
##'
##' @title Bisquare (Tukey biweight) psi function
##' @param c tuning cutoff. The default \code{c = 4.685} gives about 95\%
##'   asymptotic efficiency for the location problem.
##' @return psi_func_rcpp object usable in \code{rho.e}, \code{rho.b}.
##' @examples
##'   pf <- bisquarePsi
##'   pf@@psi(c(-6, -3, 0, 3, 6))
##' @export
makeBisquarePsi <- function(c = 4.685) {
    if (!is.numeric(c) || length(c) != 1L || !is.finite(c) || c <= 0)
        stop("'c' must be a finite positive scalar")

    ## c = Inf is the unbounded (classical) limit; c = 0 the all-rejecting limit.
    ## Handle them explicitly so robustbase's functionXal validity check
    ## (which evaluates the function at c = 0 and c = Inf) does not produce
    ## NaN from Inf * 0 in rho.
    psi_fn  <- function(x, c) {
        if (!is.finite(c)) return(x)
        if (c <= 0) return(rep(0, length(x)))
        u <- x / c
        ifelse(abs(u) <= 1, x * (1 - u^2)^2, 0)
    }
    wgt_fn  <- function(x, c) {
        if (!is.finite(c)) return(rep(1, length(x)))
        if (c <= 0) return(rep(0, length(x)))
        u <- x / c
        ifelse(abs(u) <= 1, (1 - u^2)^2, 0)
    }
    Dpsi_fn <- function(x, c) {
        if (!is.finite(c)) return(rep(1, length(x)))
        if (c <= 0) return(rep(0, length(x)))
        u <- x / c
        ifelse(abs(u) <= 1, (1 - u^2) * (1 - 5 * u^2), 0)
    }
    Dwgt_fn <- function(x, c) {
        if (!is.finite(c) || c <= 0) return(rep(0, length(x)))
        u <- x / c
        ifelse(abs(u) <= 1, -4 * x / c^2 * (1 - u^2), 0)
    }
    rho_fn  <- function(x, c) {
        if (!is.finite(c)) return(x^2 / 2)
        if (c <= 0) return(rep(0, length(x)))
        u <- x / c
        ifelse(abs(u) <= 1, c^2 / 6 * (1 - (1 - u^2)^3), c^2 / 6)
    }

    U <- stats::qnorm(1 - 1e-10)
    Erho <- function(c) {
        vapply(c, function(cc)
            2 * stats::integrate(function(z) rho_fn(z, cc) * stats::dnorm(z),
                                 0, U, rel.tol = 1e-8)$value,
            numeric(1))
    }
    Epsi2 <- function(c) {
        vapply(c, function(cc)
            2 * stats::integrate(function(z) psi_fn(z, cc)^2 * stats::dnorm(z),
                                 0, U, rel.tol = 1e-8)$value,
            numeric(1))
    }
    EDpsi <- function(c) {
        vapply(c, function(cc)
            2 * stats::integrate(function(z) Dpsi_fn(z, cc) * stats::dnorm(z),
                                 0, U, rel.tol = 1e-8)$value,
            numeric(1))
    }

    pf <- psiFunc(rho = rho_fn, psi = psi_fn, wgt = wgt_fn,
                  Dpsi = Dpsi_fn, Dwgt = Dwgt_fn,
                  Erho = Erho, Epsi2 = Epsi2, EDpsi = EDpsi,
                  name = sprintf("Tukey bisquare (c = %.4g)", c),
                  c = c)

    args <- list("psi_func_rcpp")
    for (sn in methods::slotNames("psi_func")) args[[sn]] <- methods::slot(pf, sn)
    args[["getRcppClass"]] <- function() character(0)
    args[["getInstanceWithOriginalDefaults"]] <- function() NULL
    do.call("new", args)
}

##' Default Tukey bisquare psi-function with \code{c = 4.685}.
##'
##' Pre-built \code{psi_func_rcpp} from \code{\link{makeBisquarePsi}} with
##' the conventional tuning constant \code{c = 4.685} (\eqn{\approx 95\%}
##' efficiency at the normal model).
##' @name bisquarePsi
##' @aliases bisquarePsi
##' @export bisquarePsi
NULL

setLoadAction(function(ns)
    assign("bisquarePsi", makeBisquarePsi(c = 4.685), envir = ns))
