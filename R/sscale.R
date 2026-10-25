## S-type (bounded-rho) scale psi-function for rlmer's scale equations.
##
## rlmer estimates the residual scale (and one-dimensional random-effect
## scales) from equations of the form
##   sum_i eta_i tau_i^2 w(t_i) (t_i^2 - kappa) = 0,   t_i = r_i / (tau_i sigma),
## with kappa = E[w(Z) Z^2] / E[w(Z)] (calcKappaTau). psi2propII() supplies
## w = squared Huber-type weights, for which chi(t) = w(t)(t^2 - kappa)
## increases to a finite limit chi(inf): the equation's breakdown point is
## kappa / (kappa + chi(inf)) -- 0.295 for the default (k = 1.345) and 0.156
## for k = 2.28 -- and the scale is inflated well before that.
##
## With the weight
##   W(t) = (rho(t) - b) / (t^2 - kappa0),  b = E rho(Z),  rho(sqrt(kappa0)) = b,
## (rho the bisquare rho normalised to sup 1; the singularity at
## t^2 = kappa0 is removable) the same equation becomes
##   sum_i eta_i tau_i^2 (rho(t_i) - b) = 0,
## the M-scale with bounded rho of S-estimation, whose breakdown point is
## min(b, 1 - b); calcKappaTau() returns kappa0 exactly because
## E[W(Z)(Z^2 - kappa0)] = E rho(Z) - b = 0.

##' S-type scale psi-function for high-breakdown scale estimation.
##'
##' Returns a \code{psi_func_rcpp} object for the scale arguments
##' \code{rho.sigma.e} and (for one-dimensional random-effect blocks)
##' \code{rho.sigma.b} of \code{\link{rlmer}}. With it the scale equations
##' become M-scale equations with a bounded \eqn{\rho}, as in
##' S-estimation: \eqn{\sum_i \eta_i \tau_i^2 [\rho(t_i) - b] = 0} with
##' \eqn{t_i = r_i/(\tau_i\sigma)}, \eqn{\rho} the Tukey bisquare
##' \eqn{\rho} scaled to \eqn{\sup\rho = 1} with tuning constant \code{c},
##' and \eqn{b = E\rho(Z)} for standard normal \eqn{Z} (Fisher consistency).
##' The breakdown point of the scale equation is \eqn{\min(b, 1 - b)}: 1/2
##' for the default \code{c = 1.547645}, against 0.295 for the default
##' \code{psi2propII(smoothPsi)} scale and 0.156 for
##' \code{psi2propII(smoothPsi, k = 2.28)} (\code{setting = "RSEa"}).
##'
##' It is implemented as a weight function for the existing equations,
##' \eqn{W(t) = (\rho(t) - b)/(t^2 - \kappa_0)} with \eqn{\rho(\sqrt{\kappa_0})
##' = b}, so the design-adaptive (DAS) standardisation, the consistency
##' constant and everything built on the scale estimates (standard errors,
##' Satterthwaite degrees of freedom) carry over unchanged.
##'
##' A high-breakdown scale does not make the fit high-breakdown by itself:
##' the fixed effects also need a redescending \code{rho.e} (e.g.
##' \code{\link{bisquarePsi}}) and, when the start is contaminated, a
##' high-breakdown start (\code{init = "robust"}); see \code{\link{rlmer}}.
##' The S-scale is less efficient than the default scale estimate: its
##' asymptotic efficiency for the scale at the normal is 0.54 with
##' \code{c = 1.547645}; a larger \code{c} trades breakdown point for
##' efficiency (\code{c = 2.937015} gives \eqn{b = 1/4}, breakdown point
##' 1/4 and efficiency 0.87). It is meant for \code{rho.sigma.e} and simple
##' variance components; \code{rlmer} refuses it for blocks of correlated
##' random effects and as \code{rho.e} or \code{rho.b}. To change \code{c},
##' use \code{makeSScalePsi(c = ...)}, not \code{chgDefaults}.
##'
##' @title S-type scale psi-function
##' @param c tuning constant of the bisquare \eqn{\rho}; the default
##'   1.547645 gives \eqn{b = 1/2} and breakdown point 1/2.
##' @return a \code{psi_func_rcpp} object; attribute \code{"sscale"} holds
##'   \code{b}, \code{kappa0} and \code{c}.
##' @references Rousseeuw PJ, Yohai VJ (1984). Robust regression by means
##'   of S-estimators. In \emph{Robust and Nonlinear Time Series Analysis},
##'   Lecture Notes in Statistics 26, 256--272. Springer.
##' @seealso \code{\link{psi2propII}}, \code{\link{rlmer}}.
##' @examples
##'   ss <- makeSScalePsi()
##'   attr(ss, "sscale")
##'   ## the bounded chi(t) = W(t) (t^2 - kappa0) of the scale equation
##'   t <- c(0, 0.5, 1, 2, 5, 50)
##'   ss@@wgt(t) * (t^2 - attr(ss, "sscale")$kappa0)
##' @section Lifecycle:
##'   \strong{Experimental.} The interface and the results may change in a
##'   later version.
##' @export
makeSScalePsi <- function(c = 1.547645) {
    if (!is.numeric(c) || length(c) != 1L || !is.finite(c) || c <= 0)
        stop("'c' must be a finite positive scalar")
    c0 <- c
    rhoInf <- robustbase::MrhoInf(c0, "bisquare")
    rho  <- function(x) robustbase::Mchi(x, c0, "bisquare")
    drho <- function(x) robustbase::Mpsi(x, c0, "bisquare", 0L) / rhoInf
    d2rho <- function(x) robustbase::Mpsi(x, c0, "bisquare", 1L) / rhoInf
    b <- stats::integrate(function(z) rho(z) * stats::dnorm(z), -Inf, Inf,
                          rel.tol = 1e-12)$value
    t0 <- stats::uniroot(function(t) rho(t) - b, c(1e-8, c0),
                         tol = 1e-14)$root
    kappa0 <- t0^2
    eps <- 1e-6
    W <- function(x, cc) {
        x <- abs(x)
        den <- x^2 - kappa0
        out <- (rho(x) - b) / den
        near <- !is.na(den) & abs(den) < eps
        if (any(near)) out[near] <- drho(x[near]) / (2 * x[near])
        out
    }
    Dwgt <- function(x, cc) {
        s <- sign(x); x <- abs(x)
        den <- x^2 - kappa0
        out <- (drho(x) * den - 2 * x * (rho(x) - b)) / den^2
        near <- !is.na(den) & abs(den) < eps
        if (any(near))
            out[near] <- (d2rho(x[near]) * x[near] - drho(x[near])) /
                (4 * x[near]^2)
        s * out
    }
    psi  <- function(x, cc) W(x) * x
    Dpsi <- function(x, cc) W(x) + x * Dwgt(x)
    rho_fn <- function(x, cc) rho(x)
    ## E-slots vectorise over cc for robustbase's functionXal1 validity
    ## check; the values are fixed at construction.
    U <- stats::qnorm(1 - 1e-10)
    mkE <- function(fn) {
        val <- 2 * stats::integrate(function(z) fn(z) * stats::dnorm(z), 0, U,
                                    rel.tol = 1e-8)$value
        function(cc) rep(val, length(cc))
    }
    pf <- psiFunc(rho = rho_fn, psi = psi, wgt = W, Dpsi = Dpsi, Dwgt = Dwgt,
                  Erho = mkE(rho), Epsi2 = mkE(function(x) psi(x)^2),
                  EDpsi = mkE(Dpsi),
                  name = sprintf("S-scale (bisquare rho, c = %.4g, b = %.3g)",
                                 c0, b),
                  cc = c0)
    args <- list("psi_func_rcpp")
    for (sn in methods::slotNames("psi_func")) args[[sn]] <- methods::slot(pf, sn)
    args[["getRcppClass"]] <- function() character(0)
    args[["getInstanceWithOriginalDefaults"]] <- function() NULL
    obj <- do.call("new", args)
    attr(obj, "sscale") <- list(b = b, kappa0 = kappa0, c = c0)
    obj
}

##' Default S-type scale psi-function
##'
##' Pre-built \code{psi_func_rcpp} from \code{\link{makeSScalePsi}} with
##' \code{c = 1.547645} (breakdown point 1/2), for \code{rho.sigma.e} and
##' one-dimensional \code{rho.sigma.b}.
##' @section Lifecycle:
##'   \strong{Experimental.} The interface and the results may change in a
##'   later version.
##' @name sScalePsi
##' @aliases sScalePsi
##' @seealso \code{\link{makeSScalePsi}}.
##' @export sScalePsi
NULL

setLoadAction(function(ns)
    assign("sScalePsi", makeSScalePsi(), envir = ns))
