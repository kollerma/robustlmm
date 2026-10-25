## Shipped with robustlmm: a helper of the breakdown study behind
## setting = "RSEhb" in vignette("simulationStudies"), sourced by
## bd23_map_sim.R. The S-type scale that robustlmm ships is sScalePsi().
##
## Prototype: an S-type (bounded-rho) scale for rlmer's scale equations.
##
## rlmer solves sum_i eta_i tau_i^2 W(t_i) (t_i^2 - kappa) = 0 for sigma,
## t_i = r_i / (tau_i sigma), with kappa = E[W(Z) Z^2] / E[W(Z)]. With
##   W(t) = (rho(t) - b) / (t^2 - kappa0),  rho(sqrt(kappa0)) = b = E rho(Z),
## (removable singularity at t^2 = kappa0, where W = rho'(t) / (2 t)) the
## fixed point is sum eta tau^2 (rho(t) - b) = 0: the M-scale with bounded
## rho, breakdown point min(b, 1 - b) (b = 1/2: 50%), and kappa = kappa0.
## rho is the bisquare Mchi normalised to sup 1; c = 1.547645 gives b = 1/2.

sScalePsi_proto <- function(cc = 1.547645) {
    c0 <- cc
    rho <- function(x) robustbase::Mchi(x, c0, "bisquare")
    b <- stats::integrate(function(z) rho(z) * stats::dnorm(z), -Inf, Inf,
                          rel.tol = 1e-12)$value
    t0 <- stats::uniroot(function(t) rho(t) - b, c(1e-6, c0), tol = 1e-14)$root
    kappa0 <- t0^2
    drho <- function(x) robustbase::Mpsi(x, c0, "bisquare", 0L)   # rho' (Mchi scale)
    W <- function(x, cc) {
        x <- abs(x)
        den <- x^2 - kappa0
        out <- (rho(x) - b) / den
        near <- abs(den) < 1e-6
        if (any(near)) {
            ## L'Hopital: rho'(x) / (2 x); robustbase Mpsi is the derivative of
            ## Mchi up to the normalising constant MrhoInf
            out[near] <- drho(x[near]) / robustbase::MrhoInf(c0, "bisquare") /
                (2 * x[near])
        }
        out
    }
    psi <- function(x, cc) W(x) * x
    Dwgt <- function(x, cc) {
        h <- 1e-5
        (W(x + h) - W(x - h)) / (2 * h)
    }
    Dpsi <- function(x, cc) W(x) + x * Dwgt(x)
    rho_fn <- function(x, cc) rho(x)
    U <- stats::qnorm(1 - 1e-10)
    ## E-slots vectorise over cc (robustbase's functionXal1 validity check);
    ## the value does not depend on cc, the shape is fixed at construction.
    mkE <- function(fn) {
        val <- 2 * stats::integrate(function(z) fn(z) * stats::dnorm(z), 0, U,
                                    rel.tol = 1e-8)$value
        function(cc) rep(val, length(cc))
    }
    pf <- robustbase::psiFunc(rho = rho_fn, psi = psi, wgt = W, Dpsi = Dpsi,
                              Dwgt = Dwgt,
                              Erho = mkE(rho), Epsi2 = mkE(function(x) psi(x)^2),
                              EDpsi = mkE(Dpsi),
                              name = sprintf("S-scale (bisquare, c = %.4g, b = %.3g)", c0, b),
                              cc = c0)
    args <- list("psi_func_rcpp")
    for (sn in methods::slotNames("psi_func")) args[[sn]] <- methods::slot(pf, sn)
    args[["getRcppClass"]] <- function() character(0)
    args[["getInstanceWithOriginalDefaults"]] <- function() NULL
    obj <- do.call("new", args)
    attr(obj, "sscale") <- list(b = b, kappa0 = kappa0, cc = c0)
    obj
}
