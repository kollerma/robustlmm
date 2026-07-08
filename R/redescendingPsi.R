## Redescending psi-functions for rlmer, delegating the psi / rho / weight
## / derivative evaluations to robustbase's compiled Mpsi() family so that
## the shipped functions match lmrob() exactly and additional redescenders
## (lqq, optimal, hampel, ggw) come for free.

## Internal engine shared by makeRobustbasePsi() and makeBisquarePsi().
## `family` is a robustbase Mpsi() family name, `cc` its tuning constant
## (a scalar for bisquare/optimal, a length-3/4 specification for
## hampel/ggw/lqq), `name` the display name. The rho slot is normalised as
## Mchi(x, cc, family) * MrhoInf(cc, family) so that rho(Inf) equals the
## sup of rho -- the psi_func rho convention robustlmm uses. The Fisher
## consistency expectations (Erho, Epsi2, EDpsi) are computed by numerical
## integration against the standard normal, matching the psi_func
## conventions of the other psi-functions in the package.
.buildRobustbasePsi <- function(family, cc, name) {
    ## psi / rho / wgt / Dpsi / Dwgt as functions of (x, cc). The scalar
    ## cc = Inf (unbounded / classical limit) and cc = 0 (all-rejecting
    ## limit) branches are only reached by robustbase's functionXal1
    ## validity check, which evaluates the E-slots at cc = 0 and cc = Inf;
    ## handling them explicitly avoids Inf * 0 = NaN in the normalised rho.
    psi_fn  <- function(x, cc) {
        if (length(cc) == 1L && !is.finite(cc)) return(x)
        if (length(cc) == 1L && cc <= 0)        return(rep(0, length(x)))
        robustbase::Mpsi(x, cc, family, 0L)
    }
    wgt_fn  <- function(x, cc) {
        if (length(cc) == 1L && !is.finite(cc)) return(rep(1, length(x)))
        if (length(cc) == 1L && cc <= 0)        return(rep(0, length(x)))
        robustbase::Mwgt(x, cc, family)
    }
    Dpsi_fn <- function(x, cc) {
        if (length(cc) == 1L && !is.finite(cc)) return(rep(1, length(x)))
        if (length(cc) == 1L && cc <= 0)        return(rep(0, length(x)))
        robustbase::Mpsi(x, cc, family, 1L)
    }
    rho_fn  <- function(x, cc) {
        if (length(cc) == 1L && !is.finite(cc)) return(x^2 / 2)
        if (length(cc) == 1L && cc <= 0)        return(rep(0, length(x)))
        robustbase::Mchi(x, cc, family) * robustbase::MrhoInf(cc, family)
    }
    Dwgt_fn <- function(x, cc) {
        if (length(cc) == 1L && (!is.finite(cc) || cc <= 0))
            return(rep(0, length(x)))
        p  <- robustbase::Mpsi(x, cc, family, 0L)
        dp <- robustbase::Mpsi(x, cc, family, 1L)
        ## d/dx (psi/x) = (x psi' - psi)/x^2, with the finite limit 0 at 0.
        ifelse(x == 0, 0, (dp * x - p) / x^2)
    }

    ## E-slots. When cc is a single tuning constant (bisquare, optimal)
    ## the E-functions must vectorise over cc for robustbase's
    ## functionXal1 validity check; when cc is a length-3/4 specification
    ## (hampel, ggw, lqq) it is a single tuning spec and the E-functions
    ## are plain functionXal (no cc-vectorisation / validity check).
    scalar_tuning <- length(cc) == 1L
    U <- stats::qnorm(1 - 1e-10)
    mkE <- function(fn) {
        one <- function(cci) 2 * stats::integrate(function(z)
            fn(z, cci) * stats::dnorm(z), 0, U, rel.tol = 1e-8)$value
        if (scalar_tuning) function(cc) vapply(cc, one, numeric(1))
        else               function(cc) one(cc)
    }
    Erho  <- mkE(rho_fn)
    Epsi2 <- mkE(function(x, cc) psi_fn(x, cc)^2)
    EDpsi <- mkE(Dpsi_fn)

    pf <- psiFunc(rho = rho_fn, psi = psi_fn, wgt = wgt_fn,
                  Dpsi = Dpsi_fn, Dwgt = Dwgt_fn,
                  Erho = Erho, Epsi2 = Epsi2, EDpsi = EDpsi,
                  name = name, cc = cc)

    args <- list("psi_func_rcpp")
    for (sn in methods::slotNames("psi_func")) args[[sn]] <- methods::slot(pf, sn)
    args[["getRcppClass"]] <- function() character(0)
    args[["getInstanceWithOriginalDefaults"]] <- function() NULL
    do.call("new", args)
}

.robustbasePsiName <- function(family)
    switch(family,
           bisquare = "Tukey bisquare",
           lqq      = "lqq",
           optimal  = "optimal",
           hampel   = "Hampel",
           ggw      = "Gauss weight (ggw)",
           family)

##' Construct a redescending psi-function from one of robustbase's
##' compiled \code{\link[robustbase]{Mpsi}} families.
##'
##' Returns a \code{psi_func_rcpp} object whose psi, rho, weight and
##' derivative slots delegate to robustbase's compiled psi-function
##' family, so the shipped functions match \code{\link[robustbase]{lmrob}}
##' exactly. The psi, its derivative and the weight are
##' \code{\link[robustbase]{Mpsi}(x, cc, family, 0)},
##' \code{Mpsi(x, cc, family, 1)} and
##' \code{\link[robustbase]{Mwgt}(x, cc, family)}; the rho slot is the
##' normalised \code{\link[robustbase]{Mchi}(x, cc, family) *
##' \link[robustbase]{MrhoInf}(cc, family)}, so that \code{rho(Inf)}
##' equals the supremum of rho (the \code{psi_func} rho convention used by
##' \code{\link{rlmer}}). The Fisher consistency expectations
##' (\code{Erho}, \code{Epsi2}, \code{EDpsi}) are computed by numerical
##' integration against the standard normal.
##'
##' The available families are the redescenders \code{"bisquare"} (Tukey
##' biweight), \code{"lqq"} (linear-quadratic-quadratic),
##' \code{"optimal"}, \code{"hampel"} and \code{"ggw"} (generalised
##' Gauss-weight). When \code{cc} is \code{NULL} the family's default
##' tuning (\eqn{\approx 95\%} efficiency at the normal) is taken from
##' \code{robustbase::.Mpsi.tuning.default(family)}; for \code{"lqq"},
##' \code{"ggw"} and \code{"hampel"} this is a short numeric
##' \emph{specification} vector rather than a single cutoff (see
##' \code{\link[robustbase]{lmrob.control}}).
##'
##' Of the redescenders, \code{\link{bisquarePsi}} redescends comparatively
##' fast; the \code{"lqq"} psi of Koller and Stahel (2011), used by
##' robustbase's \code{lmrob.control(setting = "KS2014")}, is the
##' recommended redescender and is pre-built as \code{\link{lqqPsi}}.
##'
##' @title Redescending psi-functions from robustbase (bisquare, lqq,
##'   optimal, hampel, ggw)
##' @param family robustbase \code{\link[robustbase]{Mpsi}} family name,
##'   one of \code{"bisquare"}, \code{"lqq"}, \code{"optimal"},
##'   \code{"hampel"}, \code{"ggw"}.
##' @param cc tuning constant (a scalar for \code{"bisquare"} /
##'   \code{"optimal"}, a specification vector for \code{"lqq"} /
##'   \code{"ggw"} / \code{"hampel"}). \code{NULL} uses the family's
##'   robustbase default (about 95\% efficiency at the normal).
##' @return psi_func_rcpp object usable as \code{rho.e} (or \code{rho.b},
##'   \code{rho.sigma.e}, \code{rho.sigma.b}) in \code{\link{rlmer}}.
##' @references
##' Koller, M. and Stahel, W. A. (2011) Sharpening Wald-type inference in
##' robust regression for small samples. \emph{Computational Statistics &
##' Data Analysis} \bold{55}(8), 2504--2515.
##' @seealso \code{\link{lqqPsi}}, \code{\link{bisquarePsi}},
##'   \code{\link[robustbase]{Mpsi}}, \code{\link[robustbase]{lmrob.control}}.
##' @examples
##'   pf <- makeRobustbasePsi("lqq")
##'   pf@@psi(c(-6, -3, 0, 3, 6))
##' @export
makeRobustbasePsi <- function(family = c("bisquare", "lqq", "optimal",
                                         "hampel", "ggw"),
                              cc = NULL) {
    family <- match.arg(family)
    if (is.null(cc)) cc <- robustbase::.Mpsi.tuning.default(family)
    if (!is.numeric(cc) || !length(cc))
        stop("'cc' must be a non-empty numeric tuning vector")
    .buildRobustbasePsi(family, cc,
                        sprintf("%s (robustbase)", .robustbasePsiName(family)))
}

##' Construct a redescending Tukey bisquare psi function.
##'
##' Returns a \code{psi_func_rcpp} object whose psi, rho, weight and
##' derivative slots implement Tukey's bisquare (biweight) function
##' \deqn{\psi(x) = x (1 - (x/c)^2)^2 \quad \text{for } |x| \le c,\quad 0 \text{ otherwise}}
##' suitable for use as \code{rho.e} (or \code{rho.b}, \code{rho.sigma.e},
##' \code{rho.sigma.b}) in \code{\link{rlmer}}.
##'
##' Since robustlmm 3.5.0 the psi, rho, weight and derivative evaluations
##' delegate to robustbase's compiled bisquare family
##' (\code{\link[robustbase]{Mpsi}}, \code{\link[robustbase]{Mwgt}},
##' \code{\link[robustbase]{Mchi}}), so they match
##' \code{\link[robustbase]{lmrob}} exactly; the returned values are
##' identical (to numerical tolerance) to the previous hand-coded
##' implementation. \code{makeBisquarePsi} is the scalar-cutoff special
##' case of the general \code{\link{makeRobustbasePsi}}.
##'
##' The Fisher consistency expectations (\code{Erho}, \code{Epsi2},
##' \code{EDpsi}) are computed by numerical integration against the
##' standard normal.
##'
##' The bisquare redescends comparatively fast; for a redescending fit the
##' \code{"lqq"} psi (\code{\link{lqqPsi}}), the recommended redescender of
##' robustbase's \code{lmrob.control(setting = "KS2014")}, is generally
##' preferable (Koller and Stahel 2011).
##'
##' @title Bisquare (Tukey biweight) psi function
##' @param c tuning cutoff. The default \code{c = 4.685} gives about 95\%
##'   asymptotic efficiency for the location problem.
##' @return psi_func_rcpp object usable in \code{rho.e}, \code{rho.b}.
##' @references
##' Koller, M. and Stahel, W. A. (2011) Sharpening Wald-type inference in
##' robust regression for small samples. \emph{Computational Statistics &
##' Data Analysis} \bold{55}(8), 2504--2515.
##' @seealso \code{\link{makeRobustbasePsi}}, \code{\link{lqqPsi}}.
##' @examples
##'   pf <- bisquarePsi
##'   pf@@psi(c(-6, -3, 0, 3, 6))
##' @export
makeBisquarePsi <- function(c = 4.685) {
    if (!is.numeric(c) || length(c) != 1L || !is.finite(c) || c <= 0)
        stop("'c' must be a finite positive scalar")
    .buildRobustbasePsi("bisquare", c,
                        sprintf("Tukey bisquare (c = %.4g)", c))
}

##' Default Tukey bisquare psi-function with \code{c = 4.685}.
##'
##' Pre-built \code{psi_func_rcpp} from \code{\link{makeBisquarePsi}} with
##' the conventional tuning constant \code{c = 4.685} (\eqn{\approx 95\%}
##' efficiency at the normal model). The bisquare redescends comparatively
##' fast; \code{\link{lqqPsi}} is the recommended redescender.
##' @name bisquarePsi
##' @aliases bisquarePsi
##' @seealso \code{\link{lqqPsi}}, \code{\link{makeRobustbasePsi}}.
##' @export bisquarePsi
NULL

##' Default lqq (linear-quadratic-quadratic) redescending psi-function.
##'
##' Pre-built \code{psi_func_rcpp} from
##' \code{\link{makeRobustbasePsi}("lqq")} with robustbase's default lqq
##' tuning (about 95\% efficiency at the normal). The lqq psi of Koller and
##' Stahel (2011) redescends more gradually than the bisquare and is the
##' recommended redescender, as used by robustbase's
##' \code{\link[robustbase]{lmrob.control}(setting = "KS2014")}. Use it as
##' \code{rho.e} in \code{\link{rlmer}}; a redescender needs a good initial
##' estimate, so pair it with \code{init = "ransac"} or
##' \code{\link{rlmer_ransac}}.
##' @name lqqPsi
##' @aliases lqqPsi
##' @references
##' Koller, M. and Stahel, W. A. (2011) Sharpening Wald-type inference in
##' robust regression for small samples. \emph{Computational Statistics &
##' Data Analysis} \bold{55}(8), 2504--2515.
##' @seealso \code{\link{makeRobustbasePsi}}, \code{\link{bisquarePsi}}.
##' @export lqqPsi
NULL

setLoadAction(function(ns)
    assign("bisquarePsi", makeBisquarePsi(c = 4.685), envir = ns))

setLoadAction(function(ns)
    assign("lqqPsi", makeRobustbasePsi("lqq"), envir = ns))
