E <- function(fun,
              lower = -Inf,
              upper = Inf,
              ...) {
    int <- integrate(function(x)
        fun(x, ...) * dnorm(x), lower, upper)
    return(int$value)
}

avarSigma <- function(psi) {
    ## Rousseeuw, P. J., Hampel, F. R., Ronchetti, E. M., & Stahel, W. A.
    ## (2011). Robust statistics: the approach based on influence functions.
    ## John Wiley & Sons. Section 2.5e, page 139
    ## chi(x) = w(x) * x ^ 2 - E[w(e) * e ^ 2] = psi(x) * x - EDpsi
    numerator <- E(function(x)
        (psi@psi(x) * x - psi@EDpsi()) ^ 2)
    denominator <-
        E(function(x)
            x * (psi@Dpsi(x) * x + psi@psi(x))) ^ 2
    asymptoticVariance <- numerator / denominator
    return(asymptoticVariance)
}

avarBeta <- function(psi) {
    ## Maronna et al, 2019, Robust Statistics, eq (2.25)
    asymptoticVariance <- psi@Epsi2() / psi@EDpsi() ^ 2
    return(asymptoticVariance)
}

dEta <- function(psi, s) {
    ## Rousseeuw, P. J., Hampel, F. R., Ronchetti, E. M., & Stahel, W. A.
    ## (2011). Robust statistics: the approach based on influence functions.
    ## John Wiley & Sons. 5.3c Paragraph 2 (Page 286)
    wgt <- psi@wgt

    funQ <- function(v)
        (v / s) ^ 2 * wgt(v) ^ 2 * dchisq(v, s)
    numerator <-
        integrate(funQ, 0, Inf)$value / (1 + 2 / s)

    funM <- function(v)
        (v / s) ^ 2 * wgt(v) * dchisq(v, s)
    denominator <- (integrate(funM, 0, Inf)$value / (1 + 2 / s)) ^ 2

    asymptoticVariance <- numerator / denominator
    return(asymptoticVariance)
}

dTau <- function(psi, s) {
    ## Rousseeuw, P. J., Hampel, F. R., Ronchetti, E. M., & Stahel, W. A.
    ## (2011). Robust statistics: the approach based on influence functions.
    ## John Wiley & Sons. 5.3c Paragraph 2 (Page 286)
    ## u^a_\tau(v) = v w^a_\eta(v) - m w^a_\delta(v)
    wgt <- psi@wgt

    funK <- function(v, kappa) {
        arg <- v - s * kappa
        return(arg * wgt(arg) * dchisq(v, s))
    }
    sol <- uniroot(function(kappa)
        integrate(funK, 0, Inf, kappa = kappa)$value,
        c(0, 2))
    kappa <- sol$root

    funQ <- function(v) {
        arg <- v - s * kappa
        return((arg * wgt(arg)) ^ 2 * dchisq(v, s))
    }
    numerator <-
        integrate(funQ, 0, Inf)$value / 2 / s

    funM <- function(v) {
        arg <- v - s * kappa
        return((v - s) * arg * wgt(arg) * dchisq(v, s))
    }
    denominator <- (integrate(funM, 0, Inf)$value / 2 / s) ^ 2

    asymptoticVariance <- numerator / denominator
    return(asymptoticVariance)
}

dMu <- function(psi, s) {
    ## Rousseeuw, P. J., Hampel, F. R., Ronchetti, E. M., & Stahel, W. A.
    ## (2011). Robust statistics: the approach based on influence functions.
    ## John Wiley & Sons. 5.3c Paragraph 2 (Page 286)
    wgt <- psi@wgt

    funQ <- function(v)
        v / s * wgt(v) ^ 2 * dchisq(v, s)
    numerator <- integrate(funQ, 0, Inf)$value

    funM <- function(v)
        v / s * wgt(v) * dchisq(v, s)
    denominator <- (integrate(funM, 0, Inf)$value) ^ 2

    asymptoticVariance <- numerator / denominator
    return(asymptoticVariance)
}

##' \code{asymptoticEfficiency} computes the theoretical asymptotic efficiency
##' for an M-estimator for various types of equations.
##'
##' The asymptotic efficiency is defined as the ratio between the asymptotic
##' variance of the maximum likelihood estimator and the asymptotic variance
##' of the (M-)estimator in question.
##'
##' The computations are only approximate, using numerical integration in the
##' general case. Depending on the regularity of the psi-function, these
##' approximations can be quite crude.
##' @title Compute Asymptotic Efficiencies
##' @param psi object of class psi_func
##' @param equation equation to base computations on. \code{"location"} and
##'   \code{"scale"} are for the univariate case. The others are for a
##'   multivariate location and scale problem. \code{"eta"} is for the shape of
##'   the covariance matrix, \code{"tau"} for the size of the covariance matrix
##'   and \code{"mu"} for the location.
##' @param dimension dimension for the multivariate location and scale problem.
##' @references Maronna, R. A., Martin, R. D., Yohai, V. J., & Salibián-Barrera,
##'   M. (2019). Robust statistics: theory and methods (with R). John Wiley &
##'   Sons., equation (2.25)
##'
##'   Rousseeuw, P. J., Hampel, F. R., Ronchetti, E. M., & Stahel, W. A. (2011).
##'   Robust statistics: the approach based on influence functions. John Wiley &
##'   Sons., Section 5.3c, Paragraph 2 (Page 286)
##' @rdname asymptoticEfficiency
##' @export
asymptoticVariance <-
    function(psi,
             equation = c("location", "scale", "eta", "tau", "mu"),
             dimension = 1) {
        name <- psi@name
        if (name == cPsi@name) {
            return(asymptoticVariance_classic(equation, dimension))
        } else if (name == "Huber, Proposal 2") {
            return(asymptoticVariance_huber_proposal2(psi, equation, dimension))
        } else {
            return(asymptoticVariance_generic(psi, equation, dimension))
        }
    }

asymptoticVariance_generic <-
    function(psi,
             equation = c("location", "scale", "eta", "tau", "mu"),
             dimension = 1) {
        equation <- match.arg(equation)
        checkDimension(dimension, equation)
        asymptoticVariance <- switch(
            equation,
            location = avarBeta(psi),
            scale = avarSigma(psi),
            eta = dEta(psi, dimension),
            tau = dTau(psi, dimension),
            mu = dMu(psi, dimension)
        )
        return(asymptoticVariance)
    }

checkDimension <- function(dimension, equation) {
    checkIsScalar(dimension, "dimension")
    if (equation %in% c("location", "scale")) {
        if (dimension != 1)
            stop(
                "dimension for argument for equation \"",
                equation,
                "\" can only be 1, but was ",
                dimension,
                "."
            )
    } else if (dimension < 2) {
        stop(
            "dimension for argument for equation \"",
            equation,
            "\" needs to be > 1, but was ",
            dimension,
            "."
        )
    }
}

asymptoticVariance_classic <-
    function(equation = c("location", "scale", "eta", "tau", "mu"),
             dimension = 1) {
        equation <- match.arg(equation)
        checkDimension(dimension, equation)
        asymptoticVariance <- switch(
            equation,
            location = 1,
            scale = 0.5,
            eta = 1,
            tau = 1,
            mu = 1
        )
        return(asymptoticVariance)
    }

asymptoticVariance_huber_proposal2 <-
    function(psi,
             equation = c("location", "scale", "eta", "tau", "mu"),
             dimension = 1) {
        equation <- match.arg(equation)
        checkDimension(dimension, equation)
        asymptoticVariance <- switch(
            equation,
            location = avarBeta(psi),
            scale = avarSigma_huber_proposal2(psi),
            eta = dEta(psi, dimension),
            tau = dTau(psi, dimension),
            mu = dMu(psi, dimension)
        )
        return(asymptoticVariance)
    }

##' Computes a partial moment for the standard normal distribution. This is the
##' expectation taken not from -Infinity to Infinity but just to \code{z}.
##' @title Compute Partial Moments
##' @param z partial moment boundary, the expectation is taken from -Inf to z.
##' @param n which moment to compute, needs to be >= 2.
##' @references Winkler, R. L., Roodman, G. M., & Britney, R. R. (1972). The
##' Determination of Partial Moments. Management Science, 19(3), 290–296.
##' http://www.jstor.org/stable/2629511, equation (2.5)
##' @examples
##'   partialMoment_standardNormal(0, 2)
##' @export
partialMoment_standardNormal <- function(z, n) {
    if (n < 2) {
        stop("Can only compute the partial moment for n >= 2")
    }
    product <- Vectorize(function(i)
        prod(n - 2 * (1:i) + 1))
    R_n <- -1 * z ^ (n - 1)
    if (n > 2) {
        lim <- (n - 2 + n %% 2) / 2
        range <- 1:lim
        R_n <- R_n - sum(product(range) * z ^ (n - 1 - 2 * range))
    }
    S_n <- 0
    if (n %% 2 == 0) {
        S_n <- product(n / 2)
    }
    partialMoment <- R_n * dnorm(z) + S_n * pnorm(z)
    return(partialMoment)
}

avarSigma_huber_proposal2 <- function(psi) {
    k <- psi@tDefs[[1]]
    EDpsi <- psi@EDpsi()
    numerator <-
        partialMoment_standardNormal(k, 4) - partialMoment_standardNormal(-k, 4) -
        2 * EDpsi * (partialMoment_standardNormal(k, 2) - partialMoment_standardNormal(-k, 2)) +
        EDpsi ^ 2 * (pnorm(k) - pnorm(-k)) + 2 * (k ^ 2 - EDpsi) ^ 2 * (1 - pnorm(k))
    denominator <-
        (2 * (
            partialMoment_standardNormal(k, 2) - partialMoment_standardNormal(-k, 2)
        )) ^ 2
    asymptoticVariance <- numerator / denominator
    return(asymptoticVariance)
}

##' @rdname asymptoticEfficiency
##' @export
asymptoticEfficiency <-
    function(psi,
             equation = c("location", "scale", "eta", "tau", "mu"),
             dimension = 1) {
        av_classic <- asymptoticVariance_classic(equation, dimension)
        av_psi <- asymptoticVariance(psi, equation, dimension)
        asymptoticEfficiency <- av_classic / av_psi
        return(asymptoticEfficiency)
    }

##' \code{findTuningParameter} finds the value for a tuning parameter \code{k}
##' that for which the desired asymptotic efficiency is achieved.
##' @param desiredEfficiency scalar, specifying the desired asymptotic
##'   efficiency, needs to be between 0 and 1.
##' @param interval interval in which to do the root search, passed on to
##'   \code{\link{uniroot}}.
##' @param ... passed on to \code{\link{uniroot}}.
##' @rdname asymptoticEfficiency
##' @export
findTuningParameter <-
    function(desiredEfficiency,
             psi,
             equation = c("location", "scale", "eta", "tau", "mu"),
             dimension = 1,
             interval = c(0.15, 50),
             ...) {
        checkIsScalar(desiredEfficiency, "desiredEfficiency")
        rootFunction <- function(k) {
            lpsi <- chgDefaults(psi, k = k)
            asymptoticEfficiency <-
                asymptoticEfficiency(lpsi, equation, dimension)
            return(desiredEfficiency - asymptoticEfficiency)
        }
        root <- uniroot(rootFunction, interval, ...)
        return(root$root)
    }
