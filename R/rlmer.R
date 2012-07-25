#######################################################
## Main methods                                      ##
#######################################################

####
## laplace approximated log likelihood
## use the following "algorithm":
##
## 1. find beta, b for given values of theta
## 2. find theta for given values of beta, b
## 3. iterate
####


##' The method is based on the robustification of the scoring
##' equations and an application of the Design Adaptive Scale
##' approach.
##'
##' Available fitting methods for theta and sigma.e:
##' \itemize{
##' \item DASexp:
##'         find theta using DAS-averaged approach, for each canditate
##'         theta, estimate beta and b and estimate sigma_e using
##'         DAS-averaged approach.
##' \item DAStau:
##'         Analogue to the DAS-estimate in robust linear regression.
##' \item Opt:
##'         optimize for (log(sigma), theta), for given theta and sigma_e,
##'         estimate beta and b
##' }
##'
##' Available fitting methods for fixed and random effects:
##' \itemize{
##' \item IRWLS:
##'            Perform iterative reweighted least squares using the
##'            estimating equations (eeq)
##' \item Rcgmin:
##'            minimize the norm part of the logarithm of the density
##'            using the Rcgmin optimizer.
##' }
##' 
##' @title Robust linear mixed models
##' @param formula a two-sided linear formula object describing the
##'   fixed-effects part of the model, with the response on the left of
##'   a \code{~} operator and the terms, separated by \code{+}
##'   operators, on the right.  The vertical bar character \code{"|"}
##'   separates an expression for a model matrix and a grouping factor.
##' @param data an optional data frame containing the variables named
##'   in \code{formula}.  By default the variables are taken from the
##'   environment from which \code{lmer} is called.
##' @param REML logical - Should the estimates be chosen to optimize
##'   the REML criterion (as opposed to the log-likelihood)?  Defaults
##'   to \code{TRUE}.
##' @param ... Additional parameters passed to lmer to find the
##'   initial estimates. See \code{\link[lme4]{lmer}}.
##' @param method method to be used for estimation of theta and sigma,
##'   see Details.
##' @param method.effects method to be used for estimation of the fixed
##'   and random effects.
##' @param rho.e object of class psi_func2, specifying the functions to
##'   use for robustifying the residuals part of the likelihood.
##' @param rho.b object of class psi_func2, specifying the functions to
##'   use for robustifying the random effects part of the likelihood.
##' @param c.sigma.e tuning parameters for rho.e when used for estimating
##'   residual scale.
##' @param c.sigma.b tuning parameters for rho.b when used for estimating
##'   variance components.
##' @param wExp.e numeric specifying the exponent used for the
##'   robustness weights for the residuals. wExp == 0 means that the
##'   standard rho function is used for estimating theta and sigma,
##'   otherwise the weights are used with the specified exponent. (wExp
##'   == 2 corresponds to Proposal II.)
##' @param wExp.b numeric specifying the exponent used for the
##'   robustness weights for the random effects. wExp == 0 means that the
##'   standard rho function is used for estimating theta and sigma,
##'   otherwise the weights are used with the specified exponent. (wExp
##'   == 2 corresponds to Proposal II.)
##' @param use.laplace whether to use the the Laplace approximation for
##'   calculating the deviance.
##' @param rel.tol relative tolerance used as criteria in the fitting
##'   process.
##' @param max.iter maximum number of iterations allowed.
##' @param verbose verbosity of output. Ranges from 0 (none) to 3
##'   (a lot of output)
##' @param use.grad list of two elements, specifying whether to use
##'   the analytical gradients in the optimization.
##' @param doFit logical scalar. When \code{doFit = FALSE} the model
##'   is not fit but instead a structure with the model matrices for the
##'   random-effects terms is returned, so they can be modified for
##'   special model forms. When \code{doFit = TRUE}, the default, the
##'   model is fit immediately.
##' @param clearCache clear the cache before returning the fitted object?
##' @param init optional lmerMod- or rlmerMod-object to use for starting
##'   values, or function producing an lmerMod object.
##' @return object of class rlmerMod.
##' @seealso \code{\link[lme4]{lmer}}
##' @examples
##'   ## dropping of VC
##'   system.time(rlmer(Yield ~ (1|Batch), Dyestuff2))
##'
##'   ## compare classic fits of various methods,
##'   ## without using lmer fit as initial estimate.
##'   ## default method
##'   system.time(rfm.DAStau <- rlmer(Yield ~ (1|Batch), Dyestuff, method="DAStau",
##'                                   rho.e = cPsi, rho.b = cPsi, init=lmerNoFit))
##'   ## DASexp method
##'   system.time(rfm.DASexp <- rlmer(Yield ~ (1|Batch), Dyestuff, method="DASexp",
##'                                   rho.e = cPsi, rho.b = cPsi, init=lmerNoFit))
##'
##'   lf.comp <- function(rfm1, rfm2, tolerance = 1e-5) {
##'       ## fix calls and other slots
##'       rfm2@@call <- rfm1@@call
##'       rfm2@@method <- rfm1@@method
##'       all.equal(rfm1, rfm2, tolerance = tolerance)
##'   }
##'   
##'   stopifnot(lf.comp(rfm.DAStau, rfm.DASexp)
##' 
##' @export
rlmer <- function(formula, data, REML = TRUE, ..., method = "DAStau",
                  method.effects = "IRWLS",
                  rho.e = smoothPsi, rho.b = smoothPsi,
                  c.sigma.e = NULL, c.sigma.b = NULL,
                  wExp.e = 2, wExp.b = 2,
                  use.laplace = REML || method %in% c("Opt"),
                  rel.tol = 1e-8, max.iter = 40*(r+1)^2, verbose = 0,
                  use.grad = c(effects = TRUE, varcomp = FALSE),
                  doFit = TRUE, clearCache = TRUE, init)
{
    lcall <- match.call()
    if (missing(init)) {
        init <- lmer(formula, data, REML, ...)
    } else if (is.function(init)) {
        init <- do.call(init,list(formula=formula, data=data, REML=REML, ...))
    }
    lobj <- as(init, "rlmerMod")
    lobj@call <- lcall
    ## clearCache before starting
    lobj@cache <- new.env()

    ## give a warning if weights or offset are used
    if (any(lobj@resp$weights != 1))
        warning("Argument weights is untested.")
    if (any(lobj@resp$offset != 0))
        warning("Argument offset is untested.")

    ## change default method for non diagonal U_b
    ## only DASexp can deal with non diagonal U_b at the moment.
    if (missing(method) && !isDiagonal(lobj@pp$U_b)) method <- "DASexp"

    ## DASexp and Opt methods can only deal with diagonal V_b(theta) at the moment
    if (!(method %in% c("Opt", "DASexp")) && !isDiagonal(lobj@pp$U_b))
        stop("Method ", method, " can only deal with diagonal V_b at the moment")
    
    ## if wExp.b > 0 warn if sigma is not estimated via DAS
    if (wExp.b > 0 && !method %in% c("DASexp", "DAStau")) {
        if (!missing(wExp.b))
            warning("wExp.b > 0 is only suitable for method DASexp and DAStau")
        wExp.b <- 0
    }
    ## if wExp.e > 0 warn if method is not DAS or DASe
    if (wExp.e > 0 && !method %in% c("DASexp", "DAStau")) {
        if (!missing(wExp.e))
            warning("wExp.e > 0 is only suitable for method DAS*")
        wExp.e <- 0
    }
    ## set arguments only relevant to rlmerMod
    lobj@rho.b <- rho.b
    lobj@rho.sigma.b <- if (is.null(c.sigma.b)) rho.b else {
        if (method %in% c("Opt")) {
            warning("method Opt does not support different c.sigma.b, ignoring argument")
            rho.e
        } else {
            tmp <- rho.b@tDefs
            if (is.null(names(c.sigma.b))) {
                tmp[1:length(c.sigma.b)] <- c.sigma.b
            } else if (all(nchar(names(c.sigma.b)) > 0)) {
                tmp[names(c.sigma.b)] <- c.sigma.b
            } else stop("c.sigma.b must either not contain any names or be fully named.")
            do.call(chgDefaults, c(list(rho.b), tmp))
        }
    }
    lobj@wExp.b <- wExp.b
    lobj@rho.e <- rho.e
    lobj@rho.sigma.e <- if (is.null(c.sigma.e)) rho.e else {
        if (method %in% c("Opt")) {
            warning("method Opt does not support different c.sigma.b, ignoring argument")
            rho.e
        } else {
            tmp <- rho.e@tDefs
            if (is.null(names(c.sigma.e))) {
                tmp[1:length(c.sigma.e)] <- c.sigma.e
            } else if (all(nchar(names(c.sigma.e)) > 0)) {
                tmp[names(c.sigma.e)] <- c.sigma.e
            } else stop("c.sigma.e must either not contain any names or be fully named.")
            do.call(chgDefaults, c(list(rho.e), tmp))
        }
    }
    lobj@wExp.e <- wExp.e
    lobj@use.laplace <- use.laplace
    lobj@method <- method
    lobj@method.effects <- method.effects
    if (substr(method, 1, 3) == "DAS") {
        lobj@pp <- as(lobj@pp, "rlmerPredD_DAS")
    }
    lobj@pp$initRho(lobj)
    lobj@pp$initMatrices(lobj)
    lobj@pp$updateMatrices()
   
    
    if (!doFit) return(lobj)
    
    ## do not start with theta == 0
    if (any(theta(lobj)[lobj@lower == 0] == 0)) {
        if (verbose > 0)
            cat("Setting variance components from 0 to 1\n")
        theta0 <- theta(lobj)
        theta0[lobj@lower == 0 & theta0 == 0] <- 1
        setTheta(lobj, theta0, fit.effects = TRUE, update.sigma = method != "Opt")
    } else {
        ## set theta at least once
        setTheta(lobj, theta(lobj), fit.effects = FALSE)
    }
    
    ## recalculate deviance
    updateDeviance(lobj)
    
    if (verbose > 0) {
        cat("\nrlmer starting values:\n")
        cat("deviance: ", deviance(lobj) , "\n")
        cat("sigma, theta: ", lobj@pp$sigma, ", ", theta(lobj), "\n")
        cat("coef: ", lobj@pp$beta, "\n")
    }
    if (verbose > 1)
        cat("b.s: ", b.s(lobj), "\n")
    
    if (use.grad['varcomp'])
        stop("analytic gradient for varcomp not implemented")
    if (!use.grad['effects'])
        stop("no gradient for effects not implemented")

    ## required for max.iter:
    r <- len(lobj, "theta")
    
    ## do fit
    lobj <- switch(method,
                   `DAStau` = rlmer.fit.DAS(lobj, verbose, max.iter, rel.tol),
                   `DASexp`= rlmer.fit.DASexp(lobj, verbose, max.iter, rel.tol),
                   `Opt`= rlmer.fit.Opt(lobj, verbose, max.iter),
                   stop("unknown fitting method"))
    
    if (verbose > 0) {
        cat("deviance: ", deviance(lobj), "\n")
        cat("sigma, theta: ", lobj@pp$sigma, ", ", theta(lobj), "\n")
        cat("coef: ", lobj@pp$beta, "\n")
    }
    if (verbose > 1) 
        cat("b.s: ", b.s(lobj), "\n")
    
    ## clear the cache
    if (clearCache) lobj@cache <- new.env()
    
    return(updateWeights(lobj))
}

## DAS method
rlmer.fit.DAS <- function(lobj, verbose, max.iter, rel.tol) {
    if (!.isREML(lobj))
        stop("can only do REML when using averaged DAS-estimate for sigma")

    lupdateTheta <- switch(lobj@method,
                           DAStau=updateThetaTau,
                           stop("method not supported by rlmer.fit.DAS:", lobj@method))
    
    ## fit
    converged <- FALSE
    theta0 <- theta(lobj)
    if (verbose > 1) {
        coef0 <- lobj@pp$beta
        b.s0 <- b.s(lobj)
        sigma0 <- lobj@pp$sigma
    }
    iter <- 0
    while (!converged && iter < max.iter) {
        iter <- iter + 1
        if (verbose > 0) cat("Iteration", iter, "\n")
        ## fit theta
        lupdateTheta(lobj, max.iter, rel.tol/10, verbose)
        theta1 <- theta(lobj)

        if (verbose > 0) {
            cat(sprintf("delta theta: %.12f\n", sum(abs(theta0 - theta1))))
            if (verbose > 1) {
                cat(sprintf("delta coef:  %.12f\n", sum(abs(coef0 - lobj@pp$beta))))
                cat(sprintf("delta u:     %.12f\n", sum(abs(b.s0 - b.s(lobj)))))
                cat(sprintf("delta sigma: %.12f\n", abs(sigma0 - lobj@pp$sigma)))
                coef0 <- lobj@pp$beta
                b.s0 <- b.s(lobj)
                sigma0 <- lobj@pp$sigma
                if (verbose > 2) {
                    cat("theta:  ", theta(lobj),"\n")
                    cat("coef:   ", lobj@pp$beta,"\n")
                    cat("b.s:    ", b.s(lobj), "\n")
                    cat("sigmae: ", lobj@pp$sigma, "\n")
                }
            }
        }

        ## all zero or change smaller than relative tolerance
        ## all zero: we can't get out of this anyway, so we have to stop.
        converged <- all(theta1 == 0) || sum(abs(theta0 - theta1)) < 200*rel.tol*sum(abs(theta0)) 
        if (verbose > 1)
            cat(sprintf("Criterion: %.12f, %.12f", sum(abs(theta0 - theta1)),
                sum(abs(theta0 - theta1)) / rel.tol / sum(abs(theta0)) / 200), "\n")
        theta0 <- theta1
    }

    ## update deviance
    updateDeviance(lobj)

    if (iter == max.iter) 
        warning("iterations did not converge, returning unconverged estimate.")
        
    lobj
}

## DAS method
rlmer.fit.DASexp <- function(lobj, verbose, max.iter, rel.tol) {
    if (!.isREML(lobj))
        stop("can only do REML when using averaged DAS-estimate for sigma")

    ## test if we are in the diagonal Lambda case (wExp.b == 0 allowed)
    if (lobj@wExp.b == 0 && !all(lobj@dim == 1)) 
        stop("For wExp.b = 0 only diagonal Lambda allowed")
    
    .jacfunc <- function(time,y,parms, ...)
        zero.theta.DASexp(y, ..., update.sigma=FALSE)
    jacfunc <- function(x, ...)
        jacobian.full(x, .jacfunc, pert=1e-3, ...)
     
    ## find root
    res <- multiroot(zero.theta.DASexp, theta(lobj), object = lobj, maxiter = max.iter,
                     rtol = rel.tol, ctol = 100*rel.tol, atol = 100*rel.tol,
                     lower = lobj@lower, method="DAS", jacfunc = jacfunc,
                     verbose = verbose > 2, VERBOSE=verbose)
    if (verbose > 0) str(res)
    
    setTheta(lobj, res$root, fit.effects = TRUE)
        
    lobj
}

## Opt method
rlmer.fit.Opt <- function(lobj, verbose, max.iter) {
    ## sigma on log scale, other thetas relative to (real) sigma
    opt.theta <- function(par, object) {
        setSigma(object, exp(par[1]))
        setTheta(object, par[-1], update.sigma = FALSE, fit.effects = TRUE)
        dev <- .deviance(object)
        if (is.infinite(dev)) {
            dev <- sign(dev) * .Machine$double.xmax
            setDeviance(object, dev) 
            warning("deviance infinite, setting dev = ", dev)
        }
        if (verbose > 2)
            cat("opt.theta(", par, ") =", dev, "\n")
        dev
    }
    
    ## fit
    opt <- bobyqa(c(log(lobj@pp$sigma), theta(lobj)), opt.theta, object = lobj,
                  lower = c(log(1e-7), lobj@lower), 
                  control = list(iprint = max(verbose - 2, 0),
                  maxfun = max.iter))
    if (opt$ierr > 0)
        warning("bobyqa did not converge (exit code ", opt$ierr,")")
    dev1 <- opt$fval
    setSigma(lobj, exp(opt$par[1]))
    setTheta(lobj, opt$par[-1], update.sigma = FALSE, fit.effects = TRUE)
    if (verbose > 0)
        cat("\nbobyqa required", opt$feval, "iterations\n")

    lobj
}
