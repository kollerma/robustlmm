##' Robust estimation of linear mixed effects models, for hierarchical nested
##' and non-nested, e.g., crossed, datasets.
##'
##' \describe{ \item{Overview:}{
##'
##' This function implements the Robust Scoring Equations estimator for linear
##' mixed effect models. It can be used much like the function
##' \code{\link[lme4]{lmer}} in the package \code{lme4}. The supported models
##' are the same as for \code{\link[lme4]{lmer}} (gaussian family only). The
##' robust approach used is based on the robustification of the scoring
##' equations and an application of the Design Adaptive Scale approach.
##'
##' Example analyses and theoretical details on the method are available in the
##' vignette (see \code{vignette("rlmer")}).
##'
##' Models are specified using the \code{formula} argument, using the same
##' syntax as for \code{\link[lme4]{lmer}}. Additionally, one also needs to
##' specify what robust scoring or weight functions are to be used (arguments
##' starting with \code{rho.}). By default a smoothed version of the Huber
##' function is used. Furthermore, the \code{method} argument can be used to
##' speed up computations at the expense of accuracy of the results. }
##'
##' \item{Computation methods:}{
##'
##' Currently, there are two different methods available for fitting models.
##' They only differ in how the consistency factors for the Design Adaptive
##' Scale estimates are computed. Available fitting methods for theta and
##' sigma.e: \itemize{
##'
##' \item \code{DAStau} (default): For this method, the consistency factors are
##' computed using numerical quadrature. This is slower but yields more accurate
##' results. This is the direct analogue to the DAS-estimate in robust linear
##' regression.
##'
##' \item \code{DASvar}: This method computes the consistency factors using a
##' direct approximation which is faster but less accurate. For complex models
##' with correlated random effects with more than one correlation term, this is
##' the only method available.
##'
##' } }
##'
##' \item{Weight functions:}{
##'
##' The tuning parameters of the weight functions \dQuote{rho} can be used to
##' adjust robustness and efficiency of the resulting estimates (arguments
##' \code{rho.e}, \code{rho.b}, \code{rho.sigma.e} and \code{rho.sigma.b}).
##' Better robustness will lead to a decrease of the efficiency. With the default
##' setting, \code{setting = "RSEn"}, the tuning parameters are set to yield
##' estimates with approximately 95\% efficiency for the fixed effects. The
##' variance components are estimated with a lower efficiency but better
##' robustness properties.
##'
##' One has to use different weight functions and tuning parameters for simple
##' variance components and for such including correlation parameters. By
##' default, they are chosen appropriately to the model at hand. However, when
##' using the \code{rho.sigma.e} and \code{rho.sigma.b} arguments, it is up to
##' the user to specify the appropriate function. See
##' \code{\link{asymptoticEfficiency}} for methods to find tuning parameters
##' that yield a given asymptotic efficiency. \itemize{
##'
##' \item For simple variance components and the residual error scale use the
##' function \code{\link{psi2propII}} to change the tuning parameters. This is
##' similar to Proposal 2 in the location-scale problem (i.e., using the
##' squared robustness weights of the location estimate for the scale estimate;
##' otherwise the scale estimate is not robust).
##'
##' \item For multi-dimensional blocks of random effects modeled, e.g.,
##' a model with correlated random intercept and slope, (referred to as
##' block diagonal case below), use the \code{\link{chgDefaults}} function to
##' change the tuning parameters. The parameter estimation problem is
##' multivariate, unlike the case without correlation where the problem was
##' univariate. For the employed estimator, this amounts to switching from
##' simple scale estimates to estimating correlation matrices. Therefore
##' different weight functions have to be used. Squaring of the weights (using
##' the function \code{\link{psi2propII}}) is no longer necessary. To yield
##' estimates with the same efficiency, the tuning parameters for the
##' block diagonal are larger than for the simple case. Tables of tuning parameters
##' are given in Table 2 and 3 of the vignette (\code{vignette("rlmer")}).
##'
##' } }
##'
##' \item{Recommended tuning parameters:}{
##'
##' For a more robust estimate, use \code{setting = "RSEn"} (the default). For
##' higher efficiency, use \code{setting = "RSEa"}. The settings described in
##' the following paragraph are used when \code{setting = "RSEa"} is specified.
##'
##' For the smoothed Huber function the tuning parameters to get approximately
##' 95\% efficiency are \eqn{k=1.345}{k=1.345} for \code{rho.e} and
##' \eqn{k=2.28}{k=2.28} for \code{rho.sigma.e} (using the squared version). For
##' simple variance components, the same can be used for \code{rho.b} and
##' \code{rho.sigma.b}. For variance components including correlation
##' parameters, use \eqn{k=5.14}{k=5.14} for both \code{rho.b} and
##' \code{rho.sigma.b}. Tables of tuning parameter are given in Table 2 and 3 of
##' the vignette (\code{vignette("rlmer")}). }
##'
##' \item{Specifying (multiple) weight functions:}{
##'
##' If custom weight functions are specified using the argument \code{rho.b}
##' (\code{rho.e}) but the argument \code{rho.sigma.b} (\code{rho.sigma.e}) is
##' missing, then the squared weights are used for simple variance components
##' and the regular weights are used for variance components including
##' correlation parameters. The same tuning parameters will be used when
##' \code{setting = "RSEn"} is used. To get
##' higher efficiency either use \code{setting = "RSEa"} (and only set arguments
##' \code{rho.e} and \code{rho.b}). Or specify the tuning parameters by hand
##' using the \code{\link{psi2propII}} and \code{\link{chgDefaults}} functions.
##'
##' To specify separate weight functions \code{rho.b} and \code{rho.sigma.b} for
##' different variance components, it is possible to pass a list instead of a
##' psi_func object. The list entries correspond to the groups as shown by
##' \code{VarCorr(.)} when applied to the model fitted with \code{lmer}. A set
##' of correlated random effects count as just one group. }
##'
##' \item{\code{lmerNoFit}:}{
##'
##' The \code{lmerNoFit} function can be used to get trivial starting values.
##' This is mainly used to verify the algorithms to reproduce the fit by
##' \code{\link{lmer}} when starting from trivial initial values. } }
##'
##' @title Robust Scoring Equations Estimator for Linear Mixed Models
##' @param formula a two-sided linear formula object describing the
##'   fixed-effects part of the model, with the response on the left of a
##'   \code{~} operator and the terms, separated by \code{+} operators, on the
##'   right.  The vertical bar character \code{"|"} separates an expression for
##'   a model matrix and a grouping factor.
##' @param data an optional data frame containing the variables named in
##'   \code{formula}.  By default the variables are taken from the environment
##'   from which \code{lmer} is called.
##' @param ... Additional parameters passed to lmer to find the initial
##'   estimates. See \code{\link[lme4]{lmer}}.
##' @param method method to be used for estimation of theta and sigma, see
##'   Details.
##' @param setting a string specifying suggested choices for the arguments
##'   \code{rho.e}, \code{rho.sigma.e}, \code{rho.b} and \code{rho.sigma.b}.
##'   Use \code{"RSEn"} (the default) or \code{"RSEa"}. Both use
##'   \code{\link{smoothPsi}} for all the \dQuote{rho} arguments. For
##'   \code{rho.sigma.e}, squared robustness weights are used (see
##'   \code{\link{psi2propII}}). \code{"RSEn"} uses the same tuning parameter as
##'   for \code{rho.e}, which leads to higher robustness but lower efficiency.
##'   \code{"RSEa"} adjusts the tuning parameter for higher asymptotic efficiency
##'   which results in lower robustness (\code{k = 2.28} for default \code{rho.e}).
##'   For diagonal random effects covariance matrices, \code{rho.sigma.b} is
##'   treated exactly as \code{rho.sigma.e}. For block diagonal random effects
##'   covariance matrices (with correlation terms), regular robustness weights
##'   are used for \code{rho.sigma.b}, not squared ones, as they're not needed.
##'   But the tuning parameters are adjusted for both \code{rho.b} and
##'   \code{rho.sigma.b} according to the dimensions of the blocks (for both
##'   \code{"RSEn"} or \code{"RSEa"}). For a block of dimension 2 (e.g.,
##'   correlated random intercept and slope) \code{k = 5.14} is used.
##' @param rho.e object of class psi_func, specifying the functions to use for
##'   the huberization of the residuals.
##' @param rho.b object of class psi_func or list of such objects (see Details),
##'   specifying the functions to use for the huberization of the random
##'   effects.
##' @param rho.sigma.e object of class psi_func, specifying the weight functions
##'   to use for the huberization of the residuals when estimating the variance
##'   components, use the \code{\link{psi2propII}} function to specify squared
##'   weights and custom tuning parameters.
##' @param rho.sigma.b (optional) object of class psi_func or list of such
##'   objects, specifying the weight functions to use for the huberization of
##'   the random effects when estimating the variance components (see Details).
##'   Use \code{\link{psi2propII}} to specify squared weights and custom tuning
##'   parameters or \code{\link{chgDefaults}} for regular weights for variance
##'   components including correlation parameters.
##' @param rel.tol relative tolerance used as criteria in the fitting process.
##' @param max.iter maximum number of iterations allowed.
##' @param verbose verbosity of output. Ranges from 0 (none) to 3 (a lot of
##'   output)
##' @param doFit logical scalar. When \code{doFit = FALSE} the model is not fit
##'   but instead a structure with the model matrices for the random-effects
##'   terms is returned (used to speed up tests). When \code{doFit = TRUE}, the
##'   default, the model is fit immediately.
##' @param init optional lmerMod- or rlmerMod-object to use for starting values,
##'   a list with elements \sQuote{fixef}, \sQuote{u}, \sQuote{sigma},
##'   \sQuote{theta}, or a function producing an lmerMod object.
##' @return object of class rlmerMod.
##' @seealso \code{\link[lme4]{lmer}}, \code{vignette("rlmer")}
##' @author Manuel Koller, with thanks to Vanda Lourenço for improvements.
##' @keywords models
##' @examples
##' ## dropping of VC
##' system.time(print(rlmer(Yield ~ (1|Batch), Dyestuff2, method="DASvar")))
##'
##' \dontrun{
##'   ## Default method "DAStau"
##'   system.time(rfm.DAStau <- rlmer(Yield ~ (1|Batch), Dyestuff))
##'   summary(rfm.DAStau)
##'   ## DASvar method (faster, less accurate)
##'   system.time(rfm.DASvar <- rlmer(Yield ~ (1|Batch), Dyestuff,
##'                                   method="DASvar"))
##'   ## compare the two
##'   compare(rfm.DAStau, rfm.DASvar)
##'
##'   ## Fit variance components with higher efficiency
##'   ## psi2propII yields squared weights to get robust estimates
##'   ## this is the same as using rlmer's argument `setting = "RSEa"`
##'   rlmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin,
##'         rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
##'         rho.sigma.b = psi2propII(smoothPsi, k = 2.28))
##'
##'   ## use chgDefaults for variance components including
##'   ## correlation terms (regular, non squared weights suffice)
##'   ## this is the same as using rlmer's argument `setting = "RSEa"`
##'   rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
##'         rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
##'         rho.b = chgDefaults(smoothPsi, k = 5.14, s=10),
##'         rho.sigma.b = chgDefaults(smoothPsi, k = 5.14, s=10))
##' }
##'
##' @importFrom lme4 lmer
##' @importFrom stats getCall
##' @export
rlmer <- function(formula, data, ..., method = c("DAStau", "DASvar"),
                  setting, rho.e, rho.b, rho.sigma.e, rho.sigma.b,
                  rel.tol = 1e-8, max.iter = 40 * (r + 1)^2, verbose = 0,
                  doFit = TRUE, init)
{
    lcall <- match.call()
    pf <- parent.frame()
    method <- match.arg(method)
    lobj <- .rlmerInit(lcall, pf, formula, data, method, rho.e, rho.b, rho.sigma.e,
                       rho.sigma.b, rel.tol, max.iter, verbose, init, setting, ...)$obj
    if (substr(lobj@method, 1, 3) == "DAS") {
        lobj@pp <- as(lobj@pp, "rlmerPredD_DAS")
        lobj@pp$method <- lobj@method
    }
    lobj@pp$initRho(lobj)
    lobj@pp$initMatrices(lobj)
    lobj@pp$updateMatrices()

    ## required for max.iter:
    r <- len(lobj, "theta")

    return(.rlmer(lobj, rel.tol, max.iter, verbose, doFit))
}

.rlmerInit <- function(lcall, pf, formula, data, method, rho.e, rho.b, rho.sigma.e,
                       rho.sigma.b, rel.tol, max.iter, verbose, init,
                       setting = c("RSEn", "RSEa"), ...) {
    if (missing(init) || is.null(init) || is.list(init)) {
        lcall2 <- lcall
        lcall2[setdiff(names(formals(rlmer)), names(formals(lmer)))] <- NULL
        lcall2$doFit <- NULL
        lcall2$REML <- TRUE
        lcall2[[1]] <- as.name("lmer")
        linit <- eval(lcall2, pf)
        if (!missing(init) && is.list(init)) {
            ## check sanity of input
            stopifnot(length(init$fixef) == length(fixef(linit)),
                      length(init$u) == length(getME(linit, "u")),
                      length(init$sigma) == length(sigma(linit)),
                      length(init$theta) == length(getME(linit, "theta")))
            ## convert object to rlmerMod
            linit <- as(linit, "rlmerMod")
            ## set all of the initial parameters, but do not fit yet
            setFixef(linit, unname(init$fixef))
            setSigma(linit, init$sigma)
            setTheta(linit, init$theta, fit.effects=FALSE, update.sigma=FALSE)
            setU(linit, init$u)
        }
        init <- linit
    } else if (is.function(init)) {
        init <- do.call(init,list(formula=formula, data=data, REML=TRUE, ...))
    } else if (is(init, "merMod") || is(init, "rlmerMod")) {
        ## check whether formula and data match with
        ## the ones in the provided in init
        if (is.null(icall <- getCall(init)))
            stop("Object 'init' should contain a 'call' component")
        if (!identical(as.character(lcall$formula), as.character(icall$formula)) |
            !identical(lcall$data, icall$data))
            warning("Arguments 'data' and 'formula' do not match with 'init': ",
                    "using model specification from 'init'")
    } else stop("Unsuitable init object, aborting.",
                "Expecting no, list (see ?rlmer), rlmerMod or merMod object")
    lobj <- as(init, "rlmerMod")
    lobj@call <- lcall

    ## give a warning if weights or offset are used
    if (any(lobj@resp$weights != 1))
        stop("Argument weights is unsave to use at the moment.")
    if (any(lobj@resp$offset != 0))
        warning("Argument offset is untested.")

    if (!missing(setting)) {
        if (!missing(rho.e) || !missing(rho.sigma.e) ||
            !missing(rho.b) || !missing(rho.sigma.b)) {
            overridden <- c()
            if (!missing(rho.e)) {
                overridden <- c(overridden, "'rho.e'")
            }
            if (!missing(rho.sigma.e)) {
                overridden <- c(overridden, "'rho.sigma.e'")
            }
            if (!missing(rho.b)) {
                overridden <- c(overridden, "'rho.b'")
            }
            if (!missing(rho.sigma.b)) {
                overridden <- c(overridden, "'rho.sigma.b'")
            }
            if (length(overridden) > 1) {
                args <- "Arguments "
            } else {
                args <- "Argument "
            }
            overridden <- paste(overridden, collapse = ", ")
            setting <- match.arg(setting)
            warning("Argument 'setting' specified together with ",
                    overridden, ". ", args, overridden,
                    " will override defaults of \"", setting, "\".")
        }
        setting <- match.arg(setting)
    } else {
        setting <- "RSEn"
    }
    if (missing(rho.e)) {
        rho.e <- smoothPsi
    }
    adjustForEfficiency <- setting == "RSEa"
    if (missing(rho.sigma.e)) {
        rho.sigma.e <- psi2propII(rho.e, adjust = adjustForEfficiency)
    }
    if (missing(rho.b)) {
        rho.b <- lapply(lobj@dim, getDefaultRhoB, rho = rho.e)
    } else if (!is.list(rho.b)) {
        rho.b <- rep.int(list(rho.b), length(lobj@dim))
        ## TODO warn if low asymptotic efficiency?
    }
    if (missing(rho.sigma.b)) {
        rho.sigma.b <- list()
        for (bt in seq_along(lobj@dim)) {
            if (lobj@dim[bt] == 1) {
                rho.sigma.b[[bt]] <- psi2propII(rho.b[[bt]], adjust = adjustForEfficiency)
            } else {
                rho.sigma.b[[bt]] <- rho.b[[bt]]
            }
        }
    } else if (!is.list(rho.sigma.b)) {
        rho.sigma.b <- rep.int(list(rho.sigma.b), length(lobj@dim))
        ## TODO warn if low asymptotic efficiency?
    }
    ## set arguments only relevant to rlmerMod
    lobj@rho.b <- rho.b
    lobj@rho.sigma.b <- rho.sigma.b
    if (!isTRUE(chk <- validObject(lobj))) stop(chk)
    lobj@rho.e <- rho.e
    lobj@rho.sigma.e <- rho.sigma.e
    if (method == "DAStau" & any(sapply(lobj@idx, nrow) > 2)) {
        warning("Method 'DAStau' does not support blocks of size larger than 2. ",
                "Falling back to method 'DASvar'.")
        method <- "DASvar"
    }
    lobj@method <- method
    return(list(obj = lobj, init = init))
}

getDefaultRhoB <- function(dimension, rho) {
    if (dimension == 1) {
        return(rho)
    }
    if (isDefaultHuberOrSmoothPsi(rho) && dimension < 8) {
        k <- switch(dimension - 1,
                    5.14, 5.55, 5.91, 6.25, 6.55, 6.84)
    } else {
        k <-
            findTuningParameter(asymptoticEfficiency(rho, "location"),
                                rho,
                                "tau",
                                dimension)
    }
    return(chgDefaults(rho, k = k))
}

.rlmer <- function(lobj, rel.tol, max.iter, verbose, doFit) {
    if (!doFit) return(updateWeights(lobj))

    ## do not start with theta == 0
    if (any(theta(lobj)[lobj@lower == 0] == 0)) {
        if (verbose > 0)
            cat("Setting variance components from 0 to 1\n")
        theta0 <- theta(lobj)
        theta0[lobj@lower == 0 & theta0 == 0] <- 1
        setTheta(lobj, theta0, fit.effects = TRUE, update.sigma = TRUE)
    } else {
        ## set theta at least once
        setTheta(lobj, theta(lobj), fit.effects = FALSE)
    }

    if (verbose > 0) {
        cat("\nrlmer starting values:\n")
        cat("sigma, theta: ", .sigma(lobj), ", ", theta(lobj), "\n")
        cat("coef: ", .fixef(lobj), "\n")
        if (verbose > 1)
            cat("b.s: ", b.s(lobj), "\n")
    }

    curWarnings <- list()
    lobj <- withCallingHandlers(.rlmer.fit(lobj, rel.tol, max.iter, verbose),
                                warning = function(w) {
                                    curWarnings <<- append(curWarnings,list(conditionMessage(w)))
                                    invokeRestart("muffleWarning")
                                })
    lobj@optinfo$warnings <- curWarnings

    if (verbose > 0) {
        cat("sigma, theta: ", .sigma(lobj), ", ", theta(lobj), "\n")
        cat("coef: ", .fixef(lobj), "\n")
        if (verbose > 1)
            cat("b.s: ", b.s(lobj), "\n")
    }

    return(updateWeights(lobj))
}

.rlmer.fit <- function(lobj, rel.tol, max.iter, verbose) {
    ## do fit: non diagonal case differently
    if (!isDiagonal(.U_b(lobj))) {
        if (lobj@method %in%  c("DASvar", "DAStau")) {
            lobj <- rlmer.fit.DAS.nondiag(lobj, verbose, max.iter, rel.tol)
        } else
            stop("Non-diagonal case only supported by DAStau and DASvar")
    } else {
        lobj <- rlmer.fit.DAS(lobj, verbose, max.iter, rel.tol)
    }
    return(lobj)
}

## DAS method
rlmer.fit.DAS.nondiag <- function(lobj, verbose, max.iter, rel.tol, method=lobj@method,
                                  checkFalseConvergence = TRUE) {
    if (!.isREML(lobj))
        stop("can only do REML when using averaged DAS-estimate for sigma")

    ## Prepare for DAStau
    if (method == "DAStau") {
        ## 4d int
        ## vectorize it!
        ghZ <- as.matrix(expand.grid(lobj@pp$ghz, lobj@pp$ghz, lobj@pp$ghz, lobj@pp$ghz))
        ghw <- apply(as.matrix(expand.grid(lobj@pp$ghw, lobj@pp$ghw, lobj@pp$ghw, lobj@pp$ghw)), 1, prod)
    } else {
        ghZ <- ghw <- c()
    }

    ## fit model using EM algorithm
    conv <- FALSE
    convBlks <- rep(FALSE, length(lobj@blocks))
    iter <- 0
    rel.tol <- sqrt(rel.tol)
    ## compute kappa
    kappas <- .kappa_b(lobj)
    ## zero pattern for T matrix
    nzT <- as.matrix(crossprod(bdiag(lobj@blocks[lobj@ind])) == 0)
    q <- lobj@pp$q
    ## false convergence indicator
    fc <- rep(FALSE, length(lobj@blocks))
    if (verbose > 0) {
        theta0 <- theta(lobj)
        if (verbose > 2) {
            coef0 <- .fixef(lobj)
            b.s0 <- b.s(lobj)
            sigma0 <- .sigma(lobj)
        }
    }
    ## iterate
    while(!conv && (iter <- iter + 1) < max.iter) {
        if (verbose > 0) cat("---- Iteration", iter, " ----\n")
        thetatilde <- theta(lobj)
        sigma <- .sigma(lobj)

        ## get expected value of cov(\tvbs)
        q <- len(lobj, "b")
        T <- switch(method,
                    DASvar=lobj@pp$Tb(),
                    DAStau=calcTau.nondiag(lobj, ghZ, ghw, .S(lobj), kappas, max.iter,
                                           rel.tol = rel.tol, verbose = verbose),
                    stop("Non-diagonal case only implemented for DASvar"))
        ## compute robustness weights and add to t and bs
        T[nzT] <- 0
        ## symmetrize T to avoid non symmetric warning, then apply chol
        T <- symmpart(T)
        ## save to cache
        lobj@pp$setT(T)
        ## apply chol to non-zero part only
        idx <- !.zeroB(lobj)
        ## stop if all are zero
        if (!any(idx)) break
        Tidx <- T[idx, idx]
        if (any(diag(Tidx) == 0.0)) {
            ## partially dropped block: can't handle yet.
            idxZeroes <- which(idx)[which(diag(Tidx) == 0.0)]
            blockAffected <- which(sapply(lobj@idx, function(cand) any(cand == idxZeroes)))
            stop("Covariance matrix for random effects block ", blockAffected,
                 " with grouping factor ", names(lobj@cnms)[blockAffected],
                 " is singular. Please simplify the model and run again.")
        }
        L <- t(chol(Tidx))
        T.bs <- numeric(q) ## set the others to zero
        T.bs[idx] <- forwardsolve(L, b.s(lobj)[idx])
        ## compute weights
        db <- .dk(lobj, sigma, FALSE, T.bs)[lobj@k]
        wbsEta <- wbsDelta <- numeric(q)
        for (type in seq_along(lobj@blocks)) {
            s <- lobj@dim[type]
            lidx <- as.vector(lobj@idx[[type]])
            if (s > 1) {
                ## for eta, we would actually need a little smaller
                ## tuning constants than for delta to get the same efficiency
                wbsEta[lidx] <- lobj@rho.sigma.b[[type]]@wgt(db[lidx])
                wbsDelta[lidx] <- (lobj@rho.sigma.b[[type]]@psi(db[lidx]) -
                                       lobj@rho.sigma.b[[type]]@psi(db[lidx] - s*kappas[type]))/s
            } else {
                lw <- lobj@rho.sigma.b[[type]]@wgt(db[lidx])
                wbsEta[lidx] <- lw
                wbsDelta[lidx] <- lw*kappas[type] ## adding kappa to wbsDelta in 1d case
            }
        }
        WbDelta <- Diagonal(x=wbsDelta)
        T <- WbDelta %*% T
        bs <- sqrt(wbsEta) * b.s(lobj)

        ## cycle block types
        for(type in seq_along(lobj@blocks)) {
            if (convBlks[type]) next
            bidx <- lobj@idx[[type]]
            if (verbose > 5) {
                cat("Tau for blocktype ", type, ":", as.vector(T[bidx[,1],bidx[,1]]), "\n")
            }
            ## catch dropped vc
            if (all(abs(bs[bidx]) < 1e-7)) {
                if (verbose > 1)
                    cat("Block", type, "dropped (all = 0), stopping iterations.\n")
                Ubtilde <- lobj@blocks[[type]]
                pat <- Ubtilde != 0
                Lind <- Ubtilde[pat]
                thetatilde[Lind] <- 0
                convBlks[type] <- TRUE
                next
            }
            s <- nrow(bidx)
            K <- ncol(bidx)
            ## right hand side
            ## sum over blocks
            rhs <- matrix(0, s, s)
            for (k in 1:K) {
                ## add weights
                rhs <- rhs + T[bidx[,k],bidx[,k]]
            }
            rhs <- rhs / K
            ## left hand side
            lbs <- matrix(bs[bidx], ncol(bidx), nrow(bidx), byrow=TRUE)
            ## add left hand side
            lhs <- crossprod(lbs / sigma) / K
            if (verbose > 2) {
                cat("LHS:", as.vector(lhs), "\n")
                cat("RHS:", as.vector(rhs), "\n")
                cat("sum(abs(LHS - RHS)):", sum(abs(lhs - rhs)), "\n")
            }
            ## if (isTRUE(all.equal(rhs, lhs, check.attributes=FALSE, tolerance = rel.tol))) {
            diff <- abs(rhs - lhs)
            if (all(diff < rel.tol * max(diff, rel.tol))) {
                if (verbose > 1)
                    cat("Estimating equations satisfied for block", type,
                        ", stopping iterations.\n")
                convBlks[type] <- TRUE
                next
            }
            deltaT <- backsolve(lchol(rhs), lchol(lhs))
            if (verbose > 1) cat("deltaT:", c(deltaT), "\n")
            ## get old parameter estimates for this block
            Ubtilde <- lobj@blocks[[type]]
            pat <- Ubtilde != 0
            Lind <- Ubtilde[pat]
            diagLind <- diag(Ubtilde)
            Ubtilde[pat] <- thetatilde[Lind]
            ## update Ubtilde by deltaT
            thetatilde[Lind] <- tcrossprod(Ubtilde, deltaT)[pat]
            ## FIXME: check boundary conditions?
            ## check if varcomp is dropped
            if (all(thetatilde[diagLind] < 1e-7)) {
                thetatilde[Lind] <- 0
                convBlks[type] <- TRUE
                next
            }
            ## check if this block is converged
            diff <- abs(thetatilde[Lind] - theta(lobj)[Lind])
            if (verbose > 3)
                cat("criterion:", sum(diff), ">=",
                    rel.tol * max(diff, rel.tol), ":",
                    sum(diff) < rel.tol * max(diff, rel.tol), "\n")
            if (sum(diff) < rel.tol * max(diff, rel.tol)) {
                convBlks[type] <- TRUE
                ## check if estimating equations are satisfied
                if (checkFalseConvergence) {
                    if (verbose > 3)
                        cat("checking estimating equations:", sum(abs(lhs - rhs)),
                            ">", sqrt(rel.tol), ":", sum(abs(lhs - rhs)) > sqrt(rel.tol), "\n")
                    if (sum(abs(lhs - rhs)) > sqrt(rel.tol))
                        fc[type] <- TRUE
                }
                next
            }
        }
        ## set theta
        setTheta(lobj, thetatilde, fit.effects = TRUE,
                 update.sigma = FALSE)
        ## update sigma without refitting effects
        updateSigma(lobj, fit.effects = FALSE)
        if (verbose > 0) {
            cat("delta theta:", format(theta0 - thetatilde, nsmall=20, scientific=FALSE),
                "\n")
            theta0 <- thetatilde
            if (verbose > 1) {
                cat(sprintf("delta coef:  %.12f\n", sum(abs(coef0 - .fixef(lobj)))))
                cat(sprintf("delta u:     %.12f\n", sum(abs(b.s0 - b.s(lobj)))))
                cat(sprintf("delta sigma: %.12f\n", abs(sigma0 - .sigma(lobj))))
                coef0 <- .fixef(lobj)
                b.s0 <- b.s(lobj)
                sigma0 <- .sigma(lobj)
                if (verbose > 2) {
                    cat("theta:", format(thetatilde, nsmall=20, scientific=FALSE), "\n")
                    cat("coef:   ", .fixef(lobj),"\n")
                    cat("b.s:    ", b.s(lobj), "\n")
                    cat("sigmae: ", .sigma(lobj), "\n")
                }
            }
        }
        if (all(convBlks)) conv <- TRUE
    }

    optinfo <- list(optimizer = "rlmer.fit.DAS.nondiag",
                    conv = list(opt = 0),
                    feval = iter,
                    warnings = list(),
                    val = diff)

    if (iter == max.iter) {
        warning(wt <- "iterations did not converge, returning unconverged estimate.")
        optinfo$warnings <- list(wt)
        optinfo$conv$opt <- 1
    }
    if (any(fc)) {
        warning(wt <- "algorithm converged, but estimating equations are not satisfied.")
        optinfo$warnings <- c(optinfo$warnings, list(wt))
        optinfo$conv$opt <- 2
    }
    lobj@optinfo <- optinfo
    return(lobj)
}

## DAS method
rlmer.fit.DAS <- function(lobj, verbose, max.iter, rel.tol) {
    if (!.isREML(lobj))
        stop("can only do REML when using averaged DAS-estimate for sigma")

    lupdateTheta <- switch(lobj@method,
                           DASvar=,
                           DAStau=updateThetaTau,
                           stop("method not supported by rlmer.fit.DAS:", lobj@method))

    ## fit
    converged <- FALSE
    theta0 <- theta(lobj)
    if (verbose > 1) {
        coef0 <- .fixef(lobj)
        b.s0 <- b.s(lobj)
        sigma0 <- .sigma(lobj)
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
                cat(sprintf("delta coef:  %.12f\n", sum(abs(coef0 - .fixef(lobj)))))
                cat(sprintf("delta u:     %.12f\n", sum(abs(b.s0 - b.s(lobj)))))
                cat(sprintf("delta sigma: %.12f\n", abs(sigma0 - .sigma(lobj))))
                coef0 <- .fixef(lobj)
                b.s0 <- b.s(lobj)
                sigma0 <- .sigma(lobj)
                if (verbose > 2) {
                    cat("theta:  ", theta(lobj),"\n")
                    cat("coef:   ", .fixef(lobj),"\n")
                    cat("b.s:    ", b.s(lobj), "\n")
                    cat("sigmae: ", .sigma(lobj), "\n")
                }
            }
        }

        ## all zero or change smaller than relative tolerance
        ## all zero: we can't get out of this anyway, so we have to stop.
        converged <- all(theta1 == 0) || sum(abs(theta0 - theta1)) <
            200*rel.tol*sum(abs(theta0))
        if (verbose > 1)
            cat(sprintf("Criterion: %.12f, %.12f", sum(abs(theta0 - theta1)),
                        sum(abs(theta0 - theta1)) / rel.tol / sum(abs(theta0)) / 200), "\n")
        theta0 <- theta1
    }

    optinfo <- list(optimizer = "rlmer.fit.DAS",
                    conv = list(opt = 0),
                    feval = iter,
                    warnings = list(),
                    val = diff)

    if (iter == max.iter) {
        warning(wt <- "iterations did not converge, returning unconverged estimate.")
        optinfo$warnings <- list(wt)
        optinfo$conv$opt <- 1
    }
    lobj@optinfo <- optinfo
    return(lobj)
}
