## code taken from rootSolve (multiroot()) and modified to suit
## this particular problem
multiroot <- function(f, start, maxiter=100, rtol=1e-6, atol=1e-8, ctol=1e-8,
                      useFortran = TRUE, lower = NULL, verbose=FALSE,
                      method = NULL, jacfunc = NULL, ...)  {
    ## revert to rootSolves implementation if no special features are requested
    if (useFortran && is.null(lower)) {
        require(rootSolve)
        if (!is.null(method)) warning("robustlmm::multiroot: ignoring method argument")
        return(rootSolve::multiroot(f, start, maxiter=maxiter, rtol=rtol, atol=atol, ctol=ctol,
                                    useFortran=useFortran, verbose=verbose, ...))
    }
    
    N        <- length(start)
    if (!is.numeric(start))
        stop("start conditions should be numeric")
    if (!is.numeric(maxiter))
        stop("`maxiter' must be numeric")
    if (as.integer(maxiter) < 1)
        stop ("maxiter must be >=1")
    if (!is.numeric(rtol))
        stop("`rtol' must be numeric")
    if (!is.numeric(atol))
        stop("`atol' must be numeric")
    if (!is.numeric(ctol))
        stop("`ctol' must be numeric")
    if (length(atol) > 1 && length(atol) != N)
        stop("`atol' must either be a scalar, or as long as `start'")
    if (length(rtol) > 1 && length(rtol) != N)
        stop("`rtol' must either be a scalar, or as long as `y'")
    if (length(ctol) > 1)
        stop("`ctol' must be a scalar")
    if (!is.null(lower) && length(lower) != length(start))
        stop("`lower' must be of the same length as start conditions")
    if (is.null(jacfunc)) {
        require(rootSolve)
        .jacfunc <- function(time,y,parms, ...) f(y, ...)
        ## (smaller pert does just destabilize the iterations)
        jacfunc <- function(x, ...)
            jacobian.full(x, .jacfunc, dy=reffx, pert=1e-3, ...)
            ## or from package numDeriv
            ## jacobian(f, x, method = "Richardson",  ...)
    }
    
    precis   <- NULL
    
    x        <- start
    jacob    <- matrix(nrow=N,ncol=N,data=0)
    reffx    <- f(x,...)     # function value,
    
    if (length (reffx) != N)
        stop("'f', function must return as many function values as elements in start")
    
    for (i in 1:maxiter) {
        refx <- x
        oldfx <- reffx 
        ## estimate jacobian (smaller pert does just destabilize the iterations)
        jacob <- jacfunc(x, ...)
        ## new estimate 
        relchange <- as.numeric(solve(jacob,-1*reffx))
        ## avoid too long jumps
        maxjump <- pmax(abs(x) / 2, 0.24)
        idx <- abs(relchange) > maxjump
        relchange[idx] <- sign(relchange[idx]) * maxjump[idx]
        x  <- x + relchange
        ## set x values to lower bound that are below the lower bound
        if (!is.null(lower)) {
            idx <- x < lower
            x[idx] <- lower[idx]
            relchange[idx] <- x[idx] - refx[idx]
        }
        reffx  <- f(x,...)     # function value,
        ## jump out of bad solution area
        if (!is.null(method)) {
            switch(method,
                   DAS={
                       ## avoid zero and near zero solution
                       ## if x was decreased, reffx is >= 0, oldfx was < 0 and
                       ## (relchange was smaller than -0.25)
                       ## or x was increased, reffx is =< 0, oldfx was > 0
                       idx <- (relchange < -0.25 & reffx >= 0 & oldfx < 0) |
                              (relchange > 0.25 & reffx <= 0 & oldfx > 0)
                       if (any(idx)) {
                           ## jump to the middle of refx and xnew
                           x[idx] <- (x[idx] + refx[idx]) / 2
                           reffx <- f(x, ...)
                           if (verbose) {
                               cat("Jumping to the middle..., idx:", idx,
                                   "\nrelchange:", relchange, "\n")
                               cat("refx:", refx, ", oldfx:", oldfx, "\n")
                               cat("xnew:", x, ", reffx:", reffx, "\n")
                           }
                       }
                   })
        }
        ## check for convergence
        if (max(abs(x - refx)) < ctol || all(abs(reffx) < atol)) break
    } # end for
    names(x) <- names(start)

    return(list(root=x,f.root=reffx,iter=i,estim.precis=precis[length(precis)]))
}  # end multiroot
