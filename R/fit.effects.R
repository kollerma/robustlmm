fitEffects <- function(par, object, rel.tol = 1e-8, max.iter = 1000) {
    n <- object@pp$n
    p <- object@pp$p
    q <- object@pp$q
    ## Mallows design weights: eta multiplies the robustness
    ## weight of each observation, so the IRLS fixed point solves the
    ## eta-weighted beta/u estimating equations
    ##   (U_e X)' H psi_e(r~) = 0,
    ##   (U_e Z U_b)' H psi_e(r~) - Lambda_b w_b u = 0,  H = diag(eta).
    ## Objects serialized before the eta field existed have length-0
    ## eta; treat as 1.
    eta <- object@resp$eta
    if (length(eta) != n) eta <- rep.int(1, n)
    idx <- c(rep(TRUE, p), !.zeroB(object))
    par[!idx] <- 0 ## set random effects corresponding to dropped vc to 0
    setU(object, par[p+1:q])
    setFixef(object, par[1:p])
    ZL <- .U_eZU_b(object)
    MAT1 <- object@pp$MAT1
    MAT1[1:n,p+1:q] <- ZL
    MAT1[n+1:q,p+1:q] <- sqrt(.Lambda_b(object))
    MAT2 <- object@pp$MAT2
    MAT2[,p+1:q] <- ZL
    ## subtract the offset: fitEffects solves for (beta, u) in the
    ## offset-corrected response; .mu() adds the offset back.
    U_ey <- .U_e(object) %*% (object@resp$y - object@resp$offset)
    converged <- FALSE
    iter <- 1
    while(!converged && iter <= max.iter) {
        ## set parameters
        ## calculate weights
        sqW <- Diagonal(x=sqrt(c(w.e <- wgt.e(object) * eta, wgt.b(object))))
        ## build matrix (a factor, resp)
        MAT <- sqW %*% MAT1
        w.y <- w.e * U_ey
        ## solve to get new parameters
        XtX <- crossprod(MAT)
        rhs <- crossprod(MAT2, w.y)
        ## warnings of the plain solve (CHOLMOD "not positive definite")
        ## are held back and dropped when the ridge below takes over
        held <- list()
        par1 <- tryCatch(withCallingHandlers(drop(solve(XtX, rhs)),
            warning = function(w) {
                held[[length(held) + 1L]] <<- w
                invokeRestart("muffleWarning")
            }), error = function(e) {
            ## With redescending psi-functions a random effect can carry
            ## no information at the current weights -- its own weight and
            ## those of all its observations are zero -- which makes the
            ## normal equations singular. A tiny ridge on the random-effect
            ## block sets such coefficients to zero (they are rejected);
            ## a fixed-effects design without information still errors.
            held <<- list()
            ridge <- 1e-10 * max(abs(diag(XtX)))
            drop(solve(XtX + Diagonal(x = c(rep(0, p), rep(ridge, q))), rhs))
        })
        for (w in held) {
            class(w) <- c("rlmerInnerWarning", class(w))
            warning(w)
        }
        ## check convergence
        converged <- sum(abs(par - par1)) < rel.tol * sum(abs(par))
        ## update parameters
        setU(object, par1[p+1:q])
        setFixef(object, par1[1:p])
        par <- par1
        iter <- iter + 1
    }
    object
}
