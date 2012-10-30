### calculate the expected value to be used in D.re
.calcE.D.re <- function(q, rho) {
    if (q == 1) return(rho@EDpsi())
    ## calculate first part
    tfun <- function(u, v) {
        uv <- u^2 + v
        dk <- sqrt(uv/q)
        ret <- (rho@Dwgt(dk))/(q*dk)*u^2*dnorm(u)*dchisq(v, q-1)
        ##ret[uv == 0 | is.infinite(uv)] <- 0
        ret
    }
    tfun2 <- function(v) integrate(tfun, -Inf, Inf, v=v)$value
    ret <- integrate(function(v) sapply(v, tfun2), 0, Inf)$value
    ## calculate the second part
    ## tfun3 <- function(s) rho@psi(s)*2^(1-q/2)*s^(q-2)*exp(-s^2/2)/gamma(q/2)
    ## the above makes problems for large q
    tfun3 <- function(s) rho@wgt(sqrt(s/q))*dchisq(s, q)
    ret + integrate(tfun3, 0, Inf)$value
}


##' Calculate the two dimensional integral G
##'
##' The \eqn{\rho}{rho}-functions can be chosen to be identical, but
##' then the scale estimate will not have bounded influence
##' (if the \eqn{\rho}{rho}-function is not redescending). As an
##' alternative, one can use \eqn{\rho_\sigma = \psi_e'^2/2}{rho_sigma = psi_e'^2}.
##' This corresponds to Huber's Proposal II.
##' 
##' @title Function \eqn{G_\sigma(a, s)}{G_sigma(a, s)}
##' @param tau denominator inside function \eqn{\rho_\sigma}{rho_sigma}
##' @param a first argument of G, leverage of observation
##' @param s second argument of G, scale of \eqn{g_{[-i]}}{g_{-i}}
##' @param rho psi_func-object
##' @param rho.sigma psi_func-object (scale estimate)
##' @param wExp exponent for robustness weights
##' @param pp rlmerPredD object
G <- function(tau = rep(1, length(a)), a, s, rho, rho.sigma, wExp, pp) {
    if (length(a) != length(s))
        stop("length of vectors a and s not equal")

    ret <- numeric(length(a))
    wgt <- .wgtTau(rho.sigma, wExp)
    for (i in 1:length(a)) {
        ## ## inner integration over z
        ## inner <- function(e) 
        ##     sum(rho@rho((e - a[i]*rho@psi(e) - s[i]*ghz)/tau[i])*ghw)
        ## ## outer integration over e
        ## ret[i] <- sum(sapply(ghz, inner)*ghw)
        x <- (pp$ghZ - a[i]*rho@psi(pp$ghZ) - s[i]*pp$ghZt)/tau[i]
        ret[i] <- sum(wgt(x)*x^2*pp$ghW)
    }
    ret
}

##' Calculate s = sqrt(v_{-i}), standard deviation of the
##' average term in the linear approximation of the residuals,
##' (in the case of diagonal Sigma_b / Lambda_\theta)
##' This scale is used in G 
##'
##' @title Calculate s
##' @param object rlmerMod object
##' @param theta logical if s is required for theta or sigma_e
##' @return vector of s
.s <- function(object, theta = FALSE, pp = object@pp) {
    if (theta) {
        M1 <- pp$K()
        M2 <- pp$L
        diag(M2) <- 0
    } else {
        M1 <- pp$A
        diag(M1) <- 0
        M2 <- pp$B()
    }
    ## setting NA to 0 (they come from 0 variance components)
    M1[is.na(M1)] <- 0
    M2[is.na(M2)] <- 0
    ## calculate s:
    ret <- pp$rho_e@Epsi2() * rowSums(M1^2)
    if (any(!pp$zeroB)) ret <- ret + drop(M2^2 %*% diag(pp$Epsi_bpsi_bt))
    sqrt(ret)
}

##' Calculate non-diagonal variant of G
##'
##' Calculate approximation of the expectation
##' E[w(d_k(\tilde u)) \tilde u \tilde u\tr]
##'
##' @title Calculate G for non-diagonal \Lambda_\theta
##' @param object rlmerMod object
##' @param S return object of .S, list of S_ks
.Gnd <- function(object, S = .S(object)) {
    wExp.b <- object@wExp.b
    psi.b <- object@rho.b@psi
    wgt.b <- object@rho.b@wgt
    wgt.sb <- .wgtTau(object@rho.sigma.b, wExp.b)

    integrand <- function(uv0, Lkk, Sk, s) {
        lidx <- 1:s
        ## variable transformation
        tmp <- 1 - uv0^2
        uv <- uv0 / tmp
        utilde <- uv[lidx] - wgt.b(sqrt(mean(uv[lidx]^2))) * Lkk %*% uv[lidx] -
            crossprod(Sk, uv[s+lidx])
        dk <- sqrt(mean(utilde^2))
        as.vector(wgt.sb(dk) * tcrossprod(utilde) *
                  prod(dnorm(uv) * (1 + uv0^2) / tmp^2))
    }
    lastSk <- lastLkk <- matrix()
    lastRet <- NA
    idx <- !object@pp$zeroB
    int <- function(k) {
        ind <- object@k == k
        ## return 0-matrix is idx[k] is FALSE
        if (!all(idx[ind])) {
            ## check if whole block is zero 
            if (!all((!idx)[ind]))
                stop("result of getZeroU ambiguous: either the whole block is zero or not")
            return(matrix(0, sum(ind), sum(ind)))
        }
        ## else calculate the integral:
        Lkk <- as.matrix(object@pp$L[ind, ind])
        Sk <- as.matrix(S[[k]])
        ## only recalculate ret if the new matrices Sk and Lkk
        ## are substantially different from the old ones
        if (isTRUE(all.equal(lastSk, Sk)) && isTRUE(all.equal(lastLkk, Lkk)))
            return(lastRet)
        s <- sum(ind)
        ## FIXME: this is just too slow
        ret <- adaptIntegrate(integrand, rep(-1, 2*s), rep(1, 2*s),
                              Lkk = Lkk, Sk = Sk, s = s,
                              fDim = sum(ind)^2, maxEval = 1000000, absError = 1e-3)
                              ## absError = 1e-5)
        if (ret$returnCode > 0) {
            stop("adaptIntegrate failed with returnCode", ret$returnCode,
                 " for k = ", k)
        }
        
        ## convert back to matrix
        ret <- matrix(ret$integral, sum(ind))
        ## update cache
        lastRet <<- ret
        lastSk <<- Sk
        lastLkk <<- Lkk

        return(ret)
    }
    lapply(seq_along(S), int)
}

##' Calculate S_k (for all ks), Cholesky decomposition of the
##' covariance matrix of the average term in the linear approximation
##' of the spherical random effects (in the case of non-diagonal
##' Sigma_b / Lambda_\theta) This is then used in Gnd
##'
##' @title Calculate all S_k
##' @param object rlmerMod object
##' @return list of matrices S_k
.S <- function(object) {
    ret <- list()
    idx <- !object@pp$zeroB
    tmpEL <- object@pp$L %*% object@pp$Epsi_bpsi_bt
    for (k in 1:max(object@k)) {
        ind <- object@k == k
        ## check if the block was dropped (u == 0)
        if (all(idx[ind])) {
            ## Ltmp <- object@pp$L[ind, !ind, drop=FALSE]
            ## use symmpart to calm chol()
            ret <- c(ret, list(drop(chol(symmpart(
                ## Matrix K_{\theta,k} consists of the
                ## rows belonging block k and all columns
                object@rho.e@Epsi2() * tcrossprod(object@pp$K()[ind, , drop=FALSE]) +
                ## Matrix L_{\theta,-k} consists of the
                ## rows belonging to block k and all but columns
                ## except the ones belonging to block k
                tcrossprod(object@pp$L[ind, !ind, drop=FALSE],
                           tmpEL[ind, !ind, drop=FALSE])
                ## tcrossprod(Ltmp, Ltmp %*%
                ##            object@pp$Epsi_bpsi_bt[!ind, !ind, drop=FALSE])
                )))))
        } else {
            ## check if whole block is zero 
            if (!all((!idx)[ind]))
                stop("result of getZeroU ambiguous: either the whole block is zero or not")
            ## return a zero block of the correct size
            ret <- c(ret, list(matrix(0, sum(ind), sum(ind))))
        }
    }
    ret
}

##' Calculate kappa for DASvar and DAStau method.
##'
##' @title Kappa for scale estimates
##' @param rho rho-function to be used
##' @param wExp exponent
##' @rdname calcKappa
calcKappaTau <- function(rho, wExp) {
    ## function to calculate the expectation
    wgt <- .wgtTau(rho, wExp)
    switch(wExp + 1,
           rho@Erho()*2  / integrate(function(x) wgt(x)*dnorm(x), -Inf, Inf)$value,
           rho@EDpsi() / integrate(function(x) wgt(x)*dnorm(x), -Inf, Inf)$value,
           rho@Epsi2() / integrate(function(x) wgt(x)*dnorm(x), -Inf, Inf)$value)
}
##' @param object rlmerMod object
##' @rdname calcKappa
calcKappaTauB <- function(object) {
    kappas <- rep(1, length(object@blocks))
    rho.sigma <- object@rho.sigma.b
    wExp <- object@wExp.b
    wgt <- .wgtTau(rho.sigma, wExp)
    for(type in seq_along(object@blocks)) {
        s <- nrow(object@idx[[type]])
        kappas[type] <- if (s==1) {
            calcKappaTau(rho.sigma, wExp)
        } else {
            tfun <- function(u, v) {
                ret <- wgt(sqrt((u^2 + v)/s))*u^2*dnorm(u)*dchisq(v, s-1)
                ret
            }
            tfun2 <- function(v) integrate(tfun, -Inf, Inf, v=v)$value
            tfun3 <- function(v) wgt(sqrt(v/s))*dchisq(v, s)
            integrate(function(v) sapply(v, tfun2), 0, Inf)$value /
                integrate(tfun3, 0, Inf)$value
        }
    }
    kappas
}

.wgtTau <- function(rho, wExp)
    switch(wExp + 1, ## wExp == 0:
           function(x) ifelse(x == 0, rho@Dpsi(0), rho@rho(x)/(x*x)*2),
           rho@wgt, ## wExp == 1
           function(x) { ## wExp == 2
               x <- rho@wgt(x)
               x*x
           })

##' Calculate tau for DAS estimate
##'
##' @title DAS-estimate
##' @param a coefficient for psi-function in approximation
##' @param s scale for rest in approximation
##' @param rho.e rho-function used for effets
##' @param rho.sigma.e rho-function used for sigma
##' @param wExp.e exponent used in scale estimation
##' @param pp rlmerPredD object
##' @param kappa kappa tuning parameter
##' @param tau initial values for tau
##' @param method method to compute tau with
##' @param rel.tol relative tolerance for calculating tau
##' @param max.it maximum number of iterations allowed
calcTau <- function(a, s, rho.e, rho.sigma.e, wExp.e, pp,
                    kappa, tau = rep(1, length(a)), method="DAStau",
                    rel.tol = 1e-6, max.it = 200) {
    stopifnot(length(a) == length(s))
    if (method == "DASvar" || !is.numeric(tau) || length(tau) != length(a)) {
        ## FIXME: this always returns tau_e irrespective of a and s...
        Tau <- with(pp, V_e - EDpsi_e * (t(A) + A) + Epsi2_e * tcrossprod(A) +
                             B() %*% tcrossprod(Epsi_bpsi_bt, B()))
        tau <- sqrt(diag(Tau))
        if (method == "DASvar") return(tau)
    }
    psi <- rho.e@psi
    psiZ <- psi(pp$ghZ)
    ## function to find root from
    wgt <- .wgtTau(rho.sigma.e, wExp.e)
    fun <- switch(wExp.e + 1, ## wExp.e == 0
                  function(tau, a, s) {
                      t <- (pp$ghZ-a*psiZ-s*pp$ghZt)/tau
                      sqrt(2*sum(rho.sigma.e@rho(t)*(tau*tau)*pp$ghW) /
                           (kappa*sum(wgt(t)*pp$ghW)))
                  }, ## wExp.e == 1:
                  function(tau, a, s) {
                      t <- (pp$ghZ-a*psiZ-s*pp$ghZt)/tau
                      sqrt(sum(rho.sigma.e@psi(t)*t*(tau*tau)*pp$ghW) /
                           (kappa*sum(wgt(t)*pp$ghW)))
                  }, ## wExp.e == 2:
                  function(tau, a, s) {
                      t <- (pp$ghZ-a*psiZ-s*pp$ghZt)/tau
                      sqrt(sum(rho.sigma.e@psi(t)^2*(tau*tau)*pp$ghW) /
                           (kappa*sum(wgt(t)*pp$ghW)))
                  })

    for (i in seq_along(a)) {
        ## check if we already calculated this
        check <- FALSE
        j <- 0
        while(!check && (j <- j + 1) < i) {
            if (isTRUE(all.equal(c(a[i], s[i]), c(a[j], s[j])))) {
                check <- TRUE
                tau[i] <- tau[j]
            }
        }
        if (!check) {
            ##cat(sprintf("calculating tau[%i]: a[%i] = %g, s[%i] = %g...\n", i, i, a[i], i, s[i]))
            it <- 0
            conv <- FALSE
            while(!conv && (it <- it + 1) < max.it && tau[i] <= 1) {
                ##print(tau[i])
                tau0 <- tau[i]
                tau[i] <- fun(tau0, a[i], s[i])
                conv <- abs(tau[i] - tau0) < rel.tol * max(rel.tol, tau0)
            }
            if (it >= max.it) {
                warning("calculation of tau[", i, "] reached max.it")
            }
            if (tau[i] > 1) {
                tau[i] <- 1
                warning("tau[", i, "] > 1, setting it to 1")
            }            
        }
    }
    
    tau
}

calcTau.nondiag <- function(object, ghZ, ghw, skbs, kappas, max.iter) {
    ## define 4d integration function
    int4d <- function(fun) drop(apply(ghZ, 1, fun) %*% ghw)
    ## initial values
    T <- object@pp$T()
    ## move into list
    TkbsI <- list()
    ## cycle blocks
    for (type in seq_along(object@blocks)) {
        bidx <- object@idx[[type]]
        for (k in 1:ncol(bidx)) ## 1:K
            TkbsI <- c(TkbsI, list(T[bidx[,k],bidx[,k]]))
    }

    ## Compute Taus
    Tbks <- list()
    wgt <- .wgtTau(object@rho.sigma.b, object@wExp.b)
    ## cycle block types
    for (type in seq_along(object@blocks)) {
        bidx <- object@idx[[type]]
        s <- nrow(bidx)
        ind <- which(object@ind == type) ## (in terms of blocks on diagonal)
        idx <- !object@pp$zeroB
        if (!all(idx[bidx])) { ## block has been dropped
            ## check if the whole block is zero
            if (!all((!idx)[bidx]))
                stop("result of getZeroU is ambiguous: either the whole block is zero or not")
            Tbks <- c(Tbks, rep(list(matrix(0, s, s)), length(ind)))
            next
        }
        ## block has not been dropped: calculate Tbk
        ## catch 1d case
        if (s == 1) {
            bidx <- drop(bidx) ## is a vector (in terms of items of b.s)
            tmp <- calcTau(diag(object@pp$L)[bidx], unlist(skbs[ind]),
                           object@rho.b, object@rho.sigma.b, object@wExp.b, object@pp, kappas[type])
            Tbks <- c(Tbks, as.list(tmp))
        } else if (s == 2) { ## 2d case
            f <- sqrt(kappas[type]) ## normalizing constant
            lastSk <- lastkk <- matrix()
            lastRet <- NA
            ## cycle blocks
            for (k in 1:ncol(bidx)) { ## 1:K
                lTbk <- as.matrix(TkbsI[[ind[k]]])
                lbidx <- bidx[,k]
                Lkk <- as.matrix(object@pp$L[lbidx, lbidx])
                Sk <- as.matrix(skbs[[ind[k]]])
                ## check for almost equal blocks
                if (!(isTRUE(all.equal(lastSk, Sk)) && isTRUE(all.equal(lastLkk, Lkk)))) {
                    ## Find Tbk for this new block...
                    lastSk <- Sk
                    lastLkk <- Lkk
                    ##cat("TbkI:", lTbk, "\n")
                    conv <- FALSE
                    iter <- 0
                    while (!conv && (iter <- iter + 1) < max.iter) {
                        qrlTbk <- qr(lTbk)
                        if (qrlTbk$rank < s) {
                            warning("Tbk not of full rank (iter=", iter, "!!\nlTbk =", as.vector(lTbk), ", setting it to Tb()\n")
                            conv <- TRUE
                            lTbk <- as.matrix(TkbsI[[ind[k]]])
                        } else {
                            funA <- function(u) {
                                btilde <- u[1:2] - wgt(sqrt(mean(u[1:2]^2))) * Lkk %*% u[1:2] - crossprod(Sk, u[3:4])
                                wgt(sqrt(mean(qr.solve(qrlTbk, btilde)^2))) ## not needed: *prod(dnorm(c(u1,u2,u3,u4))) 
                            }
                            a <- int4d(funA)
                            funB <- function(u) {
                                btilde <- u[1:2] - wgt(sqrt(mean(u[1:2]^2))) * Lkk %*% u[1:2] - crossprod(Sk, u[3:4])
                                wgt(sqrt(mean(qr.solve(qrlTbk, btilde)^2)))*tcrossprod(btilde)
                            }
                            B <- matrix(int4d(funB), s)
                            lTbk1 <- B/(f*a)
                            conv <- isTRUE(all.equal(lTbk, lTbk1)) ## rel.tol default 1e-8
                            ##cat(sprintf("k=%i, iter=%i, conv=%s\n", k, iter, conv))
                            ##cat("Tbk:", lTbk1, "\n")
                            lTbk <- lTbk1
                        }
                    }
                    lastRet <- lTbk
                }
                ## add to result
                Tbks <- c(Tbks, list(lastRet))
            }
        } else {
            warning("DAStau for blocks of dimension > 2 not defined, falling back to DASvar")
            stop("yes, do as promised")
        }    
    }

    ## combine into Matrix
    bdiag(Tbks)    
}
    
##' Update sigma: calculate averaged DAS scale
##'
##' @title Calculate averaged DAS scale
##' @param object rlmerMod-object
##' @param max.iter maximum iterations allowed in fitting
##' @param rel.tol stopping criterion
##' @param fit.effects fit effects for each new sigma?
##' @return rlmerMod-object
updateSigma <- function(object, max.iter = 100, rel.tol = 1e-6, fit.effects = TRUE) {
    rho.sigma.e <- object@rho.sigma.e
    wExp.e <- object@wExp.e
    
    switch(object@method,
           DASvar =,
           DAStau = { ## DAStau and DASvar variant (depending on class type)
               kappa <- object@pp$kappa_e
               tau <- object@pp$tau_e()
           }, { ## otherwise
               ## get matrices for r part
               a <- diag(object@pp$A)
               s <- .s(object, theta = FALSE)
               if (any(is.na(a))) {
                   cat("a, s: ", NA, "theta:", theta(object), "\n")
                   setSigma(object, as.numeric(NA))
                   return(object)
               }
               ## calculate kappas
               kappas <- G(,a,s,object@rho.e,rho.sigma.e,wExp.e,object@pp)
               sk <- sum(kappas)
               tau <- rep(1, object@pp$n)
           })
    ## cat("  tau:", tau, "\n")
    ## cat("theta:", theta(object), "\n")
    ## cat("   sk:", sk, "\n")
    
    ## ## old uniroot variant
    ## fun <- function(logscale, r) {
    ##     lr <- r/exp(logscale)
    ##     sum(.wgtxy2(rho.sigma.e, lr, lr, wExp.e)) - sk
    ## }
    ## zero.fun <- if (fit.effects) function(logscale, r) {
    ##     setSigma(object, exp(logscale))
    ##     fit.effects(c(object@pp$beta, object@pp$b.s), object)
    ##     fun(logscale, object@resp$wtres)
    ## } else fun
    ## res <- uniroot(zero.fun, c(log(1e-7), 100), tol = rel.tol,
    ##                r = object@resp$wtres)
    ## scale <- exp(res$root)
    
    ## new iterative reweighting
    fun <- if (object@method %in% c("DAStau", "DASvar")) {
        wgt <- .wgtTau(rho.sigma.e, wExp.e)
        tau2 <- tau*tau
        function(scale, r) {
            lr <- r/tau/scale
            sqrt(sum(wgt(lr)*r*r)/sum(wgt(lr)*tau2)/kappa)
        }
    } else
        function(scale, r) {
            lr <- r/tau
            sqrt(sum(.wgtxy2(rho.sigma.e, lr/scale, lr, wExp.e))/sk)
        }
    fun2 <- if (fit.effects) function(scale, r) {
        setSigma(object, scale)
        fit.effects(c(object@pp$beta, object@pp$b.s), object)
        fun(scale, object@resp$wtres)
    } else fun
    scale0 <- .sigma(object)
    ## logscale0 <- log(scale0)
    converged <- FALSE
    it <- 0
    while(!converged && (it <- it + 1) < max.iter) {
        scale <- fun2(scale0, object@resp$wtres)
        ## converged <- abs(scale/scale0 - 1) < rel.tol
        converged <- abs(scale - scale0) < rel.tol * max(rel.tol, scale0)
        scale0 <- scale
        ## logscale <- log(scale)
        ## converged <- abs(logscale - logscale0) < rel.tol * max(rel.tol, logscale0)
        ## logscale0 <- logscale
    }

    if (it >= max.iter)
        warning("Sigma iterations did not converge. Returning unconverged estimates")

    ## cat("sigma:", scale, "\n")
    setSigma(object, scale)
    invisible(object)
}

##' Calculate sum(kappa) for calculation of theta
##'
##' @title sum of kappas
##' @param object rlmerMod-object
.sumKappa <- function(object) {
    if (all(object@dim == 1)) {
        ## diagonal case, use G function to calculate kappas
        ## this allows also wExp.b == 0
        a <- diag(object@pp$L)
        s <- .s(object, theta = TRUE)
        kappas <- G(,a,s,object@rho.b,object@rho.sigma.b,object@wExp.b,object@pp)
    } else {
        ## non-diag case
        kappas <- .Gnd(object)
    }
    lapply(seq_along(object@blocks),
           function(block) Reduce("+", kappas[object@ind == block]))
}

## define the objective function used in uniroot
## this is only external since we want to use it in tests
.defineFun <- function(wExp, rho, dk, us, skbsk) {
    function(thetak) sum(.wgtxy2(rho, dk/thetak, us/thetak, wExp)) - skbsk
}

.defineFunNonDiag <- function(wExp, rho, us, skbsk, par, parIdx) {
    ## FIXME: remove this condition
    if (wExp == 0) stop("wExp > 0 required")
    function(thetak) {
        par@x[] <- thetak
        ## no need to clear @factor here, since
        ## par is of class dtCMatrix, which does not have this slot.
        us <- solve(par, us)
        sqrtdk <- rep(sqrt(rho@wgt(sqrt(colMeans(us^2)))^wExp), nrow(us))
        ## TODO: do we need to take care of infinite values?
        (tcrossprod(sqrtdk * us) - skbsk)[parIdx]
    }       
}

##' Update theta: calculate DAS scale for random effects
##' from estimated random effects (diagonal covariance
##' matrices only). Expectation correction variant.
##'
##' @title Calculate variance components
##' @param object rlmerMod-object
##' @param max.iter maximum number of iterations to fit theta
##' @param rel.tol stopping criterion
##' @param verbose level of verbosity
##' @return rlmerMod-object
updateTheta <- function(object, max.iter = 100, rel.tol = 1e-6, verbose = 0) {
    skbs <- .sumKappa(object)

    deltatheta <- rep(0, length(theta(object)))
    sigmae <- object@pp$sigma
    ds <- .dk(object, sigmae)[object@k]
    idxTheta <- 0
    ## fit theta by block
    for (block in seq_along(object@blocks)) {
        ldim <- object@dim[block]
        lq <- object@q[block]
        lind <- object@ind[object@k] == block
        idxTheta <- max(idxTheta) + 1:(ldim*(ldim+1)/2)
        us <- object@pp$b.s[lind] / sigmae
        if (all(abs(us) < 1e-7)) {
            deltatheta[idxTheta] <- 0
        } else {
            if (ldim > 1) {
                ## non-diag case (use multiroot)
                us <- Matrix(us, ldim)
                par0 <- par <- diag(ldim)
                parIdx <- lower.tri(par0, TRUE)
                par[lower.tri(par)] <- 1 ## we will need the full lower triangular matrix
                par <- as(par, "dtCMatrix") ## convert to lower triangular matrix
                fun <- .defineFunNonDiag(object@wExp.b, object@rho.sigma.b, us, skbs[[block]],
                                         par, parIdx)
                res <- multiroot(fun, par0[parIdx], rtol = rel.tol,
                                 maxiter = max.iter, verbose = verbose > 5)
            } else {
                ## simple 1d case (use uniroot)
                fun <- .defineFun(object@wExp.b, object@rho.sigma.b, ds[lind], us,
                                  drop(skbs[[block]]))
                limit <- 50
                repeat {
                    res <- try(uniroot(fun, c(0, limit), tol = rel.tol, maxiter = max.iter),
                               silent=TRUE)
                    if (!is(res, "try-error")) break
                    if (limit > 10000) stop("cannot find reasonable root for deltatheta, stopping")
                    limit <- limit*2
                }
            }
            if (verbose > 4) str(res)
            deltatheta[idxTheta] <- res$root
        }
    }

    ## check for boundary conditions
    if (any(idx <- deltatheta < object@lower)) deltatheta[idx] <- 0
    
    setTheta(object, theta(object) * deltatheta, fit.effects = TRUE,
             update.sigma = FALSE, update.deviance = FALSE)
    ## update sigma without refitting effects
    updateSigma(object, fit.effects = FALSE)

    invisible(object)
}

##' Update theta: calculate DAS scale for random effects
##' from estimated random effects (diagonal covariance
##' matrices only). Tau variant.
##'
##' @title Calculate variance components
##' @param object rlmerMod-object
##' @param max.iter maximum number of iterations to fit theta
##' @param rel.tol stopping criterion
##' @param verbose level of verbosity
##' @return rlmerMod-object
updateThetaTau <- function(object, max.iter = 100, rel.tol = 1e-6, verbose = 0) {
    ## FIXME: non-diag case?
    kappa <- object@pp$kappa_b[1]
    tau <- switch(object@method,
                  DAStau = {
                      Tbk <- object@pp$T()
                      calcTau(diag(object@pp$L), .s(object, theta = TRUE),
                              object@rho.b, object@rho.sigma.b, object@wExp.b,
                              object@pp, kappa, sqrt(diag(Tbk)))
                  },
                  DASvar = sqrt(diag(object@pp$Tb())),
                  stop("method not supported by updateThetaTau:", object@method))

    deltatheta <- rep(0, length(theta(object)))
    sigmae <- object@pp$sigma
    ds <- .dk(object, sigmae)[object@k]
    ds <- ds / tau
    tau2 <- tau*tau
    idxTheta <- 0
    ## fit theta by block
    for (block in seq_along(object@blocks)) {
        ldim <- object@dim[block]
        lq <- object@q[block]
        lind <- object@ind[object@k] == block
        idxTheta <- max(idxTheta) + 1:(ldim*(ldim+1)/2)
        us <- object@pp$b.s[lind] / sigmae
        lds <- ds[lind]
        if (all(abs(us) < 1e-7)) {
            deltatheta[idxTheta] <- 0
        } else {
            if (ldim > 1) {
                ## non-diag case (use multiroot)
                stop("non diagonal case for DAStau not implemented yet")
            } else {
                ## simple 1d case
                delta0 <- 1
                converged <- FALSE
                it <- 0
                wgt <- .wgtTau(object@rho.sigma.b, object@wExp.b)
                while(!converged && (it <- it + 1) < max.iter) {
                    w <- wgt(lds/delta0)
                    delta <- sqrt(sum(w*us*us)/sum(w*tau2[lind])/kappa)
                    converged <- abs(delta - delta0) < rel.tol * max(rel.tol, delta0)
                    delta0 <- delta
                }
                
                if (it >= max.iter)
                    warning("theta iterations did not converge. Returning unconverged estimates")
            }
            if (verbose > 4) cat("delta: ", delta, "\n")
            deltatheta[idxTheta] <- delta
        }
    }

    ## check for boundary conditions
    if (any(idx <- deltatheta < object@lower)) deltatheta[idx] <- 0
    
    setTheta(object, theta(object) * deltatheta, fit.effects = TRUE,
             update.sigma = FALSE, update.deviance = FALSE)
    ## update sigma without refitting effects
    updateSigma(object, fit.effects = FALSE)
    ## set Tbk cache
    object@pp$setT(Diagonal(x=tau2))
    
    invisible(object)
}
