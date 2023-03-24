##' @importFrom stats dchisq integrate cov2cor as.formula asOneSidedFormula deviance
##' @importFrom stats dnorm getCall model.frame model.offset model.response
##' @importFrom stats model.weights napredict naresid ppoints printCoefmat
##' @importFrom stats qnorm quantile reorder setNames symnum terms.formula
##' @importFrom stats uniroot update.formula
##' @importFrom utils packageVersion stack
##' @importFrom methods as is slot validObject
##' @importFrom grDevices devAskNewPage
## calculate the expected value to be used in D.re
.calcE.D.re <- function(s, rho) {
    if (s == 1) return(rho@EDpsi())
    tfun <- function(v) rho@Dwgt(.d2(v,s))*.Dd2(v,s)*v*dchisq(v, s)
    tfun2 <- function(v) rho@wgt(.d2(v,s))*dchisq(v, s)
    integrate(tfun, 0, Inf)$value / s + integrate(tfun2, 0, Inf)$value
}
.calcE.D.re2 <- function(rho, s) .calcE.D.re(s, rho) ## switch arguments

##' @importFrom Matrix Matrix Diagonal symmpart

.calcE.psi_bbt <- function(rho, s) {
    if (s == 1) {
        Matrix(rho@EDpsi())
    } else {
        wgt <- rho@wgt
        fun <- function(v) wgt(.d2(v,s))*v*dchisq(v,s)
        tmp <- integrate(fun, 0, Inf)$value / s
        Diagonal(x=rep(tmp, s))
    }
}

.calcE.psi_bpsi_bt <- function(rho, s) {
    if (s == 1) {
        Matrix(rho@Epsi2())
    } else {
        wgt <- rho@wgt
        fun <- function(v) wgt(.d2(v,s))^2*v*dchisq(v, s)
        tmp <- integrate(fun, 0, Inf)$value / s
        Diagonal(x=rep(tmp, s))
    }
}

## Calculate the two dimensional integral G
##
## The \eqn{\rho}{rho}-functions can be chosen to be identical, but then the
## scale estimate will not have bounded influence (if the
## \eqn{\rho}{rho}-function is not redescending). As an alternative, one can use
## \eqn{\rho_\sigma = \psi_e'^2/2}{rho_sigma = psi_e'^2}. This corresponds to
## Huber's Proposal 2.
##
## @title Function \eqn{G_\sigma(a, s)}{G_sigma(a, s)}
## @param tau denominator inside function \eqn{\rho_\sigma}{rho_sigma}
## @param a first argument of G, leverage of observation
## @param s second argument of G, scale of \eqn{g_{[-i]}}{g_{-i}}
## @param rho psi_func-object
## @param rho.sigma psi_func-object (scale estimate)
## @param pp rlmerPredD object
G <- function(tau = rep(1, length(a)), a, s, rho, rho.sigma, pp) {
    if (length(a) != length(s))
        stop("length of vectors a and s not equal")

    ret <- numeric(length(a))
    psi <- rho.sigma@psi
    for (i in seq_along(a)) {
        x <- (pp$ghZ - a[i]*rho@psi(pp$ghZ) - s[i]*pp$ghZt)/tau[i]
        ret[i] <- sum(psi(x)*x*pp$ghW)
    }
    ret
}

## Calculate s = sqrt(v_{-i}), standard deviation of the average term in the
## linear approximation of the residuals, (in the case of diagonal Sigma_b /
## Lambda_\theta) This scale is used in G
##
## @title Calculate s
## @param object rlmerMod object
## @param theta logical if s is required for theta or sigma_e
## @return vector of s
.s <- function(object, theta = FALSE, pp = object@pp, B = pp$B()) {
    if (theta) {
        M1 <- pp$Kt
        M2 <- pp$L
        diag(M2) <- 0
    } else {
        M2 <- B
    }
    ## setting NA to 0 (they come from 0 variance components)
    if (any(naIdx <- is.na(M2))) {
        M2[naIdx] <- 0
    }
    ## calculate s:
    if (theta) {
        ret <- computeDiagonalOfCrossproductMatrix(M1)
    } else {
        ret <- pp$diagAAt - pp$diagA^2
    }
    ret <- pp$rho_e@Epsi2() * ret
    if (any(!.zeroB(pp=pp))) ret <- ret + drop(M2^2 %*% diag(pp$Epsi_bpsi_bt))
    sqrt(ret)
}

## Calculate S_k (for all ks), Cholesky decomposition of the covariance matrix
## of the average term in the linear approximation of the spherical random
## effects (in the case of non-diagonal Sigma_b / Lambda_\theta).
##
## @title Calculate all S_k
## @param object rlmerMod object
## @return list of matrices S_k
.S <- function(object) {
    ret <- vector("list", max(object@k))
    idx <- !.zeroB(object)
    tmpEL <- object@pp$L %*% object@pp$Epsi_bpsi_bt
    for (k in seq_along(ret)) {
        ind <- object@k == k
        notInd <- which(!ind)
        ind <- which(ind)
        ## check if the block was dropped (u == 0)
        if (any(idx[ind])) {
            ret[[k]] <- drop(lchol(
                ## Matrix K_{\theta,k} consists of the
                ## rows belonging block k and all columns
                object@rho.e@Epsi2() *
                    crossproductColumnSubMatrix(object@pp$Kt, ind) +
                    ## Matrix L_{\theta,-k} consists of the
                    ## rows belonging to block k and all but columns
                    ## except the ones belonging to block k
                    tCrossproductColumnRowSubMatrices(object@pp$L, tmpEL, ind, notInd)
            ))
        } else {
            ## return a zero block of the correct size
            ret[[k]] <- matrix(0, sum(ind), sum(ind))
        }
    }
    ret
}

## Calculate kappa for DASvar and DAStau method.
##
## @title Kappa for scale estimates
## @param rho rho function to be used
## @param s dimension of block
## @rdname calcKappa
calcKappaTau <- function(rho, s=1) {
    if (s > 38) return(1) ## otherwise there might be numerical problems
    psi <- rho@psi
    if (s == 1) {
        wgt <- rho@wgt
        tfun <- function(e) psi(e)*e*dnorm(e) ## do not subtract 1 here
        tfun2 <- function(e) wgt(e)*dnorm(e)  ## or use .d2, resp.
        integrate(tfun, -Inf, Inf)$value / integrate(tfun2, -Inf, Inf)$value
    } else {
        tfun3 <- function(v, kappa) psi(v-s*kappa)*dchisq(v, s)
        tfun4 <- function(kappa) integrate(tfun3, 0, Inf, kappa=kappa)$value
        tolerance <- sqrt(.Machine$double.eps)
        if (tfun4(1) > -tolerance) {
            return(1)
        }
        uniroot(tfun4, c(0, 1), tol = tolerance)$root
    }
}

## @param object rlmerMod object
## @rdname calcKappa
calcKappaTauB <- function(object) {
    object@pp$btapply(object@rho.sigma.b, calcKappaTau)
}

## Calculate tau for DAS estimate
##
## @title DAS-estimate
## @param a coefficient for psi-function in approximation
## @param s scale for rest in approximation
## @param rho.e rho-function used for effets
## @param rho.sigma.e rho-function used for sigma
## @param pp rlmerPredD object
## @param kappa kappa tuning parameter
## @param tau initial values for tau
## @param method method to compute tau with
## @param rel.tol relative tolerance for calculating tau
## @param max.it maximum number of iterations allowed

calcTau <- function(a, s, rho.e, rho.sigma.e, pp,
                    kappa, tau = rep(1, length(a)),
                    rel.tol = 1e-6, max.it = 200) {
    stopifnot(length(a) == length(s))
    psi <- rho.e@psi
    psiZ <- psi(pp$ghZ)
    ## function to find root from
    fun <- function(tau, a, s) {
        t <- (pp$ghZ-a*psiZ-s*pp$ghZt)/tau
        sqrt(sum(rho.sigma.e@psi(t)*t*(tau*tau)*pp$ghW) /
             (kappa*sum(rho.sigma.e@wgt(t)*pp$ghW)))
    }

    for (i in seq_along(a)) {
        ## check if we already calculated this
        check <- FALSE
        a_i <- a[i]
        s_i <- s[i]
        j <- 0
        while ((j <- j + 1) < i) {
            if( abs(a_i - a[j]) < 1e-8 && abs(s_i - s[j]) < 1e-8 ) {
                check <- TRUE
                tau[i] <- tau[j]
                break
            }
        }

        if (!check) {
            it <- 0
            conv <- FALSE
            while(!conv && (it <- it + 1) < max.it && tau[i] <= 1) {
                tau0 <- tau[i]
                tau[i] <- fun(tau0, a_i, s_i)
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

calcTau.nondiag <- function(object, ghZ12, ghZ34, ghw, skbs, kappas, max.iter,
                            rel.tol = 1e-4, verbose = 0) {
    ## initial values
    TkbsI <- object@pp$TList()

    ## Compute Taus
    Tbks <- list()
    ## cycle block types
    for (type in seq_along(object@blocks)) {
        if (verbose > 5)
            cat("computing tau for blocktype", type, "...\n")
        bidx <- object@idx[[type]]
        s <- nrow(bidx)
        ind <- which(object@ind == type) ## (in terms of blocks on diagonal)
        idx <- !.zeroB(object)
        if (!any(idx[bidx])) { ## block has been dropped
            Tbks <- c(Tbks, rep(list(matrix(0, s, s)), length(ind)))
            next
        }
        ## block has not been dropped: calculate Tbk
        ## catch 1d case
        if (s == 1) {
            bidx <- drop(bidx) ## is a vector (in terms of items of b.s)
            tmp <- calcTau(diag(object@pp$L)[bidx], unlist(skbs[ind]),
                           object@rho.b[[type]], object@rho.sigma.b[[type]],
                           object@pp, kappas[type])
            Tbks <- c(Tbks, as.list(tmp*tmp))
        } else if (s == 2) { ## 2d case
            wgt <- object@rho.b[[type]]@wgt
            wgt.sigma <- object@rho.sigma.b[[type]]@wgt
            psi.sigma <- object@rho.sigma.b[[type]]@psi
            skappa <- s*object@pp$kappa_b[type]
            wgtDelta <- function(u) (psi.sigma(u) - psi.sigma(u-skappa))/s
            lastSk <- lastLkk <- matrix()
            lastRet <- NA
            tmp0 <- wgt(.d(ghZ12,2)) * ghZ12
            ## cycle blocks
            for (k in 1:ncol(bidx)) { ## 1:K
                lTbk <- TkbsI[[ind[k]]]
                lbidx <- bidx[,k]
                Lkk <- as.matrix(object@pp$L[lbidx, lbidx])
                Sk <- as.matrix(skbs[[ind[k]]])
                ## check for almost equal blocks
                diff <- if (any(dim(lastSk) != dim(Sk))) 1 else
                    abs(c(lastSk - Sk, lastLkk - Lkk))
                if (any(diff >= rel.tol * max(diff, rel.tol))) {
                    ## Find Tbk for this new block...
                    lastSk <- Sk
                    lastLkk <- Lkk
                    if (verbose > 5)
                        cat("TbkI for k =", k, ":", lTbk, "\n")
                    btilde <- ghZ12 - tmp0 %*% Lkk - ghZ34 %*% Sk
                    btilde12 <- btilde[,1] * btilde[,2]
                    btildeSq <- btilde * btilde
                    conv <- FALSE
                    iter <- 0
                    while (!conv && (iter <- iter + 1) < max.iter) {
                        lLTbk <- try(chol(lTbk), silent=TRUE)
                        if (is(lLTbk, "try-error")) {
                            warning("chol(lTbk) failed: ", lLTbk, "\nlTbk was: ", paste(lTbk),
                                    "\nSetting it to Tb()\n")
                            conv <- TRUE
                            lTbk <- TkbsI[[ind[k]]]
                        } else {
                            tmp1 <- computeDiagonalOfCrossproductNumericMatrix(backsolve(lLTbk, t(btilde)))
                            a <- sum(wgtDelta(tmp1) * ghw)
                            if (abs(a) < 1e-7) {
                                a <- lasta
                            } else {
                                lasta <- a
                            }
                            tmp2 <- wgt.sigma(tmp1) * ghw
                            tmp3 <- sum(tmp2 * btilde12)
                            B <- matrix(c(sum(tmp2 * btildeSq[,1]),
                                          tmp3, tmp3,
                                          sum(tmp2 * btildeSq[,2])), 2)
                            lTbk1 <- B/a
                            diff <- abs(c(lTbk - lTbk1))
                            conv <- all(diff < rel.tol * max(diff, rel.tol))
                            if (verbose > 5) {
                                cat(sprintf("k=%i, iter=%i, conv=%s\n", k, iter, conv))
                                cat("Tbk:", lTbk1, "\n")
                                if (verbose > 6) {
                                    cat("B:", B, "\n")
                                    cat("a:", a, "\n")
                                    cat("lLTbk:", lLTbk, "\n")
                                    cat("btilde[1,]:", btilde[1,], "\n")
                                    cat("backsolve()[,1]:", backsolve(lLTbk, t(btilde))[,1], "\n")
                                    cat("fivenum(tmp1):", fivenum(tmp1), "\n")
                                    cat("wgtDelta(fivenum(tmp1)):", wgtDelta(fivenum(tmp1)), "\n")
                                    cat("fivenum(wgtDelta(tmp1)):", fivenum(wgtDelta(tmp1)), "\n")
                                }
                            }
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

    ## update cache
    object@pp$setTList(Tbks)

    ## combine into Matrix
    out <- .bdiag(Tbks)
    return(out)
}

## Update sigma: calculate averaged DAS scale
##
## @title Calculate averaged DAS scale
## @param object rlmerMod-object
## @param max.iter maximum iterations allowed in fitting
## @param rel.tol stopping criterion
## @param fit.effects fit effects for each new sigma?
## @return rlmerMod-object
updateSigma <- function(object, max.iter = 100, rel.tol = 1e-6, fit.effects = TRUE) {
    rho.sigma.e <- object@rho.sigma.e

    kappa <- .kappa_e(object)
    tau <- object@pp$tau_e()

    ## use iterative reweighting
    fun <- if (object@method %in% c("DAStau", "DASvar")) {
        wgt <- rho.sigma.e@wgt
        tau2 <- tau*tau
        function(scale, r) {
            lw <- wgt(r/tau/scale)
            sqrt(sum(lw*r*r)/sum(lw*tau2)/kappa)
        }
    }
    fun2 <- if (fit.effects) function(scale, r) {
        setSigma(object, scale)
        fitEffects(c(.fixef(object), .u(object)), object)
        fun(scale, object@resp$wtres)
    } else fun
    scale0 <- .sigma(object)
    converged <- FALSE
    it <- 0
    while(!converged && (it <- it + 1) < max.iter) {
        scale <- fun2(scale0, object@resp$wtres)
        converged <- abs(scale - scale0) < rel.tol * max(rel.tol, scale0)
        scale0 <- scale
    }

    if (it >= max.iter)
        warning("Sigma iterations did not converge. Returning unconverged estimates")

    ## cat("sigma:", scale, "\n")
    setSigma(object, scale)
    invisible(object)
}

## Calculate sum(kappa) for calculation of theta
##
## @title sum of kappas
## @param object rlmerMod-object
.sumKappa <- function(object) {
    if (all(object@dim == 1)) {
        ## diagonal case, use G function to calculate kappas
        a <- diag(object@pp$L)
        s <- .s(object, theta = TRUE)
        kappas <- numeric(0)
        for (bt in seq_along(object@blocks)) {
            bind <- as.vector(object@idx[[bt]])
            kappas <- c(kappas, G(,a[bind],s[bind],object@rho.b[[bt]],object@rho.sigma.b[[bt]],object@pp))
        }
    } else {
        ## non-diag case
        stop("Non diag case not supported for method", object@method)
    }
    lapply(seq_along(object@blocks),
           function(block) Reduce("+", kappas[object@ind == block]))
}

## Update theta: calculate DAS scale for random effects from estimated random
## effects (diagonal covariance matrices only). Tau variant.
##
## @title Calculate variance components
## @param object rlmerMod-object
## @param max.iter maximum number of iterations to fit theta
## @param rel.tol stopping criterion
## @param verbose level of verbosity
## @return rlmerMod-object
updateThetaTau <- function(object, max.iter = 100, rel.tol = 1e-6, verbose = 0) {
    kappas <- .kappa_b(object)
    tau <- switch(object@method,
                  DAStau = sqrt(diag(calcTau.nondiag(object,skbs=.s(object, theta=TRUE),
                                                     kappas=kappas, max.iter=max.iter))),
                  DASvar = sqrt(diag(object@pp$Tb())),
                  stop("method not supported by updateThetaTau:", object@method))

    deltatheta <- rep(0, length(theta(object)))
    sigmae <- .sigma(object)
    ds <- .dk(object, sigmae)[object@k]
    ## FIXME remove diag above and do stuff as matrix.
    ds <- ds / tau
    tau2 <- tau*tau
    idxTheta <- 0
    ## fit theta by block
    for (block in seq_along(object@blocks)) {
        ldim <- object@dim[block]
        lq <- object@q[block]
        lind <- object@ind[object@k] == block
        idxTheta <- max(idxTheta) + 1:(ldim*(ldim+1)/2)
        us <- b.s(object)[lind] / sigmae
        lds <- ds[lind]
        if (all(abs(us) < 1e-7)) {
            deltatheta[idxTheta] <- 0
        } else {
            if (ldim > 1) {
                ## non-diag case
                stop("non diagonal case for DAStau not implemented yet")
            } else {
                ## simple 1d case
                delta0 <- 1
                converged <- FALSE
                it <- 0
                wgt <- object@rho.sigma.b[[block]]@wgt
                while(!converged && (it <- it + 1) <= max.iter) {
                    w <- wgt(lds/delta0)
                    delta <- sqrt(sum(w*us*us)/sum(w*tau2[lind])/kappas[block])
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
             update.sigma = FALSE)
    ## update sigma without refitting effects
    updateSigma(object, fit.effects = FALSE)

    invisible(object)
}
