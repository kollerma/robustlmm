.d <- function(bs, s=length(bs)) {
  if (s == 1) return(bs)
  if (is.matrix(bs)) rowSums(bs*bs) else sum(bs*bs)
}
.d2 <- function(sbs2, s) {
  if (s == 1) stop("s must be larger than 1") ## disable this test?
  sbs2
}
.Dd2 <- function(sbs2, s) {
  2.
}
.zeroB <- function(object, pp = object@pp)
  if (inherits(pp, "rlmerPredD_test")) pp$zeroB else pp$zeroB()

lchol <- function(x) {
  r <- try(chol.default(x), silent=TRUE)
  ## if chol fails, return sqrt of diagonal
  if (is(r, "try-error")) {
    Diagonal(x = sqrt(diag(x)))
  } else r
}

## calculate the expected value to be used in D.re
calcE.D.re <- function(s, rho, rel.tol = .Machine$double.eps^0.5) {
  if (s == 1) return(rho@EDpsi())
  tfun <- function(v) rho@Dwgt(.d2(v,s))*.Dd2(v,s)*v*dchisq(v, s)
  tfun2 <- function(v) rho@wgt(.d2(v,s))*dchisq(v, s)
  integrate(tfun, 0, Inf, rel.tol = rel.tol)$value / s +
      integrate(tfun2, 0, Inf, rel.tol = rel.tol)$value
}
calcE.D.re2 <- function(rho, s) calcE.D.re(s, rho)

calcE.psi_bbt <- function(rho, s, rel.tol = .Machine$double.eps^0.5) {
  if (s == 1) {
    Matrix(rho@EDpsi())
  } else {
    wgt <- rho@wgt
    fun <- function(v) wgt(.d2(v,s))*v*dchisq(v,s)
    tmp <- integrate(fun, 0, Inf, rel.tol = rel.tol)$value / s
    Diagonal(x=rep(tmp, s))
  }
}

calcE.psi_bpsi_bt <- function(rho, s, rel.tol = .Machine$double.eps^0.5) {
  if (s == 1) {
    Matrix(rho@Epsi2())
  } else {
    wgt <- rho@wgt
    fun <- function(v) wgt(.d2(v,s))^2*v*dchisq(v, s)
    tmp <- integrate(fun, 0, Inf, rel.tol = rel.tol)$value / s
    Diagonal(x=rep(tmp, s))
  }
}

## Calculate s = sqrt(v_{-i}), standard deviation of the
## average term in the linear approximation of the residuals,
## (in the case of diagonal Sigma_b / Lambda_\theta)
## This scale is used in G
##
## @title Calculate s
## @param object rlmerMod object
## @param theta logical if s is required for theta or sigma_e
## @return vector of s
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
  if (any(!.zeroB(pp = pp))) ret <- ret + drop(M2^2 %*% diag(pp$Epsi_bpsi_bt))
  sqrt(ret)
}

## Calculate S_k (for all ks), Cholesky decomposition of the
## covariance matrix of the average term in the linear approximation
## of the spherical random effects (in the case of non-diagonal
## Sigma_b / Lambda_\theta).
##
## @title Calculate all S_k
## @param object rlmerMod object
## @return list of matrices S_k
.S <- function(object) {
  ret <- list()
  idx <- !.zeroB(object)
  tmpEL <- object@pp$L() %*% object@pp$Epsi_bpsi_bt()
  for (k in 1:max(object@k)) {
    ind <- object@k == k
    ## check if the block was dropped (u == 0)
    if (any(idx[ind])) {
      ## Ltmp <- object@pp$L[ind, !ind, drop=FALSE]
      ## use symmpart to calm chol()
      ret <- c(ret, list(drop(lchol(symmpart(
        ## Matrix K_{\theta,k} consists of the
        ## rows belonging block k and all columns
        object@rho.e@Epsi2() * crossprod(object@pp$Kt()[, ind, drop=FALSE]) +
          ## Matrix L_{\theta,-k} consists of the
          ## rows belonging to block k and all but columns
          ## except the ones belonging to block k
          tcrossprod(object@pp$L()[ind, !ind, drop=FALSE],
                     tmpEL[ind, !ind, drop=FALSE])
        ## tcrossprod(Ltmp, Ltmp %*%
        ##            object@pp$Epsi_bpsi_bt[!ind, !ind, drop=FALSE])
      )))))
    } else {
      ## return a zero block of the correct size
      ret <- c(ret, list(matrix(0, sum(ind), sum(ind))))
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
    uniroot(tfun4, c(0, 1), tol = .Machine$double.eps^0.5)$root
  }
}

## @param object rlmerMod object
## @rdname calcKappa
calcKappaTauB <- function(object, pp = object@pp, rhoSigmaB = object@rho.sigma.b) {
  pp$btapply(rhoSigmaB, calcKappaTau)
}

##############################################
## other functions used for testing

.sigma <- robustlmm:::.sigma
uArranged <- robustlmm:::uArranged
rho.b <- robustlmm:::rho.b
rho.e <- robustlmm:::rho.e
b.s <- robustlmm:::b.s
.kappa_b <- robustlmm:::.kappa_b

## std.b: Return the spherical random effects or
##   "Standardize" the Matrix matrix: \eqn{\Lambda^{-1} matrix / \sigma}{Lambda^-1 matrix / sigma}
##
## @title Standardized values
## @param object rlmerMod object
## @param sigma to use for standardization
## @param matrix matrix to standardize
## @param drop apply drop to result?
## @param t transpose result
## @rdname std
std.b <- function(object, sigma = .sigma(object), matrix, drop=TRUE, t=FALSE)
  object@pp$stdB(sigma, matrix, drop, t)

## std.e: Calculate the standardized residuals or
##   "Standardize" the Matrix sigma: \eqn{R^{-1} matrix / \sigma}{R^-1 matrix / sigma}
##
## @rdname std
std.e <- function(object, sigma = .sigma(object), matrix, drop=TRUE) {
  if (missing(matrix)) return(object@resp$wtres / sigma)
  ## for the moment: just divide by sigma
  if (drop) matrix <- drop(matrix)
  matrix/sigma
}

## Calculate scaled squared Mahalanobis distances per group
##
## @title Distances per group
## @param object rlmerMod object
## @param sigma to use for standardization
## @param bs spherical random effects
## @param center whether to return the centered distances.
## @param ... ignored
.dk <- function(object, sigma, center, bs = b.s(object), ...) {
  ua <- uArranged(object, bs/sigma)
  unlist(lapply(seq_along(object@blocks), function(bt) {
    us <- ua[[bt]]
    s <- ncol(us)
    if (s == 1) return(us)
    ## else: square, sum and subtract s
    ret <- rowSums(us*us)
    if (center) ret <- ret - .kappa_b(object)[bt]*s
    ret
  }))
}

## dist.b: Calculate the distance from 0 standardized by sigma.
##   This is just value divided by sigma for uncorrelated
##   observations. For correlated items, this is the Mahalanobis
##   distance from 0. If \code{shifted=TRUE} and correlated items,
##   the squared distances are centered by -kappa_b*s. This is
##   required to compute the weights used for the size of the
##   covariance matrix of the random effects.
##
## @title Calculate distance
## @param object object to use
## @param sigma scale for standardization
## @param center whether to use the centered distances
## @param ... passed on to internal functions.
## @rdname dist
dist.b <- function(object, sigma = .sigma(object), center=FALSE, ...) {
  db <- .dk(object, sigma, center, ...)
  db[object@k]
}

## dist.e: Calculate dist for residuals
##   always assume they are uncorrelated
##
## @rdname dist
dist.e <- function(object, sigma = .sigma(object)) {
  std.e(object, sigma) ## just the usual rescaled residuals
}

## wgt.b: Calculate the robustness weights psi(d) / d,
##   standardized by sigma. The robustness weights are calculated
##   with d the Mahalanobis distance. Each group of correlated items
##   then gets a constant weight.
##   The robustness weights for the random effects themselves are
##   different than the ones used for estimating the size of the
##   covariance matrix of the random effects. Those are additionally
##   centered. That way, inlier can also be downweighted.
##   If \code{center=TRUE}, then the centered distances are used to
##   compute the robustness weights and the weight function given
##   by rho.sigma.b is used.
##
## @title Calculate robustness weights
## @param object object to use
## @param sigma scale for standardization
## @param center whether return the centered robustness weights, see Details.
## @rdname wgt
## @export
wgt.b <- function(object, sigma = .sigma(object), center = FALSE) {
  db <- dist.b(object, sigma, center)
  rho <- rho.b(object, if (center) "sigma" else "default")
  ret <- numeric()
  for (bt in seq_along(object@blocks)) {
    bind <- as.vector(object@idx[[bt]])
    ret <- c(ret, rho[[bt]]@wgt(db[bind]))
  }
  ret
}

## wgt.e: robustness weights of residuals
##
## @param use.rho.sigma return the weights computed using rho.sigma.
## @rdname wgt
## @export
wgt.e <- function(object, sigma = .sigma(object), use.rho.sigma = FALSE)
  if (use.rho.sigma) object@rho.sigma.e@wgt(dist.e(object, sigma)) else
    object@rho.e@wgt(dist.e(object, sigma))

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
  psiZ <- rho.e@psi(pp$ghZ)
  ## function to find root from
  fun <- function(tau, a, s) {
    t <- (pp$ghZ-a*psiZ-s*pp$ghZt)/tau
    numerator <- sum(rho.sigma.e@psi(t)*t*(tau*tau)*pp$ghW)
    denominator <- sum(rho.sigma.e@wgt(t)*pp$ghW)*kappa
    newValue = sqrt(numerator /denominator)
    ## cat(" Update sqrt(", numerator, "/", denominator, ") = ", newValue, "\n")
    newValue
  }

  for (i in seq_along(a)) {
    ## check if we already calculated this
    check <- FALSE
    a_i <- a[i]
    s_i <- s[i]
    j <- 0
    while ((j <- j + 1) < i) {
      ## if( (a_i == a[j]) && (s_i == s[j]) ) {
      if( abs(a_i - a[j]) < 1e-8 && abs(s_i - s[j]) < 1e-8 ) {
        check <- TRUE
        tau[i] <- tau[j]
        break
      }
    }

    if (!check) {
      ## cat(sprintf("calculating tau[%i]: a[%i] = %.16g, s[%i] = %.16g, init = %g...\n", i, i, a[i], i, s[i], tau[i]))
      it <- 0
      conv <- FALSE
      while(!conv && (it <- it + 1) < max.it && tau[i] <= 1) {
        ##print(tau[i])
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

calcTau.nondiag <- function(object, ghZ, ghw, skbs, kappas, max.iter,
                            rel.tol = 1e-4, verbose = 0, pp = object@pp) {
  ## define 4d integration function
  ## int4d <- function(fun) drop(apply(ghZ, 1, fun) %*% ghw)
  ## initial values
  T <- pp$T()
  ## move into list
  TkbsI <- list()
  ## cycle blocks
  for (type in seq_along(object@blocks)) {
    bidx <- object@idx[[type]]
    for (k in 1:ncol(bidx)) ## 1:K
      TkbsI <- c(TkbsI, list(as.matrix(T[bidx[,k],bidx[,k]])))
  }

  ## Compute Taus
  Tbks <- list()
  ## cycle block types
  for (type in seq_along(object@blocks)) {
    if (verbose > 5)
      cat("computing tau for blocktype", type, "...\n")
    bidx <- object@idx[[type]]
    s <- nrow(bidx)
    ind <- which(object@ind == type) ## (in terms of blocks on diagonal)
    idx <- !.zeroB(pp=pp)
    if (!any(idx[bidx])) { ## block has been dropped
      Tbks <- c(Tbks, rep(list(matrix(0, s, s)), length(ind)))
      next
    }
    ## block has not been dropped: calculate Tbk
    ## catch 1d case
    if (s == 1) {
      bidx <- drop(bidx) ## is a vector (in terms of items of b.s)
      tmp <- calcTau(diag(pp$L)[bidx], unlist(skbs[ind]),
                     object@rho.b[[type]], object@rho.sigma.b[[type]],
                     pp, kappas[type])
      Tbks <- c(Tbks, as.list(tmp*tmp))
    } else if (s == 2) { ## 2d case
      wgt <- object@rho.b[[type]]@wgt
      wgt.sigma <- object@rho.sigma.b[[type]]@wgt
      psi.sigma <- object@rho.sigma.b[[type]]@psi
      skappa <- s*pp$kappa_b[type]
      if (verbose > 6) {
        cat("type = ", type, "\n")
        cat("skappa = ", format(skappa, digits = 15), "\n")
      }
      wgtDelta <- function(u) (psi.sigma(u) - psi.sigma(u-skappa))/s
      lastSk <- lastLkk <- matrix()
      lastRet <- NA
      ## cycle blocks
      for (k in 1:ncol(bidx)) { ## 1:K
        lTbk <- as.matrix(TkbsI[[ind[k]]])
        lbidx <- bidx[,k]
        Lkk <- as.matrix(pp$L[lbidx, lbidx])
        Sk <- as.matrix(skbs[[ind[k]]])
        ## check for almost equal blocks
        ## if (!(isTRUE(all.equal(lastSk, Sk)) && isTRUE(all.equal(lastLkk, Lkk)))) {
        diff <- if (any(dim(lastSk) != dim(Sk))) 1 else
          abs(c(lastSk - Sk, lastLkk - Lkk))
        if (any(diff >= rel.tol * max(diff, rel.tol))) {
          ## Find Tbk for this new block...
          lastSk <- Sk
          lastLkk <- Lkk
          if (verbose > 5)
            cat("TbkI for k =", k, ":", lTbk, "\n")
          if (verbose > 6) {
            cat("Lkk = ", format(Lkk, digits = 15), "\n")
            cat("Sk = ", format(Sk, digits = 15), "\n")
          }
          conv <- FALSE
          iter <- 0
          while (!conv && (iter <- iter + 1) < max.iter) {
            lUTbk <- try(chol(lTbk), silent=TRUE)
            if (is(lUTbk, "try-error")) {
              warning("chol(lTbk) failed: ", lUTbk, "\nlTbk was: ", paste(lTbk),
                      "\nSetting it to Tb()\n")
              conv <- TRUE
              lTbk <- as.matrix(TkbsI[[ind[k]]])
            } else {
              ## browser()
              if (FALSE) {

                ## old:
                funA <- function(u) {
                  btilde <- u[1:2] - wgt(.d(u[1:2],2)) * Lkk %*% u[1:2] - crossprod(Sk, u[3:4])
                  wgtDelta(drop(crossprod(backsolve(lUTbk, btilde))))
                  ## not needed: *prod(dnorm(u))
                }
                a <- int4d(funA)
                funB <- function(u) {
                  btilde <- u[1:2] - wgt(.d(u[1:2],2)) * Lkk %*% u[1:2] -
                    crossprod(Sk, u[3:4])
                  wgt.sigma(drop(crossprod(backsolve(lUTbk, btilde))))*tcrossprod(btilde)
                }
                B <- matrix(int4d(funB), s)

              }
              btilde <- ghZ[,1:2] - wgt(.d(ghZ[,1:2],2)) * ghZ[, 1:2] %*% Lkk -
                ghZ[, 3:4] %*% Sk
              tmp1 <- colSums(backsolve(lUTbk, t(btilde))^2)
              tmp2 <- btilde[,1] * btilde[,2]
              a <- sum(wgtDelta(tmp1) * ghw)
              if (abs(a) < 1e-7) {
                a <- lasta
              } else {
                lasta <- a
              }
              B <- matrix(colSums(wgt.sigma(tmp1) * ghw *
                                    matrix(c(btilde[,1]*btilde[,1], tmp2, tmp2,
                                             btilde[,2]*btilde[,2]), length(ghw))),2)
              lTbk1 <- B/a
              ## conv <- isTRUE(all.equal(lTbk, lTbk1)) ## rel.tol default 1e-8
              diff <- abs(c(lTbk - lTbk1))
              conv <- all(diff < rel.tol * max(diff, rel.tol))
              if (verbose > 5) {
                cat(sprintf("k=%i, iter=%i, conv=%s\n", k, iter, conv))
                if (verbose > 6) {
                  cat("B:", B, "\n")
                  cat("a:", a, "\n")
                }
                cat("Tbk:", lTbk1, "\n")
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

  ## combine into Matrix
  out <- bdiag(Tbks)
  ## print(pp$Tb()[1:4, 1:4])
  ## print(out[1:4, 1:4])
  return(out)
}
