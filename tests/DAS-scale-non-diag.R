## Test DAS-scale
## non diagonal variant

require(robustlmm)

## import functions
.s <- robustlmm:::.s
G <- robustlmm:::G
.S <- robustlmm:::.S
.Gnd <- robustlmm:::.Gnd
getZeroU <- robustlmm:::getZeroU

## test functions for classical case
.S.classic <- function(object) {
    ret <- list()
    idx <- !object@pp$zeroB
    for (k in 1:max(object@k)) {
        ind <- object@k == k
        ## check if the block was dropped (u == 0)
        if (all(idx[ind])) {
            ret <- c(ret, list(chol(tcrossprod(object@pp$K()[ind, , drop=FALSE]) +
                                    tcrossprod(object@pp$L[ind, !ind, drop=FALSE])
                                    )))
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
.Gnd.classic <- function(object, S) {
    Ekappa <- tcrossprod(Diagonal(object@pp$q) - object@pp$L) + tcrossprod(object@pp$K())
    lapply(1:max(object@k), function(k) {
        ind <- object@k == k
        Ekappa[ind, ind]
    })
}
.Gnd.classic.simplified <- function(object, S = .S(object)) {
    lapply(1:max(object@k), function(k) {
        ind <- object@k == k
        ## Omegak <- tcrossprod(object@pp$K()[ind, , drop=FALSE]) +
        ##     tcrossprod(object@pp$L[ind, !ind, drop=FALSE])
        Omegak <- crossprod(as.matrix(S[[k]]))
        tcrossprod(Diagonal(sum(ind)) - object@pp$L[ind, ind]) + Omegak
    })
}

## test function
testFun <- function(formula, data, seed) {
    lcall<- match.call()
    cat("\nTesting", lcall$data, "...\n")
    
    ## test classical case
    rfm <- rlmer(formula, data, method="DASexp", doFit=FALSE,
                 rho.e = cPsi, rho.b = cPsi)
    
    ## test .S function
    Sk <- .S(rfm)
    ESk <- .S.classic(rfm)
    stopifnot(all.equal(Sk, ESk,check.attr=FALSE))
    
    ## test kappa
    ## calculate expected value of tcrossprod(\tilde b.s) directly
    Ekappas <- lapply(.Gnd.classic(rfm), as.matrix)
    ## using the simplifed version
    Ekappas.simplified <- lapply(.Gnd.classic.simplified(rfm, Sk), as.matrix)
    stopifnot(all.equal(Ekappas, Ekappas.simplified, check.attr = FALSE, tol=1e-5))
    ## numerical version
    kappas <- .Gnd(rfm)
    ## FIXME: this should not fail, then add stopifnot
    stopifnot(all.equal(kappas, Ekappas.simplified, check.attr=FALSE, tol=1e-3),
              all.equal(kappas, Ekappas, check.attr=FALSE, tol=1e-3))
    
    ## test robust case
    rfm <- rlmer(formula, data, rho.e = smoothPsi, rho.b = smoothPsi,
                  wExp.b = 2, wExp.e = 2, method="DASexp", doFit = FALSE)
    
    ## build functions
    rho.e <- rho.e(rfm)
    rho.b <- rho.b(rfm)
    
    d <- function(u, k) sqrt(as.matrix(aggregate(u^2, list(k), mean)[,-1]))[k,]
    wgt <- function(u, k, rho) rho@wgt(d(u, k))
    utilde <- function(e, u, K, L, k, rho.e, rho.b)
        u - K %*% rho.e@psi(e) - L %*% (wgt(u, k, rho.b) * u)
    f.intfun <- function(lutilde, k, rho.b, wExp.b)
        tcrossprod(wgt(lutilde, k, rho.b)^wExp.b * lutilde, lutilde) / ncol(lutilde)
    f.integrand <- function(e, u, K, L, rho.e, rho.b, wExp.b, k) {
        lutilde <- utilde(e, u, K, L, k, rho.e, rho.b)
        f.intfun(lutilde, k, rho.b, wExp.b)
    }
    
    ## test functions
    q <- len(rfm, "b")
    n <- nobs(rfm)
    K <- as.matrix(rfm@pp$K())
    L <- as.matrix(rfm@pp$L)
    k <- rfm@k
    wExp.b <- wExp.b(rfm)
    rep <- 50
    set.seed(seed)
    e <- matrix(rnorm(rep*n), n)
    u <- matrix(rnorm(rep*q), q)
    ## vectorized utilde
    stopifnot(dim(utilde(e, u, K, L, k, rho.e, rho.b)) == c(q, rep))
    ## vectorized d
    stopifnot(dim(d(u, k)) == c(q, rep))
    ## vectorized wgt
    stopifnot(dim(wgt(u, k, rho.e)) == dim(d(u, k)))
    ## f.integrand should always be q x q
    stopifnot(dim(f.integrand(e, u, K, L, rho.e, rho.b, wExp.b, k)) == c(q, q))

    ## the integration function itself
    f.int <- function(object, rep) {
        q <- len(rfm, "b")
        n <- nobs(rfm)
        e <- matrix(rnorm(rep*n), n)
        u <- matrix(rnorm(rep*q), q)
        K <- as.matrix(rfm@pp$K())
        L <- as.matrix(rfm@pp$L)
        ret <- f.integrand(e, u, K, L, rho.e(object), rho.b(object), wExp.b(object), object@k)
        attr(ret, "e") <- e
        attr(ret, "u") <- u
        ret
    }

    ## Disabling the following tests to save time in R CMD check.
    
    ## ## --- calculate expectation using importance sampling
    ## test.1 <- f.int(rfm, 10000)
    
    ## ## --- calculate expectation using .S and .Gnd from rlmer
    ## res <- .Gnd(rfm)
    ## test.4 <- as.matrix(bdiag(res))

    ## ## --- compare blocks by block type
    ## ## drop other observations not on the "block diagonal"
    ## test.1[test.4 == 0] <- 0
    ## ## FIXME: increase accuracy
    ## print(all.equal(test.1, test.4, check.attr = FALSE))

    ## ## --- if diagonal: compare with G function
    ## if (all(rfm@dim) == 1) {
    ##     idx <- !getZeroU(rfm)
    ##     kappas <- rep(0, length(idx))
    ##     a <-  diag(rfm@pp$L)
    ##     s <- .s(rfm, theta = TRUE)
    ##     kappas[idx] <- G(,a,s,rho.b(rfm),rho.b(rfm, "sigma"),wExp.b(rfm),rfm@pp)
    ##     ## FIXME: increase accuracy and use stopifnot
    ##     print(all.equal(diag(test.4), kappas))
    ## }
    
    ## invisible(list(test.1, test.4))

    invisible(TRUE)
}


## ## test diagonal example first
## system.time(res1 <- testFun(Yield ~ (1 | Batch), Dyestuff, 1))

## ## body weight example
## data(BodyWeight, package='nlme')
## str(BodyWeight, give.attr = FALSE)
## system.time(res2 <- testFun(weight ~ Time * Diet + (1 + Time|Rat), BodyWeight, 2))

## sleepstudy example
system.time(res3 <- testFun(Reaction ~ Days + (Days|Subject), sleepstudy, 3))

## ## sleepstudy 2 example
## sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])
## system.time(res4 <- testFun(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2, 4))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
