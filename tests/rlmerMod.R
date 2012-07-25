require(robustlmm)
require(numDeriv) ## for jacobian()

## errorHandler <- function() {
##     dump.frames("rlmerMod.dump", TRUE)
##     cat("Dumped to rlmerMod.dump\n")
##     print(sessionInfo())
##     quit(status=1)
## }
## options(error=errorHandler)

std.b <- robustlmm:::std.b
std.e <- robustlmm:::std.e
sumRho.b <- robustlmm:::sumRho.b
sumRho.e <- robustlmm:::sumRho.e
len <- robustlmm:::len
Lambda <- robustlmm:::Lambda
u <- robustlmm:::u
b <- robustlmm:::b
u.lmerMod <- robustlmm:::u.lmerMod
u.rlmerMod <- robustlmm:::u.rlmerMod
b.lmerMod <- robustlmm:::b.lmerMod
b.rlmerMod <- robustlmm:::b.rlmerMod

t.solve.omit.zeroes <- function(A, B) {
    ## omit zero rows of A
    A <- as.matrix(A)
    idx <- colSums(abs(A)) != 0
    if (any(idx) && !Matrix:::isSymmetric(A) && !Matrix:::isTriangular(A))
        stop("If there are zero columns in A, it must be symmetric or triangular")
    if (missing(B)) {
        ret <- A * 0
        if (any(idx))
            ret[idx,idx] <- solve(A[idx,idx])
    } else {
        B <- as.matrix(B)
        ret <- B * 0
        if (any(idx))
            ret[idx, ] <- solve(A[idx,idx], B[idx, ])
    }
    ## set infinite, NA, NaN terms to 0 (they correspond to dropped vc)
    if (any(idx <- !is.finite(ret)))
        ret[idx] <- 0
    ret
}

t.fun <- function(rfm, fm) {
    if (!missing(fm)) {
        ## test conversion features
        ## merMod -> rlmerMod -> merMod
        stopifnot(all.equal(fm, as(as(fm, 'rlmerMod'), 'lmerMod')))
        ## rlmerMod -> merMod -> rlmerMod
        ## we have to cheat a little
        conv <- as(as(rfm, 'lmerMod'), 'rlmerMod')
        conv@pp <- rfm@pp
        conv@method <- rfm@method
        conv@use.laplace <- rfm@use.laplace
        conv@wExp.b <- rfm@wExp.b
        conv@wExp.e <- rfm@wExp.e
        stopifnot(all.equal(rfm, conv))
    }
    
    rfm@pp$initMatrices(rfm)
    
    ## test standardization functions
    tmp1 <- as.matrix(std.b(rfm, sigma(rfm), rfm@pp$Zt))
    tmp2 <- as.matrix(std.b(rfm, sigma(rfm), Diagonal(len(rfm, "b"))))
    tmp1[is.infinite(tmp1)] <- 0
    tmp2[is.infinite(tmp2)] <- 0
    stopifnot(all.equal(std.b(rfm, sigma(rfm)),u(rfm)/sigma(rfm)),
              all.equal(std.e(rfm, sigma(rfm)), resid(rfm)/sigma(rfm)),
              all.equal(std.e(rfm, sigma(rfm), t(rfm@pp$X), drop=FALSE),
                        t(rfm@pp$X) / sigma(rfm)),
              all.equal(tmp1,
                        t.solve.omit.zeroes(sigma(rfm)*Lambda(rfm), rfm@pp$Zt)),
              all.equal(tmp2,
                        t.solve.omit.zeroes(sigma(rfm)*Lambda(rfm)),
                        check.attributes = FALSE))
    
    ## test sumRho functions: classical case
    ## FIXME: take correlation into account
    stopifnot(all.equal(sumRho.b(rfm, sigma(rfm), lambda = FALSE),
                        sum((u(rfm)/sigma(rfm))^2/2)),
              all.equal(sumRho.b(rfm, sigma(rfm), lambda = TRUE),
                        sum((u(rfm)/sigma(rfm))^2/2)),
              all.equal(sumRho.e(rfm, sigma(rfm), lambda = FALSE),
                        sum((resid(rfm)/sigma(rfm))^2/2)),
              all.equal(sumRho.e(rfm, sigma(rfm), lambda = TRUE),
                        sum((resid(rfm)/sigma(rfm))^2/2)))
    
    ## test determinant functions
    idx <- colSums(abs(Lambda(rfm))) != 0
    stopifnot(all.equal(robustlmm:::detCov.e(rfm),
                        robustlmm:::detCov.e(rfm, TRUE, sigma(rfm)),
                        2*length(rfm@resp$y)*log(sigma(rfm))),
              all.equal(robustlmm:::detCov.e(rfm, FALSE),
                        sigma(rfm)^(2*length(rfm@resp$y))),
              all.equal(robustlmm:::detCov.b(rfm),
                        if (any(idx)) 2*log(det(sigma(rfm)*
                                                Lambda(rfm)[idx,idx]))
                        else 0),
              all.equal(robustlmm:::detCov.b(rfm, FALSE),
                        if (any(idx)) det(sigma(rfm)*
                                          Lambda(rfm)[idx,idx])^2
                        else 0))

    ## test Q function (determinant part)
    Q <- robustlmm:::Q.rlmerMod(rfm, FALSE)
    if (!isREML(rfm)) {
        ## for ML case: only take determinant from subset
        p <- len(rfm, "beta")
        q <- sum(!robustlmm:::getZeroU(rfm))
        Q <- if (q > 0) Q[p+1:q, p+1:q] else Matrix(0, 1, 1)
    }
    if (isREML(rfm) || q > 0)
        stopifnot(all.equal(determinant(Q)$modulus, robustlmm:::Q.rlmerMod(rfm, TRUE, TRUE)),
                  all.equal(determinant(Q, FALSE)$modulus, robustlmm:::Q.rlmerMod(rfm, TRUE, FALSE)))
    
    ## test robll function:
    stopifnot(all.equal(robustlmm:::robll(rfm, norm.only = TRUE),
                        sumRho.e(rfm, sigma(rfm), lambda = TRUE) +
                        sumRho.b(rfm, sigma(rfm), lambda = TRUE)))

    ## test Laplace function:
    stopifnot(all.equal(robustlmm:::laplace.rlmerMod(rfm),
                        robustlmm:::laplace.rlmerMod(as(fm, "rlmerMod"))),
              all.equal(robustlmm:::laplace.rlmerMod(rfm),
                        unname(fm@devcomp$cmp[ifelse(isREML(rfm), "REML", "dev")])))
    
    ## check gradient function
    p <- len(rfm, "beta")
    q <- len(rfm, "b")
    rfm.tmp <- eval(rfm@call)
    u.tmp <- rep(1, q)
    u.tmp[robustlmm:::getZeroU(rfm.tmp)] <- 0
    robustlmm:::u(rfm.tmp) <- u.tmp
    robustlmm:::fixef(rfm.tmp) <- rep(1, p)
    opt.effects <- function(par) {
        robustlmm:::u(rfm.tmp) <- par[p+1:q]
        robustlmm:::fixef(rfm.tmp) <- par[1:p]
        robustlmm:::robll(rfm.tmp, norm.only=TRUE) ## calculate norm only
    }
    stopifnot(all.equal(unname(robustlmm:::gradient.rlmerMod(rfm.tmp)),
                        grad(opt.effects, c(rep(1, p), u.tmp))))
    
    return("Test passed")
}

update.merMod <- function(...) lme4:::update.merMod(...)

## one-way anova
fm1 <- lmer(Yield ~ (1 | Batch), Dyestuff)
rfm1 <- rlmer(Yield ~ (1 | Batch), Dyestuff, doFit = FALSE,
              rho.e = cPsi, rho.b = cPsi)
t.fun(rfm1, fm1)
t.fun(update.merMod(rfm1, REML=FALSE), update.merMod(fm1, REML=FALSE))

## model with 0 variance components
fm2 <- lmer(Yield ~ (1 | Batch), Dyestuff2)
rfm2 <- rlmer(Yield ~ (1 | Batch), Dyestuff2, doFit = FALSE,
              rho.e = cPsi, rho.b = cPsi)
t.fun(rfm2, fm2)
t.fun(update.merMod(rfm2, REML=FALSE), update.merMod(fm2, REML=FALSE))

## crossed random effects
fm3 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
rfm3 <- rlmer(diameter ~ (1|plate) + (1|sample), Penicillin, doFit = FALSE,
              rho.e = cPsi, rho.b = cPsi)
t.fun(rfm3, fm3)
t.fun(update.merMod(rfm3, REML=FALSE), update.merMod(fm3, REML=FALSE))

## correlated random effects
fm4 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
rfm4 <- rlmer(Reaction ~ Days + (Days|Subject), sleepstudy, doFit = FALSE,
              rho.e = cPsi, rho.b = cPsi)
t.fun(rfm4, fm4)
t.fun(update.merMod(rfm4, REML=FALSE), update.merMod(fm4, REML=FALSE))

sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])

## correlated random effects and zeroes in the variance components
fm5 <- lmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2)
rfm5 <- rlmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2, doFit = FALSE,
              rho.e = cPsi, rho.b = cPsi)
t.fun(rfm5, fm5)
t.fun(update.merMod(rfm5, REML=FALSE), update.merMod(fm5, REML=FALSE))


## test updating fixef
robustlmm:::fixef(rfm1) <- 1500
stopifnot(fixef(rfm1) == 1500,
          all(rfm1@resp$mu == 1500 + robustlmm:::getZ(rfm1) %*% b(rfm1)),
          all(rfm1@resp$wtres == rfm1@resp$y - rfm1@resp$mu)## ,
          ## all(resid(rfm1) == rfm1@resp$res * rfm1@resp$weights)
          )
robustlmm:::fixef(rfm4) <- c(250, 10)
stopifnot(fixef(rfm4) == c(250, 10),
          all(rfm4@resp$mu == robustlmm:::getX(rfm4) %*% c(250, 10) +
              robustlmm:::getZ(rfm4) %*% b(rfm4)),
          all(rfm4@resp$wtres == rfm4@resp$y - rfm4@resp$mu)## ,
          ## all(resid(rfm4) == rfm4@resp$res * rfm4@resp$weights)
          )

## test updating theta
robustlmm:::theta(rfm1, fit.effects = TRUE, update.sigma = FALSE) <- 2
stopifnot(theta(rfm1) == 2,
          diag(Lambda(rfm1)) == 2)
robustlmm:::theta(rfm4, fit.effects = TRUE, update.sigma = FALSE) <- c(1, 0, 2)
stopifnot(theta(rfm4) == c(1, 0, 2),
          diag(Lambda(rfm4)) == c(1, 2))
robustlmm:::theta(rfm5, fit.effects = TRUE, update.sigma = FALSE) <- c(1, 0, 0, 0)
stopifnot(theta(rfm5) == c(1, 0, 0, 0))
robustlmm:::theta(rfm5, fit.effects = TRUE, update.sigma = FALSE) <- c(1, 1, 0, 0)
stopifnot(theta(rfm5) == c(1, 1, 1e-7, 0))
robustlmm:::theta(rfm5, fit.effects = TRUE, update.sigma = FALSE) <- c(0, 1, 0, 1)
stopifnot(theta(rfm5) == c(0, 0, 0, 1))

## test updating u
## test integrity of rfm1:
stopifnot(all(Lambda(rfm1) %*% u(rfm1) == b(rfm1)),
          all(std.b(rfm1, 1, Matrix(b(rfm1))) == u(rfm1)))
rfm1a <- rfm1
rfm1b <- rfm1
robustlmm:::theta(rfm1a, fit.effects = TRUE, update.sigma = FALSE) <- 3
robustlmm:::theta(rfm1b, fit.effects = TRUE, update.sigma = FALSE) <- 3
## before setting u or b
stopifnot(all.equal(u(rfm1a), u(rfm1b)),
          all.equal(drop(Lambda(rfm1a) %*% u(rfm1a)), b(rfm1a)),
          all.equal(std.b(rfm1a, 1, Matrix(b(rfm1a))), u(rfm1a)),
          all.equal(b(rfm1a), b(rfm1b)),
          all.equal(drop(Lambda(rfm1b) %*% u(rfm1b)), b(rfm1b)),
          all.equal(std.b(rfm1b, 1, Matrix(b(rfm1b))), u(rfm1b)))
robustlmm:::u(rfm1a) <- u(rfm1)
robustlmm:::b(rfm1b) <- drop(Lambda(rfm1b) %*% u(rfm1))
## after:
stopifnot(all.equal(u(rfm1a), u(rfm1)),
          all.equal(b(rfm1b), drop(Lambda(rfm1b) %*% u(rfm1))),
          all.equal(u(rfm1a), u(rfm1b)),
          all.equal(drop(Lambda(rfm1a) %*% u(rfm1a)), b(rfm1a)),
          all.equal(std.b(rfm1a, 1, Matrix(b(rfm1a))), u(rfm1a)),
          all.equal(b(rfm1a), b(rfm1b)),
          all.equal(drop(Lambda(rfm1b) %*% u(rfm1b)), b(rfm1b)),
          all.equal(std.b(rfm1b, 1, Matrix(b(rfm1b))), u(rfm1b)))

## Check Namespace (otherwise R will fail to find them
## and there will be an "attempt to apply non-function" error)
print(selectMethod("t", "Cholesky"))
print(selectMethod("update", "dCHMsimpl"))

## check calculation of deviance when we let theta_i go to 0
rfm <- update.merMod(rfm1, REML=FALSE)
fm <- update.merMod(fm1, REML=FALSE)
i <- 1
theta <- theta(fm)
### compare with lme4:
devfun <- update(fm, devFunOnly=TRUE) 

tfun <- Vectorize(function(val) {
    theta[i] <- exp(-val)
    lmerval <- devfun(theta)
    robustlmm:::theta(rfm, fit.effects = TRUE) <- theta
    ## print(all.equal(environment(devfun)$pp$delb, unname(fixef(rfm))))
    ## print(all.equal(environment(devfun)$pp$delu, unname(u(rfm))))
    ## print(sumRho.b(rfm, sigma(rfm), lambda = TRUE))
    ## print(sumRho.e(rfm, sigma(rfm), lambda = TRUE))
    ## nmp <- nobs(rfm) - if(isREML(rfm)) len(rfm, "beta") else 0
    ## wrss <- environment(devfun)$resp$wrss()
    ## ussq <- environment(devfun)$pp$sqrL(1)
    ## sigma0 <- sqrt((wrss+ussq)/nmp)
    #print(all.equal(sigma0, sigma(rfm), tolerance=1e-5))
    c(robustlmm:::updateDeviance(rfm), lmerval)
})
(vals <- tfun(c(40:60, Inf)))
stopifnot(all.equal(vals[1,], vals[2,], tolerance=1e-3))
(vals <- tfun(c(1e-8, 0.001, 0.01, 0.1, 1, Inf)))
stopifnot(all.equal(vals[1,], vals[2,], tolerance=1e-3))

rfm <- rfm1
fm <- fm1
theta <- theta(fm)
devfun <- update(fm, devFunOnly=TRUE)
(vals <- tfun(c(40:60, Inf)))
stopifnot(all.equal(vals[1,], vals[2,]))
(vals <- tfun(c(1e-8, 0.001, 0.01, 0.1, 1, Inf)))
stopifnot(all.equal(vals[1,], vals[2,], tolerance=1e-4))

rfm <- update.merMod(rfm4, REML=FALSE)
fm <- update.merMod(fm4, REML=FALSE)
theta <- theta(fm)
i <- 3
### compare with lme4:
devfun <- update(fm, devFunOnly=TRUE)
(vals <- tfun(c(40:60, Inf)))
stopifnot(all.equal(vals[1,], vals[2,], tolerance=1e-3))
(vals <- tfun(c(1:20, Inf)))
stopifnot(all.equal(vals[1,], vals[2,], tolerance=1e-3))
(vals <- tfun(c(1e-8, 0.001, 0.01, 0.1, 1, Inf)))
## FIXME: this tolerance is really large...
stopifnot(all.equal(vals[1,], vals[2,], tolerance=1e-2))

rfm <- update.merMod(rfm5) #, REML=FALSE)
fm <- update.merMod(fm5) #, REML=FALSE)
theta <- theta(fm)
i <- 4
### compare with lme4:
devfun <- update(fm, devFunOnly=TRUE)
(vals <- tfun(c(40:60, Inf)))
stopifnot(all.equal(vals[1,], vals[2,], tolerance=1e-3))
(vals <- tfun(c(1e-8, 0.001, 0.01, 0.1, 1, Inf)))
stopifnot(all.equal(vals[1,], vals[2,], tolerance=1e-3))


## test dependency on initial values of theta
testInit <- function(formula, data, ...) {
    fm_1     <- lmerNoFit(formula, data)
    fm_10000 <- lmerNoFit(formula, data, initTheta = 10000)
    fm_0     <- lmerNoFit(formula, data, initTheta = 0)

    o1 <- capture.output(print(rlmer(formula, data, ..., init = fm_1)))
    o2 <- capture.output(print(rlmer(formula, data, ..., init = fm_10000)))
    o3 <- capture.output(print(rlmer(formula, data, ..., init = fm_0)))

    stopifnot(all.equal(o1, o2),
              all.equal(o1, o3))
}

testInit(Yield ~ (1 | Batch), Dyestuff)
testInit(Yield ~ (1 | Batch), Dyestuff, rho.e = smoothPsi, rho.b = smoothPsi)


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
