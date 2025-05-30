## Test disabled
quit()

## Test DAS-scale

require(robustlmm)

## Import functions
## rlmer
G <- robustlmm:::G
updateSigma <- robustlmm:::updateSigma
calcTau <- robustlmm:::calcTau
`theta<-` <- robustlmm:::`theta<-`
getX <- robustlmm:::getX
.s <- robustlmm:::.s
.fixef <- robustlmm:::.fixef
Lambda <- robustlmm:::Lambda
u <- robustlmm:::u
b <- robustlmm:::b
u.lmerMod <- robustlmm:::u.lmerMod
u.rlmerMod <- robustlmm:::u.rlmerMod
b.lmerMod <- robustlmm:::b.lmerMod
b.rlmerMod <- robustlmm:::b.rlmerMod
rho.e <- robustlmm:::rho.e
## lme4
nobars <- reformulas:::nobars
## robustbase
lmrob.hatmatrix <- robustbase:::lmrob.hatmatrix

## Test cases
rfm1 <- rlmer(Yield ~ 1 | Batch, Dyestuff, doFit=FALSE,
              rho.e = cPsi, rho.b = cPsi)
rfm2 <- rlmer(Yield ~ 1 | Batch, Dyestuff2, doFit=FALSE,
              rho.e = cPsi, rho.b = cPsi)
rfm3 <- rlmer(Reaction ~ Days + (Days|Subject), sleepstudy, doFit=FALSE,
              rho.e = cPsi, rho.b = cPsi)
sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])
rfm4 <- rlmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2, doFit=FALSE,
              rho.e = cPsi, rho.b = cPsi)

####
## Test row-wise computation of A
####

testCalculateA <- function(rfm) {
    U_eX <- rfm@pp$U_eX
    U_eZU_b <- rfm@pp$U_eZU_b
    M <- rfm@pp$M()
    groupsA <- 0:(nobs(rfm)-1)

    tmp1 <- U_eZU_b %*% M$M_bb
    tmp2 <- as.matrix(U_eZU_b %*% M$M_bB)
    tmp3 <- as.matrix(U_eX %*% M$M_BB)

    A <- tcrossprod(U_eX, tmp2) +
        tcrossprod(tmp2, U_eX) +
        tcrossprod(U_eX, tmp3) +
        tcrossprod(tmp1, U_eZU_b)

    expectedDiagA <- diag(A)
    expectedDiagAAt <- diag(tcrossprod(A))

    res <- robustlmm:::calculateA(as(as(U_eX, "generalMatrix"), "unpackedMatrix"),
                                  U_eZU_b,
                                  as(as(tmp1, "generalMatrix"), "unpackedMatrix"),
                                  as(as(M$M_bB, "generalMatrix"), "unpackedMatrix"),
                                  as(as(M$M_BB, "generalMatrix"), "unpackedMatrix"),
                                  as.integer(groupsA))

    stopifnot(all.equal(expectedDiagA, res$diagA, check.attributes = FALSE),
              all.equal(expectedDiagAAt, res$diagAAt, check.attributes = FALSE))
}
testCalculateA(rfm1)
testCalculateA(rfm2)
testCalculateA(rfm3)
testCalculateA(rfm4)

testA <- function(rfm) {
    A <- rfm@pp$A()
    rfm@pp$updateMatrices()
    stopifnot(all.equal(diag(A), rfm@pp$diagA, check.attributes = FALSE),
              all.equal(rowSums(A^2), rfm@pp$diagAAt, check.attributes = FALSE))
}
testA(rfm1)
testA(rfm2)
testA(rfm3)
testA(rfm4)

####
## Test row-wise computation of tau
####

testComputeDiagonalOfProduct <- function() {
    A <- Matrix(1:4, 2, 2)
    B <- Matrix(11:14, 2, 2)
    C <- Matrix(1:6, 2, 3)
    D <- Matrix(11:16, 3, 2)
    E <- Matrix(1:6, 2, 3)

    stopifnot(
        all.equal(diag(A %*% B), robustlmm:::computeDiagonalOfProduct(A, B)),
        all.equal(diag(C %*% D), robustlmm:::computeDiagonalOfProduct(C, D)),
        all.equal(diag(A %*% E), robustlmm:::computeDiagonalOfProduct(A, E)),
        all.equal(diag(D %*% E), robustlmm:::computeDiagonalOfProduct(D, E)),
        all.equal(diag(t(E) %*% A), robustlmm:::computeDiagonalOfProduct(t(E), A)))

    testthat::expect_error(robustlmm:::computeDiagonalOfProduct(C, B))
    testthat::expect_error(robustlmm:::computeDiagonalOfProduct(NULL, B))
    testthat::expect_error(robustlmm:::computeDiagonalOfProduct(A, NULL))
    testthat::expect_error(robustlmm:::computeDiagonalOfProduct(NULL, NULL))
    ## the following crashes R:
    ## testthat::expect_error(robustlmm:::computeDiagonalOfProduct())

    cat("testComputeDiagonalOfProduct passed\n")
}
testComputeDiagonalOfProduct()

testTau <- function(rfm) {
    fullA <- rfm@pp$A()
    B <- rfm@pp$B()
    Tau <- diag(rfm@pp$v_e) - rfm@pp$EDpsi_e * (t(fullA) + fullA) +
        rfm@pp$Epsi2_e * tcrossprod(fullA) +
        B %*% tcrossprod(rfm@pp$Epsi_bpsi_bt, B)
    target <- sqrt(diag(Tau))
    rfm@pp$method <- "DASvar"
    current <- rfm@pp$tau_e()
    rfm@pp$method <- rfm@method
    stopifnot(all.equal(target, current, check.attributes = FALSE))
}
testTau(rfm1)
testTau(rfm2)
testTau(rfm3)
testTau(rfm4)

####
## Test function G
####
G2 <- Vectorize(function(tau=1, a, s, rho, rho.sigma,
                        numpoints = 13) {
    gh <- robustbase:::ghq(numpoints)
    ghz <- gh$nodes
    ghw <- gh$weights*dnorm(gh$nodes)

    ## inner integration over z
    inner <- function(e) {
        x <- (e - a*rho@psi(e) - s*ghz)/tau
        sum(rho.sigma@wgt(x)*x^2*ghw)
    }
    ## outer integration over e
    sum(sapply(ghz, inner)*ghw)
}, vectorize.args = c("tau", "a", "s"))

G.int <- Vectorize(function(tau=1, a, s, rho, rho.sigma) {
    fun <- function(z, e) {
        x <- (e - a*rho@psi(e) - s*z)/tau
        rho.sigma@wgt(x)*x^2*dnorm(e)*dnorm(z)
    }
    inner <- Vectorize(function(e) integrate(fun, -Inf, Inf, e = e)$value)
    integrate(inner, -Inf, Inf)$value
}, vectorize.args = c("tau", "a", "s"))

testG <- function(object, theta=FALSE) {
    object@pp$updateMatrices()
    if (theta) {
        a <- diag(object@pp$L)
        s <- .s(object, theta = TRUE)
    } else {
        a <- diag(object@pp$A())
        s <- .s(object)
    }

    stopifnot(all.equal(G(,a,s,rho.e(object),rho.e(object, "sigma"),object@pp),
                        unname(G.int(,a,s,rho.e(object),rho.e(object, "sigma"))),
                        tolerance = 1e-4))
}

testG(rfm1)
testG(rfm2)
testG(rfm3)
testG(rfm4)
for (val in c(0, 0.01, 0.1, 0.5, 1, 2, 10, 100)) {
    theta(rfm1, fit.effects = TRUE) <- val
    testG(rfm1, theta=TRUE)
}

a.test <- seq(0, 1, length.out = 10)
s.test <- rep(a.test, each = 10)
a.test <- rep(a.test, 10)

tau.test <- seq(0, 1, length.out = 100)
A <- G(tau.test, a.test, s.test, smoothPsi, chgDefaults(smoothPsi, k=2, s=10), rfm3@pp)
B <- G2(tau.test, a.test, s.test, smoothPsi, chgDefaults(smoothPsi, k=2, s=10))
all.equal(A, B)

######
### DAS-theta
######

rfm1 <- update(rfm1, rho.e = smoothPsi, rho.b = smoothPsi)
rfm1@pp$updateMatrices()

.dk <- robustlmm:::.dk
dist.b <- robustlmm:::dist.b
## here dk and abs(b / theta) and abs(u) should be identical
stopifnot(all.equal(.dk(rfm1, 1, center=FALSE), unname(b(rfm1)) / diag(Lambda(rfm1))),
          all.equal(.dk(rfm1, 1, center=FALSE), unname(u(rfm1))))

### test DAS-theta for rfm4
rfm4 <- update(rfm4, rho.e = smoothPsi, rho.b = smoothPsi)

tmp <- rep(.dk(rfm4, 1, center=FALSE), each=2)[c(1:37, 40, 42, 44)]
stopifnot(all.equal(tmp, dist.b(rfm4, 1)))

####
## test optimality of sigma
## testing this only for the new lme4 version.
if (packageVersion("lme4") >= "0.99999911.0") {
    rfm <- rlmer(Yield ~ (1 | Batch), Dyestuff)
    fm <- lmer(Yield ~ (1 | Batch), Dyestuff)
    devfun <- lmer(Yield ~ (1 | Batch), Dyestuff, devFunOnly=TRUE)

    runTests <- function(rfm, ltheta) {
        robustlmm:::theta(rfm, fit.effects = TRUE) <- ltheta
        res <- devfun(ltheta)
        wrss <- environment(devfun)$resp$wrss()
        ussq <- environment(devfun)$pp$sqrL(1)
        sigma0 <- unname(sqrt((wrss + ussq) / fm@devcomp$dims['nmp']))
        print(c(sigma(rfm), sigma0))
        beta <- environment(devfun)$pp$beta(1)
        print(data.frame(.fixef(rfm), beta))
        u <- environment(devfun)$pp$u(1)
        print(data.frame(u(rfm), u))
        ## test for equality
        cat("effects difference:", all.equal(c(beta, u), c(.fixef(rfm), unname(u(rfm)))), "\n")
        stopifnot(all.equal(c(beta, u), c(.fixef(rfm), unname(u(rfm))),
                            tolerance = 1e-7),
                  all.equal(theta(rfm), ltheta),
                  all.equal(sigma(rfm), sigma0, tolerance = 1e-5))
    }


    ## test final theta
    rfm <- rlmer(Yield ~ (1 | Batch), Dyestuff, rho.e = cPsi, rho.b = cPsi)
    runTests(rfm, theta(rfm))

    ## test final theta, without init
    rfm <- rlmer(Yield ~ (1 | Batch), Dyestuff, rho.e = cPsi, rho.b = cPsi,
                 init = lmerNoFit(Yield ~ (1 | Batch), Dyestuff))
    runTests(rfm, theta(rfm))
}

####
## test find blocks function

findBlocksOld <- function(obj) {
    LambdaInd <- t(obj$Lambdat())
    LambdaInd@x[] <- (1:length(obj$theta))[obj$Lind]
    LambdaInd <- as(LambdaInd, "matrix") ## to avoid attempt to apply non function error
    bend <- unique(apply(LambdaInd != 0, 2, function(x) max(which(x))))
    nblocks <- length(bend)
    bstart <- c(1, bend[-nblocks]+1)
    bidx <- lapply(1:nblocks, function(i) seq.int(bstart[i],bend[i]))
    blocks <- lapply(bidx, function(idx) LambdaInd[idx,idx])
    bind <- match(blocks, ublocks <- unique(blocks))
    bdim <- sapply(bidx, length)
    q <- aggregate(bdim, list(bind), sum)$x
    k <- unlist(lapply(1:nblocks, function(i) rep(i, bdim[i])))
    list(blocks = ublocks, ind = bind, idx = bidx, dim = bdim, q = q, k = k)
}

findBlocks <- robustlmm:::findBlocks

testBlocks <- function(rfm) {
    blocks <- findBlocks(rfm@pp)
    blocksOld <- findBlocksOld(rfm@pp)

    blocks[["idx"]] <- blocksOld[["idx"]] <- NULL
    blocksOld[["dim"]] <- unique(blocksOld[["dim"]])

    stopifnot(all.equal(blocks, blocksOld))
}

testBlocks(rfm1)
testBlocks(rfm2)
testBlocks(rfm3)
testBlocks(rfm4)

####
## test .wgtxy2

.wgtxy2 <- robustlmm:::.wgtxy2
set.seed(0)
r1 <- rnorm(20)
r2 <- rnorm(20)

smoothProp2 <- psi2propII(smoothPsi)

stopifnot(all.equal(.wgtxy2(smoothPsi, r1, r1), smoothPsi@psi(r1)*r1),
          all.equal(.wgtxy2(smoothProp2, r1, r1), smoothPsi@psi(r1)^2),
          all.equal(.wgtxy2(smoothPsi, r1, r2), smoothPsi@wgt(r1)*r2*r2),
          all.equal(.wgtxy2(smoothProp2, r1, r2), smoothPsi@wgt(r1)^2*r2*r2))

####
## Test computeDiagonalOf(T)Crossproduct(Numeric)Matrix
####

testComputeDiagonalOfCrossproduct <- function() {
    mats <-
        list(
            Matrix(1:4, 2, 200000),
            Matrix(1:6, 200000, 3),
            Matrix(11:16, 3, 2),
            Matrix(1:4, 2, 2),
            Matrix(1:6, 2, 3),
            Matrix(11:16, 3, 2)
        )
    mats[[4]][1, 1] <- NA
    mats[[5]][1, 1] <- NaN
    mats[[6]][1, 1] <- Inf

    for (mat in mats) {
        numMat <- as.matrix(mat)
        stopifnot(
            all.equal(
                colSums(mat ^ 2, na.rm = TRUE),
                robustlmm:::computeDiagonalOfCrossproductMatrix(mat)
            ),
            all.equal(
                rowSums(mat ^ 2, na.rm = TRUE),
                robustlmm:::computeDiagonalOfTCrossproductMatrix(mat)
            ),
            all.equal(
                colSums(numMat ^ 2, na.rm = TRUE),
                robustlmm:::computeDiagonalOfCrossproductNumericMatrix(numMat)
            ),
            all.equal(
                rowSums(numMat ^ 2, na.rm = TRUE),
                robustlmm:::computeDiagonalOfTCrossproductNumericMatrix(numMat)
            )
        )
    }
}

testComputeDiagonalOfCrossproduct()

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
