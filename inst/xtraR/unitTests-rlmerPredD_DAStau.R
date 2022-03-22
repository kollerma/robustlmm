xtraR <- system.file("xtraR", package="robustlmm")
source(file.path(xtraR, "unitTestBase.R"))
source(file.path(xtraR, "unitTestObjects.R"))

########################################################
## test ScalarTauParameters

testScalarTauParameters <- function(tol = 1e-8) {
  cat("testScalarTauParameters...")
  i <- 1
  a <- diag(rPDAStaus[[i]]$A())[1]
  s <- rPDAStaus[[i]]$s_e()[1]
  kappa <- rPDAStaus[[i]]$kappa_e()
  rho <- attr(rPDAStaus[[i]], "input")$rho_e_src
  rhoSigma <- attr(rPDAStaus[[i]], "input")$rhoSigma_e_src
  tau <- 0.8

  ltest <- function(x, y) {
    t <- (x-a*rho@psi(x)-s*y)/tau
    expectedNumerator <- rhoSigma@psi(t)*t*(tau*tau)
    expectedDenominator <- rhoSigma@wgt(t)

    actual <-
      robustlmm:::test_ScalarTauParameters(a[1], s[1], kappa,
                                           rho@getInstanceWithOriginalDefaults()$.pointer,
                                           rhoSigma@getInstanceWithOriginalDefaults()$.pointer,
                                           tau, x, y)

    stopifnot(all.equal(expectedNumerator, actual[1]),
              all.equal(expectedDenominator, actual[2]))
  }

  ltest(0, 0)
  ltest(-0.5, 0.5)
  ltest(-.25, -.1)
  cat("ok\n")
}

testScalarTauParameters()

########################################################
## test calcTau

testCalcTau <- function(i, nodes, tol = 1e-8) {
  cat("testCalcTau for", names(rPDAStaus)[i], "...")
  a <- diag(rPDAStaus[[i]]$A())[1]
  s <- rPDAStaus[[i]]$s_e()[1]
  kappa <- rPDAStaus[[i]]$kappa_e()
  rho <- attr(rPDAStaus[[i]], "input")$rho_e_src
  rhoSigma <- attr(rPDAStaus[[i]], "input")$rhoSigma_e_src

  expected <- calcTau(a, s, rho, rhoSigma, rPDASTests[[i]], kappa)
  fun <- function(j) {
    robustlmm:::test_calcTau(1, 1e-6, 200, a[j], s[j], kappa,
                             rho@getInstanceWithOriginalDefaults()$.pointer,
                             rhoSigma@getInstanceWithOriginalDefaults()$.pointer,
                             nodes)
  }
  actual <- sapply(seq_along(a), fun)
  stopifnot(all.equal(expected, actual, tolerance = tol))
  cat("ok\n")
}

for (i in seq_along(rPDAStaus)) testCalcTau(i, length(rPDASTests[[i]]$ghz))

########################################################
## test calcTauVectorized

testCalcTauVectorized <- function(i, nodes, tol = 1e-8) {
  cat("testCalcTauVectorized for", names(rPDAStaus)[i], "...")
  a <- diag(rPDAStaus[[i]]$A())
  s <- rPDAStaus[[i]]$s_e()
  kappa <- rPDAStaus[[i]]$kappa_e()
  rho <- attr(rPDAStaus[[i]], "input")$rho_e_src
  rhoSigma <- attr(rPDAStaus[[i]], "input")$rhoSigma_e_src

  expected <- calcTau(a, s, rho, rhoSigma, rPDASTests[[i]], kappa,
                      rPDASTests[[i]]$tau_e(), max.it = 200)
  actual <-
    robustlmm:::test_calcTauVectorized(rPDASTests[[i]]$tau_e(),
                                       1e-6, 200, a, s, kappa,
                                       rho@getInstanceWithOriginalDefaults()$.pointer,
                                       rhoSigma@getInstanceWithOriginalDefaults()$.pointer,
                                       nodes)
  stopifnot(all.equal(unname(expected), actual, tolerance = tol))
  cat("ok\n")
}

for (i in seq_along(rPDAStaus)) testCalcTauVectorized(i, length(rPDASTests[[i]]$ghz))

########################################################
## test DAStau_e()

test.DAStau_e <- function(i, tol = 5e-6) {
  cat("test.DAStau_e for", names(rPDAStaus)[i], "...")
  ltest <- function() {
    a <- diag(rPDASTests[[i]]$A)
    s <- unname(.s(theta=FALSE, pp=rPDASTests[[i]]))
    kappa <- rPDASTests[[i]]$kappa_e
    rho <- attr(rPDASs[[i]], "input")$rho_e_src
    rhoSigma <- attr(rPDASs[[i]], "input")$rhoSigma_e_src

    ## check inputs
    stopifnot(all.equal(a, diag(rPDAStaus[[i]]$A()),
                        tolerance = 5e-6, check.attributes = FALSE),
              all.equal(s, rPDAStaus[[i]]$s_e(),
                        tolerance = 5e-6, check.attributes = FALSE),
              all.equal(kappa, rPDAStaus[[i]]$kappa_e(),
                        tolerance = 5e-6, check.attributes = FALSE))

    expected <- calcTau(a, s, rho, rhoSigma, rPDASTests[[i]],
                        kappa, rPDASs[[i]]$tau_e(), "DAStau")
    actual0 <-
      robustlmm:::test_calcTauVectorized(rPDASs[[i]]$tau_e(),
                                         1e-6, 200, a, s, kappa,
                                         rho@getInstanceWithOriginalDefaults()$.pointer,
                                         rhoSigma@getInstanceWithOriginalDefaults()$.pointer,
                                         attr(rPDASs[[i]], "input")$nodes)
    stopifnot(all.equal(actual0, expected, check.attributes = FALSE,
                        tolerance = tol))
    actual <- rPDAStaus[[i]]$tau_e()
    stopifnot(all.equal(actual, actual0, check.attributes = FALSE,
                        tolerance = tol),
              all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = tol))
  }
  lsetTheta <- function(theta) {
    rPDAStaus[[i]]$setTheta(theta)
    rPDASs[[i]]$setTheta(theta)
    rPDASTests[[i]]$setTheta(theta)
  }
  on.exit({
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)
  cat("ok\n")
}

for (i in seq_along(rPDAStaus)) test.DAStau_e(i)

########################################################
## test MatrixTauParameters

test.MatrixTauParameters <- function(i) {
  cat("test.MatrixTauParameters for", names(rPDAStaus)[i], "...")
  lobj <- objs[[i]]

  ltest <- function() {
    pp <- lobj@pp
    type <- 1
    bidx <- lobj@idx[[type]]
    lbidx <- bidx[,1]
    skbs <- .S(lobj)
    s <- nrow(bidx)
    ind <- which(lobj@ind == type)
    wgt <- lobj@rho.b[[type]]@wgt
    wgt.sigma <- lobj@rho.sigma.b[[type]]@wgt
    psi.sigma <- lobj@rho.sigma.b[[type]]@psi
    wgtDelta <- function(u) (psi.sigma(u) - psi.sigma(u-skappa))/s

    args <- list(skappa = skappa <- s*pp$kappa_b()[type],
                 Lkk = Lkk <- as.matrix(pp$L()[lbidx, lbidx]),
                 Sk = Sk <- as.matrix(skbs[[ind[1]]]),
                 rho = lobj@rho.b[[type]]@getInstanceWithOriginalDefaults()$.pointer,
                 rhoSigma = lobj@rho.sigma.b[[type]]@getInstanceWithOriginalDefaults()$.pointer,
                 Tbk = lTbk <- as.matrix(pp$Tb()[lbidx,lbidx]),
                 u = c(0, 0, 0, 0))
    lLTbk <- chol(lTbk)

    expected <- function(u) {
      btilde <- u[1:2] - wgt(.d(u[1:2],2)) * Lkk %*% u[1:2] - crossprod(Sk, u[3:4])
      crossprodSolve <- drop(crossprod(backsolve(lLTbk, btilde)))
      list(ndim = nrow(bidx) * 2,
           fdimB = sum(lower.tri(lTbk, diag = TRUE)),
           btilde = drop(btilde),
           crossprodSolve = crossprodSolve,
           wgtDelta = wgtDelta(u[1]),
           funA = wgtDelta(crossprodSolve),
           funB = wgt.sigma(crossprodSolve)*tcrossprod(btilde)
      )
    }
    res <- robustlmm:::testMatrixTauParameters(args)
    stopifnot(all.equal(expected(args$u), res))

    args$u <- 1:4
    res <- robustlmm:::testMatrixTauParameters(args)
    stopifnot(all.equal(expected(args$u), res))

    args$u <- (1:4) / 10
    res <- robustlmm:::testMatrixTauParameters(args)
    stopifnot(all.equal(expected(args$u), res))
  }
  lsetTheta <- function(theta) {
    rPDAStaus[[i]]$setTheta(theta)
    rPDASTests[[i]]$setTheta(theta)
    lobj@pp$setTheta(theta)
  }
  on.exit({
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)
  if (i == 4) {
    lsetTheta(c(1.46725, 0.5210227, 0.1587546))
    ltest()
    lsetTheta(c(0.2533617, 0.6969634, 0.5566632 ))
    ltest()
  }
  cat("ok\n")
}

for (i in grep("sleepstudy", names(rPDAStaus))) test.MatrixTauParameters(i)

########################################################
## test Tau_b (DAStau)

test.Tau_b2 <- function(i) {
  cat("test.Tau_b2 (DAStau) for", names(rPDAStaus)[i], "...")
  lobj <- objs[[i]]
  ghZ <- as.matrix(expand.grid(rPDASTests[[i]]$ghz, rPDASTests[[i]]$ghz, rPDASTests[[i]]$ghz, rPDASTests[[i]]$ghz))
  ghw <- apply(as.matrix(expand.grid(rPDASTests[[i]]$ghw, rPDASTests[[i]]$ghw, rPDASTests[[i]]$ghw, rPDASTests[[i]]$ghw)), 1, prod)

  ltest <- function() {
    s <- .S(lobj)
    kappas <- .kappa_b(lobj)
    expected <- as.matrix(calcTau.nondiag(lobj, ghZ, ghw, s, kappas, 1000,
                                          1e-5, 0, pp = rPDASTests[[i]]))
    ## cat("Expected: \n")
    ## print(expected[1:2, 1:2])
    actual <- as.matrix(rPDAStaus[[i]]$Tb())
    ## cat("Actual: \n")
    ## print(actual[1:2, 1:2])
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-4))
  }
  lsetTheta <- function(theta) {
    ## cat("Set theta to ", theta, "\n")
    rPDAStaus[[i]]$setTheta(theta)
    rPDASTests[[i]]$setTheta(theta)
    lobj@pp$setTheta(theta)
  }
  on.exit({
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)
  if (i == 4) {
    lsetTheta(c(1.46725, 0.5210227, 0.1587546))
    ltest()
    lsetTheta(c(0.2533617, 0.6969634, 0.5566632 ))
    ltest()
  }
  cat("ok\n")
}

for (i in grep("sleepstudy", names(rPDAStaus))) test.Tau_b2(i)
