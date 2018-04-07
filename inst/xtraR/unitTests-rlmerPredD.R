xtraR <- system.file("xtraR", package="robustlmm")
source(file.path(xtraR, "unitTestBase.R"))
source(file.path(xtraR, "unitTestObjects.R"))

########################################################
## test init

testInit <- function(i) {
  cat("testInit for", names(rPDs)[i],"...")
  ltest <- function() {
    Lambdat <- as.matrix(rPDs[[i]]$Lambdat())
    expected <- as.matrix(rPDTests[[i]]$Lambdat())
    stopifnot(all.equal(Lambdat, expected, check.attributes = FALSE))

    theta <- rPDs[[i]]$theta()
    expected <- rPDTests[[i]]$theta
    stopifnot(all.equal(theta, expected, check.attributes = FALSE))

    n <- rPDs[[i]]$n
    expected <- rPDTests[[i]]$n
    stopifnot(all.equal(n, expected, check.attributes = FALSE))

    p <- rPDs[[i]]$p
    expected <- rPDTests[[i]]$p
    stopifnot(all.equal(p, expected, check.attributes = FALSE))

    actual <- rPDs[[i]]$q
    expected <- rPDTests[[i]]$q
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))

    U_e <- as.matrix(rPDs[[i]]$U_e())
    expected <- as.matrix(rPDTests[[i]]$U_e)
    stopifnot(all.equal(U_e, expected, check.attributes = FALSE))

    beta <- rPDs[[i]]$beta()
    expected <- rPDTests[[i]]$beta
    stopifnot(all.equal(beta, expected, check.attributes = FALSE))

    b_s <- rPDs[[i]]$b_s()
    expected <- rPDTests[[i]]$b.s
    stopifnot(all.equal(b_s, expected, check.attributes = FALSE))

    b <- rPDs[[i]]$b()
    expected <- rPDTests[[i]]$b.r
    stopifnot(all.equal(b, expected, check.attributes = FALSE))

    sigma <- rPDs[[i]]$sigma()
    expected <- rPDTests[[i]]$sigma
    stopifnot(all.equal(sigma, expected, check.attributes = FALSE))

    zeroB <- rPDs[[i]]$zeroB()
    expected <- rPDTests[[i]]$zeroB
    stopifnot(all.equal(zeroB, expected + 0, check.attributes = FALSE))
  }

  ltest()
  cat("ok\n")
}

for (i in seq_along(rPDs)) testInit(i)

########################################################
## test setTheta, Lambdat and invU_btZtU_et

testSetTheta <- function(i) {
  cat("testSetTheta for", names(rPDs)[i],"...")
  ltest <- function() {
    actual <- rPDs[[i]]$theta()
    expected <- rPDTests[[i]]$theta
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))

    actual <- as.matrix(rPDs[[i]]$Lambdat())
    expected <- as.matrix(rPDTests[[i]]$Lambdat())
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))

    actual <- as.matrix(rPDs[[i]]$invU_btZtU_et())
    expected <- as.matrix(rPDTests[[i]]$U_btZt.U_et)
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))
  }
  lsetTheta <- function(theta) {
    rPDs[[i]]$setTheta(theta)
    rPDTests[[i]]$setTheta(theta)
  }
  on.exit({
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltheta <- getME(fms[[i]], "theta")
  for (j in seq_along(ltheta)) {
    ltheta[j] <- j
    lsetTheta(ltheta)
    ltest()
  }
  cat("ok\n")
}

for (i in seq_along(rPDs)) testSetTheta(i)

########################################################
## test rlmerPredD::M() and M class

testM <- function(i) {
  cat("testM for", names(rPDs)[i],"...")
  ltest <- function() {
    actual <- rPDs[[i]]$M()
    expected <- rPDTests[[i]]$M()

    stopifnot(all.equal(as.matrix(rPDs[[i]]$U_b()), as.matrix(rPDTests[[i]]$U_b), check.attributes = FALSE),
              all.equal(actual$XZ(), as.matrix(expected$M_XZ), check.attributes = FALSE),
              all.equal(actual$bb(), as.matrix(expected$M_bb), check.attributes = FALSE),
              all.equal(actual$bB(), as.matrix(expected$M_bB), check.attributes = FALSE),
              all.equal(actual$BB(), as.matrix(expected$M_BB), check.attributes = FALSE),
              all.equal(actual$bbinv(), solve(as.matrix(expected$M_bb)), check.attributes = FALSE))
  }
  lsetTheta <- function(theta) {
    rPDs[[i]]$setTheta(theta)
    rPDTests[[i]]$setTheta(theta)
  }
  on.exit({
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)

  cat("ok\n")
}

for (i in seq_along(rPDs)) testM(i)

########################################################
## test setZeroB

testSetZeroB <- function(i) {
  cat("testSetZeroB for", names(rPDs)[i],"...")
  ltest <- function() {
    actual <- rPDs[[i]]$zeroB()
    expected <- rPDTests[[i]]$zeroB
    stopifnot(all.equal(actual, expected + 0, check.attributes = FALSE))
  }
  lsetTheta <- function(theta) {
    rPDs[[i]]$setTheta(theta)
    rPDTests[[i]]$setTheta(theta)
  }
  on.exit({
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)
  cat("ok\n")
}

for (i in seq_along(rPDs)) testSetZeroB(i)

########################################################
## test setB_s

testSetB_s <- function(i) {
  cat("testSetB_s for", names(rPDs)[i],"...")
  ltest <- function() {
    lsetB_s(rnorm(length(getME(fms[[i]], "u"))))
    actual <- rPDs[[i]]$b_s()
    expected <- rPDTests[[i]]$b.s
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))

    actual <- rPDs[[i]]$b()
    expected <- rPDTests[[i]]$b.r
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))
  }
  lsetB_s <- function(b_s) {
    rPDs[[i]]$setB_s(b_s)
    rPDTests[[i]]$setB.s(b_s)
  }
  lsetTheta <- function(theta) {
    rPDs[[i]]$setTheta(theta)
    rPDTests[[i]]$setTheta(theta)
  }
  on.exit({
    lsetB_s(getME(fms[[i]], "u"))
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)
  cat("ok\n")
}

for (i in seq_along(rPDs)) testSetB_s(i)

########################################################
## test setB_r

testSetB_r <- function(i) {
  cat("testSetB_r for", names(rPDs)[i],"...")
  ltest <- function() {
    vec <- rnorm(length(getME(fms[[i]], "u")))
    lsetB_r(vec)
    actual <- rPDs[[i]]$b_s()
    expected <- rPDTests[[i]]$b.s
    if (!isTRUE(all.equal(actual, expected, check.attributes = FALSE))) {
      cat("U_b = ")
      print(rPDTests[[i]]$U_b[1:2,1:2])
      cat("vec = ", vec[1:4], "\n")
      cat("actual =", actual[1:4], "\nexpected =", expected[1:4], "\n")
    }
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))

    actual <- rPDs[[i]]$b()
    expected <- rPDTests[[i]]$b.r
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))
  }
  lsetB_r <- function(b_r) {
    rPDs[[i]]$setB_r(b_r)
    rPDTests[[i]]$setB(b_r)
  }
  lsetTheta <- function(theta) {
    rPDs[[i]]$setTheta(theta)
    rPDTests[[i]]$setTheta(theta)
  }
  on.exit({
    lsetB_r(drop(getME(fms[[i]], "b")))
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)
  cat("ok\n")
}

for (i in seq_along(rPDs)) testSetB_r(i)

########################################################
## test initRho

testInitRho <- function(i) {
  cat("testInitRho for", names(rPDs)[i],"...")
  actual <- rPDs[[i]]$M_XZ0()
  expected <- rPDTests[[i]]$M_XZ0
  stopifnot(all.equal(actual, as.matrix(expected), check.attributes = FALSE))

  actual <- rPDs[[i]]$M_ZZ0()
  expected <- rPDTests[[i]]$M_ZZ0
  stopifnot(all.equal(as.matrix(actual), as.matrix(expected), check.attributes = FALSE))

  actual <- rPDs[[i]]$M_ZZ0_sub_M_ZX0invM_XXMZZ0()
  expected <- rPDTests[[i]]$M_ZZ0 - rPDTests[[i]]$M_ZX0M_XX.M_ZZ0
  stopifnot(all.equal(as.matrix(actual), as.matrix(expected), check.attributes = FALSE))

  actual <- rPDs[[i]]$Lambda_b()
  expected <- rPDTests[[i]]$Lambda_b
  stopifnot(all.equal(actual, as.matrix(expected), check.attributes = FALSE,
                      tolerance = 5e-6))

  actual <- rPDs[[i]]$Lambda_bD_b()
  expected <- rPDTests[[i]]$Lambda_bD_b
  stopifnot(all.equal(actual, as.matrix(expected), check.attributes = FALSE))

  actual <- rPDs[[i]]$Epsi_bbt()
  expected <- rPDTests[[i]]$Epsi_bbt
  stopifnot(all.equal(as.matrix(actual), as.matrix(expected),
                      check.attributes = FALSE, tolerance = 5e-7))

  actual <- rPDs[[i]]$Epsi_bpsi_bt()
  expected <- rPDTests[[i]]$Epsi_bpsi_bt
  stopifnot(all.equal(as.matrix(actual), as.matrix(expected),
                      check.attributes = FALSE, tolerance = 5e-6))

  actual <- rPDs[[i]]$Epsi2_b()
  expected <- rPDTests[[i]]$Epsi2_b
  stopifnot(all.equal(as.matrix(actual), as.matrix(expected),
                      check.attributes = FALSE, tolerance = 5e-6))

  gc()
  cat("ok\n")
}

for (i in seq_along(rPDs)) testInitRho(i)

########################################################
## test unsc

testUnsc <- function(i) {
  cat("testUnsc for", names(rPDs)[i],"...")
  ltest <- function() {
    actual <- rPDs[[i]]$unsc()
    expected <- rPDTests[[i]]$unsc()
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-5))
  }
  lsetTheta <- function(theta) {
    rPDs[[i]]$setTheta(theta)
    rPDTests[[i]]$setTheta(theta)
  }
  on.exit({
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)
  cat("ok\n")
}

for (i in seq_along(rPDs)) testUnsc(i)

########################################################
## test mu

test.mu <- function(i) {
  cat("test.mu for", names(rPDs)[i], "...")
  ltest <- function() {
    expected <- with(rPDTests[[i]], drop(crossprod(Zt, b.r) + (X %*% beta)))
    actual <- rPDs[[i]]$mu()
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-14))
  }
  lsetBBeta <- function(b, beta) {
    b <- drop(b)
    beta <- drop(beta)
    rPDTests[[i]]$setB(b)
    rPDTests[[i]]$beta <- beta
    rPDs[[i]]$setB_r(b)
    rPDs[[i]]$setBeta(beta)
  }
  on.exit({
    lsetBBeta(getME(fms[[i]], "b"), getME(fms[[i]], "beta"))
    gc()
  })

  ltest()
  lsetBBeta(getME(fms[[i]], "b") * 2, getME(fms[[i]], "beta") * 2)
  ltest()
  lsetBBeta(rnorm(getME(fms[[i]], "q")), rnorm(getME(fms[[i]], "p")))
  ltest()

  cat("ok\n")
}

for (i in seq_along(rPDs)) test.mu(i)

########################################################
## test distB

test.distB <- function(i) {
  cat("test.distB for", names(rPDs)[i], "...")

  idx <- attr(rPDs[[i]], "input")[["idx"]]
  ind <- attr(rPDs[[i]], "input")[["ind"]]
  dim <- attr(rPDs[[i]], "input")[["dim"]]
  k <- attr(rPDs[[i]], "input")[["k"]]
  dk <- unlist(lapply(robustlmm:::uArranged(b.s = with(rPDTests[[i]], b.s / sigma), idx = idx),
                      function(us) if (ncol(us) == 1) us else rowSums(us * us)))
  bidx <- ind %in% which(dim > 1)
  dk[bidx] <- sqrt(dk[bidx])
  expected <- dk[k]

  actual <- rPDs[[i]]$distB()

  stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                      tolerance = 5e-14))

  cat("ok\n")
}

for (i in seq_along(rPDs)) test.distB(i)

########################################################
## test effects

testEffects <- function(i) {
  cat("testEffects for", names(rPDs)[i],"...")
  ltest <- function() {
    actual <- rPDs[[i]]$effects()
    expected <- c(rPDTests[[i]]$beta, rPDTests[[i]]$b.s)
    idx <- c(rep(TRUE, rPDTests[[i]]$p), !rPDTests[[i]]$zeroB)
    expected[!idx] <- 0
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-5))
  }
  lsetTheta <- function(theta) {
    rPDs[[i]]$setTheta(theta)
    rPDTests[[i]]$setTheta(theta)
  }
  on.exit({
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)
  cat("ok\n")
}

for (i in seq_along(rPDs)) testEffects(i)
