xtraR <- system.file("xtraR", package="robustlmm")
source(file.path(xtraR, "unitTestBase.R"))
source(file.path(xtraR, "unitTestObjects.R"))

########################################################
## test s

test.s <- function(i) {
  cat("test.s for", names(rPDs)[i],"...")
  ltest <- function() {
    actual <- rPDASs[[i]]$s_e()
    expected <- .s(theta = FALSE, pp = rPDASTests[[i]])
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-6))
  }
  lsetTheta <- function(theta) {
    rPDASs[[i]]$setTheta(theta)
    rPDASTests[[i]]$setTheta(theta)
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

for (i in seq_along(rPDASs)) test.s(i)

########################################################
## test calcKappaTau

test.kappaTau <- function(i) {
  cat("test.kappaTau for", names(rPDs)[i],"...")
  ltest <- function() {
    expected <- rPDASTests[[i]]$kappa_e
    actual <- rPDASs[[i]]$kappa_e()
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-6))

    expected <- rPDASTests[[i]]$kappa_b
    actual <- rPDASs[[i]]$kappa_b()
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-5))
  }

  ltest()

  gc()
  cat("ok\n")
}

for (i in seq_along(rPDASs)) test.kappaTau(i)

########################################################
## test updateMatrices

test.updateMatrices <- function(i) {
  cat("test.updateMatrices for", names(rPDs)[i],"...")
  ltest <- function() {
    ## A
    expected <- as.matrix(rPDASTests[[i]]$A)
    actual <- rPDASs[[i]]$A()
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))

    ## Kt
    expected <- as.matrix(rPDASTests[[i]]$Kt)
    actual <- rPDASs[[i]]$Kt()
    stopifnot(all.equal(actual, expected, check.attributes = FALSE))

    ## L
    expected <- as.matrix(rPDASTests[[i]]$L)
    actual <- rPDASs[[i]]$L()
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 1e-6))
  }
  lsetTheta <- function(theta) {
    rPDASs[[i]]$setTheta(theta)
    rPDASTests[[i]]$setTheta(theta)
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

for (i in seq_along(rPDASs)) test.updateMatrices(i)

########################################################
## test tau_e

test.tau_e <- function(i) {
  cat("test.tau_e for", names(rPDs)[i], "...")
  ltest <- function() {
    Tau <- with(rPDASTests[[i]], V_e - EDpsi_e * (t(A) + A) + Epsi2_e * tcrossprod(A) +
                  B() %*% tcrossprod(Epsi_bpsi_bt, B()))
    expected <- sqrt(diag(Tau))
    actual <- rPDASs[[i]]$tau_e()
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-6))
  }
  lsetTheta <- function(theta) {
    rPDASs[[i]]$setTheta(theta)
    rPDASTests[[i]]$setTheta(theta)
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

for (i in seq_along(rPDASs)) test.tau_e(i)

########################################################
## test Tau_b

test.Tau_b <- function(i) {
  cat("test.Tau_b for", names(rPDs)[i], "...")
  ltest <- function() {
    expected <- as.matrix(rPDASTests[[i]]$Tb())
    actual <- as.matrix(rPDASs[[i]]$Tb())
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-5))
  }
  lsetTheta <- function(theta) {
    rPDASs[[i]]$setTheta(theta)
    rPDASTests[[i]]$setTheta(theta)
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

for (i in seq_along(rPDASs)) test.Tau_b(i)
