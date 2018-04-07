xtraR <- system.file("xtraR", package="robustlmm")
source(file.path(xtraR, "unitTestBase.R"))
source(file.path(xtraR, "unitTestObjects.R"))

########################################################
## test mu
test.mu <- function(i) {
  cat("test.mu for", names(resps)[i], "...")
  actual <- resps[[i]]$mu
  expected <- getME(fms[[i]], "mu")
  all.equal(expected, actual)
  cat("ok\n")
}

for (i in seq_along(resps)) test.mu(i)

########################################################
## test offset
test.offset <- function(i) {
  cat("test.offset for", names(resps)[i], "...")
  actual <- resps[[i]]$offset
  expected <- getME(fms[[i]], "offset")
  all.equal(expected, actual)
  cat("ok\n")
}

for (i in seq_along(resps)) test.offset(i)

########################################################
## test sqrtrwt
test.sqrtrwt <- function(i) {
  cat("test.sqrtrwt for", names(resps)[i], "...")
  actual <- resps[[i]]$sqrtrwt
  expected <- fms[[i]]@resp$sqrtrwt
  all.equal(expected, actual)
  cat("ok\n")
}

for (i in seq_along(resps)) test.sqrtrwt(i)

########################################################
## test weights
test.weights <- function(i) {
  cat("test.weights for", names(resps)[i], "...")
  actual <- resps[[i]]$weights
  expected <- weights(fms[[i]])
  all.equal(expected, actual)
  cat("ok\n")
}

for (i in seq_along(resps)) test.weights(i)

########################################################
## test wtres
test.wtres <- function(i) {
  cat("test.wtres for", names(resps)[i], "...")
  actual <- resps[[i]]$wtres
  expected <- fms[[i]]@resp$wtres
  all.equal(expected, actual)
  cat("ok\n")
}

for (i in seq_along(resps)) test.wtres(i)

########################################################
## test y
test.y <- function(i) {
  cat("test.y for", names(resps)[i], "...")
  actual <- resps[[i]]$y
  expected <- getME(fms[[i]], "y")
  all.equal(expected, actual)
  cat("ok\n")
}

for (i in seq_along(resps)) test.y(i)

########################################################
## test updateMu
test.updateMu <- function(i) {
  cat("test.updateMu for", names(resps)[i], "...")
  lsetMu <- function(mu) {
    resps[[i]]$updateMu(mu)
    fms[[i]]@resp$updateMu(mu)
  }
  ltest <- function(mu) {
    lsetMu(mu)
    actual <- resps[[i]]$wtres
    expected <- fms[[i]]@resp$wtres
    all.equal(expected, actual)
  }
  oldMu <- getME(fms[[i]], "mu")
  on.exit(lsetMu(oldMu))

  ltest(getME(fms[[i]], "mu") * 2)
  ltest(rnorm(getME(fms[[i]], "n")))
  ltest(oldMu)

  cat("ok\n")
}

for (i in seq_along(resps)) test.updateMu(i)
