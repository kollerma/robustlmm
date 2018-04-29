xtraR <- system.file("xtraR", package="robustlmm")
source(file.path(xtraR, "unitTestBase.R"))
source(file.path(xtraR, "unitTestObjects.R"))

########################################################
## test wgt_e, wgt_b

test.wgt <- function(i) {
  cat("test.wgt for", names(objs)[i], "...")
  ltest <- function() {
    expected <- wgt.e(objs[[i]])
    actual <- robustlmm:::wgt_e(rPDASs[[i]]$.pointer, resps[[i]]$.pointer)
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 1e-12))

    expected <- wgt.b(objs[[i]])
    actual <- robustlmm:::wgt_b(rPDASs[[i]]$.pointer, resps[[i]]$.pointer)
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 1e-12))
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

for (i in seq_along(objs)) test.wgt(i)

########################################################
## test fit effects class

test.FitEffects <- function(i) {
  cat("test.FitEffects for", names(objs)[i], "...")
  object <- objs[[i]]

  n <- object@pp$n
  p <- object@pp$p
  q <- object@pp$q
  LtZt <- robustlmm:::.U_btZt.U_et(object)
  MAT1 <- Matrix(0, n+q, p+q)
  MAT1[1:n,1:p] <- robustlmm:::..U_eX(object)
  MAT1[1:n,p+1:q] <- t(LtZt)
  MAT1[n+1:q,p+1:q] <- sqrt(robustlmm:::.Lambda_b(object))
  MAT2 <- Matrix(0, p+q, n)
  MAT2[1:p,] <- t(robustlmm:::..U_eX(object))
  MAT2[p+1:q,] <- LtZt
  .U_ey <- solve(robustlmm:::.U_e(object), object@resp$y)
  sqW <- Diagonal(x=sqrt(c(w.e <- robustlmm:::wgt.e(object), robustlmm:::wgt.b(object))))

  actual <- new(robustlmm:::FitEffects, object@pp$.pointer,
                object@resp$.pointer, 1e-8, 50)

  stopifnot(all.equal(as.matrix(t(MAT1)), as.matrix(actual$mat1()),
                      tolerance = 1e-15, check.attributes = FALSE),
            all.equal(as.matrix(MAT2), as.matrix(actual$mat2()),
                      tolerance = 1e-15, check.attributes = FALSE),
            all.equal(as.matrix(.U_ey), as.matrix(actual$invU_ey()),
                      tolerance = 1e-15, check.attributes = FALSE),
            all.equal(diag(sqW)^2, actual$W(),
                      tolerance = 1e-15, check.attributes = FALSE)
  )

  cat(" matrices ok ")

  MAT <- sqW %*% MAT1
  w.y <- w.e * .U_ey
  expected <- drop(solve(crossprod(MAT), MAT2 %*% w.y))

  stopifnot(all.equal(expected, actual$nextValue(),
                      tolerance = 1e-12, check.attributes = FALSE))

  cat("ok\n")
}

for (i in seq_along(objs)) test.FitEffects(i)

########################################################
## test S

test.S <- function(i) {
  cat("test.S for", names(objs)[i],"...")
  ltest <- function() {
    actual <- objs[[i]]@pp$S_b()
    expected <- lapply(.S(objs[[i]]), as.matrix)
    stopifnot(all.equal(actual, expected, check.attributes = FALSE,
                        tolerance = 5e-6))
  }
  lsetTheta <- function(theta) {
    objs[[i]]@pp$setTheta(theta)
  }
  on.exit({
    lsetTheta(getME(fms[[i]], "theta"))
    gc()
  })

  ltest()

  testBlocks(fms[[i]], ltest, lsetTheta)

  cat("ok\n")
}

for (i in seq_along(objs)) test.S(i)

