xtraR <- system.file("xtraR", package="robustlmm")
source(file.path(xtraR, "unitTestBase.R"))

########################################################
## test calcED_re routine
testCalcED_re <- function(psiFunc) {
  instance <- psiFunc@getInstanceWithOriginalDefaults()
  cat("testCalcED_re for", instance$show(), "...")
  stopifnot(all.equal(calcE.D.re(1, psiFunc),
                      robustlmm:::calcED_re(instance$.pointer, 1), tolerance = 1e-8),
            all.equal(calcE.D.re(2, psiFunc),
                      robustlmm:::calcED_re(instance$.pointer, 2), tolerance = 1e-5),
            all.equal(calcE.D.re(3, psiFunc),
                      robustlmm:::calcED_re(instance$.pointer, 3), tolerance = 1e-5))
  cat("ok\n")
}

testForAllPsiFunc(testCalcED_re)

########################################################
## test calcEpsi_bbt routine
testCalcEpsi_bbt <- function(psiFunc) {
  instance <- psiFunc@getInstanceWithOriginalDefaults()
  cat("testCalcEpsi_bbt for", instance$show(), "...")
  stopifnot(all.equal(as.matrix(calcE.psi_bbt(psiFunc, 1)),
                      robustlmm:::calcEpsi_bbt(instance$.pointer, 1),
                      tolerance = 1e-8, check.attributes = FALSE),
            all.equal(as.matrix(calcE.psi_bbt(psiFunc, 2)),
                      robustlmm:::calcEpsi_bbt(instance$.pointer, 2),
                      tolerance = 1e-5, check.attributes = FALSE),
            all.equal(as.matrix(calcE.psi_bbt(psiFunc, 3)),
                      robustlmm:::calcEpsi_bbt(instance$.pointer, 3),
                      tolerance = 1e-5, check.attributes = FALSE))
  cat("ok\n")
}

testForAllPsiFunc(testCalcEpsi_bbt)

########################################################
## test calcEpsi_bpsi_bt routine
testCalcEpsi_bpsi_bt <- function(psiFunc) {
  instance <- psiFunc@getInstanceWithOriginalDefaults()
  cat("testCalcEpsi_bpsi_bt for", instance$show(), "...")
  stopifnot(all.equal(as.matrix(calcE.psi_bpsi_bt(psiFunc, 1)),
                      robustlmm:::calcEpsi_bpsi_bt(instance$.pointer, 1),
                      tolerance = 1e-8, check.attributes = FALSE),
            all.equal(as.matrix(calcE.psi_bpsi_bt(psiFunc, 2)),
                      robustlmm:::calcEpsi_bpsi_bt(instance$.pointer, 2),
                      tolerance = 1e-5, check.attributes = FALSE),
            all.equal(as.matrix(calcE.psi_bpsi_bt(psiFunc, 3)),
                      robustlmm:::calcEpsi_bpsi_bt(instance$.pointer, 3),
                      tolerance = 1e-5, check.attributes = FALSE))
  cat("ok\n")
}

testForAllPsiFunc(testCalcEpsi_bpsi_bt)

########################################################
## test GaussHermiteQuadrature

testGH <- function(func, nNodes = 20, tol = 1e-5) {
  cat("testGH for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_GaussHermiteQuadrature(nNodes, func)
  if (doBench) stop <- get_nanotime()
  rule <- fastGHQuad::gaussHermiteData(nNodes)
  expected <- fastGHQuad::ghQuad(function(x) func(x) * exp(x^2), rule)
  stopifnot(all.equal(expected, actual, tolerance = tol))
  cat(".")
  expected <- integrate(func, -Inf, Inf, rel.tol = tol / 10)$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testGH(identity)
testGH(dnorm)
testGH(function(x) x^2 * dnorm(x))
testGH(function(x) x^3 * dnorm(x))

testGH2 <- function(func)
  testGH(func, nNodes = 50, tol = 1e-12)

testGH2(identity)
testGH2(dnorm)
testGH2(function(x) x^2 * dnorm(x))
testGH2(function(x) x^3 * dnorm(x))

########################################################
## test DqagNormalExpectation

testDNExpectation <- function(func, tol = 1e-8) {
  cat("testDNExpectation for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_DqagNormalExpectation(func)
  if (doBench) stop <- get_nanotime()
  expected <- integrate(function(x) func(x)*dnorm(x), -Inf, Inf)$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testDNExpectation(function(x) rep.int(0., length(x)))
testDNExpectation(identity)
testDNExpectation(function(x) x^2)
testDNExpectation(function(x) x^3)

########################################################
## test GaussHermiteNormalExpectation

testGHExpectation <- function(func, nNodes = 20, tol = 1e-5) {
  cat("testGHExpectation for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_GaussHermiteNormalExpectation(nNodes, func)
  if (doBench) stop <- get_nanotime()
  expected <- integrate(function(x) func(x)*dnorm(x),
                        -Inf, Inf, rel.tol = tol / 100)$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testGHExpectation(function(x) rep.int(0., length(x)))
testGHExpectation(identity)
testGHExpectation(function(x) x^2)
testGHExpectation(function(x) x^3)

testGHExpectation2 <- function(func)
  testGHExpectation(func, nNodes = 50, tol = 1e-10)

testGHExpectation2(function(x) rep.int(0., length(x)))
testGHExpectation2(identity)
testGHExpectation2(function(x) x^2)
testGHExpectation2(function(x) x^3)

########################################################
## test DqagIntegration2d_ninfInf

testDqagIntegration2d_ninfInf <- function(func, tol = 1e-8) {
  cat("testDqagIntegration2d_ninfInf for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_DqagIntegration2d_ninfInf(func)
  if (doBench) stop <- get_nanotime()
  inner <- function(y) {
    sapply(y, function(yi) integrate(func, -Inf, Inf, y = yi)$value)
  }
  expected <- integrate(inner, -Inf, Inf)$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testDqagIntegration2d_ninfInf(function(x, y) rep.int(0., length(x)))
testDqagIntegration2d_ninfInf(function(x, y) dnorm(x) * dnorm(y))
testDqagIntegration2d_ninfInf(function(x, y) x^2 * dnorm(x) * dnorm(y))
testDqagIntegration2d_ninfInf(function(x, y) y^2 * dnorm(x) * dnorm(y))
testDqagIntegration2d_ninfInf(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y))

########################################################
## test DqagIntegration2d_aInf

testDqagIntegration2d_aInf <- function(func, a, tol = 1e-8) {
  cat("testDqagIntegration2d_aInf for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_DqagIntegration2d_aInf(func, a)
  if (doBench) stop <- get_nanotime()
  inner <- function(y) {
    sapply(y, function(yi) integrate(func, a[1], Inf, y = yi)$value)
  }
  expected <- integrate(inner, a[2], Inf)$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testDqagIntegration2d_aInf(function(x, y) rep.int(0., length(x)), c(0, 0))
testDqagIntegration2d_aInf(function(x, y) dnorm(x) * dnorm(y), c(0, 0))
testDqagIntegration2d_aInf(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, 0))
testDqagIntegration2d_aInf(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, 0))
testDqagIntegration2d_aInf(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(0, 0))

testDqagIntegration2d_aInf(function(x, y) rep.int(0., length(x)), c(-3, 0))
testDqagIntegration2d_aInf(function(x, y) dnorm(x) * dnorm(y), c(-3, 0))
testDqagIntegration2d_aInf(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, -3)) ## try c(-20, 0)
testDqagIntegration2d_aInf(function(x, y) y^2 * dnorm(x) * dnorm(y), c(-3, 0))
testDqagIntegration2d_aInf(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(-3, 0))

########################################################
## test DqagIntegration2d_ninfB

testDqagIntegration2d_ninfB <- function(func, b, tol = 1e-7) {
  cat("testDqagIntegration2d_ninfB for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_DqagIntegration2d_ninfB(func, b)
  if (doBench) stop <- get_nanotime()
  inner <- function(y) {
    sapply(y, function(yi) integrate(func, -Inf, b[1], y = yi)$value)
  }
  expected <- integrate(inner, -Inf, b[2])$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testDqagIntegration2d_ninfB(function(x, y) rep.int(0., length(x)), c(0, 0))
testDqagIntegration2d_ninfB(function(x, y) dnorm(x) * dnorm(y), c(0, 0))
testDqagIntegration2d_ninfB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, 0))
testDqagIntegration2d_ninfB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, 0))
testDqagIntegration2d_ninfB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(0, 0))

testDqagIntegration2d_ninfB(function(x, y) rep.int(0., length(x)), c(3, -1))
testDqagIntegration2d_ninfB(function(x, y) dnorm(x) * dnorm(y), c(3, -1))
testDqagIntegration2d_ninfB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, -1))
testDqagIntegration2d_ninfB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, -1))
testDqagIntegration2d_ninfB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(3, -1))

########################################################
## test DqagIntegration2d_aB

testDqagIntegration2d_aB <- function(func, a, b, tol = 1e-8) {
  cat("testDqagIntegration2d_aB for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_DqagIntegration2d_aB(func, a, b)
  if (doBench) stop <- get_nanotime()
  inner <- function(y) {
    sapply(y, function(yi) integrate(func, a[1], b[1], y = yi)$value)
  }
  expected <- integrate(inner, a[2], b[2])$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testDqagIntegration2d_aB(function(x, y) rep.int(0., length(x)), c(0, 0), c(10, 10))
testDqagIntegration2d_aB(function(x, y) dnorm(x) * dnorm(y), c(0, 0), c(10, 10))
testDqagIntegration2d_aB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, 0), c(10, 10))
testDqagIntegration2d_aB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, 0), c(10, 10))
testDqagIntegration2d_aB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(0, 0), c(10, 10))

testDqagIntegration2d_aB(function(x, y) rep.int(0., length(x)), c(-3, 0), c(3, 5))
testDqagIntegration2d_aB(function(x, y) dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))
testDqagIntegration2d_aB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))
testDqagIntegration2d_aB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))
testDqagIntegration2d_aB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))

########################################################
## test Hcubature_ninfInf integration

testHcubature_ninfInf <- function(func, ndim, fdim, expected, tol = 1e-5) {
  if (!("cubature" %in% installed.packages()[, 1])) {
    cat("testHcubature_ninfInf skipping as cubature package is missing\n")
    return()
  }

  cat("testHcubature_ninfInf for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_Hcubature_ninfInf(func, ndim, fdim)
  if (doBench) stop <- get_nanotime()
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testHcubature_ninfInf(dnorm, 1, 1, 1)
testHcubature_ninfInf(function(x) prod(dnorm(x)), 2, 1, 1)
testHcubature_ninfInf(function(x) c(prod(dnorm(x)), prod(x^2*dnorm(x))), 2, 2, c(1, 1))

## the following tests are disabled as pcubature somehow only
## produces NaNs here
if (FALSE) {

  ########################################################
  ## test Pcubature_ninfInf integration

  testPcubature_ninfInf <- function(func, ndim, fdim, expected, tol = 1e-5) {
    if (!("cubature" %in% installed.packages()[, 1])) {
      cat("testPcubature_ninfInf skipping as cubature package is missing\n")
      return()
    }

    cat("testPcubature_ninfInf for", deparse(func), "...")
    if (doBench) start <- get_nanotime()
    actual <- robustlmm:::test_Pcubature_ninfInf(func, ndim, fdim)
    if (doBench) stop <- get_nanotime()
    stopifnot(all.equal(expected, actual, tolerance = tol))
    timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
    cat("ok", timing, "\n")
  }

  testPcubature_ninfInf(dnorm, 1, 1, 1)
  testPcubature_ninfInf(function(x) prod(dnorm(x)), 2, 1, 1)
  testPcubature_ninfInf(function(x) c(prod(dnorm(x)), prod(x^2*dnorm(x))), 2, 2, c(1, 1))

  ########################################################
  ## test Pcubature_aInf integration

  testPcubature_aInf <- function(func, ndim, fdim, bound, expected, tol = 1e-5) {
    if (!("cubature" %in% installed.packages()[, 1])) {
      cat("testPcubature_aInf skipping as cubature package is missing\n")
      return()
    }

    cat("testPcubature_aInf for", deparse(func), "...")
    if (doBench) start <- get_nanotime()
    actual <- robustlmm:::test_Pcubature_aInf(func, ndim, fdim, bound)
    if (doBench) stop <- get_nanotime()
    stopifnot(all.equal(expected, actual, tolerance = tol))
    timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
    cat("ok", timing, "\n")
  }

  testPcubature_aInf(dnorm, 1, 1, 2, 1 - pnorm(2))
  testPcubature_aInf(function(x) prod(dnorm(x)), 2, 1, c(1, 2), (1 -  pnorm(1)) * (1 - pnorm(2)))
  testPcubature_aInf(function(x) c(prod(dnorm(x)), prod(x^2*dnorm(x))), 2, 2, c(1, 2),
                     c((1 - pnorm(1)) * (1 - pnorm(2)), 0.05237427))

  ########################################################
  ## test Pcubature_ninfB integration

  testPcubature_ninfB <- function(func, ndim, fdim, bound, expected, tol = 1e-5) {
    if (!("cubature" %in% installed.packages()[, 1])) {
      cat("testPcubature_ninfB skipping as cubature package is missing\n")
      return()
    }

    cat("testPcubature_ninfB for", deparse(func), "...")
    if (doBench) start <- get_nanotime()
    actual <- robustlmm:::test_Pcubature_ninfB(func, ndim, fdim, bound)
    if (doBench) stop <- get_nanotime()
    stopifnot(all.equal(expected, actual, tolerance = tol))
    timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
    cat("ok", timing, "\n")
  }

  testPcubature_ninfB(dnorm, 1, 1, 2, pnorm(2))
  testPcubature_ninfB(function(x) prod(dnorm(x)), 2, 1, c(1, 2), pnorm(1) * pnorm(2))
  testPcubature_ninfB(function(x) c(prod(dnorm(x)), prod(x^2*dnorm(x))), 2, 2, c(1, 2),
                      c(pnorm(1) * pnorm(2), 0.5210166))

} ## end disabling tests

########################################################
## test Pcubature_aB integration

testPcubature_aB <- function(func, ndim, fdim, a, b, expected, tol = 1e-5) {
  if (!("cubature" %in% installed.packages()[, 1])) {
    cat("testPcubature_aB skipping as cubature package is missing\n")
    return()
  }

  cat("testPcubature_aB for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_Pcubature_aB(func, ndim, fdim, a, b)
  if (doBench) stop <- get_nanotime()
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testPcubature_aB(dnorm, 1, 1, -1, 2, pnorm(2) - pnorm(-1))
testPcubature_aB(function(x) prod(dnorm(x)), 2, 1, c(-3, -4), c(1, 2),
                 (pnorm(1) - pnorm(-3)) * (pnorm(2) - pnorm(-4)))
testPcubature_aB(function(x) c(prod(dnorm(x)), prod(x^2*dnorm(x))), 2, 2,
                 c(-3, -4), c(1, 2),
                 c( (pnorm(1) - pnorm(-3)) * (pnorm(2) - pnorm(-4)), 0.5079543))

########################################################
## test Hcubature2d_ninfInf

testHcubature2d_ninfInf <- function(func, tol = 1e-5) {
  if (!("cubature" %in% installed.packages()[, 1])) {
    cat("testHcubature2d_ninfInf skipping as cubature package is missing\n")
    return()
  }

  cat("testHcubature2d_ninfInf for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_Hcubature2d_ninfInf(func)
  if (doBench) stop <- get_nanotime()
  inner <- function(y) {
    sapply(y, function(yi) integrate(func, -Inf, Inf, y = yi)$value)
  }
  expected <- integrate(inner, -Inf, Inf)$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testHcubature2d_ninfInf(function(x, y) rep.int(0., length(x)))
testHcubature2d_ninfInf(function(x, y) dnorm(x) * dnorm(y))
testHcubature2d_ninfInf(function(x, y) x^2 * dnorm(x) * dnorm(y))
testHcubature2d_ninfInf(function(x, y) y^2 * dnorm(x) * dnorm(y))
testHcubature2d_ninfInf(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y))

########################################################
## test Hcubature2d_aInf

testHcubature2d_aInf <- function(func, a, tol = 1e-5) {
  if (!("cubature" %in% installed.packages()[, 1])) {
    cat("testHcubature2d_aInf skipping as cubature package is missing\n")
    return()
  }

  cat("testHcubature2d_aInf for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_Hcubature2d_aInf(func, a)
  if (doBench) stop <- get_nanotime()
  inner <- function(y) {
    sapply(y, function(yi) integrate(func, a[1], Inf, y = yi)$value)
  }
  expected <- integrate(inner, a[2], Inf)$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testHcubature2d_aInf(function(x, y) rep.int(0., length(x)), c(0, 0))
testHcubature2d_aInf(function(x, y) dnorm(x) * dnorm(y), c(0, 0))
testHcubature2d_aInf(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, 0))
testHcubature2d_aInf(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, 0))
testHcubature2d_aInf(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(0, 0))

testHcubature2d_aInf(function(x, y) rep.int(0., length(x)), c(-3, 0))
testHcubature2d_aInf(function(x, y) dnorm(x) * dnorm(y), c(-3, 0))
testHcubature2d_aInf(function(x, y) x^2 * dnorm(x) * dnorm(y), c(-3, 0))
testHcubature2d_aInf(function(x, y) y^2 * dnorm(x) * dnorm(y), c(-3, 0))
testHcubature2d_aInf(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(-3, 0))

########################################################
## test Hcubature2d_ninfB

testHcubature2d_ninfB <- function(func, b, tol = 1e-5) {
  if (!("cubature" %in% installed.packages()[, 1])) {
    cat("testHcubature2d_ninfB skipping as cubature package is missing\n")
    return()
  }

  cat("testHcubature2d_ninfB for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_Hcubature2d_ninfB(func, b)
  if (doBench) stop <- get_nanotime()
  inner <- function(y) {
    sapply(y, function(yi) integrate(func, -Inf, b[1], y = yi)$value)
  }
  expected <- integrate(inner, -Inf, b[2])$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testHcubature2d_ninfB(function(x, y) rep.int(0., length(x)), c(0, 0))
testHcubature2d_ninfB(function(x, y) dnorm(x) * dnorm(y), c(0, 0))
testHcubature2d_ninfB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, 0))
testHcubature2d_ninfB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, 0))
testHcubature2d_ninfB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(0, 0))

testHcubature2d_ninfB(function(x, y) rep.int(0., length(x)), c(3, -1))
testHcubature2d_ninfB(function(x, y) dnorm(x) * dnorm(y), c(3, -1))
testHcubature2d_ninfB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, -1))
testHcubature2d_ninfB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, -1))
testHcubature2d_ninfB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(3, -1))

########################################################
## test Hcubature2d_aB

testHcubature2d_aB <- function(func, a, b, tol = 1e-5) {
  if (!("cubature" %in% installed.packages()[, 1])) {
    cat("testHcubature2d_aB skipping as cubature package is missing\n")
    return()
  }

  cat("testHcubature2d_aB for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_Hcubature2d_aB(func, a, b)
  if (doBench) stop <- get_nanotime()
  inner <- function(y) {
    sapply(y, function(yi) integrate(func, a[1], b[1], y = yi)$value)
  }
  expected <- integrate(inner, a[2], b[2])$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testHcubature2d_aB(function(x, y) rep.int(0., length(x)), c(0, 0), c(10, 10))
testHcubature2d_aB(function(x, y) dnorm(x) * dnorm(y), c(0, 0), c(10, 10))
testHcubature2d_aB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, 0), c(10, 10))
testHcubature2d_aB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, 0), c(10, 10))
testHcubature2d_aB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(0, 0), c(10, 10))

testHcubature2d_aB(function(x, y) rep.int(0., length(x)), c(-3, 0), c(3, 5))
testHcubature2d_aB(function(x, y) dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))
testHcubature2d_aB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))
testHcubature2d_aB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))
testHcubature2d_aB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))

## the following tests are disabled as pcubature somehow only
## produces NaNs here
if (FALSE) {

  ########################################################
  ## test Pcubature2d_ninfInf

  testPcubature2d_ninfInf <- function(func, tol = 1e-5) {
    if (!("cubature" %in% installed.packages()[, 1])) {
      cat("testPcubature2d_ninfInf skipping as cubature package is missing\n")
      return()
    }

    cat("testPcubature2d_ninfInf for", deparse(func), "...")
    if (doBench) start <- get_nanotime()
    actual <- robustlmm:::test_Pcubature2d_ninfInf(func)
    if (doBench) stop <- get_nanotime()
    inner <- function(y) {
      sapply(y, function(yi) integrate(func, -Inf, Inf, y = yi)$value)
    }
    expected <- integrate(inner, -Inf, Inf)$value
    stopifnot(all.equal(expected, actual, tolerance = tol))
    timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
    cat("ok", timing, "\n")
  }

  testPcubature2d_ninfInf(function(x, y) rep.int(0., length(x)))
  testPcubature2d_ninfInf(function(x, y) dnorm(x) * dnorm(y))
  testPcubature2d_ninfInf(function(x, y) x^2 * dnorm(x) * dnorm(y))
  testPcubature2d_ninfInf(function(x, y) y^2 * dnorm(x) * dnorm(y))
  testPcubature2d_ninfInf(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y))

  ########################################################
  ## test Pcubature2d_aInf

  testPcubature2d_aInf <- function(func, a, tol = 1e-5) {
    if (!("cubature" %in% installed.packages()[, 1])) {
      cat("testPcubature2d_aInf skipping as cubature package is missing\n")
      return()
    }

    cat("testPcubature2d_aInf for", deparse(func), "...")
    if (doBench) start <- get_nanotime()
    actual <- robustlmm:::test_Pcubature2d_aInf(func, a)
    if (doBench) stop <- get_nanotime()
    inner <- function(y) {
      sapply(y, function(yi) integrate(func, a[1], Inf, y = yi)$value)
    }
    expected <- integrate(inner, a[2], Inf)$value
    stopifnot(all.equal(expected, actual, tolerance = tol))
    timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
    cat("ok", timing, "\n")
  }

  testPcubature2d_aInf(function(x, y) rep.int(0., length(x)), c(0, 0))
  testPcubature2d_aInf(function(x, y) dnorm(x) * dnorm(y), c(0, 0))
  testPcubature2d_aInf(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, 0))
  testPcubature2d_aInf(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, 0))
  testPcubature2d_aInf(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(0, 0))

  testPcubature2d_aInf(function(x, y) rep.int(0., length(x)), c(-3, 0))
  testPcubature2d_aInf(function(x, y) dnorm(x) * dnorm(y), c(-3, 0))
  testPcubature2d_aInf(function(x, y) x^2 * dnorm(x) * dnorm(y), c(-3, 0))
  testPcubature2d_aInf(function(x, y) y^2 * dnorm(x) * dnorm(y), c(-3, 0))
  testPcubature2d_aInf(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(-3, 0))

  ########################################################
  ## test Pcubature2d_ninfB

  testPcubature2d_ninfB <- function(func, b, tol = 1e-5) {
    if (!("cubature" %in% installed.packages()[, 1])) {
      cat("testPcubature2d_ninfB skipping as cubature package is missing\n")
      return()
    }

    cat("testPcubature2d_ninfB for", deparse(func), "...")
    if (doBench) start <- get_nanotime()
    actual <- robustlmm:::test_Pcubature2d_ninfB(func, b)
    if (doBench) stop <- get_nanotime()
    inner <- function(y) {
      sapply(y, function(yi) integrate(func, -Inf, b[1], y = yi)$value)
    }
    expected <- integrate(inner, -Inf, b[2])$value
    stopifnot(all.equal(expected, actual, tolerance = tol))
    timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
    cat("ok", timing, "\n")
  }

  testPcubature2d_ninfB(function(x, y) rep.int(0., length(x)), c(0, 0))
  testPcubature2d_ninfB(function(x, y) dnorm(x) * dnorm(y), c(0, 0))
  testPcubature2d_ninfB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, 0))
  testPcubature2d_ninfB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, 0))
  testPcubature2d_ninfB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(0, 0))

  testPcubature2d_ninfB(function(x, y) rep.int(0., length(x)), c(3, -1))
  testPcubature2d_ninfB(function(x, y) dnorm(x) * dnorm(y), c(3, -1))
  testPcubature2d_ninfB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, -1))
  testPcubature2d_ninfB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, -1))
  testPcubature2d_ninfB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(3, -1))

} ## end disabling tests

########################################################
## test Pcubature2d_aB

testPcubature2d_aB <- function(func, a, b, tol = 1e-5) {
  if (!("cubature" %in% installed.packages()[, 1])) {
    cat("testPcubature2d_aB skipping as cubature package is missing\n")
    return()
  }

  cat("testPcubature2d_aB for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_Pcubature2d_aB(func, a, b)
  if (doBench) stop <- get_nanotime()
  inner <- function(y) {
    sapply(y, function(yi) integrate(func, a[1], b[1], y = yi)$value)
  }
  expected <- integrate(inner, a[2], b[2])$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testPcubature2d_aB(function(x, y) rep.int(0., length(x)), c(0, 0), c(10, 10))
testPcubature2d_aB(function(x, y) dnorm(x) * dnorm(y), c(0, 0), c(10, 10))
testPcubature2d_aB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(0, 0), c(10, 10))
testPcubature2d_aB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(0, 0), c(10, 10))
testPcubature2d_aB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(0, 0), c(10, 10))

testPcubature2d_aB(function(x, y) rep.int(0., length(x)), c(-3, 0), c(3, 5))
testPcubature2d_aB(function(x, y) dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))
testPcubature2d_aB(function(x, y) x^2 * dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))
testPcubature2d_aB(function(x, y) y^2 * dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))
testPcubature2d_aB(function(x, y) (x^2 + y^2) * dnorm(x) * dnorm(y), c(-3, 0), c(3, 5))

########################################################
## test DqagNormalExpectation2d_ninfInf

testDqagNormalExpectation2d_ninfInf <- function(func, tol = 1e-8) {
  cat("testDqgagExpectation2d_ninfInf for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_DqagNormalExpectation2d_ninfInf(func)
  if (doBench) stop <- get_nanotime()
  func1 <- function(x, y) func(x, y) * dnorm(x) * dnorm(y)
  inner <- function(y) {
    sapply(y, function(yi) integrate(func1, -Inf, Inf, y = yi)$value)
  }
  expected <- integrate(inner, -Inf, Inf)$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testDqagNormalExpectation2d_ninfInf(function(x, y) rep.int(0., length(x)))
testDqagNormalExpectation2d_ninfInf(function(x, y) x)
testDqagNormalExpectation2d_ninfInf(function(x, y) x^2)
testDqagNormalExpectation2d_ninfInf(function(x, y) x * 0 + y^2)
testDqagNormalExpectation2d_ninfInf(function(x, y) x^2 + y^2)

########################################################
## test GaussHermiteNormalExpectation2d

testGaussHermiteNormalExpectation2d <- function(func, tol = 1e-8) {
  cat("testGaussHermiteNormalExpectation2d for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_GaussHermiteNormalExpectation2d(func, 20)
  if (doBench) stop <- get_nanotime()
  func1 <- function(x, y) func(x, y) * dnorm(x) * dnorm(y)
  inner <- function(y) {
    sapply(y, function(yi) integrate(func1, -Inf, Inf, y = yi)$value)
  }
  expected <- integrate(inner, -Inf, Inf)$value
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testGaussHermiteNormalExpectation2d(function(x, y) rep.int(0., length(x)))
testGaussHermiteNormalExpectation2d(function(x, y) x)
testGaussHermiteNormalExpectation2d(function(x, y) x^2)
testGaussHermiteNormalExpectation2d(function(x, y) x * 0 + y^2)
testGaussHermiteNormalExpectation2d(function(x, y) x^2 + y^2)

testExpectation2d_numerator <- function(nodes = 20, tol = 1e-5) {
  cat("testExpectation2d_numerator...")
  a <- 0.2002626964786321
  s <- 0.3867406996647227
  tau <- 1.0
  rho.e <- smoothPsi
  rho.sigma.e <- psi2propII(smoothPsi, 0.8)

  calcTauNumerator <- function(x, y) {
    t <- (x-a*rho.e@psi(x)-s*y)/tau
    numerator <- rho.sigma.e@psi(t)*t*(tau*tau)
    return(numerator)
  }
  actual0 <- robustlmm:::test_DqagNormalExpectation2d_ninfInf(calcTauNumerator)
  func1 <- function(x, y) calcTauNumerator(x, y) * dnorm(x) * dnorm(y)
  inner <- function(y) {
    sapply(y, function(yi) integrate(func1, -Inf, Inf, y = yi)$value)
  }
  expected0 <- integrate(inner, -Inf, Inf)$value

  ## cat(sprintf("\nExpected0 = %g, actual0 = %g\n", expected0, actual0))
  stopifnot(all.equal(expected0, actual0, tolerance = tol))

  actual <- robustlmm:::test_GaussHermiteNormalExpectation2d(calcTauNumerator, nodes)

  gh <- robustlmm:::ghq(nodes, FALSE)
  ghz <- gh$nodes * sqrt(2)
  ghw <- gh$weights / sqrt(pi) #*dnorm(gh$nodes)
  ghZ <- matrix(ghz, nodes, nodes)
  ghZt <- t(ghZ)
  ghW <- ghw %o% ghw

  ## cat(sprintf("Result = %g or %g\n",
  ##             sum(calcTauNumerator(ghz, ghz[1])*ghw),
  ##             sum(calcTauNumerator(ghz[1], ghz)*ghw)
  ##             ))
  expected <- sum(calcTauNumerator(ghZ, ghZt)*ghW)

  ## cat(sprintf("Expected = %g, actual = %g\n", expected, actual))
  stopifnot(all.equal(expected, actual, tolerance = tol),
            all.equal(actual0, actual, tolerance = tol))
  cat("ok\n")
}

testExpectation2d_numerator(400)
## For 13 nodes:
## R: 0.3500566 (using dnorm)
## C++: 0.351420 (using change of variable)

########################################################
## test HcubatureNormalExpectation_ninfInf integration

testHcubatureNormalExpectation_ninfInf <- function(func, ndim, fdim, expected, tol = 1e-5) {
  if (!("cubature" %in% installed.packages()[, 1])) {
    cat("testHcubatureNormalExpectation_ninfInf skipping as cubature package is missing\n")
    return()
  }

  cat("testHcubatureNormalExpectation_ninfInf for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_HcubatureNormalExpectation_ninfInf(func, ndim, fdim)
  if (doBench) stop <- get_nanotime()
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

testHcubatureNormalExpectation_ninfInf(function(x) 1, 1, 1, 1)
testHcubatureNormalExpectation_ninfInf(function(x) 2, 2, 1, 2)
testHcubatureNormalExpectation_ninfInf(function(x) prod(x) + 2, 2, 1, 2)
## testHcubatureNormalExpectation_ninfInf(function(x) c(1, prod(x)), 2, 2, c(1, 0))
## testHcubatureNormalExpectation_ninfInf(function(x) c(1, prod(x), prod(x*x)), 2, 3, c(1, 0, 1))

########################################################
## test test_GaussianQuadratureNd_ninfInf integration

test_GaussianQuadratureNd_ninfInf <- function(func, ndim, fdim, expected, nnodes, tol = 5e-5) {
  cat("test_GaussianQuadratureNd_ninfInf for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_GaussianQuadratureNd_ninfInf(func, ndim, fdim, nnodes);
  if (doBench) stop <- get_nanotime()
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

for (nnodes in c(13, 50)) {
  test_GaussianQuadratureNd_ninfInf(function(x) dnorm(x), 1, 1, 1, nnodes)
  test_GaussianQuadratureNd_ninfInf(function(x) x*dnorm(x), 1, 1, 0, nnodes)
  test_GaussianQuadratureNd_ninfInf(function(x) x*x*dnorm(x), 1, 1, 1, nnodes)
  test_GaussianQuadratureNd_ninfInf(function(x) prod(dnorm(x)), 2, 1, 1, nnodes)
  test_GaussianQuadratureNd_ninfInf(function(x) prod(x^2*dnorm(x)), 2, 1, 1, nnodes)
  test_GaussianQuadratureNd_ninfInf(function(x) c(prod(dnorm(x)), prod(x^2*dnorm(x))), 2, 2, c(1, 1), nnodes)
}

########################################################
## test test_GaussianQuadratureNdNormalExpectation_ninfInf integration

test_GaussianQuadratureNdNormalExpectation_ninfInf <- function(func, ndim, fdim, expected, nnodes, tol = 5e-5) {
  cat("test_GaussianQuadratureNdNormalExpectation_ninfInf for", deparse(func), "...")
  if (doBench) start <- get_nanotime()
  actual <- robustlmm:::test_GaussianQuadratureNdNormalExpectation_ninfInf(func, ndim, fdim, nnodes);
  if (doBench) stop <- get_nanotime()
  stopifnot(all.equal(expected, actual, tolerance = tol))
  timing <- if (doBench) paste(((stop - start) / 1e6), "ms") else ""
  cat("ok", timing, "\n")
}

for (nnodes in c(13, 50)) {
  test_GaussianQuadratureNdNormalExpectation_ninfInf(function(x) 1, 1, 1, 1, nnodes)
  test_GaussianQuadratureNdNormalExpectation_ninfInf(function(x) x, 1, 1, 0, nnodes)
  test_GaussianQuadratureNdNormalExpectation_ninfInf(function(x) x*x, 1, 1, 1, nnodes)
  test_GaussianQuadratureNdNormalExpectation_ninfInf(function(x) 1, 2, 1, 1, nnodes)
  test_GaussianQuadratureNdNormalExpectation_ninfInf(function(x) prod(x^2), 2, 1, 1, nnodes)
  test_GaussianQuadratureNdNormalExpectation_ninfInf(function(x) c(1, prod(x^2)), 2, 2, c(1, 1), nnodes)
}
