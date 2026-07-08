## Disable test
quit()

## github-transition: SLOW TEST (~210s, the slowest in the suite) --
## when merging to the github-transition (CRAN-release) branch, disable
## it there by adding `quit()` at the top, per that branch's convention.
## WS16 step 1 regression tests: .scoreByCluster / .vcov_varpar.
##
## 1. Exact identity: the per-cluster score columns sum to
##    .scoreVec(par0, fit) -- on a diagonal-V_b DASvar fit (weighted)
##    and on the hard cell, block-V_b DAStau with weights + offset.
##    Error metric uses the global un-cancelled scale max|S| (the
##    u-rows of S are themselves ~0 at the fit: each u_j equation
##    involves only its own cluster's data, so its per-cluster sum is
##    the full, solved equation -- asserted below as a property).
## 2. .vcov_varpar returns a symmetric PSD (1+L)x(1+L) matrix with the
##    n.clusters attribute, in (sigma, theta) coordinates.
## 3. Magnitude pin on the clean classical fit (cPsi, Dyestuff):
##    Var(sigma-hat) within a generous factor of the iid Gaussian
##    sigma^2 / (2 (n - 1)) -- catches the factor-12
##    (caseweight-identity misuse) and similar gross errors while
##    tolerating the J = 6 sandwich noise and the mixed-model
##    inflation. The real df calibration evidence is the WS16 anchor
##    (IF-thread1/thread1_validation/ws16_satterthwaite_anchor.R).
## 4. Guards: crossed grouping factors still error out (WS17 Phase 2);
##    nested designs ARE supported (the score decomposes over the
##    top-level clusters); Mallows
##    design weights are now carried (WS11 step 5) -- score identity holds.

suppressMessages(require(robustlmm))

## ----------------------------------------------------------------
## 1a. Identity, diagonal V_b, weighted DASvar smoothPsi
## ----------------------------------------------------------------
set.seed(42)
w <- runif(nrow(Dyestuff), 0.5, 3)
fitD <- suppressWarnings(
    rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = w,
          method = "DASvar", init = lmerNoFit))
SD <- robustlmm:::.scoreByCluster(fitD)
parD <- c(robustlmm:::.fixef(fitD), fitD@pp$b.s,
          log(robustlmm:::.sigma(fitD)), getME(fitD, "theta"))
svD <- robustlmm:::.scoreVec(parD, fitD)
relD <- max(abs(rowSums(SD) - svD)) / max(abs(SD))
cat("identity (diag DASvar, weighted): rel err", format(relD), "\n")
stopifnot(relD < 1e-10)

## ----------------------------------------------------------------
## 1b. Identity, block V_b, weighted + offset DAStau
## ----------------------------------------------------------------
set.seed(202)
wB  <- runif(nrow(sleepstudy), 0.5, 3)
ssB <- within(sleepstudy, off <- rnorm(nrow(sleepstudy), sd = 20))
fitB <- suppressWarnings(suppressMessages(
    rlmer(Reaction ~ Days + (Days | Subject) + offset(off), ssB,
          weights = wB, method = "DAStau", init = lmerNoFit)))
SB <- robustlmm:::.scoreByCluster(fitB)
parB <- c(robustlmm:::.fixef(fitB), fitB@pp$b.s,
          log(robustlmm:::.sigma(fitB)), getME(fitB, "theta"))
svB <- suppressWarnings(robustlmm:::.scoreVec(parB, fitB))
relB <- max(abs(rowSums(SB) - svB)) / max(abs(SB))
cat("identity (block DAStau, weighted+offset): rel err",
    format(relB), "\n")
## DAStau tolerance: .scoreVec reconverges tau from the cache-cleared
## state while .scoreByCluster uses the estimator's cached tau; both
## sit at the same fixed point to ~1e-7 (see PLAN-WS7-WS9 bug 7).
stopifnot(relB < 1e-6)

## u-rows of S are the full per-cluster u-equations, ~0 at the fit
p <- fitB@pp$p; q <- fitB@pp$q
relU0 <- max(abs(SB[p + seq_len(q), ])) / max(abs(SB))
cat("u-rows of S at fit (should be ~0):", format(relU0), "\n")
stopifnot(relU0 < 1e-6)

## ----------------------------------------------------------------
## 2. A: shape, symmetry, PSD, attribute
## ----------------------------------------------------------------
IFB <- suppressWarnings(implicitIF_full(fitB))
AB  <- robustlmm:::.vcov_varpar(fitB, IFB)
LB  <- length(getME(fitB, "theta"))
stopifnot(identical(dim(AB), c(1L + LB, 1L + LB)),
          isTRUE(all.equal(AB, t(AB), check.attributes = FALSE)),
          identical(attr(AB, "boundary"), FALSE),
          min(eigen(AB, symmetric = TRUE,
                    only.values = TRUE)$values) > -1e-10 * max(abs(AB)),
          identical(attr(AB, "n.clusters"), nlevels(ssB$Subject)),
          identical(rownames(AB),
                    c("sigma", paste0("theta", seq_len(LB)))))
cat("A (block DAStau): diag", format(diag(AB)), "\n")

## ----------------------------------------------------------------
## 3. Magnitude pin, clean classical fit
## ----------------------------------------------------------------
fitC <- rlmer(Yield ~ 1 + (1 | Batch), Dyestuff,
              rho.e = cPsi, rho.b = cPsi, method = "DASvar",
              init = lmerNoFit)
AC <- robustlmm:::.vcov_varpar(fitC)
ref <- sigma(fitC)^2 / (2 * (nrow(Dyestuff) - 1))
ratio <- AC["sigma", "sigma"] / ref
cat("clean cPsi Var(sigma-hat):", format(AC["sigma", "sigma"]),
    " iid ref:", format(ref), " ratio:", format(ratio), "\n")
stopifnot(ratio > 0.2, ratio < 6)

## ----------------------------------------------------------------
## 4. Guards
## ----------------------------------------------------------------
## crossed (plate x sample): no common top-level cluster, so the one-way
## cluster builder (.scoreByCluster) rejects it -- but the multiway
## Cameron-Gelbach-Miller covariance (.vcov_varpar, via the per-
## observation .scoreByObs) handles it (WS17 Phase 2). The covariance
## itself is MC-validated in PLAN-WS17 (analytic/empirical SE ratios ~1).
fitP <- suppressMessages(suppressWarnings(
    rlmer(diameter ~ (1 | plate) + (1 | sample), Penicillin,
          rho.e = cPsi, rho.b = cPsi, method = "DASvar",
          init = lmerNoFit)))
stopifnot(is.null(robustlmm:::.top_cluster(fitP)),
          inherits(tryCatch(robustlmm:::.scoreByCluster(fitP),
                            error = function(e) e), "error"))
## per-observation score columns sum to the full score (~0 at the fit)
GP   <- robustlmm:::.scoreByObs(fitP)
parP <- c(robustlmm:::.fixef(fitP), fitP@pp$b.s,
          log(robustlmm:::.sigma(fitP)), getME(fitP, "theta"))
stopifnot(identical(ncol(GP), fitP@pp$n),
          max(abs(rowSums(GP) - robustlmm:::.scoreVec(parP, fitP))) < 1e-5)
## the multiway covariance is PSD and finite
AP <- suppressWarnings(robustlmm:::.vcov_varpar(fitP))
stopifnot(isTRUE(attr(AP, "multiway")), all(is.finite(AP)),
          min(eigen(AP, symmetric = TRUE, only.values = TRUE)$values) >=
              -1e-10)

## WS17 nested (1 | school/class): the top cluster is the coarsest
## factor (school); the score still decomposes (columns sum to
## .scoreVec ~ 0) and the df is finite. The covariance itself is
## MC-validated in PLAN-WS17 (analytic/empirical SE ratios ~1).
set.seed(317)
dN <- expand.grid(class = factor(1:3), school = factor(1:8), rep = 1:4)
csN <- interaction(dN$school, dN$class, drop = TRUE)
dN$y <- 1 + rnorm(8, 0, 1.2)[dN$school] +
            rnorm(nlevels(csN), 0, 0.7)[csN] + rnorm(nrow(dN))
fitN <- suppressWarnings(rlmer(y ~ 1 + (1 | school/class), dN,
                               method = "DAStau"))
tcN <- robustlmm:::.top_cluster(fitN)
stopifnot(!is.null(tcN), tcN$name == "school", nlevels(tcN$g) == 8L)
SN   <- robustlmm:::.scoreByCluster(fitN)
parN <- c(robustlmm:::.fixef(fitN), fitN@pp$b.s,
          log(robustlmm:::.sigma(fitN)), getME(fitN, "theta"))
stopifnot(identical(dim(SN), c(fitN@pp$p + fitN@pp$q + 1L + 2L, 8L)),
          max(abs(rowSums(SN) - robustlmm:::.scoreVec(parN, fitN))) < 1e-5,
          all(is.finite(robustlmm:::.satterthwaite_df(fitN))))

## WS11 step 5: Mallows design weights ARE now carried by the
## per-cluster score (and the whole influence layer). The per-cluster
## columns must still sum to .scoreVec (~0 at the fit), now with both
## sides eta-weighted.
set.seed(11)
fitE <- suppressWarnings(
    rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, method = "DASvar",
          init = lmerNoFit,
          design.weights = runif(nrow(Dyestuff), 0.5, 1)))
SE   <- robustlmm:::.scoreByCluster(fitE)
parE <- c(robustlmm:::.fixef(fitE), fitE@pp$b.s,
          log(robustlmm:::.sigma(fitE)), getME(fitE, "theta"))
stopifnot(max(abs(rowSums(SE) - robustlmm:::.scoreVec(parE, fitE))) < 1e-6,
          is.numeric(robustlmm:::.satterthwaite_df(fitE)))   # WS16 df now works
cat("guards: ok\n")

## ----------------------------------------------------------------
## 5. Boundary fit (WS14): a variance component estimated at 0 is a
##    REDUCIBLE boundary -- conditional on it being 0, the remaining
##    (sigma) covariance is regular, so .vcov_varpar zeroes that
##    component's row/col, sets boundary = TRUE and reducible = TRUE, and
##    does NOT warn; the conditional Satterthwaite df is finite. (For a
##    single grouping factor only sigma remains, so it is always
##    reducible.) A genuinely singular fit (reduced block degenerate)
##    still warns and sets reducible = FALSE.
## ----------------------------------------------------------------
set.seed(5)
nb <- 40L
gb <- factor(rep(1:8, each = 5L))
yb <- rnorm(nb)
yb <- yb - ave(yb, gb)   # cluster-demeaned: theta-hat = 0 by construction
db <- data.frame(y = yb, g = gb)
fitBnd <- suppressMessages(suppressWarnings(
    rlmer(y ~ 1 + (1 | g), db, method = "DASvar", init = lmerNoFit)))
if (getME(fitBnd, "theta") < 1e-8) {
    wmsg <- character(0)
    ABnd <- withCallingHandlers(
        robustlmm:::.vcov_varpar(fitBnd),
        warning = function(w) {
            wmsg <<- c(wmsg, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    stopifnot(identical(attr(ABnd, "boundary"), TRUE),
              identical(attr(ABnd, "reducible"), TRUE),
              identical(as.integer(attr(ABnd, "boundary.components")), 1L),
              !any(grepl("boundary", wmsg)),       # reducible -> no warning
              all(ABnd[2L, ] == 0), all(ABnd[, 2L] == 0),  # theta zeroed
              ABnd[1L, 1L] > 0)                    # sigma block intact
    dfBnd <- suppressWarnings(robustlmm:::.satterthwaite_df(fitBnd))
    stopifnot(isTRUE(attr(dfBnd, "boundary")),
              isTRUE(attr(dfBnd, "reducible")),
              all(is.finite(as.numeric(dfBnd))))   # usable conditional df
    cat("boundary guard: ok (reducible; theta-hat =",
        getME(fitBnd, "theta"), ")\n")
} else {
    cat("boundary guard: fit not at boundary for this seed; skipped\n")
}

## ----------------------------------------------------------------
## 5b. WS14 key identity: at a reducible boundary the conditional df
##     EQUALS the df of the reduced model (the model with the zeroed
##     component dropped). Nested (1|g1/g2) with the g1:g2 variance
##     truly 0 -> top cluster g1 is identical to the reduced (1|g1) fit,
##     a clean apples-to-apples reference.
## ----------------------------------------------------------------
set.seed(3)
J1b <- 22L; nclb <- 3L; mb <- 5L
dW <- expand.grid(rep = seq_len(mb), g2 = factor(seq_len(nclb)),
                  g1 = factor(seq_len(J1b)))
dW$x <- rnorm(nrow(dW))
dW$y <- 1 + 0.5 * dW$x + rnorm(J1b, 0, 1.3)[as.integer(dW$g1)] +
        rnorm(nrow(dW))                       # no g1:g2 effect
fW   <- suppressWarnings(rlmer(y ~ x + (1 | g1/g2), dW, method = "DASvar"))
fWr  <- suppressWarnings(rlmer(y ~ x + (1 | g1),    dW, method = "DASvar"))
if (min(getME(fW, "theta")) < 1e-8) {        # inner component on boundary
    AW <- suppressWarnings(robustlmm:::.vcov_varpar(fW))
    stopifnot(isTRUE(attr(AW, "boundary")), isTRUE(attr(AW, "reducible")))
    dfW  <- suppressWarnings(robustlmm:::.satterthwaite_df(fW))
    dfWr <- robustlmm:::.satterthwaite_df(fWr)
    ## relative tolerance: the only discrepancy is FD-gradient noise (the
    ## central FD over theta), ~1e-7 relative on each df.
    stopifnot(max(abs(as.numeric(dfW) - as.numeric(dfWr)) /
                  as.numeric(dfWr)) < 1e-4)
    cat("WS14 nested boundary == reduced-model df: ok\n")
} else {
    cat("WS14 nested boundary: inner component not at boundary; skipped\n")
}

## ----------------------------------------------------------------
## 6. .satterthwaite_df: package FD-of-unsc pipeline vs the
##    closed-form-V route on a clean cPsi fit (the WS16 anchor
##    design, one deterministic dataset). For cPsi the default vcov
##    IS the classical sigma^2 (X' Omega^-1 X)^-1, so the two df
##    computations must agree to FD accuracy.
## ----------------------------------------------------------------
set.seed(20260611 + 6024)
Jd <- 24L; n_per <- 5L; nd <- Jd * n_per
gd <- factor(rep(seq_len(Jd), each = n_per))
xd <- rnorm(nd)
yd <- 1 + 0.5 * xd + rnorm(Jd, sd = 0.5)[as.integer(gd)] + rnorm(nd)
dd <- data.frame(yd, xd, gd)
fitS <- rlmer(yd ~ xd + (1 | gd), dd, rho.e = cPsi, rho.b = cPsi,
              method = "DASvar", init = lmerNoFit)
IFS <- suppressWarnings(implicitIF_full(fitS))
dfP <- robustlmm:::.satterthwaite_df(fitS, IF = IFS)

## closed-form route (as in ws16_satterthwaite_anchor.R)
AS   <- robustlmm:::.vcov_varpar(fitS, IFS)
Xd   <- as.matrix(fitS@pp$X)
Ztd  <- fitS@pp$Zt
sigS <- sigma(fitS); thS <- getME(fitS, "theta")
vfun <- function(sig, th) {
    Om <- th^2 * crossprod(as.matrix(Ztd)) + diag(nd)
    diag(sig^2 * solve(t(Xd) %*% solve(Om, Xd)))
}
v0  <- vfun(sigS, thS)
h   <- 1e-5 * max(thS, 0.1)
gth <- (vfun(sigS, thS + h) - vfun(sigS, thS - h)) / (2 * h)
gsg <- 2 * v0 / sigS
dfC <- sapply(1:2, function(j) {
    gj <- c(gsg[j], gth[j])
    2 * v0[j]^2 / drop(t(gj) %*% AS %*% gj)
})
relDf <- max(abs(dfP - dfC) / dfC)
cat("satterthwaite_df vs closed form: df_pkg",
    format(dfP, digits = 5), " df_closed", format(dfC, digits = 5),
    " max rel err", format(relDf), "\n")
stopifnot(relDf < 1e-4, all(is.finite(dfP)), all(dfP > 0))

## pp state restored after the FD sweep
stopifnot(identical(getME(fitS, "theta"), thS),
          isTRUE(all.equal(as.matrix(vcov(fitS)),
                           sigS^2 * as.matrix(fitS@pp$unsc()),
                           check.attributes = FALSE)))
cat("test 6 (.satterthwaite_df closed-form pin + state restore): ok\n")

cat("vcov-varpar: all tests passed\n")
