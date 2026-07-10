## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Regression tests for the EXPERIMENTAL Monte-Carlo DAS-tau calibration
## ("DASmc": options robustlmm.dastau.mc / robustlmm.dastau.mc.all /
## robustlmm.dasmc.nsim / robustlmm.dasmc.seed, see ?`robustlmm-options`).
##
## The MC path replaces the s == 2 Gauss-Hermite fixed point of
## calcTau.nondiag() with a plain-Monte-Carlo fixed point that is valid
## for any block dimension s >= 2, drawn ONCE per fit with common random
## numbers (deterministically seeded, moment-matched to exact zero mean
## and identity second moment). Simulation-validated (it removes the
## DASvar sphericity calibration residual at s = 4 in the companion ar1
## study); guarded here by fast deterministic checks:
##
## 1. Determinism: repeated fits under the option are identical, and the
##    caller's .Random.seed is untouched by the CRN draw.
## 2. Classical exactness: under cPsi the moment-matched MC fixed point
##    reproduces the DASvar closed form exactly, so the full DAStau-MC
##    fit must equal the DASvar fit to solver precision (~1e-8 here;
##    the prototype gate showed 8e-11 on an NC = 4 fit).
## 3. s > 2 enablement: option ON, method = "DAStau" must NOT fall back
##    to DASvar for an NC = 3 ar1 block (no fallback warning, method
##    stays "DAStau") and the calibration must measurably differ from
##    DASvar; option OFF must preserve the pre-existing fallback warning.
## 4. s = 2 default path untouched: option ON without .mc.all leaves an
##    NC = 2 fit bit-identical to the Gauss-Hermite DAStau fit.
## 5. Option hygiene: robustlmm.dasmc.nsim is respected (tiny-N fit runs
##    and differs from a large-N fit; smoke only).
##
## All options are set locally and restored via on.exit().

suppressMessages(require(robustlmm))

## balanced ar1 repeated-measures generator (validation-study design,
## shrunk to test size)
makeData <- function(seed, nc, J, reps = 2, rho0 = -0.5,
                     sigma0 = 1.5, resid = 0.7) {
    set.seed(seed)
    R <- rho0^abs(outer(seq_len(nc), seq_len(nc), "-"))
    ch <- chol(sigma0^2 * R)
    fl <- factor(seq_len(nc))
    d <- do.call(rbind, lapply(seq_len(J), function(j) {
        b <- as.vector(crossprod(ch, rnorm(nc)))
        do.call(rbind, lapply(seq_len(reps), function(r) {
            data.frame(subj = j, f = fl, y = 10 + b + rnorm(nc, 0, resid))
        }))
    }))
    d$subj <- factor(d$subj)
    d
}
form <- y ~ 1 + ar1(0 + f | subj)
pars <- function(fit) c(fixef(fit), sigma = sigma(fit), theta(fit))
d3 <- makeData(101, nc = 3, J = 24)

## ----------------------------------------------------------------
## 1. Determinism (CRN) + caller's RNG state untouched
## ----------------------------------------------------------------
local({
    oo <- options(robustlmm.dastau.mc = TRUE, robustlmm.dasmc.nsim = 5000)
    on.exit(options(oo))
    set.seed(42)
    rngBefore <- .Random.seed
    fA <- suppressMessages(rlmer(form, data = d3, method = "DAStau"))
    rngAfter <- get(".Random.seed", envir = globalenv())
    fB <- suppressMessages(rlmer(form, data = d3, method = "DAStau"))
    stopifnot(identical(theta(fA), theta(fB)),
              identical(sigma(fA), sigma(fB)),
              identical(fixef(fA), fixef(fB)),
              identical(rngBefore, rngAfter))
    cat("1. determinism + RNG-state restoration: OK\n")
})

## ----------------------------------------------------------------
## 2. Classical exactness: cPsi => DAStau-MC == DASvar
## ----------------------------------------------------------------
## Under cPsi the radial weights are constant, so the moment-matched MC
## fixed point returns the closed-form second moment EXACTLY (independent
## of nsim); the two fits may differ only by solver tolerance.
local({
    fVar <- rlmer(form, data = d3, method = "DASvar",
                  rho.e = cPsi, rho.b = cPsi,
                  rho.sigma.e = cPsi, rho.sigma.b = cPsi)
    oo <- options(robustlmm.dastau.mc = TRUE, robustlmm.dasmc.nsim = 2000)
    on.exit(options(oo))
    fMc <- suppressMessages(
        rlmer(form, data = d3, method = "DAStau",
              rho.e = cPsi, rho.b = cPsi,
              rho.sigma.e = cPsi, rho.sigma.b = cPsi))
    dPar <- max(abs(pars(fMc) - pars(fVar)))
    cat("2. classical (cPsi) DAStau-MC vs DASvar: max |diff| =",
        format(dPar), "\n")
    stopifnot(dPar < 1e-7)
})

## ----------------------------------------------------------------
## 3. s > 2: option lifts the DASvar fallback; OFF keeps it
## ----------------------------------------------------------------
local({
    ## option OFF (default): the fallback warning must still fire
    warned <- character()
    fOff <- withCallingHandlers(
        rlmer(form, data = d3, method = "DAStau"),
        warning = function(w) {
            warned <<- c(warned, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    stopifnot(any(grepl("does not support blocks of size larger than 2",
                        warned)),
              identical(fOff@method, "DASvar"))

    ## option ON: no fallback warning, method stays DAStau, the
    ## experimental-path message fires, and the calibration measurably
    ## differs from DASvar
    oo <- options(robustlmm.dastau.mc = TRUE, robustlmm.dasmc.nsim = 5000)
    on.exit(options(oo))
    warned2 <- msgs <- character()
    fOn <- withCallingHandlers(
        rlmer(form, data = d3, method = "DAStau"),
        warning = function(w) {
            warned2 <<- c(warned2, conditionMessage(w))
            invokeRestart("muffleWarning")
        },
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    stopifnot(!any(grepl("does not support blocks of size larger than 2",
                         warned2)),
              identical(fOn@method, "DAStau"),
              any(grepl("EXPERIMENTAL Monte-Carlo DAS-tau", msgs)),
              all(is.finite(pars(fOn))))
    dTheta <- max(abs(theta(fOn) - theta(fOff)))
    cat("3. s > 2: fallback preserved OFF; ON stays DAStau,",
        "max |theta diff| vs DASvar =", format(dTheta), "\n")
    stopifnot(dTheta > 1e-6)
})

## ----------------------------------------------------------------
## 4. s = 2 unchanged by default: option ON (without .mc.all) must
##    leave the Gauss-Hermite DAStau fit bit-identical
## ----------------------------------------------------------------
local({
    d2 <- makeData(202, nc = 2, J = 30, rho0 = 0.6)
    fGH <- rlmer(form, data = d2, method = "DAStau")
    oo <- options(robustlmm.dastau.mc = TRUE)
    on.exit(options(oo))
    msgs <- character()
    fOn <- withCallingHandlers(
        rlmer(form, data = d2, method = "DAStau"),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    stopifnot(!any(grepl("Monte-Carlo", msgs)),
              identical(theta(fGH), theta(fOn)),
              identical(sigma(fGH), sigma(fOn)),
              identical(fixef(fGH), fixef(fOn)))
    cat("4. s = 2 without .mc.all: bit-identical to GH DAStau: OK\n")
})

## ----------------------------------------------------------------
## 5. Option hygiene: robustlmm.dasmc.nsim is respected (smoke)
## ----------------------------------------------------------------
local({
    oo <- options(robustlmm.dastau.mc = TRUE, robustlmm.dasmc.nsim = 500)
    on.exit(options(oo))
    fSmall <- suppressWarnings(suppressMessages(
        rlmer(form, data = d3, method = "DAStau")))
    options(robustlmm.dasmc.nsim = 20000)
    fLarge <- suppressMessages(rlmer(form, data = d3, method = "DAStau"))
    stopifnot(!identical(theta(fSmall), theta(fLarge)),
              all(is.finite(pars(fSmall))),
              all(is.finite(pars(fLarge))))
    cat("5. nsim respected (500 vs 20000 differ, both finite): OK\n")
})

## options restored by the local() blocks above
stopifnot(is.null(getOption("robustlmm.dastau.mc")),
          is.null(getOption("robustlmm.dastau.mc.all")),
          is.null(getOption("robustlmm.dasmc.nsim")),
          is.null(getOption("robustlmm.dasmc.seed")))
cat("dastau-mc: all tests passed\n")
