## Disable test
quit()

## github-transition: MONTE-CARLO TEST (a few seconds) -- when merging to
## the github-transition (CRAN-release) branch, disable it there by adding
## `quit()` at the top, per that branch's convention.
##
## Seeded leverage-breakdown smoke test for the Mallows design weights
## (PLAN-3.5.0-stabilization.md item 5.4). The redescending psi of the RSE
## bounds the influence of response outliers but NOT of high-leverage
## design points: a "good-leverage-bad" point (extreme covariate, small
## residual under a spurious slope) drags the plain RSE slope just as it
## drags least squares. The Mallows weights eta downweight such points by
## their design leverage, so the slope stays near the truth.
##
## This is a regression guard against gross breakage of that protection,
## not a precise efficiency study (that lives in the build-ignored
## inst/simulationStudy/leverageBreakdown.R). It asserts, over a small
## seeded Monte Carlo with planted leverage contamination, that the
## Mallows slope bias is small AND much smaller than the plain RSE's.

suppressMessages(require(robustlmm))

nsim   <- 25L
b0     <- 1; b1 <- 1                 # truth: intercept 1, slope 1
bc     <- -1                         # contaminant slope (drags the fit)
mc.cores <- if (.Platform$OS.type == "windows") 1L else 2L

## one replicate: single grouping factor, ~7% high-leverage contamination
## following the spurious slope bc; returns the plain-RSE and Mallows slope.
oneRep <- function(i) {
    set.seed(5000L + i)
    J <- 15L; m <- 6L; n <- J * m
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(n)
    re <- rnorm(J, 0, 0.8)[as.integer(g)]
    y <- b0 + b1 * x + re + rnorm(n)
    ## planted leverage points: extreme x, response on the spurious slope
    nc <- max(2L, round(0.07 * n))
    ci <- seq_len(nc)                # first nc obs (deterministic given seed)
    x[ci] <- runif(nc, 8, 12)
    y[ci] <- b0 + bc * x[ci] + rnorm(nc, 0, 0.3)
    d <- data.frame(y, x, g)
    sl <- function(f) tryCatch(unname(fixef(f)["x"]),
                               error = function(e) NA_real_)
    fR <- tryCatch(suppressWarnings(rlmer(y ~ x + (1 | g), d,
                   method = "DASvar")), error = function(e) NULL)
    fM <- tryCatch(suppressWarnings(rlmer(y ~ x + (1 | g), d,
                   method = "DASvar", design.weights = "mcd")),
                   error = function(e) NULL)
    c(rse = if (is.null(fR)) NA_real_ else sl(fR),
      mallows = if (is.null(fM)) NA_real_ else sl(fM))
}

M <- do.call(rbind, parallel::mclapply(seq_len(nsim), oneRep,
                                       mc.cores = mc.cores))
ok <- is.finite(M[, "rse"]) & is.finite(M[, "mallows"])
stopifnot(sum(ok) >= 0.9 * nsim)
biasRse     <- mean(M[ok, "rse"])     - b1
biasMallows <- mean(M[ok, "mallows"]) - b1
cat(sprintf("leverage breakdown: slope bias  RSE=%.2f  Mallows=%.2f\n",
            biasRse, biasMallows))

## The two absolute bounds carry margin around the seeded values and
## together capture the claim -- Mallows stays near truth while the plain
## RSE breaks down -- so by construction Mallows is the better of the two.
## biasRse uses no design weights, hence no covMcd: it is fully
## deterministic across platforms; only biasMallows depends on covMcd's
## subsampling, which the leverage points (x in 8-12) make robust.
stopifnot(abs(biasMallows) < 0.45,                 # Mallows stays near truth
          abs(biasRse)     > 0.5)                  # plain RSE breaks down
cat("design-weights-breakdown: all checks passed\n")
