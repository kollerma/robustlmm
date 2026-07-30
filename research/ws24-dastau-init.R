## WS24 / Q4a: does the (buggy) DASvar tau_e initial value change the
## converged DAStau tau?  The fix commit asserted "DAStau unchanged
## (its calcTau fixed point never depended on the DASvar init)".
## Run against the PRE-FIX library.
##   Rscript ws24-dastau-init.R <libpath>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages({library(lme4); library(robustlmm)})
cat("robustlmm", as.character(packageVersion("robustlmm")), "\n\n")

MASTER_SEED <- 20260611L
J <- 50L; n_per <- 5L; n <- J * n_per
set.seed(MASTER_SEED)
w_true <- exp(runif(n, log(0.25), log(4)))
x <- rnorm(n); g <- factor(rep(seq_len(J), each = n_per))
set.seed(MASTER_SEED + 1L)
b <- rnorm(J, sd = 0.5); e <- rnorm(n, sd = 1 / sqrt(w_true))
dat <- data.frame(y = 1 + 0.5 * x + b[as.integer(g)] + e, x = x, g = g)

fit <- suppressWarnings(suppressMessages(
    rlmer(y ~ x + (1 | g), dat, weights = w_true,
          rho.e = cPsi, rho.b = cPsi, method = "DAStau")))
pp <- fit@pp
cat("converged sigma  :", sigma(fit), "\n")

## the two candidate initial values, at THIS pp (identical a, s)
Btmp   <- pp$B()
tmp    <- Matrix::tcrossprod(pp$Epsi_bpsi_bt, Btmp)
bterm  <- robustlmm:::computeDiagonalOfProduct(as(Btmp, "unpackedMatrix"),
                                               as(tmp,  "unpackedMatrix"))
lin    <- pp$EDpsi_e * 2 * pp$diagA
quad   <- pp$Epsi2_e * pp$diagAAt
tau2_old <- pp$v_e - lin + quad + bterm     ## pre-fix leading term
tau2_new <- 1       - lin + quad + bterm     ## post-fix leading term

cat("\n-- initial tau^2 --\n")
cat(sprintf("old: range [%.4f, %.4f]; > 1: %d of %d; < 0.01 (floor): %d\n",
            min(tau2_old), max(tau2_old), sum(tau2_old > 1), n, sum(tau2_old < 0.01)))
cat(sprintf("new: range [%.4f, %.4f]; > 1: %d of %d; < 0.01 (floor): %d\n",
            min(tau2_new), max(tau2_new), sum(tau2_new > 1), n, sum(tau2_new < 0.01)))
clamp <- function(t2) {t2[t2 < 0.01] <- 0.01; sqrt(t2)}
init_old <- clamp(tau2_old); init_new <- clamp(tau2_new)
cat(sprintf("old init tau > 1: %d of %d  (calcTau's `while (... && tau[i] <= 1)` guard skips these)\n",
            sum(init_old > 1), n))

## same (a, s), same rho, only the starting value differs
s  <- robustlmm:::.s(fit, theta = FALSE, pp = pp, B = Btmp)
a  <- pp$diagA
cT <- function(init) robustlmm:::calcTau(a, s, fit@rho.e, fit@rho.sigma.e,
                                         pp, pp$kappa_e, init)
tau_from_old <- suppressWarnings(cT(init_old))
tau_from_new <- suppressWarnings(cT(init_new))
tau_from_1   <- suppressWarnings(cT(rep(1, n)))

cat("\n-- converged tau after ONE calcTau pass, same (a, s), different start --\n")
cat(sprintf("max |tau(old init) - tau(new init)| = %.3e\n",
            max(abs(tau_from_old - tau_from_new))))
cat(sprintf("max |tau(new init) - tau(1)|        = %.3e\n",
            max(abs(tau_from_new - tau_from_1))))
cat(sprintf("n with tau(old init) == 1 exactly (clamped): %d\n", sum(tau_from_old == 1)))
cat(sprintf("tau(old init) range [%.4f, %.4f]; tau(new init) range [%.4f, %.4f]\n",
            min(tau_from_old), max(tau_from_old), min(tau_from_new), max(tau_from_new)))
i <- which.max(abs(tau_from_old - tau_from_new))
cat(sprintf("worst obs %d: w = %.3f, init old %.4f -> %.4f, init new %.4f -> %.4f\n",
            i, w_true[i], init_old[i], tau_from_old[i], init_new[i], tau_from_new[i]))

## does the fit's own converged tau agree with the (a, s) fixed point?
cat(sprintf("\nfit's stored tau_e vs tau(new init): max abs diff %.3e\n",
            max(abs(pp$tau_e() - tau_from_new))))
