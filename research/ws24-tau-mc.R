## WS24 / Q1, decisive check: validate the LEADING CONSTANT of the
## DASvar tau^2 against tau's own definition, without using lmer.
##
## tau_{e,i} is the standard deviation of the standardised residual
##   r_i = sqrt(w_i) (y_i - x_i'beta-hat - z_i'b-hat) / sigma
## under the linearisation of the estimating equations.  At
## rho.e = rho.b = cPsi the linearisation is exact, so with theta and
## sigma held at their true values the Monte-Carlo sd of r_i must
## equal tau_i to MC error -- for ANY prior weights.
##
## Compares:
##   tau_new = 1   - 2 lambda diagA + Epsi2 diagAAt + B-term   (post-fix)
##   tau_old = v_e - 2 lambda diagA + Epsi2 diagAAt + B-term   (pre-fix)
## against the empirical sd.  Run against the FIXED library.
##   Rscript ws24-tau-mc.R <libpath> <outdir> [nrep]
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
outdir <- args[2]; NREP <- if (length(args) > 2) as.integer(args[3]) else 4000L
suppressPackageStartupMessages({library(lme4); library(robustlmm); library(parallel)})

SEED <- 20260726L
J <- 30L; m <- 6L; n <- J * m
set.seed(SEED)
g      <- factor(rep(seq_len(J), each = m))
x      <- rnorm(n)
w_true <- exp(runif(n, log(0.05), log(20)))     # 400-fold weight range
sigma_true <- 1.3; sb_true <- 0.8
theta_true <- sb_true / sigma_true               # relative sd, lme4 parametrisation
beta_true  <- c(1, 0.5)
WGT <- w_true

sim <- function(seed) {
    set.seed(seed)
    b <- rnorm(J, sd = sb_true)
    e <- rnorm(n, sd = sigma_true / sqrt(w_true))
    beta_true[1] + beta_true[2] * x + b[as.integer(g)] + e
}
dat <- data.frame(y = sim(SEED + 1L), x = x, g = g)
form <- y ~ x + (1 | g)

## one fit, used only to build the predictor module and to read tau
fit0 <- suppressMessages(suppressWarnings(
    rlmer(form, dat, weights = WGT, rho.e = cPsi, rho.b = cPsi,
          method = "DASvar", init = lmerNoFit)))

## analytic tau at the TRUE theta -----------------------------------
robustlmm:::setSigma(fit0, sigma_true)
robustlmm:::setTheta(fit0, theta_true, fit.effects = TRUE, update.sigma = FALSE)
pp   <- fit0@pp
Btmp <- pp$B()
tmp  <- Matrix::tcrossprod(pp$Epsi_bpsi_bt, Btmp)
bt   <- robustlmm:::computeDiagonalOfProduct(as(Btmp, "unpackedMatrix"),
                                             as(tmp, "unpackedMatrix"))
corr <- pp$EDpsi_e * 2 * pp$diagA - pp$Epsi2_e * pp$diagAAt - bt
tau_new <- sqrt(1       - corr)
tau_old <- sqrt(pmax(pp$v_e - corr, 0.01))
tau_pkg <- pp$tau_e()

## Monte-Carlo sd of the standardised residual ----------------------
one <- function(rep) {
    d <- dat; d$y <- sim(SEED + 1000L + rep)
    f <- suppressMessages(suppressWarnings(
        rlmer(form, d, weights = WGT, rho.e = cPsi, rho.b = cPsi,
              method = "DASvar", init = lmerNoFit)))
    robustlmm:::setSigma(f, sigma_true)
    robustlmm:::setTheta(f, theta_true, fit.effects = TRUE, update.sigma = FALSE)
    f@resp$wtres / sigma_true
}
R <- do.call(rbind, mclapply(seq_len(NREP), one, mc.cores = 5L))
tau_mc  <- apply(R, 2, sd)
mcse    <- tau_mc / sqrt(2 * (NREP - 1))          # sd of an sd estimate
z_new   <- (tau_new - tau_mc) / mcse
z_old   <- (tau_old - tau_mc) / mcse

res <- list(seed = SEED, nrep = NREP, n = n, weightRange = range(w_true),
            theta_true = theta_true, sigma_true = sigma_true,
            tau_mc = tau_mc, mcse = mcse, tau_new = tau_new, tau_old = tau_old,
            tau_pkg = tau_pkg, z_new = z_new, z_old = z_old,
            stamp = c(robustlmm = as.character(packageVersion("robustlmm")),
                      lme4 = as.character(packageVersion("lme4")),
                      R = R.version.string))
saveRDS(res, file.path(outdir, "tau-mc.rds"))

cat(sprintf("MC reps %d, n = %d, weights in [%.3f, %.2f]\n\n", NREP, n,
            min(w_true), max(w_true)))
cat(sprintf("package tau_e == tau_new to %.2e (sanity)\n\n",
            max(abs(tau_pkg - tau_new))))
f <- function(z, lab) cat(sprintf(
    "%-22s max|z| %8.2f   mean z %8.2f   |z|>3 in %3d of %d obs   max rel dev %7.3f%%\n",
    lab, max(abs(z)), mean(z), sum(abs(z) > 3), n,
    100 * max(abs(if (lab == "post-fix (leading 1)") tau_new else tau_old) / tau_mc - 1)))
f(z_new, "post-fix (leading 1)")
f(z_old, "pre-fix  (leading v_e)")
