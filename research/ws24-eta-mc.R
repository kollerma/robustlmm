## WS24 / Q4b: is the Mallows design weight eta double-counted in the
## same way the prior weight v_e was?
##
## Structurally it cannot be: U_e = diag(sqrt(v_e)) whitens with the
## PRIOR weights only; eta never enters the whitening, it enters the
## bread through D_e = diag(eta * lambda_e) and the meat through the
## explicit H / H^2 factors.  This checks it numerically against tau's
## definition, in the hardest case -- prior weights AND design weights
## active at the same time.
##
## At rho.e = rho.b = cPsi with theta, sigma held at the truth, the MC
## sd of the standardised residual must equal tau_e.  Arms:
##   eta = 1,   w varying     (baseline, = the Q1 check)
##   eta varying, w = 1        (design weights alone)
##   eta varying, w varying    (both -- the double-count candidate)
## Each is also compared against the two naive alternatives that WOULD
## show a double count: leading term eta (or eta * v_e) instead of 1.
##   Rscript ws24-eta-mc.R <libpath> <outdir> [nrep]
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
outdir <- args[2]; NREP <- if (length(args) > 2) as.integer(args[3]) else 4000L
suppressPackageStartupMessages({library(lme4); library(robustlmm); library(parallel)})

SEED <- 20260726L
J <- 30L; m <- 6L; n <- J * m
set.seed(SEED)
g <- factor(rep(seq_len(J), each = m)); x <- rnorm(n)
w_var  <- exp(runif(n, log(0.05), log(20)))
eta_var <- pmin(1, pmax(0.15, runif(n, 0.15, 1)))
sigma_true <- 1.3; sb_true <- 0.8
theta_true <- sb_true / sigma_true; beta_true <- c(1, 0.5)

WGT <- NULL; ETA <- NULL
form <- y ~ x + (1 | g)
simy <- function(seed, w) {
    set.seed(seed)
    beta_true[1] + beta_true[2] * x + rnorm(J, sd = sb_true)[as.integer(g)] +
        rnorm(n, sd = sigma_true / sqrt(w))
}
mkfit <- function(d) {
    f <- suppressMessages(suppressWarnings(
        rlmer(form, d, weights = WGT, design.weights = ETA,
              rho.e = cPsi, rho.b = cPsi, method = "DASvar", init = lmerNoFit)))
    robustlmm:::setSigma(f, sigma_true)
    robustlmm:::setTheta(f, theta_true, fit.effects = TRUE, update.sigma = FALSE)
    f
}
arms <- list(
    "eta=1, w varying"       = list(w = w_var,      eta = NULL),
    "eta varying, w=1"       = list(w = rep(1, n),  eta = eta_var),
    "eta varying, w varying" = list(w = w_var,      eta = eta_var))

out <- list()
for (an in names(arms)) {
    WGT <<- arms[[an]]$w; ETA <<- arms[[an]]$eta
    etaVec <- if (is.null(ETA)) rep(1, n) else ETA
    d0 <- data.frame(y = simy(SEED + 1L, WGT), x = x, g = g)
    f0 <- mkfit(d0); pp <- f0@pp
    tau_pkg <- pp$tau_e()
    ## the three candidate leading terms, correction term held fixed
    Amat  <- as.matrix(pp$A())
    Btmp  <- pp$B()
    bt    <- robustlmm:::computeDiagonalOfProduct(
                as(Btmp, "unpackedMatrix"),
                as(Matrix::tcrossprod(pp$Epsi_bpsi_bt, Btmp), "unpackedMatrix"))
    corr  <- etaVec * pp$EDpsi_e * 2 * pp$diagA -
             pp$Epsi2_e * as.numeric(Amat^2 %*% (etaVec^2)) - bt
    cand  <- cbind(one       = 1        - corr,
                   eta       = etaVec   - corr,
                   v_e       = pp$v_e   - corr,
                   eta_x_v_e = etaVec * pp$v_e - corr)
    one <- function(rep) {
        d <- d0; d$y <- simy(SEED + 1000L + rep, WGT)
        mkfit(d)@resp$wtres / sigma_true
    }
    R <- do.call(rbind, mclapply(seq_len(NREP), one, mc.cores = 5L))
    tau_mc <- apply(R, 2, sd); mcse <- tau_mc / sqrt(2 * (NREP - 1))
    z <- apply(cand, 2, function(t2) (sqrt(pmax(t2, 1e-8)) - tau_mc) / mcse)
    out[[an]] <- list(tau_pkg = tau_pkg, tau_mc = tau_mc, mcse = mcse,
                      cand = cand, z = z, eta = etaVec, w = WGT)
    cat(sprintf("\n== %s ==  (package tau_e matches candidate `one` to %.2e)\n",
                an, max(abs(tau_pkg - sqrt(pmax(cand[, "one"], 1e-8))))))
    for (cn in colnames(z))
        cat(sprintf("   leading %-10s max|z| %8.2f  mean z %8.2f  |z|>3: %3d/%d\n",
                    cn, max(abs(z[, cn])), mean(z[, cn]), sum(abs(z[, cn]) > 3), n))
}
saveRDS(list(arms = out, nrep = NREP, seed = SEED,
             stamp = c(robustlmm = as.character(packageVersion("robustlmm")),
                       lme4 = as.character(packageVersion("lme4")),
                       R = R.version.string)),
        file.path(outdir, "eta-tau-mc.rds"))
