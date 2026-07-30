## WS24 / bug 6 -- version sweep.
##
## One fixed weighted dataset, one fixed model, fitted with
##   lmer(weights = w)                                    -- reference
##   rlmer(weights = w, rho.e = rho.b = cPsi, DASvar)     -- must equal lmer
##   rlmer(weights = w, rho.e = rho.b = cPsi, DAStau)     -- must equal lmer
## At cPsi the RSE reduces to REML, so any deviation from lmer is bias.
##
## Run once per installed robustlmm version:
##   Rscript ws24-version-sweep.R <libpath> <label> <outdir>
## Results are written to <outdir>/sweep-<label>.rds.

args <- commandArgs(trailingOnly = TRUE)
libpath <- args[1]; label <- args[2]; outdir <- args[3]
.libPaths(c(libpath, .libPaths()))

suppressPackageStartupMessages({library(lme4); library(robustlmm)})

## ---- the check dataset (same generator as IF-thread1
## thread1_validation/diagnose_weighted_default_vcov.R) ----------------
MASTER_SEED <- 20260611L
J <- 50L; n_per <- 5L; n <- J * n_per
beta <- c(1, 0.5); sigma <- 1; sigma_b <- 0.5
set.seed(MASTER_SEED)
w_true <- exp(runif(n, log(0.25), log(4)))
x <- rnorm(n)
g <- factor(rep(seq_len(J), each = n_per))
set.seed(MASTER_SEED + 1L)
b <- rnorm(J, sd = sigma_b)
e <- rnorm(n, sd = sigma / sqrt(w_true))
dat <- data.frame(y = beta[1] + beta[2] * x + b[as.integer(g)] + e, x = x, g = g)

form <- y ~ x + (1 | g)

grab <- function(fit) {
    if (is.null(fit)) return(c(sigma = NA_real_, slope = NA_real_, se = NA_real_,
                               sigma_b = NA_real_))
    vc <- tryCatch(as.matrix(vcov(fit)), error = function(e) matrix(NA_real_, 2, 2))
    c(sigma   = as.numeric(sigma(fit)),
      slope   = as.numeric(fixef(fit)[2]),
      se      = sqrt(vc[2, 2]),
      sigma_b = as.numeric(sqrt(unlist(VarCorr(fit))[1])))
}

fitRlmer <- function(method) {
    w <- NULL
    out <- withCallingHandlers(
        tryCatch(rlmer(form, dat, weights = w_true, rho.e = cPsi, rho.b = cPsi,
                       method = method),
                 error = function(e) {w <<- c(w, paste("ERROR:", conditionMessage(e))); NULL}),
        warning = function(cnd) {w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")},
        message = function(cnd) invokeRestart("muffleMessage"))
    list(fit = out, msgs = w)
}

fl <- tryCatch(lmer(form, dat, weights = w_true, REML = TRUE), error = function(e) NULL)
rv <- fitRlmer("DASvar")
rt <- fitRlmer("DAStau")

res <- list(
    label      = label,
    robustlmm  = as.character(packageVersion("robustlmm")),
    lme4       = as.character(packageVersion("lme4")),
    Matrix     = as.character(packageVersion("Matrix")),
    Rversion   = R.version.string,
    seed       = MASTER_SEED,
    weightRange= range(w_true),
    lmer       = grab(fl),
    DASvar     = grab(rv$fit),
    DAStau     = grab(rt$fit),
    DASvar_msgs= rv$msgs,
    DAStau_msgs= rt$msgs)

saveRDS(res, file.path(outdir, paste0("sweep-", label, ".rds")))
cat(sprintf("%-10s robustlmm %-8s lmer sigma %.6f | DASvar %.6f (rel %+.4f) | DAStau %.6f (rel %+.4f)\n",
            label, res$robustlmm, res$lmer["sigma"],
            res$DASvar["sigma"], res$DASvar["sigma"]/res$lmer["sigma"] - 1,
            res$DAStau["sigma"], res$DAStau["sigma"]/res$lmer["sigma"] - 1))
