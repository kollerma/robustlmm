## WS11 PLAN section 5.5: efficiency of Mallows-RSE at the clean model.
##
## Mallows weighting buys leverage robustness (section 5.4) at the cost
## of efficiency: it downweights the natural tail of the covariate
## distribution even when the data are clean. This study quantifies that
## cost across the tuning grid (gamma in {1, 2}, cutoff in {0.95, 0.975,
## 0.99}) so a default can be chosen and the loss documented for users.
##
## Relative efficiency is reported against lmer (the MLE, fully efficient
## at the clean Gaussian model): RE = Var(lmer beta-hat) / Var(. beta-hat),
## so RE < 1 means less efficient. We also report the marginal cost
## versus the plain RSE (Var(RSE)/Var(Mallows)) -- the price of adding
## leverage protection on top of the existing robustness.

suppressMessages({
    library(robustlmm); library(lme4); library(parallel)
})

## ---- configuration (clean model; same design as leverageBreakdown) -
J        <- 30L; m <- 6L; n <- J * m
beta0    <- 0; beta1 <- 1
sigma_b  <- 1; sigma_e <- 1
gammas   <- c(1, 2)
cutoffs  <- c(0.95, 0.975, 0.99)
nreps    <- 2000L
ncores   <- 10L
seed0    <- 20260618L

settings <- expand.grid(gamma = gammas, cutoff = cutoffs,
                        stringsAsFactors = FALSE)

makeClean <- function() {
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(n)
    b <- rnorm(J, 0, sigma_b)[as.integer(g)]
    y <- beta0 + beta1 * x + b + rnorm(n, 0, sigma_e)
    data.frame(y = y, x = x, g = g)
}

## mcd squared Mahalanobis distances of the non-constant design columns
## (here: x), then eta = min(1, (qchisq(cutoff, p*) / d2)^(gamma/2)).
mcdEta <- function(x, gamma, cutoff, pstar = 1L) {
    mcd <- robustbase::covMcd(cbind(x))
    d2  <- stats::mahalanobis(cbind(x), mcd$center, mcd$cov)
    pmin(1, (stats::qchisq(cutoff, df = pstar) /
             pmax(d2, .Machine$double.eps))^(gamma / 2))
}

b1 <- function(f, robust) {
    if (is.null(f)) return(NA_real_)
    unname((if (robust) robustlmm:::.fixef(f) else lme4::fixef(f))[2])
}

oneRep <- function(r) {
    set.seed(seed0 + r)
    d  <- makeClean()
    fl <- tryCatch(suppressMessages(suppressWarnings(
              lmer(y ~ x + (1 | g), d, REML = TRUE))), error = function(e) NULL)
    fr <- tryCatch(suppressMessages(suppressWarnings(
              rlmer(y ~ x + (1 | g), d))), error = function(e) NULL)
    out <- c(lmer = b1(fl, FALSE), RSE = b1(fr, TRUE))
    ndown <- numeric(nrow(settings))
    mall  <- numeric(nrow(settings))
    for (i in seq_len(nrow(settings))) {
        eta <- mcdEta(d$x, settings$gamma[i], settings$cutoff[i])
        ndown[i] <- mean(eta < 1)
        fm <- tryCatch(suppressMessages(suppressWarnings(
                  rlmer(y ~ x + (1 | g), d, design.weights = eta))),
                  error = function(e) NULL)
        mall[i] <- b1(fm, TRUE)
    }
    list(b1 = c(out, setNames(mall, paste0("M", seq_len(nrow(settings))))),
         ndown = ndown)
}

t0   <- Sys.time()
reps <- mclapply(seq_len(nreps), oneRep, mc.cores = ncores)
runtime <- as.numeric(Sys.time() - t0, units = "mins")
cat(sprintf("ran %d reps on %d cores in %.1f min\n", nreps, ncores, runtime))

B    <- do.call(rbind, lapply(reps, function(z) z$b1))
ND   <- do.call(rbind, lapply(reps, function(z) z$ndown))
varV <- apply(B, 2L, var, na.rm = TRUE)
mseV <- apply(B, 2L, function(v) mean((v - beta1)^2, na.rm = TRUE))

tab <- data.frame(
    method  = c("lmer", "RSE", sprintf("Mallows g%d c%.3f",
                                       settings$gamma, settings$cutoff)),
    pct_down = c(NA, NA, round(100 * colMeans(ND), 1)),
    RE_vs_lmer = round(varV["lmer"] / varV, 3),
    RE_vs_RSE  = round(varV["RSE"]  / varV, 3),
    rmse_b1    = round(sqrt(mseV), 4),
    row.names = NULL)
print(tab, row.names = FALSE)

saveRDS(list(tab = tab, settings = settings, varV = varV, mseV = mseV,
             config = list(J = J, m = m, nreps = nreps, seed0 = seed0,
                           runtime_min = runtime)),
        file.path("inst", "simulationStudy", "mallowsEfficiency_results.rds"))
cat("saved inst/simulationStudy/mallowsEfficiency_results.rds\n")
