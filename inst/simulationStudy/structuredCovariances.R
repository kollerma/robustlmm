## 3.5.0 structured random-effect covariances: consistency + robustness.
##
## Two structures are studied, each fitted classically (lmer) and
## robustly (rlmer, DASvar):
##   * cs:  3-dim within-subject factor, heterogeneous compound symmetry
##          (marginal sds 2, 1.5, 1.2; one common correlation rho = 0.5).
##   * ar1: 4-dim within-subject factor, homogeneous AR(1)
##          (single sd 1.5; Cor(i, j) = rho^|i - j|, rho = 0.6).
##
## For each structure we report the estimated structural correlation and
## marginal sd under (a) clean data and (b) random-effects-level
## contamination (the contamination model robustlmm targets): a fraction
## eps of the subjects get an outlying random-effect vector drawn
## independently from N(0, (cfac * mean_sd)^2 I) -- large in magnitude
## and uncorrelated across the within-subject dimensions. The story:
##   * clean: lmer and rlmer both recover the truth (the projection keeps
##     the fit on the structured manifold, so it is consistent);
##   * contaminated: the outlying subjects inflate the classical marginal
##     sds and attenuate the common correlation toward 0, while the robust
##     fit downweights them and stays near the truth.
##
## Build-ignored research script; writes structuredCovariances_results.rds.

suppressMessages({
    library(robustlmm); library(lme4); library(parallel)
})

## ---- configuration -------------------------------------------------
nsubj   <- 40L
nrep    <- 5L            # replicates per subject x factor-level cell
sigma_e <- 0.8
eps     <- 0.10          # contamination fraction
cfac    <- 6             # gross-error inflation factor
beta0   <- 10
nreps   <- 500L
ncores  <- 8L
seed0   <- 20260623L

## cs truth: 3 levels, heterogeneous sds, common correlation
cs_sd   <- c(2, 1.5, 1.2)
cs_rho  <- 0.5
## ar1 truth: 4 levels, homogeneous sd, lag-decaying correlation
ar1_sd  <- 1.5
ar1_rho <- 0.6

cholOf <- function(sds, R) t(chol(diag(sds, length(sds)) %*% R %*%
                                  diag(sds, length(sds))))
csR  <- { R <- matrix(cs_rho, 3, 3); diag(R) <- 1; R }
ar1R <- ar1_rho^abs(outer(1:4, 1:4, "-"))
csL  <- cholOf(cs_sd, csR)
ar1L <- cholOf(rep(ar1_sd, 4), ar1R)

## ---- one dataset of a given structure ------------------------------
makeData <- function(L, nlev) {
    f <- factor(letters[seq_len(nlev)])
    d <- expand.grid(f = f, rep = seq_len(nrep), subj = factor(seq_len(nsubj)))
    b <- t(L %*% matrix(rnorm(nlev * nsubj), nlev, nsubj))  # nsubj x nlev
    e <- rnorm(nrow(d), 0, sigma_e)
    idx <- cbind(as.integer(d$subj), as.integer(d$f))
    d$y_clean <- beta0 + b[idx] + e
    ## contaminate the random effects of a fraction of subjects: an
    ## outlying, uncorrelated, large-magnitude RE vector.
    bc <- b
    meanSd <- mean(sqrt(rowSums(L^2)))            # mean marginal RE sd
    ko <- sample.int(nsubj, max(1L, round(eps * nsubj)))
    bc[ko, ] <- matrix(rnorm(length(ko) * nlev, 0, cfac * meanSd),
                       length(ko), nlev)
    d$y_cont <- beta0 + bc[idx] + e
    d
}

## ---- extract (rho, mean marginal sd) from a fit --------------------
## A robust fit that failed to converge within max.it (rare: ~3% of
## contaminated ar1 fits hit a slow non-converging path) is treated as
## missing so it neither pollutes the means nor dominates the runtime.
fitStats <- function(fit) {
    if (is.null(fit)) return(c(rho = NA_real_, sd = NA_real_))
    if (is(fit, "rlmerMod") && !isTRUE(fit@optinfo$conv$opt == 0))
        return(c(rho = NA_real_, sd = NA_real_))
    vc <- VarCorr(fit)[[1L]]
    cm <- attr(vc, "correlation")
    sd <- attr(vc, "stddev")
    c(rho = unname(cm[2, 1]), sd = unname(mean(sd)))
}

fitBoth <- function(form, data, structure) {
    f <- as.formula(form)
    fl <- tryCatch(suppressMessages(suppressWarnings(
        lmer(f, data, REML = TRUE))), error = function(e) NULL)
    ## max.iter caps the rare pathological projection loop (default is
    ## much larger); healthy fits converge in well under this.
    fr <- tryCatch(suppressMessages(suppressWarnings(
        rlmer(f, data, method = "DASvar", max.iter = 100L))),
        error = function(e) NULL)
    c(setNames(fitStats(fl), paste0("lmer.", c("rho", "sd"))),
      setNames(fitStats(fr), paste0("rlmer.", c("rho", "sd"))))
}

oneRep <- function(r) {
    set.seed(seed0 + r)
    dcs  <- makeData(csL, 3L)
    dar1 <- makeData(ar1L, 4L)
    rbind(
        cs_clean  = fitBoth("y_clean ~ 1 + cs(0 + f | subj)",  dcs,  "cs"),
        cs_cont   = fitBoth("y_cont ~ 1 + cs(0 + f | subj)",   dcs,  "cs"),
        ar1_clean = fitBoth("y_clean ~ 1 + ar1(0 + f | subj)", dar1, "ar1"),
        ar1_cont  = fitBoth("y_cont ~ 1 + ar1(0 + f | subj)",  dar1, "ar1"))
}

## ---- run -----------------------------------------------------------
t0   <- Sys.time()
reps <- mclapply(seq_len(nreps), oneRep, mc.cores = ncores)
runtime <- as.numeric(Sys.time() - t0, units = "mins")
ok <- vapply(reps, function(z) is.matrix(z) && all(dim(z) == c(4L, 4L)), NA)
reps <- reps[ok]
cat(sprintf("ran %d/%d reps on %d cores in %.1f min\n",
            sum(ok), nreps, ncores, runtime))

## stack into one array: rep x scenario x column
scen <- rownames(reps[[1]])
cols <- colnames(reps[[1]])
arr  <- array(NA_real_, c(length(reps), 4L, 4L),
              dimnames = list(NULL, scen, cols))
for (i in seq_along(reps)) arr[i, , ] <- reps[[i]]

## ---- truth + summary table -----------------------------------------
truth <- list(
    cs_clean  = c(rho = cs_rho,  sd = mean(cs_sd)),
    cs_cont   = c(rho = cs_rho,  sd = mean(cs_sd)),
    ar1_clean = c(rho = ar1_rho, sd = ar1_sd),
    ar1_cont  = c(rho = ar1_rho, sd = ar1_sd))

summ <- do.call(rbind, lapply(scen, function(s) {
    tr <- truth[[s]]
    m  <- colMeans(arr[, s, ], na.rm = TRUE)
    data.frame(
        scenario   = s,
        true_rho   = tr["rho"],
        lmer_rho   = round(m["lmer.rho"], 3),
        rlmer_rho  = round(m["rlmer.rho"], 3),
        true_sd    = round(tr["sd"], 3),
        lmer_sd    = round(m["lmer.sd"], 3),
        rlmer_sd   = round(m["rlmer.sd"], 3),
        row.names  = NULL)
}))
print(summ, row.names = FALSE)

saveRDS(list(summary = summ, array = arr, truth = truth,
             config = list(nsubj = nsubj, nrep = nrep, sigma_e = sigma_e,
                           eps = eps, cfac = cfac, nreps = sum(ok),
                           seed0 = seed0, runtime_min = runtime,
                           cs_sd = cs_sd, cs_rho = cs_rho,
                           ar1_sd = ar1_sd, ar1_rho = ar1_rho)),
        file.path("inst", "simulationStudy", "structuredCovariances_results.rds"))
cat("saved inst/simulationStudy/structuredCovariances_results.rds\n")
