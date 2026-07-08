## WS11 PLAN section 5.5 (companion): how the Mallows efficiency cost
## scales with the number of CONTINUOUS covariates p*. The single-x
## table (mallowsEfficiency.R) shows a mild ~1-4% cost; that understates
## the picture, because Mallows-type weighting is "notoriously"
## inefficient mainly as the covariate dimension grows -- more leverage
## directions, a noisier covMcd, and information removed in p* dimensions
## per downweighted point. This run holds the tuning at the default
## (gamma = 1, cutoff = 0.975) and sweeps p* so users can see the cost
## they pay for leverage robustness in their own design's dimension.

suppressMessages({
    library(robustlmm); library(lme4); library(parallel)
})

J <- 30L; m <- 6L; n <- J * m
sigma_b <- 1; sigma_e <- 1
gamma <- 1; cutoff <- 0.975
pstar_grid <- c(1L, 2L, 4L, 6L)
nreps  <- 1500L
ncores <- 10L
seed0  <- 20260618L

mkForm <- function(pstar)
    as.formula(paste0("y ~ ", paste0("x", seq_len(pstar), collapse = " + "),
                      " + (1 | g)"))

oneRep <- function(r, pstar) {
    set.seed(seed0 * 7L + r * 11L + pstar)
    g  <- factor(rep(seq_len(J), each = m))
    Xc <- matrix(rnorm(n * pstar), n, pstar)
    colnames(Xc) <- paste0("x", seq_len(pstar))
    beta <- rep(1, pstar)
    b  <- rnorm(J, 0, sigma_b)[as.integer(g)]
    y  <- as.numeric(Xc %*% beta) + b + rnorm(n, 0, sigma_e)
    d  <- data.frame(y = y, g = g); d[colnames(Xc)] <- as.data.frame(Xc)
    fo <- mkForm(pstar)
    mcd <- robustbase::covMcd(Xc)
    d2  <- stats::mahalanobis(Xc, mcd$center, mcd$cov)
    eta <- pmin(1, (stats::qchisq(cutoff, df = pstar) /
                    pmax(d2, .Machine$double.eps))^(gamma / 2))
    grab <- function(f, rob) if (is.null(f)) rep(NA_real_, pstar) else
        unname((if (rob) robustlmm:::.fixef(f) else lme4::fixef(f))[-1])
    fr <- tryCatch(suppressMessages(suppressWarnings(rlmer(fo, d))),
                   error = function(e) NULL)
    fm <- tryCatch(suppressMessages(suppressWarnings(
              rlmer(fo, d, design.weights = eta))), error = function(e) NULL)
    c(RSE  = sum((grab(fr, TRUE) - beta)^2),
      Mall = sum((grab(fm, TRUE) - beta)^2),
      pdown = mean(eta < 1))
}

t0 <- Sys.time()
rows <- list()
for (pstar in pstar_grid) {
    res <- mclapply(seq_len(nreps), function(r) oneRep(r, pstar),
                    mc.cores = ncores)
    M <- do.call(rbind, res)
    ## sum-of-squared-error over the pstar slopes; RE = mean(SSE_RSE)/mean(SSE_M)
    rows[[length(rows) + 1L]] <- data.frame(
        pstar    = pstar,
        pct_down = round(100 * mean(M[, "pdown"], na.rm = TRUE), 1),
        RE_vs_RSE = round(mean(M[, "RSE"], na.rm = TRUE) /
                          mean(M[, "Mall"], na.rm = TRUE), 3))
}
tab <- do.call(rbind, rows)
runtime <- as.numeric(Sys.time() - t0, units = "mins")
cat(sprintf("ran %d reps/pstar on %d cores in %.1f min (gamma=1, cutoff=0.975)\n",
            nreps, ncores, runtime))
print(tab, row.names = FALSE)
saveRDS(list(tab = tab, gamma = gamma, cutoff = cutoff, nreps = nreps,
             runtime_min = runtime),
        file.path("inst", "simulationStudy", "mallowsEfficiencyDim_results.rds"))
cat("saved inst/simulationStudy/mallowsEfficiencyDim_results.rds\n")
