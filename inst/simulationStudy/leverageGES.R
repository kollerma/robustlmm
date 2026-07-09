## WS11 PLAN section 5.6: the bounded-vs-unbounded influence picture.
##
## Gross-error sensitivity in the design space: the influence of a single
## contaminating observation at (x0, y0) on beta-hat, as x0 moves out
## along the covariate axis. For the plain RSE the psi-functions bound
## the response direction but NOT the design direction, so the slope
## influence grows ~linearly in x0 (unbounded GES in x). Mallows weights
## eta(x0) -> 0 as x0 leaves the bulk, so eta(x0) * x0 stays bounded --
## a bounded GES in x. This is the empirical sensitivity curve
## (refit-based) averaged over base datasets; it produces the
## bounded/unbounded figure for the paper / vignette.

suppressMessages({
    library(robustlmm); library(lme4); library(parallel)
})

J <- 30L; m <- 6L; n <- J * m
beta0 <- 0; beta1 <- 1; sigma_b <- 1; sigma_e <- 1
x0_grid <- c(0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15)
y0_off  <- 8          # contaminating residual offset (saturates psi_e)
nbase   <- 200L       # base datasets to average the sensitivity curve over
ncores  <- 10L
seed0   <- 20260619L

oneBase <- function(r) {
    set.seed(seed0 + r)
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(n)
    b <- rnorm(J, 0, sigma_b)[as.integer(g)]
    y <- beta0 + beta1 * x + b + rnorm(n, 0, sigma_e)
    d <- data.frame(y = y, x = x, g = g)
    fitR <- tryCatch(suppressMessages(suppressWarnings(rlmer(y ~ x + (1|g), d))),
                     error = function(e) NULL)
    fitM <- tryCatch(suppressMessages(suppressWarnings(
                rlmer(y ~ x + (1|g), d, design.weights = "mcd"))),
                error = function(e) NULL)
    if (is.null(fitR) || is.null(fitM)) return(NULL)
    b1R0 <- robustlmm:::.fixef(fitR)[2]
    b1M0 <- robustlmm:::.fixef(fitM)[2]
    ## add one contaminating point to subject 1; y0 far above the line
    scR <- scM <- numeric(length(x0_grid))
    for (j in seq_along(x0_grid)) {
        x0 <- x0_grid[j]
        dc <- rbind(d, data.frame(y = beta0 + beta1 * x0 + y0_off,
                                  x = x0, g = factor(1L, levels = levels(g))))
        fR <- tryCatch(suppressMessages(suppressWarnings(rlmer(y ~ x + (1|g), dc))),
                       error = function(e) NULL)
        fM <- tryCatch(suppressMessages(suppressWarnings(
                  rlmer(y ~ x + (1|g), dc, design.weights = "mcd"))),
                  error = function(e) NULL)
        ## sensitivity curve = (n+1) * change in slope from one added point
        scR[j] <- if (is.null(fR)) NA else (n + 1) * (robustlmm:::.fixef(fR)[2] - b1R0)
        scM[j] <- if (is.null(fM)) NA else (n + 1) * (robustlmm:::.fixef(fM)[2] - b1M0)
    }
    rbind(RSE = scR, Mallows = scM)
}

t0 <- Sys.time()
res <- mclapply(seq_len(nbase), oneBase, mc.cores = ncores)
res <- res[!vapply(res, is.null, logical(1))]
runtime <- as.numeric(Sys.time() - t0, units = "mins")

scR <- rowMeans(vapply(res, function(z) z["RSE", ],     numeric(length(x0_grid))), na.rm = TRUE)
scM <- rowMeans(vapply(res, function(z) z["Mallows", ], numeric(length(x0_grid))), na.rm = TRUE)
curve <- data.frame(x0 = x0_grid, SC_RSE = scR, SC_Mallows = scM)
cat(sprintf("averaged over %d base datasets, %d cores, %.1f min\n",
            length(res), ncores, runtime))
print(curve, row.names = FALSE, digits = 3)
cat(sprintf("\nsup |SC| over the grid:  RSE = %.2f (at x0=%g),  Mallows = %.2f (at x0=%g)\n",
            max(abs(scR)), x0_grid[which.max(abs(scR))],
            max(abs(scM)), x0_grid[which.max(abs(scM))]))

## ---- figure ------------------------------------------------------
pdf(file.path("inst", "simulationStudy", "leverageGES.pdf"),
    width = 6.5, height = 4.5)
op <- par(mar = c(4.2, 4.2, 2.5, 1))
ylim <- range(0, scR, scM, na.rm = TRUE)
plot(x0_grid, scR, type = "b", pch = 19, col = "firebrick",
     ylim = ylim, xlab = expression(leverage~location~x[0]),
     ylab = expression(slope~sensitivity~curve~SC(x[0])),
     main = "Influence of a leverage point: RSE vs Mallows-RSE")
lines(x0_grid, scM, type = "b", pch = 17, col = "steelblue")
abline(h = 0, col = "grey70", lty = 3)
legend("topleft", bty = "n",
       legend = c("RSE (unbounded in x)", "Mallows-RSE (bounded in x)"),
       col = c("firebrick", "steelblue"), pch = c(19, 17), lwd = 1)
par(op); dev.off()

saveRDS(list(curve = curve, config = list(J = J, m = m, x0_grid = x0_grid,
            y0_off = y0_off, nbase = length(res), runtime_min = runtime)),
        file.path("inst", "simulationStudy", "leverageGES_results.rds"))
cat("saved inst/simulationStudy/leverageGES.pdf and leverageGES_results.rds\n")
