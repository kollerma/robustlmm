## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## T2 validation: thresholded vs raw bootstrap generator for
## anova(test = "boot") (warp-speed Monte Carlo, Giacomini, Politis and
## White 2013). The decision rule was fixed before the run.
##
## DGP: y = 1 + 0.5 x + b0_j + b1_j x + e, b0 ~ N(0, 1), b1 ~ N(0, s_x^2),
## e ~ N(0, 1); x, z ~ N(0, 1) per observation; m = 8 per group.
## H0: the z-slope variance is zero. Nuisance: the x-slope variance
## (s_x = 0 puts it on the boundary).
## fit0: y ~ x + (1|g) + (0+x|g); fit1 adds (0+z|g).
## Per replicate: D_obs, and one bootstrap draw from each generator.
## Size at 5%: mean(D_obs > q_0.95 of the pooled draws of that generator).

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 2000L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "t2_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

one_rep <- function(seed, J, sx, m = 8L) {
    set.seed(seed)
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(J * m); z <- rnorm(J * m)
    b0 <- rnorm(J); b1 <- rnorm(J, sd = sx)
    d <- data.frame(g = g, x = x, z = z,
                    y = 1 + 0.5 * x + b0[g] + b1[g] * x + rnorm(J * m))
    q <- function(expr) tryCatch(suppressMessages(suppressWarnings(expr)),
                                 error = function(e) NULL)
    f0 <- q(rlmer(y ~ x + (1 | g) + (0 + x | g), d, method = "DASvar"))
    f1 <- q(rlmer(y ~ x + (1 | g) + (0 + x | g) + (0 + z | g), d,
                  method = "DASvar"))
    if (is.null(f0) || is.null(f1))
        return(c(D_obs = NA, D_raw = NA, D_thr = NA, trig = NA))
    qd <- robustlmm:::.quasi_deviance
    s0 <- robustlmm:::.sigma(f0)
    D_obs <- qd(f0, scale = s0) - qd(f1, scale = s0)
    pp0 <- f0@pp
    th0 <- getME(f0, "theta"); th1 <- getME(f1, "theta")
    Xb0 <- as.numeric(pp0$X %*% robustlmm:::.fixef(f0))
    Z0 <- t(pp0$Zt)
    dat <- robustlmm:::.boot_data_frame(f0, f1)
    thr <- robustlmm:::.boot_threshold_theta(f0, th0)
    Ub_thr <- pp0$U_b
    if (length(thr$zeroed)) Ub_thr@x <- thr$theta[pp0$Lind]
    draw <- function(Ub) {
        ds <- robustlmm:::.boot_sim_data(dat, "y", Xb0, Z0, Ub, s0,
                                         pp0$n, pp0$q)
        g0 <- robustlmm:::.boot_refit(f0, ds, th0)
        g1 <- robustlmm:::.boot_refit(f1, ds, th1)
        if (is.null(g0) || is.null(g1)) return(NA_real_)
        sg <- robustlmm:::.sigma(g0)
        qd(g0, scale = sg) - qd(g1, scale = sg)
    }
    c(D_obs = D_obs, D_raw = draw(pp0$U_b), D_thr = draw(Ub_thr),
      trig = as.numeric(length(thr$zeroed) > 0))
}

cells <- expand.grid(J = c(18L, 50L), sx = c(0, 0.5))
library(parallel)
cl <- makeCluster(ncores)
clusterCall(cl, function() {
    suppressMessages(library(robustlmm))
    NULL
})
clusterExport(cl, "one_rep")
for (i in seq_len(nrow(cells))) {
    J <- cells$J[i]; sx <- cells$sx[i]
    f <- file.path(outdir, sprintf("cell_J%d_sx%.1f.rds", J, sx))
    if (file.exists(f)) next
    t0 <- proc.time()
    seeds <- 100000L * i + seq_len(R_rep)
    res <- parLapplyLB(cl, seeds, function(s, J, sx) one_rep(s, J, sx),
                       J = J, sx = sx)
    res <- do.call(rbind, res)
    el <- proc.time() - t0
    saveRDS(list(J = J, sx = sx, res = res, wall = el[["elapsed"]],
                 ncores = ncores), f)
    cat(sprintf("cell J=%d sx=%.1f done in %.0f s\n", J, sx, el[["elapsed"]]))
}
stopCluster(cl)

## summary
summ <- do.call(rbind, lapply(list.files(outdir, "^cell_.*\\.rds$",
                                         full.names = TRUE), function(f) {
    o <- readRDS(f); r <- o$res
    ok <- stats::complete.cases(r)
    r <- r[ok, , drop = FALSE]
    size <- function(Db) mean(r[, "D_obs"] > stats::quantile(Db, 0.95))
    data.frame(J = o$J, sx = o$sx, n_ok = nrow(r), n_fail = sum(!ok),
               size_raw = size(r[, "D_raw"]), size_thr = size(r[, "D_thr"]),
               mcse = sqrt(0.05 * 0.95 / nrow(r)),
               trigger_rate = mean(r[, "trig"]),
               wall_s = o$wall, cpu_min = o$wall * o$ncores / 60)
}))
print(summ, digits = 3)
saveRDS(summ, file.path(outdir, "summary.rds"))
