## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## T1 validation: anova(test = "score") with two added components
## (cone statistic), warp-speed Monte Carlo. The decision rule was fixed
## before the run.
##
## DGP: y = 1 + 0.5 x + b0_j + b1_j x + b2_j z + e, m = 8 per group,
## b0 ~ N(0, 1), e ~ N(0, 1), x, z ~ N(0, 1) per observation.
##   null:   b1 = b2 = 0
##   contam: null, plus +5 on every observation of 10% of the groups
##   power:  var(b1) = var(b2) = 0.25
## Null model y ~ x + (1|g); tested directions x and z.
## Per replicate: T_obs and one bootstrap draw T* (null refit only).

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 2000L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "t1_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

one_rep <- function(seed, J, scen, m = 8L) {
    set.seed(seed)
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(J * m); z <- rnorm(J * m)
    sd_s <- if (scen == "power") 0.5 else 0
    y <- 1 + 0.5 * x + rnorm(J)[g] + rnorm(J, sd = sd_s)[g] * x +
        rnorm(J, sd = sd_s)[g] * z + rnorm(J * m)
    if (scen == "contam") {
        bad <- sample.int(J, max(1L, round(0.1 * J)))
        y[g %in% bad] <- y[g %in% bad] + 5
    }
    d <- data.frame(g = g, x = x, z = z, y = y)
    f0 <- tryCatch(suppressMessages(suppressWarnings(
        rlmer(y ~ x + (1 | g), d, method = "DASvar"))),
        error = function(e) NULL)
    if (is.null(f0)) return(c(T_obs = NA, T_star = NA, mec = NA))
    Zadd <- cbind(x = x, z = z)
    consts <- robustlmm:::.score_multi_consts(f0@rho.e)
    st <- robustlmm:::.anova_score_stat_multi(f0, Zadd, consts)
    pp0 <- f0@pp
    ds <- robustlmm:::.boot_sim_data(f0@frame, "y",
                                     as.numeric(pp0$X %*% robustlmm:::.fixef(f0)),
                                     t(pp0$Zt), pp0$U_b, robustlmm:::.sigma(f0),
                                     pp0$n, pp0$q)
    fs <- robustlmm:::.boot_refit(f0, ds, getME(f0, "theta"))
    Ts <- if (is.null(fs)) NA_real_ else
        tryCatch(robustlmm:::.anova_score_stat_multi(fs, Zadd, consts)$T,
                 error = function(e) NA_real_)
    c(T_obs = st$T, T_star = Ts, mec = st$min_eig_cor)
}

cells <- expand.grid(J = c(18L, 50L), scen = c("null", "contam", "power"),
                     stringsAsFactors = FALSE)
library(parallel)
cl <- makeCluster(ncores)
clusterCall(cl, function() {
    suppressMessages(library(robustlmm))
    NULL
})
clusterExport(cl, "one_rep")
for (i in seq_len(nrow(cells))) {
    J <- cells$J[i]; scen <- cells$scen[i]
    f <- file.path(outdir, sprintf("cell_J%d_%s.rds", J, scen))
    if (file.exists(f)) next
    t0 <- proc.time()
    seeds <- 700000L + 10000L * i + seq_len(R_rep)
    res <- parLapplyLB(cl, seeds, function(s, J, scen) one_rep(s, J, scen),
                       J = J, scen = scen)
    res <- do.call(rbind, res)
    el <- proc.time() - t0
    saveRDS(list(J = J, scen = scen, res = res, wall = el[["elapsed"]],
                 ncores = ncores), f)
    cat(sprintf("cell J=%d %s done in %.0f s\n", J, scen, el[["elapsed"]]))
}
stopCluster(cl)

## Warp-speed sizes: each cell's critical value is the 0.95 quantile of
## its own pooled bootstrap draws (reject_own). (An earlier version of this
## comment said the null cell's critical value was used for the
## contaminated and power cells; the code has always used each cell's own
## draws.)
summ <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    o <- readRDS(file.path(outdir, sprintf("cell_J%d_%s.rds", cells$J[i],
                                           cells$scen[i])))
    r <- o$res[stats::complete.cases(o$res), , drop = FALSE]
    crit_own <- stats::quantile(r[, "T_star"], 0.95)
    data.frame(J = o$J, scen = o$scen, n_ok = nrow(r),
               n_fail = nrow(o$res) - nrow(r),
               reject_own = mean(r[, "T_obs"] > crit_own),
               mcse = sqrt(0.05 * 0.95 / nrow(r)),
               atom0 = mean(r[, "T_obs"] <= 1e-10),
               near_singular = mean(r[, "mec"] < 0.01),
               wall_s = o$wall, cpu_min = o$wall * o$ncores / 60)
}))
print(summ, digits = 3)
saveRDS(summ, file.path(outdir, "summary.rds"))
