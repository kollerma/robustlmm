## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## T4 validation: centred M-scale vs Qn for RANSAC candidate selection.
## The decision rule was fixed before the run.
##
## DGP: y = 2 + x + b_j + e, J = 20 groups, m = 10, b ~ N(0, 1),
## e ~ N(0, 1), x ~ N(0, 1). Contaminated cell: 20% of the observations
## (at random) get e + 8. Clean cell: no contamination.
## Per replicate and criterion: ransac_lme4 (K = 200, adaptive) and the
## final rlmer_ransac fit (bisquare rho.e, same seed); record the
## intercept and slope of the selected start and of the final fit.

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 200L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "t4_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

one_rep <- function(seed, contam, J = 20L, m = 10L) {
    set.seed(seed)
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(J * m)
    e <- rnorm(J * m)
    if (contam) {
        out <- sample.int(J * m, round(0.2 * J * m))
        e[out] <- e[out] + 8
    }
    d <- data.frame(g = g, x = x, y = 2 + x + rnorm(J)[g] + e)
    res <- list()
    for (crit in c("centred", "Qn")) {
        sf <- if (crit == "Qn") robustbase::Qn else NULL
        st <- tryCatch(suppressWarnings(
            ransac_lme4(y ~ x + (1 | g), d, K = 200L, seed = seed,
                        scale_fn = sf)), error = function(e) NULL)
        b_start <- if (is.null(st)) c(NA, NA) else unname(lme4::fixef(st$fit))
        ## final fit from that start (rlmer with the RANSAC init)
        fin <- if (is.null(st)) NULL else tryCatch(suppressWarnings(
            rlmer(y ~ x + (1 | g), d, method = "DASvar",
                  rho.e = bisquarePsi,
                  rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                  rho.b = chgDefaults(smoothPsi, k = 5.14),
                  rho.sigma.b = chgDefaults(smoothPsi, k = 5.14),
                  init = robustlmm:::.ransacInitList(st$fit))),
            error = function(e) NULL)
        b_fin <- if (is.null(fin)) c(NA, NA) else unname(fixef(fin))
        res[[crit]] <- c(b_start, b_fin)
    }
    c(centred = res$centred, Qn = res$Qn)
}

library(parallel)
cl <- makeCluster(ncores)
clusterCall(cl, function() {
    suppressMessages(library(robustlmm))
    NULL
})
clusterExport(cl, "one_rep")
for (contam in c(TRUE, FALSE)) {
    f <- file.path(outdir, sprintf("cell_contam%d.rds", contam))
    if (file.exists(f)) next
    t0 <- proc.time()
    seeds <- (if (contam) 500000L else 600000L) + seq_len(R_rep)
    res <- parLapplyLB(cl, seeds, function(s, contam) one_rep(s, contam),
                       contam = contam)
    res <- do.call(rbind, res)
    colnames(res) <- c(paste0("centred_", c("int_start", "slope_start",
                                            "int_final", "slope_final")),
                       paste0("Qn_", c("int_start", "slope_start",
                                       "int_final", "slope_final")))
    el <- proc.time() - t0
    saveRDS(list(contam = contam, res = res, wall = el[["elapsed"]],
                 ncores = ncores), f)
    cat(sprintf("cell contam=%d done in %.0f s\n", contam, el[["elapsed"]]))
}
stopCluster(cl)

summ <- do.call(rbind, lapply(list.files(outdir, "^cell_.*\\.rds$",
                                         full.names = TRUE), function(f) {
    o <- readRDS(f); r <- o$res
    err <- function(col, truth) abs(r[, col] - truth)
    rmse <- function(col, truth) sqrt(mean((r[, col] - truth)^2, na.rm = TRUE))
    data.frame(contam = o$contam,
               n = nrow(r),
               med_int_err_start_centred = median(err("centred_int_start", 2), na.rm = TRUE),
               med_int_err_start_Qn = median(err("Qn_int_start", 2), na.rm = TRUE),
               med_int_err_final_centred = median(err("centred_int_final", 2), na.rm = TRUE),
               med_int_err_final_Qn = median(err("Qn_int_final", 2), na.rm = TRUE),
               rmse_int_start_centred = rmse("centred_int_start", 2),
               rmse_int_start_Qn = rmse("Qn_int_start", 2),
               rmse_slope_start_centred = rmse("centred_slope_start", 1),
               rmse_slope_start_Qn = rmse("Qn_slope_start", 1),
               wall_s = o$wall, cpu_min = o$wall * o$ncores / 60)
}))
print(t(summ), digits = 3)
saveRDS(summ, file.path(outdir, "summary.rds"))
