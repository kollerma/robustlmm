## Shipped with robustlmm: a script of the breakdown study behind
## setting = "RSEhb" in vignette("simulationStudies"). Run it with Rscript
## from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## BD6b: breakdown map of the shipped configuration, rlmer(setting = "RSEhb")
## (S-screen start, bisquare c = 3.444, bisquare rho.b, S-type scales),
## against the default, on the map's cells with fresh seeds, both with the
## default method = "DAStau" (DASvar under-estimates the residual scale by
## about 4.5% for a redescending rho.e). The numbers quoted in ?rlmer come
## from this run.
## make_data() and cells from the BD5 harness; sourced first because it also
## defines R_rep, ncores and outdir, which are set below for this run
source(system.file("simulationStudy", "bd5_starts_sim.R", package = "robustlmm"), local = TRUE)
args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 100L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "bd6b_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
one_rep <- function(seed, type, eps) {
    d <- make_data(seed, type, eps)
    q <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)),
                                 error = function(e) NULL)
    s <- function(f) if (is.null(f)) rep(NA_real_, 4) else
        c(unname(fixef(f)), sigma(f), sigma(f) * getME(f, "theta"))
    c(s(q(rlmer(y ~ x + (1 | g), d))),
      s(q(rlmer(y ~ x + (1 | g), d, setting = "RSEhb"))))
}
if (sys.nframe() == 0L) {
    library(parallel)
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {
        suppressMessages(library(robustlmm)); NULL
    })
    clusterExport(cl, c("one_rep", "make_data"))
    for (i in seq_len(nrow(cells))) {
        f <- file.path(outdir, sprintf("cell_%s_%03d.rds", cells$type[i],
                                       round(100 * cells$eps[i])))
        if (file.exists(f)) next
        t0 <- proc.time()
        res <- parLapplyLB(cl, 5000000L + 10000L * i + seq_len(R_rep),
                           function(s, type, eps) one_rep(s, type, eps),
                           type = cells$type[i], eps = cells$eps[i])
        res <- do.call(rbind, res)
        colnames(res) <- paste(rep(c("RSEn", "RSEhb"), each = 4),
                               c("int", "slope", "sigma", "sd_b"), sep = "_")
        el <- proc.time() - t0
        saveRDS(list(type = cells$type[i], eps = cells$eps[i], res = res,
                     wall = el[["elapsed"]], ncores = ncores), f)
        cat(sprintf("cell %s eps=%.2f done in %.0f s\n", cells$type[i], cells$eps[i],
                    el[["elapsed"]]))
    }
    stopCluster(cl)
}
