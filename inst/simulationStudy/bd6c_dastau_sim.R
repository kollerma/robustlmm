## Shipped with robustlmm: a script of the breakdown study behind
## setting = "RSEhb" in vignette("simulationStudies"). Run it with Rscript
## from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## BD6c: BD6 (bd6_efficiency_sim.R) re-run with the default method =
## "DAStau" on the final build: the numbers in
## ?rlmer and NEWS come from this run. Same DGP, cells and seeds as BD6; HB
## is setting = "RSEhb" (the shipped bundle, including init = "robust" with
## its fixed-seed S-screen), and the sandwich intervals request
## df = "cluster" explicitly.
##
## DGP: y = 2 + x + b_j + e, m = 10, b ~ N(0, 1), e ~ N(0, 1), x ~ N(0, 1);
## model y ~ x + (1 | g), method DAStau. Cells: clean J = 20 and J = 50
## (500 replicates), and at J = 50 with 10% vertical outliers (+8) and with
## 10% outlying clusters (+8) (300 replicates).
## Fits:
##   lmer   lme4::lmer (REML)
##   RSEn   rlmer defaults
##   HB     rho.e = bisquare c = 3.444 (85%), rho.b = bisquarePsi,
##          rho.sigma.e = rho.sigma.b = sScalePsi, init = "robust"
##   HB95   HB with rho.e = bisquarePsi (c = 4.685)
## Per fit: intercept, slope, sigma, sd_b, and 95% CI coverage of intercept
## and slope for Wald-z (default vcov), Wald-t with Satterthwaite df, and
## the CR2 sandwich (confint(vcov_type = "sandwich")).

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 500L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "bd6c_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fits <- c("lmer", "RSEn", "HB", "HB95")

one_rep <- function(seed, J, type) {
    set.seed(seed)
    m <- 10L
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(J * m); b <- rnorm(J); e <- rnorm(J * m)
    if (type == "C8") { bad <- sample.int(J, round(0.1 * J)); b[bad] <- b[bad] + 8 }
    y <- 2 + x + b[g] + e
    if (type == "V8") { bad <- sample.int(J * m, round(0.1 * J * m)); y[bad] <- y[bad] + 8 }
    d <- data.frame(g = g, x = x, y = y)
    truth <- c(2, 1)
    q <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)),
                                 error = function(e) NULL)
    cover <- function(ci) if (is.null(ci)) c(NA, NA) else
        as.numeric(ci[1:2, 1] <= truth & truth <= ci[1:2, 2])
    rl <- function(...) q(rlmer(y ~ x + (1 | g), d, ...))
    ff <- list(lmer = q(lme4::lmer(y ~ x + (1 | g), d)), RSEn = rl(),
               HB = rl(setting = "RSEhb"),
               HB95 = rl(rho.e = bisquarePsi, rho.b = bisquarePsi,
                         rho.sigma.e = sScalePsi, rho.sigma.b = sScalePsi,
                         init = "robust"))
    unlist(lapply(fits, function(nm) {
        f <- ff[[nm]]
        if (is.null(f)) return(rep(NA_real_, 10))
        isr <- is(f, "rlmerMod")
        est <- c(unname(fixef(f)), sigma(f),
                 sigma(f) * (if (isr) getME(f, "theta") else lme4::getME(f, "theta")))
        if (isr) {
            cz <- cover(q(confint(f, method = "Wald")))
            ct <- cover(q(confint(f, method = "Wald", df = "satterthwaite")))
            cs <- cover(q(confint(f, method = "Wald", vcov_type = "sandwich",
                                  df = "cluster")))
        } else {
            se <- sqrt(diag(as.matrix(vcov(f))))
            cz <- as.numeric(abs(fixef(f) - truth) <= qnorm(0.975) * se)
            ct <- cs <- c(NA, NA)
        }
        c(est, cz, ct, cs)
    }))
}

cells <- data.frame(J = c(20L, 50L, 50L, 50L), type = c("clean", "clean", "V8", "C8"),
                    reps = c(R_rep, R_rep, round(0.6 * R_rep), round(0.6 * R_rep)))
cn <- as.vector(outer(c("int", "slope", "sigma", "sd_b", "covz_int", "covz_slope",
                        "covt_int", "covt_slope", "covs_int", "covs_slope"), fits,
                      function(a, b) paste(b, a, sep = "_")))
if (sys.nframe() == 0L && !isTRUE(getOption("bd6c.source_only"))) {
    library(parallel)
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {
        suppressMessages(library(robustlmm))
        NULL
    })
    clusterExport(cl, c("one_rep", "fits"))
    for (i in seq_len(nrow(cells))) {
        f <- file.path(outdir, sprintf("cell_J%d_%s.rds", cells$J[i], cells$type[i]))
        if (file.exists(f)) next
        t0 <- proc.time()
        res <- parLapplyLB(cl, 4000000L + 10000L * i + seq_len(cells$reps[i]),
                           function(s, J, type) one_rep(s, J, type),
                           J = cells$J[i], type = cells$type[i])
        res <- do.call(rbind, res); colnames(res) <- cn
        el <- proc.time() - t0
        saveRDS(list(J = cells$J[i], type = cells$type[i], res = res,
                     wall = el[["elapsed"]], ncores = ncores), f)
        cat(sprintf("cell J=%d %s done in %.0f s\n", cells$J[i], cells$type[i],
                    el[["elapsed"]]))
    }
    stopCluster(cl)
}
