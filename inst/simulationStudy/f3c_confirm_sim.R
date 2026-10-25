## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## F3c: confirmatory study for the sandwich default (CR2 + approximate
## Hotelling) on a design independent of F3/F3b. The decision rule was fixed before
## the run. End-to-end through the package: anova(vcov_type = "sandwich", ddf = "none" / "cluster") and
## confint(vcov_type = "sandwich", df = "cluster").
##
## Design: J in {10, 20, 40}, unbalanced cluster sizes m_j ~ U{4, ..., 16};
## y = 1 + 0.5 x + b0_j + b1_j x + e, sd(b0) = 1, sd(b1) = 0.5,
## e ~ t_4 / sqrt(2) (unit variance, heavy tails); null terms:
##   z  between-cluster binary, exactly half the clusters     (k = 1)
##   u  between-cluster log-normal exp(N(0, 1)) (skewed)       (k = 1)
##   v  within-cluster N(0, 1)                                  (k = 1)
##   zu joint restriction z = u = 0, anova(fit0, fit1)          (k = 2)
## Model y ~ x + z + u + v + (x | g), DASvar, default psi-functions.

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 1000L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "f3c_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

one_rep <- function(seed, J) {
    set.seed(seed)
    m <- sample(4:16, J, replace = TRUE)
    g <- factor(rep(seq_len(J), m)); n <- length(g)
    x <- rnorm(n); v <- rnorm(n)
    z <- sample(rep(0:1, J / 2))[g]
    u <- exp(rnorm(J))[g]
    y <- 1 + 0.5 * x + rnorm(J)[g] + rnorm(J, sd = 0.5)[g] * x +
        rt(n, 4) / sqrt(2)
    d <- data.frame(g = g, x = x, v = v, z = z, u = u, y = y)
    q <- function(e) tryCatch(suppressMessages(suppressWarnings(e)),
                              error = function(e) NULL)
    na <- rep(NA_real_, 16)
    fit <- q(rlmer(y ~ x + z + u + v + (x | g), d, method = "DASvar"))
    fit0 <- q(rlmer(y ~ x + v + (x | g), d, method = "DASvar"))
    if (is.null(fit) || is.null(fit0)) return(na)
    a_old <- q(anova(fit, vcov_type = "sandwich", ddf = "none"))
    a_new <- q(anova(fit, vcov_type = "sandwich", ddf = "cluster"))
    p_old <- q(anova(fit0, fit, vcov_type = "sandwich", ddf = "none"))
    p_new <- q(anova(fit0, fit, vcov_type = "sandwich", ddf = "cluster"))
    pv <- function(a, col, rows) if (is.null(a)) rep(NA_real_, length(rows)) else
        a[rows, col]
    ## coverage of intercept and x slope
    truth <- c(1, 0.5)
    b <- fixef(fit)[1:2]
    ci_new <- q(confint(fit, vcov_type = "sandwich", df = "cluster"))
    Vg <- q(as.matrix(vcov(fit, type = "sandwich")))
    cov_new <- if (is.null(ci_new)) c(NA, NA) else
        as.numeric(ci_new[1:2, 1] <= truth & truth <= ci_new[1:2, 2])
    se_g <- if (is.null(Vg)) c(NA, NA) else sqrt(diag(Vg))[1:2]
    cov_gz <- as.numeric(abs(b - truth) <= qnorm(0.975) * se_g)
    cov_gt <- as.numeric(abs(b - truth) <= qt(0.975, J - 1) * se_g)
    c(pv(a_old, "Pr(>Chisq)", c("z", "u", "v")), pv(p_old, "Pr(>Chisq)", "fit1"),
      pv(a_new, "Pr(>F)", c("z", "u", "v")), pv(p_new, "Pr(>F)", "fit1"),
      cov_new, cov_gz, cov_gt)
}

cn <- c(paste0("old_", c("z", "u", "v", "zu")), paste0("new_", c("z", "u", "v", "zu")),
        "covCR2_int", "covCR2_x", "covG1z_int", "covG1z_x", "covG1t_int", "covG1t_x")
Js <- c(10L, 20L, 40L)
if (sys.nframe() == 0L && !isTRUE(getOption("f3c.source_only"))) {
    library(parallel)
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {
        suppressMessages(library(robustlmm))
        NULL
    })
    clusterExport(cl, "one_rep")
    for (J in Js) {
        f <- file.path(outdir, sprintf("cell_J%d.rds", J))
        if (file.exists(f)) next
        t0 <- proc.time()
        res <- parLapplyLB(cl, 1500000L + 10000L * J + seq_len(R_rep),
                           function(s, J) one_rep(s, J), J = J)
        res <- do.call(rbind, res); colnames(res) <- cn
        el <- proc.time() - t0
        saveRDS(list(J = J, res = res, wall = el[["elapsed"]], ncores = ncores), f)
        cat(sprintf("cell J=%d done in %.0f s\n", J, el[["elapsed"]]))
    }
    stopCluster(cl)
    summ <- do.call(rbind, lapply(Js, function(J) {
        o <- readRDS(file.path(outdir, sprintf("cell_J%d.rds", J))); r <- o$res
        rej <- function(v) mean(v < 0.05, na.rm = TRUE)
        do.call(rbind, lapply(c("z", "u", "v", "zu"), function(tt)
            data.frame(J = J, term = tt,
                       n_ok = sum(is.finite(r[, paste0("new_", tt)])),
                       old = rej(r[, paste0("old_", tt)]),
                       new = rej(r[, paste0("new_", tt)]),
                       covCR2_int = mean(r[, "covCR2_int"], na.rm = TRUE),
                       covCR2_x = mean(r[, "covCR2_x"], na.rm = TRUE),
                       covG1z_int = mean(r[, "covG1z_int"], na.rm = TRUE),
                       covG1t_int = mean(r[, "covG1t_int"], na.rm = TRUE),
                       covG1t_x = mean(r[, "covG1t_x"], na.rm = TRUE),
                       cpu_min = o$wall * o$ncores / 60)))
    }))
    print(summ, digits = 3)
    saveRDS(summ, file.path(outdir, "summary.rds"))
}
