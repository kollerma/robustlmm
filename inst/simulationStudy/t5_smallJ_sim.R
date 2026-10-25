## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## T5: small-J coverage of sandwich Wald intervals. The decision rule was fixed
## before the run.
##
## y = 250 + 10 t + b0_j + b1_j t + e, t = 0..9, sd(b0) = 25, sd(b1) = 6,
## sd(e) = 25; model y ~ t + (t | g), DASvar.
## Candidates (95% intervals for intercept and slope):
##   a: vcov_sandwich (G1) + normal quantile      (current)
##   b: vcov_sandwich (G1) + t(J - 1)
##   c: full-IF sandwich (G1; sigma/theta estimation included) + t(J - 1)
##   d: leave-one-cluster-out jackknife + t(J - 1)
##   ref: default vcov + Satterthwaite t (confint(df = "satterthwaite"))

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 1000L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "t5_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

one_rep <- function(seed, J) {
    set.seed(seed)
    m <- 10L
    g <- factor(rep(seq_len(J), each = m))
    t <- rep(0:9, J)
    y <- 250 + 10 * t + rnorm(J, sd = 25)[g] + rnorm(J, sd = 6)[g] * t +
        rnorm(J * m, sd = 25)
    d <- data.frame(g = g, t = t, y = y)
    q <- function(e) tryCatch(suppressMessages(suppressWarnings(e)),
                              error = function(e) NULL)
    fit <- q(rlmer(y ~ t + (t | g), d, method = "DASvar"))
    na <- rep(NA, 10)
    if (is.null(fit)) return(na)
    b <- fixef(fit); truth <- c(250, 10)
    cover <- function(se, crit) as.numeric(abs(b - truth) <= crit * se)
    tq <- qt(0.975, J - 1); zq <- qnorm(0.975)
    ## a, b
    Vs <- q(vcov_sandwich(fit))
    se_s <- if (is.null(Vs)) c(NA, NA) else sqrt(diag(Vs))
    ## c: full IF (beta block), G1
    se_f <- q({
        IF <- robustlmm:::implicitIF_full(fit)
        S  <- robustlmm:::.scoreByCluster(fit)
        Xi <- robustlmm:::.applyJparInv(IF, S)
        A  <- tcrossprod(Xi[1:2, , drop = FALSE]) * J / (J - 1)
        sqrt(diag(A))
    })
    if (is.null(se_f)) se_f <- c(NA, NA)
    ## d: jackknife over clusters
    bj <- vapply(levels(g), function(lv) {
        f <- q(rlmer(y ~ t + (t | g), droplevels(d[d$g != lv, ]),
                     method = "DASvar"))
        if (is.null(f)) c(NA_real_, NA_real_) else unname(fixef(f))
    }, numeric(2))
    se_j <- if (anyNA(bj)) c(NA, NA) else
        sqrt(diag((J - 1) / J * tcrossprod(bj - rowMeans(bj))))
    ## ref: default vcov + Satterthwaite
    ci_ref <- q(confint(fit, method = "Wald", df = "satterthwaite"))
    cov_ref <- if (is.null(ci_ref)) c(NA, NA) else
        as.numeric(ci_ref[, 1] <= truth & truth <= ci_ref[, 2])
    c(a = cover(se_s, zq), b = cover(se_s, tq), c = cover(se_f, tq),
      d = cover(se_j, tq), ref = cov_ref)
}

Js <- c(8L, 12L, 18L, 30L)
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
    res <- parLapplyLB(cl, 800000L + 1000L * J + seq_len(R_rep),
                       function(s, J) one_rep(s, J), J = J)
    res <- do.call(rbind, res)
    colnames(res) <- paste(rep(c("a", "b", "c", "d", "ref"), each = 2),
                           c("int", "slope"), sep = "_")
    el <- proc.time() - t0
    saveRDS(list(J = J, res = res, wall = el[["elapsed"]], ncores = ncores), f)
    cat(sprintf("cell J=%d done in %.0f s\n", J, el[["elapsed"]]))
}
stopCluster(cl)

summ <- do.call(rbind, lapply(Js, function(J) {
    o <- readRDS(file.path(outdir, sprintf("cell_J%d.rds", J)))
    cv <- colMeans(o$res, na.rm = TRUE)
    data.frame(J = J, t(round(cv, 3)), n_na = sum(!stats::complete.cases(o$res)),
               cpu_min = o$wall * o$ncores / 60)
}))
print(summ)
saveRDS(summ, file.path(outdir, "summary.rds"))
