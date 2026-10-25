## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## F3d: replication of the sandwich default with fresh seeds. The CR2
## sandwich with approximate-Hotelling df became the default reference for
## vcov_type = "sandwich" through the F3c rule, which is looser than F3b's
## rule and was written after F3b's results. F3d re-tests the default on
## BOTH earlier designs, end-to-end through the package. The decision rule
## (summarise_f3d() below) was fixed before the run.
##
## Design A (F3b's DGP): m = 10, t = 0..9, y = 250 + 10 t + b0 + b1 t + e,
##   sd(b0) = 25, sd(b1) = 6, sd(e) = 25; null terms w (between N(0, 1),
##   k = 1), f (between 3-level factor, levels cycled over clusters, k = 2),
##   x (within N(0, 1), k = 1); model y ~ t + w + f + x + (t | g);
##   J in {8, 12, 18, 30}; coverage of the intercept and t.
## Design B (F3c's DGP): m_j ~ U{4..16}, y = 1 + 0.5 x + b0 + b1 x + e,
##   sd(b0) = 1, sd(b1) = 0.5, e ~ t_4 / sqrt(2); null terms z (between
##   binary, half the clusters), u (between exp(N(0, 1))), v (within N(0, 1)),
##   joint z = u = 0 (anova(fit0, fit1), k = 2); model
##   y ~ x + z + u + v + (x | g); J in {10, 20, 40}; coverage of the
##   intercept and x.
## Both: DASvar, default psi-functions, 1000 replicates per cell; per term
## the p-values of anova(vcov_type = "sandwich", ddf = "none") (G1 +
## chi-square) and ddf = "cluster" (CR2 + approximate Hotelling), and the
## 95% coverage of confint(vcov_type = "sandwich", df = "cluster").
## Seeds: design A 3100000 + 10000 J + r, design B 3600000 + 10000 J + r
## (disjoint from F3 / F3b / F3c / T5 / F2).

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 1000L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "f3d_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

q <- function(e) tryCatch(suppressMessages(suppressWarnings(e)),
                          error = function(e) NULL)
pv <- function(a, col, rows)
    if (is.null(a)) rep(NA_real_, length(rows)) else a[rows, col]
cover <- function(ci, truth)
    if (is.null(ci)) c(NA_real_, NA_real_) else
        as.numeric(ci[1:2, 1] <= truth & truth <= ci[1:2, 2])

rep_A <- function(seed, J) {
    set.seed(seed)
    m <- 10L
    g <- factor(rep(seq_len(J), each = m))
    t <- rep(0:9, J)
    w <- rnorm(J)[g]
    f <- factor(rep_len(c("a", "b", "c"), J))[g]
    x <- rnorm(J * m)
    y <- 250 + 10 * t + rnorm(J, sd = 25)[g] + rnorm(J, sd = 6)[g] * t +
        rnorm(J * m, sd = 25)
    d <- data.frame(g = g, t = t, w = w, f = f, x = x, y = y)
    fit <- q(rlmer(y ~ t + w + f + x + (t | g), d, method = "DASvar"))
    if (is.null(fit)) return(rep(NA_real_, 8))
    a_old <- q(anova(fit, vcov_type = "sandwich", ddf = "none"))
    a_new <- q(anova(fit, vcov_type = "sandwich", ddf = "cluster"))
    ci <- q(confint(fit, vcov_type = "sandwich", df = "cluster"))
    c(pv(a_old, "Pr(>Chisq)", c("w", "f", "x")),
      pv(a_new, "Pr(>F)", c("w", "f", "x")), cover(ci, c(250, 10)))
}

rep_B <- function(seed, J) {
    set.seed(seed)
    m <- sample(4:16, J, replace = TRUE)
    g <- factor(rep(seq_len(J), m)); n <- length(g)
    x <- rnorm(n); v <- rnorm(n)
    z <- sample(rep(0:1, J / 2))[g]
    u <- exp(rnorm(J))[g]
    y <- 1 + 0.5 * x + rnorm(J)[g] + rnorm(J, sd = 0.5)[g] * x +
        rt(n, 4) / sqrt(2)
    d <- data.frame(g = g, x = x, v = v, z = z, u = u, y = y)
    fit <- q(rlmer(y ~ x + z + u + v + (x | g), d, method = "DASvar"))
    fit0 <- q(rlmer(y ~ x + v + (x | g), d, method = "DASvar"))
    if (is.null(fit) || is.null(fit0)) return(rep(NA_real_, 10))
    a_old <- q(anova(fit, vcov_type = "sandwich", ddf = "none"))
    a_new <- q(anova(fit, vcov_type = "sandwich", ddf = "cluster"))
    p_old <- q(anova(fit0, fit, vcov_type = "sandwich", ddf = "none"))
    p_new <- q(anova(fit0, fit, vcov_type = "sandwich", ddf = "cluster"))
    ci <- q(confint(fit, vcov_type = "sandwich", df = "cluster"))
    c(pv(a_old, "Pr(>Chisq)", c("z", "u", "v")), pv(p_old, "Pr(>Chisq)", "fit1"),
      pv(a_new, "Pr(>F)", c("z", "u", "v")), pv(p_new, "Pr(>F)", "fit1"),
      cover(ci, c(1, 0.5)))
}

cells <- rbind(data.frame(design = "A", J = c(8L, 12L, 18L, 30L)),
               data.frame(design = "B", J = c(10L, 20L, 40L)))
terms <- list(A = c("w", "f", "x"), B = c("z", "u", "v", "zu"))
colnames_for <- function(des)
    c(paste0("old_", terms[[des]]), paste0("new_", terms[[des]]),
      "cov_int", "cov_slope")

## The decision rule, evaluated by summarise_f3d()
summarise_f3d <- function(outdir) {
    rows <- list()
    for (i in seq_len(nrow(cells))) {
        des <- cells$design[i]; J <- cells$J[i]
        o <- readRDS(file.path(outdir, sprintf("cell_%s_J%d.rds", des, J)))
        r <- o$res
        for (tt in terms[[des]]) {
            ok <- is.finite(r[, paste0("new_", tt)]) & is.finite(r[, paste0("old_", tt)])
            rows[[length(rows) + 1L]] <- data.frame(
                design = des, J = J, term = tt, n = nrow(r), n_ok = sum(ok),
                old = mean(r[ok, paste0("old_", tt)] < 0.05),
                new = mean(r[ok, paste0("new_", tt)] < 0.05),
                cov_int = mean(r[, "cov_int"], na.rm = TRUE),
                cov_slope = mean(r[, "cov_slope"], na.rm = TRUE),
                cpu_min = o$wall * o$ncores / 60)
        }
    }
    s <- do.call(rbind, rows)
    s$c1 <- abs(s$new - 0.05) < abs(s$old - 0.05) |
        (abs(s$new - 0.05) <= 0.01 & abs(s$old - 0.05) <= 0.01)
    s$c2 <- s$J < 20 | (s$new >= 0.02 & s$new <= 0.08)
    s$c3 <- s$J < 20 | (s$cov_int >= 0.93 & s$cov_slope >= 0.93)
    s$c4 <- s$n_ok >= 0.98 * s$n
    ## F3b's stricter band, reported only (not gating)
    s$f3b_band <- ifelse(s$J >= 12, s$new >= 0.025 & s$new <= 0.065, s$new <= 0.08)
    s$f3b_cov <- s$J < 12 | (s$cov_int >= 0.935 & s$cov_slope >= 0.935)
    pass <- all(s$c1, s$c2, s$c3, s$c4)
    list(table = s, pass = pass)
}

if (sys.nframe() == 0L && !isTRUE(getOption("f3d.source_only"))) {
    library(parallel)
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {
        suppressMessages(library(robustlmm))
        NULL
    })
    clusterExport(cl, c("rep_A", "rep_B", "q", "pv", "cover"))
    for (i in seq_len(nrow(cells))) {
        des <- cells$design[i]; J <- cells$J[i]
        f <- file.path(outdir, sprintf("cell_%s_J%d.rds", des, J))
        if (file.exists(f)) next
        base <- if (des == "A") 3100000L else 3600000L
        repfun <- if (des == "A") rep_A else rep_B
        t0 <- proc.time()
        ## (repfun, not fun: a `fun =` argument would bind to parLapplyLB's
        ## own `fun`; the first launch failed on that before producing data)
        res <- parLapplyLB(cl, base + 10000L * J + seq_len(R_rep),
                           function(s, J, repfun) repfun(s, J), J = J,
                           repfun = repfun)
        res <- do.call(rbind, res); colnames(res) <- colnames_for(des)
        el <- proc.time() - t0
        saveRDS(list(design = des, J = J, res = res, wall = el[["elapsed"]],
                     ncores = ncores,
                     version = as.character(packageVersion("robustlmm"))), f)
        cat(sprintf("cell %s J=%d done in %.0f s\n", des, J, el[["elapsed"]]))
    }
    stopCluster(cl)
    sm <- summarise_f3d(outdir)
    print(sm$table, digits = 3)
    cat("\nDecision rule:", if (sm$pass) "PASSED" else "FAILED", "\n")
    saveRDS(sm, file.path(outdir, "summary.rds"))
}
