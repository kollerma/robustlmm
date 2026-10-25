## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## F1: design-adaptive breakdown constant b for the RANSAC selection scale.
## The decision rule was fixed before the run.
##
## Design: balanced binary between-cluster factor, J = 20 clusters (10 per
## arm), m = 10: y = 1 + 2 arm + b_j + e, b ~ N(0, 1), e ~ N(0, 1);
## model y ~ arm + (1 | g). All rows of arm 0 lie on one hyperplane through
## the origin, k/n = 1/2: design_breakdown() gives b_opt = 1/4, a ceiling of
## about 1/4, and an empty Yohai guarantee at b = 1/2.
## Cells (contamination confined to arm 1):
##   clean
##   clus3   3 of the 10 arm-1 clusters shifted by +8 (15% of all rows)
##   clus4   4 of the 10 arm-1 clusters shifted by +8 (20%)
##   obs30   30% of the arm-1 observations shifted by +8 (15%)
## Criteria: the bisquare M-scale about zero with b = 1/2 (c = 1.547645,
## current) and with b = b_opt = 1/4 (c solving E_Phi rho_c = 1/4).
## Per criterion: the RANSAC start's arm coefficient and intercept, and the
## final rlmer fit (bisquare rho.e, as T4) from that start; plus the final
## fit from rlmer's default start.
## Exploratory (not part of the rule): "marg" scores each candidate by the
## b = 1/2 M-scale of its MARGINAL residuals y - X beta_cand. ransac_lme4()
## scores conditional residuals (in-sample clusters keep their predicted
## random effects), under which a contaminated cluster inside the subsample
## is absorbed by its own random effect.

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 200L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "f1_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

c_for_b <- function(b) stats::uniroot(function(cc)
    stats::integrate(function(z) robustbase::Mchi(z, cc, "bisquare") *
                         stats::dnorm(z), -Inf, Inf)$value - b,
    c(0.5, 10), tol = 1e-10)$root
C_HALF <- 1.547645
C_QUARTER <- c_for_b(0.25)

## Same subsample draws as ransac_lme4(seed = ) (fixed K, no early stop),
## scored on marginal residuals.
ransac_marginal <- function(formula, d, K, seed, cc = C_HALF, b = 0.5) {
    set.seed(seed)
    info <- robustlmm:::.ransac_design_info(formula, d)
    strata <- robustlmm:::.ransac_strata(formula, d)
    X <- stats::model.matrix(reformulas::nobars(formula), d)
    best <- Inf; best_fit <- NULL
    for (k in seq_len(K)) {
        draw <- robustlmm:::.ransac_nonsingular_subsample(formula, d, 0.5,
                                                          strata, info)
        if (!draw$ok) next
        f <- tryCatch(suppressMessages(suppressWarnings(
            lme4::lmer(formula, d[draw$idx, , drop = FALSE], REML = TRUE))),
            error = function(e) NULL)
        if (is.null(f) || length(lme4::fixef(f)) != ncol(X)) next
        r <- d$y - as.numeric(X %*% lme4::fixef(f))
        s <- robustlmm:::.ransac_mscale0(r, cc = cc, b = b)
        if (is.finite(s) && s < best) { best <- s; best_fit <- f }
    }
    best_fit
}

one_rep <- function(seed, cell, J = 20L, m = 10L) {
    set.seed(seed)
    g <- factor(rep(seq_len(J), each = m))
    arm_cl <- rep(0:1, each = J / 2)
    arm <- arm_cl[g]
    y <- 1 + 2 * arm + rnorm(J)[g] + rnorm(J * m)
    if (cell %in% c("clus3", "clus4")) {
        nb <- if (cell == "clus3") 3L else 4L
        bad <- sample(which(arm_cl == 1L), nb)
        y[g %in% bad] <- y[g %in% bad] + 8
    } else if (cell == "obs30") {
        a1 <- which(arm == 1L)
        bad <- sample(a1, round(0.3 * length(a1)))
        y[bad] <- y[bad] + 8
    }
    d <- data.frame(g = g, arm = arm, y = y)
    fin <- function(init) tryCatch(suppressWarnings(suppressMessages(
        rlmer(y ~ arm + (1 | g), d, method = "DASvar", rho.e = bisquarePsi,
              rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
              rho.b = chgDefaults(smoothPsi, k = 5.14),
              rho.sigma.b = chgDefaults(smoothPsi, k = 5.14),
              init = init))), error = function(e) NULL)
    f0 <- fin(NULL)
    out <- if (is.null(f0)) c(NA, NA) else unname(fixef(f0))
    for (cc in c(C_HALF, C_QUARTER)) {
        b <- if (cc == C_HALF) 0.5 else 0.25
        sf <- function(r) robustlmm:::.ransac_mscale0(r, cc = cc, b = b)
        st <- tryCatch(suppressWarnings(ransac_lme4(y ~ arm + (1 | g), d,
                K = 200L, seed = seed, scale_fn = sf)), error = function(e) NULL)
        if (is.null(st) || is.null(st$fit)) { out <- c(out, rep(NA, 4)); next }
        init <- robustlmm:::.ransacInitList(st$fit)
        f <- fin(init)
        out <- c(out, init$fixef, if (is.null(f)) c(NA, NA) else unname(fixef(f)))
    }
    fm <- tryCatch(ransac_marginal(y ~ arm + (1 | g), d, K = 200L, seed = seed),
                   error = function(e) NULL)
    if (is.null(fm)) return(c(out, rep(NA, 4)))
    init <- robustlmm:::.ransacInitList(fm)
    f <- fin(init)
    c(out, init$fixef, if (is.null(f)) c(NA, NA) else unname(fixef(f)))
}

cells <- c("clean", "clus3", "clus4", "obs30")
cn <- c("def_int_fin", "def_arm_fin",
        paste0("b50_", c("int_start", "arm_start", "int_fin", "arm_fin")),
        paste0("b25_", c("int_start", "arm_start", "int_fin", "arm_fin")),
        paste0("marg_", c("int_start", "arm_start", "int_fin", "arm_fin")))
if (sys.nframe() == 0L && !isTRUE(getOption("f1.source_only"))) {
    library(parallel)
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {
        suppressMessages(library(robustlmm))
        NULL
    })
    clusterExport(cl, c("one_rep", "ransac_marginal", "C_HALF", "C_QUARTER"))
    for (i in seq_along(cells)) {
        f <- file.path(outdir, sprintf("cell_%s.rds", cells[i]))
        if (file.exists(f)) next
        t0 <- proc.time()
        res <- parLapplyLB(cl, 1400000L + 10000L * i + seq_len(R_rep),
                           function(s, cell) one_rep(s, cell), cell = cells[i])
        res <- do.call(rbind, res); colnames(res) <- cn
        el <- proc.time() - t0
        saveRDS(list(cell = cells[i], res = res, wall = el[["elapsed"]],
                     ncores = ncores, c_quarter = C_QUARTER), f)
        cat(sprintf("cell %s done in %.0f s\n", cells[i], el[["elapsed"]]))
    }
    stopCluster(cl)
}
