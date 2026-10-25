## Shipped with robustlmm: a script of the breakdown study behind
## setting = "RSEhb" in vignette("simulationStudies"). Run it with Rscript
## from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## BD5: high-breakdown starts for (beta, sigma, theta). The decision rule
## was fixed before the run.
##
## Same DGP and contamination cells as bd23_map_sim.R. Starts:
##   lmer      the internal lmer start (rlmer's default)
##   ransac    ransac_lme4() defaults: fixef, sigma, theta of the best candidate
##   mm_scr    lmrob() MM on the marginal model, screen |r| <= 2.576 scale,
##             lmer on the kept rows (the current init = "robust")
##   s_scr     the same screen with the S-estimate (lmrob(method = "S"))
##   s_dec     S-estimate coefficients, sigma_e0 / sigma_b0 from a robust
##             within/between decomposition of its residuals
## For each start: the start's intercept, slope, sigma, sd_b, and the final
## fit of the full S-type estimator (rho.e = rho.b = bisquarePsi, both scales
## makeSScalePsi()) from it.

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 100L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "bd5_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
start_names <- c("lmer", "ransac", "mm_scr", "s_scr", "s_dec")

make_data <- function(seed, type, eps, J = 20L, m = 10L) {
    set.seed(seed)
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(J * m); b <- rnorm(J); e <- rnorm(J * m)
    if (type == "C8" && eps > 0) {
        bad <- sample.int(J, round(eps * J)); b[bad] <- b[bad] + 8
    }
    y <- 2 + x + b[g] + e
    if (type == "V8" && eps > 0) {
        bad <- sample.int(J * m, round(eps * J * m)); y[bad] <- y[bad] + 8
    } else if (type %in% c("L4", "L8") && eps > 0) {
        bad <- sample.int(J * m, round(eps * J * m))
        x[bad] <- rnorm(length(bad), if (type == "L4") 4 else 8, 0.3)
        y[bad] <- rnorm(length(bad))
    }
    data.frame(g = g, x = x, y = y)
}

starts <- function(d, seed) {
    q <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)),
                                 error = function(e) NULL)
    lmer_list <- function(f) list(fixef = unname(lme4::fixef(f)), sigma = sigma(f),
                                  theta = max(unname(lme4::getME(f, "theta")), 0.05))
    screen <- function(res, scale) {
        keep <- abs(res) <= stats::qnorm(0.995) * scale
        q(lmer_list(lme4::lmer(y ~ x + (1 | g), d[keep, , drop = FALSE])))
    }
    out <- list()
    out$lmer <- q(lmer_list(lme4::lmer(y ~ x + (1 | g), d)))
    rr <- q(ransac_lme4(y ~ x + (1 | g), d, K = 200L, seed = seed))
    out$ransac <- if (is.null(rr) || is.null(rr$fit)) NULL else {
        l <- robustlmm:::.ransacInitList(rr$fit); l$theta <- max(l$theta, 0.05); l }
    mm <- q(robustbase::lmrob(y ~ x, d))
    out$mm_scr <- if (is.null(mm)) NULL else screen(stats::residuals(mm), mm$scale)
    s <- q(robustbase::lmrob(y ~ x, d, method = "S"))
    out$s_scr <- if (is.null(s)) NULL else screen(stats::residuals(s), s$scale)
    out$s_dec <- if (is.null(s)) NULL else {
        r <- stats::residuals(s)
        med <- tapply(r, d$g, stats::median)
        m <- mean(table(d$g))
        se <- stats::mad(r - med[as.character(d$g)]) * sqrt(m / (m - 1))
        sb2 <- stats::mad(med)^2 - (pi / 2) * se^2 / m
        list(fixef = unname(stats::coef(s)), sigma = se,
             theta = sqrt(max(sb2, 0.0025 * se^2)) / se)
    }
    out
}

one_rep <- function(seed, type, eps) {
    d <- make_data(seed, type, eps)
    st <- starts(d, seed)
    ss <- makeSScalePsi()
    vapply(start_names, function(s) {
        init <- st[[s]]
        if (is.null(init)) return(rep(NA_real_, 8))
        f <- tryCatch(suppressWarnings(suppressMessages(
            rlmer(y ~ x + (1 | g), d, method = "DASvar", rho.e = bisquarePsi,
                  rho.b = bisquarePsi, rho.sigma.e = ss, rho.sigma.b = ss,
                  init = init))), error = function(e) NULL)
        fin <- if (is.null(f)) rep(NA_real_, 4) else
            c(unname(fixef(f)), sigma(f), sigma(f) * getME(f, "theta"))
        c(init$fixef, init$sigma, init$sigma * init$theta, fin)
    }, numeric(8))
}

cells <- rbind(data.frame(type = "V8", eps = c(0, 0.1, 0.2, 0.25, 0.3, 0.4)),
               data.frame(type = "C8", eps = c(0.1, 0.2, 0.3, 0.4)),
               data.frame(type = "L4", eps = c(0.05, 0.1, 0.2)),
               data.frame(type = "L8", eps = c(0.05, 0.1, 0.2, 0.3)))
if (sys.nframe() == 0L && !isTRUE(getOption("bd5.source_only"))) {
    library(parallel)
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {
        suppressMessages(library(robustlmm))
        NULL
    })
    clusterExport(cl, c("one_rep", "make_data", "starts", "start_names"))
    for (i in seq_len(nrow(cells))) {
        f <- file.path(outdir, sprintf("cell_%s_%03d.rds", cells$type[i],
                                       round(100 * cells$eps[i])))
        if (file.exists(f)) next
        t0 <- proc.time()
        ## fresh seeds (independent of bd23_map_sim.R)
        res <- parLapplyLB(cl, 3000000L + 10000L * i + seq_len(R_rep),
                           function(s, type, eps) one_rep(s, type, eps),
                           type = cells$type[i], eps = cells$eps[i])
        el <- proc.time() - t0
        saveRDS(list(type = cells$type[i], eps = cells$eps[i], res = res,
                     wall = el[["elapsed"]], ncores = ncores), f)
        cat(sprintf("cell %s eps=%.2f done in %.0f s\n", cells$type[i], cells$eps[i],
                    el[["elapsed"]]))
    }
    stopCluster(cl)
}
