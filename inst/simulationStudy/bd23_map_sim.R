## Shipped with robustlmm: a script of the breakdown study behind
## setting = "RSEhb" in vignette("simulationStudies"). Run it with Rscript
## from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## BD2 + BD3: breakdown map of the current rlmer settings and of the
## high-breakdown scale candidates. The decision rule was fixed before the run.
##
## DGP: y = 2 + x + b_j + e, J = 20, m = 10, b ~ N(0, 1), e ~ N(0, 1),
## x ~ N(0, 1); model y ~ x + (1 | g), method DASvar. Contamination:
##   V8   fraction eps of observations: e + 8            (vertical)
##   C8   fraction eps of clusters: b_j + 8               (outlying clusters)
##   L4   fraction eps of observations: x ~ N(4, 0.3^2), y ~ N(0, 1)
##   L8   fraction eps of observations: x ~ N(8, 0.3^2), y ~ N(0, 1)
## Estimators (rho.b = smoothPsi and the RSEn b-scale unless stated;
## "S-scale" = the bounded-rho scale of sscale_proto.R, bisquare b = 1/2):
##   lmer          lme4::lmer, REML
##   RSEn          rlmer defaults (smoothPsi, psi2propII(smoothPsi))
##   RSEa          rlmer setting = "RSEa" (scale k = 2.28)
##   bisq_RSEn     rho.e = bisquarePsi, scales psi2propII(smoothPsi, k = 1.345)
##   bisq_RSEa     rho.e = bisquarePsi, scales psi2propII(smoothPsi, k = 2.28)
##   bisq_S        rho.e = bisquarePsi, rho.sigma.e = S-scale
##   full_S        rho.e = rho.b = bisquarePsi, both scales S-scale
##   bisq_S_rs     bisq_S from the robust start
##   full_S_rs     full_S from the robust start
##   full_S_rs85   full_S_rs with rho.e = bisquare c = 3.444 (85% efficiency)
##   MM_oracle     rho.e = bisquarePsi, sigma and theta fixed at the truth
##                 (init = truth, max.iter = 0): benchmark, not an estimator
##   MM_ransac     rho.e = bisquarePsi, start = ransac_lme4() fixed effects,
##                 sigma0 = M-scale about zero of its full-data residuals,
##                 theta0 from the candidate; sigma, theta held (max.iter = 0)
## Robust start (random-intercept models): beta0 = lmrob(y ~ x) (default
## bisquare MM from a 50%-breakdown fast-S estimate; not setting = "KS2014",
## whose initial estimate broke down on a far-leverage smoke case) on the
## marginal model; with its residuals r, sigma_e0 = MAD of the within-cluster
## deviations r - median_j(r) (times sqrt(m/(m-1))),
## sigma_b0^2 = max(0, MAD(median_j r)^2 - (pi/2) sigma_e0^2/m),
## theta0 = sigma_b0 / sigma_e0.
##   full_S_rsl    full_S from the screened-lmer start: the same lmrob fit,
##                 then lmer on the observations with |r| <= 2.576 scale
##                 (any random-effects structure)

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 100L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "bd23_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

est_names <- c("lmer", "RSEn", "RSEa", "bisq_RSEn", "bisq_RSEa", "bisq_S",
               "full_S", "bisq_S_rs", "full_S_rs", "full_S_rs85", "full_S_rsl",
               "MM_oracle", "MM_ransac")

robust_start <- function(d) {
    lr <- robustbase::lmrob(y ~ x, d)
    b0 <- if (lr$converged) coef(lr) else lr$init.S$coefficients
    r <- d$y - as.numeric(cbind(1, d$x) %*% b0)
    med <- tapply(r, d$g, stats::median)
    m <- mean(table(d$g))
    se <- stats::mad(r - med[as.character(d$g)]) * sqrt(m / (m - 1))
    sb2 <- stats::mad(med)^2 - (pi / 2) * se^2 / m
    list(fixef = unname(b0), sigma = se, theta = sqrt(max(sb2, 0.01 * se^2)) / se)
}
screened_lmer_start <- function(d) {
    lr <- robustbase::lmrob(y ~ x, d)
    keep <- abs(stats::residuals(lr)) <= stats::qnorm(0.995) * lr$scale
    f <- suppressMessages(suppressWarnings(
        lme4::lmer(y ~ x + (1 | g), d[keep, , drop = FALSE])))
    list(fixef = unname(lme4::fixef(f)), sigma = sigma(f),
         theta = max(unname(lme4::getME(f, "theta")), 0.1))
}

one_rep <- function(seed, type, eps, J = 20L, m = 10L) {
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
    d <- data.frame(g = g, x = x, y = y)
    q <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)),
                                 error = function(e) NULL)
    summ <- function(f, isr = TRUE) {
        if (is.null(f)) return(rep(NA_real_, 5))
        sig <- sigma(f)
        th <- if (isr) getME(f, "theta") else lme4::getME(f, "theta")
        conv <- if (isr) as.numeric(length(f@optinfo$warnings) == 0) else 1
        c(unname(fixef(f)), sig, sig * th, conv)
    }
    ss <- sScalePsi_proto()
    rl <- function(...) q(rlmer(y ~ x + (1 | g), d, method = "DASvar", ...))
    sm <- psi2propII(smoothPsi, k = 1.345); sa <- psi2propII(smoothPsi, k = 2.28)
    b85 <- makeBisquarePsi(3.443689)
    rs <- q(robust_start(d))
    from_rs <- function(...) if (is.null(rs)) rep(NA_real_, 5) else summ(rl(init = rs, ...))
    out <- list(
        lmer = summ(q(lme4::lmer(y ~ x + (1 | g), d)), isr = FALSE),
        RSEn = summ(rl()),
        RSEa = summ(rl(setting = "RSEa")),
        bisq_RSEn = summ(rl(rho.e = bisquarePsi, rho.sigma.e = sm, rho.b = smoothPsi,
                            rho.sigma.b = sm)),
        bisq_RSEa = summ(rl(rho.e = bisquarePsi, rho.sigma.e = sa, rho.b = smoothPsi,
                            rho.sigma.b = sa)),
        bisq_S = summ(rl(rho.e = bisquarePsi, rho.sigma.e = ss, rho.b = smoothPsi,
                         rho.sigma.b = sm)),
        full_S = summ(rl(rho.e = bisquarePsi, rho.sigma.e = ss, rho.b = bisquarePsi,
                         rho.sigma.b = ss)),
        bisq_S_rs = from_rs(rho.e = bisquarePsi, rho.sigma.e = ss, rho.b = smoothPsi,
                            rho.sigma.b = sm),
        full_S_rs = from_rs(rho.e = bisquarePsi, rho.sigma.e = ss, rho.b = bisquarePsi,
                            rho.sigma.b = ss),
        full_S_rs85 = from_rs(rho.e = b85, rho.sigma.e = ss, rho.b = bisquarePsi,
                              rho.sigma.b = ss),
        full_S_rsl = {
            rsl <- q(screened_lmer_start(d))
            if (is.null(rsl)) rep(NA_real_, 5) else
                summ(rl(init = rsl, rho.e = bisquarePsi, rho.sigma.e = ss,
                        rho.b = bisquarePsi, rho.sigma.b = ss))
        },
        MM_oracle = summ(rl(rho.e = bisquarePsi, rho.sigma.e = sm, rho.b = smoothPsi,
                            rho.sigma.b = sm,
                            init = list(fixef = c(2, 1), sigma = 1, theta = 1),
                            max.iter = 0)))
    st <- q(ransac_lme4(y ~ x + (1 | g), d, K = 200L, seed = seed))
    out$MM_ransac <- if (is.null(st) || is.null(st$fit)) rep(NA_real_, 5) else {
        init <- robustlmm:::.ransacInitList(st$fit)
        r <- d$y - stats::predict(st$fit, newdata = d, allow.new.levels = TRUE)
        init$sigma <- robustlmm:::.ransac_mscale0(r)
        summ(rl(rho.e = bisquarePsi, rho.sigma.e = sm, rho.b = smoothPsi,
                rho.sigma.b = sm, init = init, max.iter = 0))
    }
    unlist(out[est_names])
}

cells <- rbind(data.frame(type = "V8", eps = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4)),
               data.frame(type = "C8", eps = c(0.1, 0.2, 0.3, 0.4)),
               data.frame(type = "L4", eps = c(0.05, 0.1, 0.2)),
               data.frame(type = "L8", eps = c(0.05, 0.1, 0.2, 0.3)))
cn <- as.vector(outer(c("int", "slope", "sigma", "sd_b", "conv"), est_names,
                      function(a, b) paste(b, a, sep = "_")))
if (sys.nframe() == 0L && !isTRUE(getOption("bd23.source_only"))) {
    library(parallel)
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {
        suppressMessages(library(robustlmm))
        source(system.file("simulationStudy", "sscale_proto.R", package = "robustlmm"))
        NULL
    })
    clusterExport(cl, c("one_rep", "est_names", "robust_start", "screened_lmer_start"))
    for (i in seq_len(nrow(cells))) {
        f <- file.path(outdir, sprintf("cell_%s_%03d.rds", cells$type[i],
                                       round(100 * cells$eps[i])))
        if (file.exists(f)) next
        t0 <- proc.time()
        res <- parLapplyLB(cl, 2000000L + 10000L * i + seq_len(R_rep),
                           function(s, type, eps) one_rep(s, type, eps),
                           type = cells$type[i], eps = cells$eps[i])
        res <- do.call(rbind, res); colnames(res) <- cn
        el <- proc.time() - t0
        saveRDS(list(type = cells$type[i], eps = cells$eps[i], res = res,
                     wall = el[["elapsed"]], ncores = ncores), f)
        cat(sprintf("cell %s eps=%.2f done in %.0f s\n", cells$type[i], cells$eps[i],
                    el[["elapsed"]]))
    }
    stopCluster(cl)
}
