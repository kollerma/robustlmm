## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## Re-run of F2 after a fix to how rlmer applies list starts (they were
## applied with the classical psi before). Same script, seeds and designs
## as the first run.
## F2: RANSAC start-selection criteria, judged on the final rlmer fit.
## The decision rule was fixed before the run.
##
## Selection criteria (scale of a candidate's full-data residuals r):
##   A  bisquare M-scale about zero (current default, T4)
##   B  Qn(r) (robustlmm 3.5.0)
##   C  B, then the start's intercept re-centred by median(r)
##   D  profiled M-scale min_mu s(r - mu); the start's intercept is moved by
##      the minimising mu (the S-location of the residuals)
##   E  tau-scale about zero, robustbase::scaleTau2(r, mu0 = 0)
## plus A and D with the start's sigma set to the criterion's robust scale
## instead of the subsample lmer sigma ("A_rs", "D_rs").
##
## Part 1 ("basin", primary): the six Koller-Stahel generators of
## inst/simulationStudy/ransacConsensus.R (sleepstudy design, bisquare rho.e,
## K = 80, single start), where a start can land in the phony
## |rho-hat| -> 1 attractor. Per variant: phony indicator, final slope and
## intercept, start slope and intercept.
## Part 2 ("invariance"): T4's DGP with 20% vertical outliers (+8) and a
## bad-leverage cell (10% of observations at x ~ N(4, 0.3^2), y ~ N(0, 1)),
## 50 replicates: does the final fit depend on the start at all? Each
## variant's final fit is compared with the fit from rlmer's default start.

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 200L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "f2_rerun_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

mscale0 <- function(r) robustlmm:::.ransac_mscale0(r)
s_loc <- function(r) {
    ## S-location of r for the bisquare M-scale: argmin_mu s(r - mu)
    m <- stats::median(r); s <- stats::mad(r)
    if (!is.finite(s) || s <= 0) return(list(mu = m, s = 0))
    o <- stats::optimize(function(mu) mscale0(r - mu), m + c(-2, 2) * s)
    list(mu = o$minimum, s = o$objective)
}
crit_fun <- list(
    A = mscale0,
    B = robustbase::Qn,
    D = function(r) s_loc(r)$s,
    E = function(r) robustbase::scaleTau2(r, mu0 = 0))
variants <- c("A", "B", "C", "D", "E", "A_rs", "D_rs")

## starts for every variant from one RANSAC run per criterion
starts_for <- function(formula, d, yname, K, seed) {
    out <- list()
    for (cr in c("A", "B", "D", "E")) {
        st <- tryCatch(suppressWarnings(
            ransac_lme4(formula, d, K = K, seed = seed,
                        scale_fn = crit_fun[[cr]])), error = function(e) NULL)
        if (is.null(st) || is.null(st$fit)) next
        init <- robustlmm:::.ransacInitList(st$fit)
        r <- d[[yname]] - stats::predict(st$fit, newdata = d,
                                         allow.new.levels = TRUE)
        if (cr == "A") {
            out$A <- init
            out$A_rs <- init; out$A_rs$sigma <- mscale0(r)
        } else if (cr == "B") {
            out$B <- init
            out$C <- init; out$C$fixef[1] <- init$fixef[1] + stats::median(r)
        } else if (cr == "D") {
            sl <- s_loc(r)
            out$D <- init; out$D$fixef[1] <- init$fixef[1] + sl$mu
            out$D_rs <- out$D; out$D_rs$sigma <- sl$s
        } else out$E <- init
    }
    out
}

## ---- Part 1: phony-basin generators (ransacConsensus.R) ----
setup_basin <- function() {
    lm0 <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
    L_b <- t(chol(matrix(lme4::VarCorr(lm0)$Subject[1:4], 2, 2)))
    huber_propII <- function(x, k = 1.345, max.iter = 200, rel.tol = 1e-9) {
        mu <- median(x); s <- mad(x)
        kappa <- integrate(function(z) pmin(1, k/abs(z))^2 * z^2 * dnorm(z),
                           -Inf, Inf)$value
        for (it in seq_len(max.iter)) {
            w <- pmin(1, k / pmax(abs((x - mu) / s), 1e-12))
            mu_new <- sum(w * x) / sum(w)
            s_new <- sqrt(mean(pmin(1, k / abs((x - mu_new) / s))^2 *
                               (x - mu_new)^2) / kappa)
            if (abs(mu_new - mu) + abs(s_new - s) <
                rel.tol * (abs(mu) + abs(s))) break
            mu <- mu_new; s <- s_new
        }
        c(location = mu, scale = s)
    }
    set.seed(20260524)
    STD <- lapply(c("N", "CN", "t3", "skt3"), function(dd) {
        hp <- huber_propII(raw_draw(dd, 2e5))
        list(loc = hp[["location"]], scale = hp[["scale"]])
    })
    names(STD) <- c("N", "CN", "t3", "skt3")
    list(beta = lme4::fixef(lm0), sigma = sigma(lm0), L_b = L_b, STD = STD)
}
raw_draw <- function(dist, n) switch(dist,
    "N" = rnorm(n), "CN" = ifelse(runif(n) < 0.1, rnorm(n, 4, 1), rnorm(n)),
    "t3" = rt(n, 3), "skt3" = skewt::rskt(n, df = 3, gamma = 2))
GENS <- c("N/N", "N/CN", "CN/N", "CN/CN", "t3/t3", "skt3/skt3")

one_basin <- function(r, gen, B) {
    samp <- function(dd, n) (raw_draw(dd, n) - B$STD[[dd]]$loc) / B$STD[[dd]]$scale
    p <- strsplit(gen, "/", fixed = TRUE)[[1]]
    seed <- 20260927L + r
    set.seed(seed)
    d <- lme4::sleepstudy
    subj <- unique(d$Subject)
    for (i in seq_along(subj)) {
        b <- as.numeric(B$L_b %*% samp(p[1], 2L)); rows <- d$Subject == subj[i]
        d$Reaction[rows] <- B$beta[1] + B$beta[2] * d$Days[rows] +
            b[1] + b[2] * d$Days[rows] + B$sigma * samp(p[2], sum(rows))
    }
    fml <- Reaction ~ Days + (Days | Subject)
    st <- starts_for(fml, d, "Reaction", K = 80L, seed = seed)
    vapply(variants, function(v) {
        init <- st[[v]]
        if (is.null(init)) return(rep(NA_real_, 5))
        f <- tryCatch(suppressWarnings(suppressMessages(rlmer(fml, d,
            method = "DASvar", rho.e = bisquarePsi,
            rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
            rho.b = chgDefaults(smoothPsi, k = 5.14, s = 10),
            rho.sigma.b = chgDefaults(smoothPsi, k = 5.14, s = 10),
            init = init))), error = function(e) NULL)
        if (is.null(f)) return(c(NA, NA, NA, init$fixef))
        rho <- tryCatch(robustlmm:::.re_max_abscor(f), error = function(e) NA_real_)
        c(as.numeric(rho > 0.9999), unname(fixef(f)), init$fixef)
    }, numeric(5))
}

## ---- Part 2: start-invariance in T4's DGP ----
one_invar <- function(seed, cell, J = 20L, m = 10L) {
    set.seed(seed)
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(J * m); e <- rnorm(J * m)
    y <- 2 + x + rnorm(J)[g] + e
    if (cell == "vertical") {
        out <- sample.int(J * m, round(0.2 * J * m)); y[out] <- y[out] + 8
    } else {
        out <- sample.int(J * m, round(0.1 * J * m))
        x[out] <- rnorm(length(out), 4, 0.3); y[out] <- rnorm(length(out))
    }
    d <- data.frame(g = g, x = x, y = y)
    fin <- function(init) tryCatch(suppressWarnings(suppressMessages(
        rlmer(y ~ x + (1 | g), d, method = "DASvar", rho.e = bisquarePsi,
              rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
              rho.b = chgDefaults(smoothPsi, k = 5.14),
              rho.sigma.b = chgDefaults(smoothPsi, k = 5.14),
              init = init))), error = function(e) NULL)
    f0 <- fin(NULL)
    b0 <- if (is.null(f0)) c(NA, NA) else unname(fixef(f0))
    st <- starts_for(y ~ x + (1 | g), d, "y", K = 200L, seed = seed)
    c(b0, vapply(variants, function(v) {
        init <- st[[v]]
        if (is.null(init)) return(rep(NA_real_, 4))
        f <- fin(init)
        c(if (is.null(f)) c(NA, NA) else unname(fixef(f)), init$fixef)
    }, numeric(4)))
}

if (sys.nframe() == 0L && !isTRUE(getOption("f2.source_only"))) {
    library(parallel)
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {
        suppressMessages(library(robustlmm))
        NULL
    })
    B <- setup_basin()
    clusterExport(cl, c("one_basin", "one_invar", "starts_for", "mscale0",
                        "s_loc", "crit_fun", "variants", "raw_draw", "B"))
    for (gi in seq_along(GENS)) {
        f <- file.path(outdir, sprintf("basin_%d.rds", gi))
        if (file.exists(f)) next
        t0 <- proc.time()
        res <- parLapplyLB(cl, seq_len(R_rep), function(r, gen) one_basin(r, gen, B),
                           gen = GENS[gi])
        el <- proc.time() - t0
        saveRDS(list(gen = GENS[gi], res = res, wall = el[["elapsed"]],
                     ncores = ncores), f)
        cat(sprintf("basin %s done in %.0f s\n", GENS[gi], el[["elapsed"]]))
    }
    for (cell in c("vertical", "leverage")) {
        f <- file.path(outdir, sprintf("invar_%s.rds", cell))
        if (file.exists(f)) next
        t0 <- proc.time()
        res <- parLapplyLB(cl, 1300000L + seq_len(50L),
                           function(s, cell) one_invar(s, cell), cell = cell)
        el <- proc.time() - t0
        saveRDS(list(cell = cell, res = do.call(rbind, res),
                     wall = el[["elapsed"]], ncores = ncores), f)
        cat(sprintf("invariance %s done in %.0f s\n", cell, el[["elapsed"]]))
    }
    stopCluster(cl)
}
