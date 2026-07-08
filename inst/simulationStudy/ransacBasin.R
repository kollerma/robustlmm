## WS12 validation study: (1) how well the RANSAC + bisquare estimator
## does, and (2) whether the basin-radius warning works as a diagnostic.
##
## Design (Koller-Stahel 2022 / IF-thread1 Phase-C generators). Sleepstudy
## template (18 subjects x 10 days, random intercept + slope, so the
## phony |rho| -> 1 random-correlation failure is possible). True params =
## classical lme4 fit on sleepstudy. Per replicate and generator we fit:
##   - default : rlmer, monotone smoothPsi (the RSEn baseline);
##   - ransac  : rlmer_ransac + bisquare rho.e (redescending) -- the WS12
##               path, with the basin-radius gate active.
## Generators are RE-dist / error-dist with each in {N, CN, t3, skt3},
## all standardised to Huber-Proposal-2 location 0 / scale 1 (per
## ks_2022_setup.R); CN = 0.9 N(0,1) + 0.1 N(4,1).
##
## A fit is "phony" if |rho-hat| > 0.9999 (the Phase-C threshold). For the
## basin gate we record, per ransac+bisquare fit, whether the warning
## fired (the converged fixed effects left r*(c) of the RANSAC start) and
## cross-tabulate it against phony / large-slope-error, to measure the
## gate's sensitivity (fires when the fit is bad) and specificity (silent
## when it is fine).
##
## ---- FINDINGS (500 reps/gen, K<=100, 2026-06-21, robustlmm WS12) ------
## 1. RANSAC + bisquare is NOT a uniform improvement. Phony rate vs the
##    default monotone RSE: better under RE-contamination / clean errors
##    (N/N 0.116->0.080, CN/N 0.052->0.032) but WORSE under error
##    contamination / heavy tails (N/CN 0.208->0.274, CN/CN 0.116->0.188,
##    skt3/skt3 0.118->0.174); slope RMSE is slightly worse in every cell.
##    => bisquare is not ready to be a default even with the RANSAC start;
##    it earns its keep only for RE-side contamination. (Matches the
##    cautious Phase-C reading and paper5's OBR deprioritisation.)
## 2. The basin-radius warning is a WEAK detector of the phony-rho
##    failure: pooled P(warn | phony) = 0.166 vs P(warn | not phony) =
##    0.113 -- it misses ~83% of phony fits and fires on ~11% of good
##    ones. The reason is structural: r*(c) bounds the FIXED-EFFECT
##    support-preservation radius (beta), while the phony failure is a
##    RANDOM-EFFECT rho -> +-1 boundary attractor (a different basin).
##    The gate correctly guards its STATED property (beta support
##    preservation) but is the wrong instrument for the rho-phony failure.
##    => the warning text overclaims "possibly phony basin"; reword to the
##    beta-support guarantee, and/or add a direct rho-phony check
##    (warn when |rho-hat| -> 1 / the RE covariance is near-singular).
## ----------------------------------------------------------------------

suppressMessages({
    library(robustlmm); library(lme4); library(parallel)
})

NREPS   <- 500L
NCORES  <- 10L
RANSACK <- 100L            # max RANSAC draws (adaptive early-stop on top)
SEED0   <- 20260621L
GENS    <- c("N/N", "N/CN", "CN/N", "CN/CN", "t3/t3", "skt3/skt3")

## ---- truth from the classical fit -------------------------------------
lm0       <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
beta_true <- lme4::fixef(lm0)
sigma_true<- sigma(lm0)
V_b_true  <- matrix(VarCorr(lm0)$Subject[1:4], 2, 2)
L_b       <- t(chol(V_b_true))
subj      <- unique(sleepstudy$Subject)
nS        <- length(subj)

## ---- Huber Proposal-2 location/scale (to standardise each dist) -------
huber_propII <- function(x, k = 1.345, max.iter = 200, rel.tol = 1e-9) {
    mu <- median(x); s <- mad(x)
    kappa <- integrate(function(z) pmin(1, k/abs(z))^2 * z^2 * dnorm(z),
                       -Inf, Inf)$value
    for (it in seq_len(max.iter)) {
        z <- (x - mu) / s
        w <- pmin(1, k / pmax(abs(z), 1e-12))
        mu_new <- sum(w * x) / sum(w)
        z2 <- (x - mu_new) / s
        s_new <- sqrt(mean(pmin(1, k/abs(z2))^2 * (x - mu_new)^2) / kappa)
        if (abs(mu_new - mu) + abs(s_new - s) <
            rel.tol * (abs(mu) + abs(s))) break
        mu <- mu_new; s <- s_new
    }
    c(location = mu, scale = s)
}

## raw (unstandardised) draw per distribution
raw_draw <- function(dist, n) {
    switch(dist,
        "N"    = rnorm(n),
        "CN"   = ifelse(runif(n) < 0.1, rnorm(n, 4, 1), rnorm(n)),
        "t3"   = rt(n, df = 3),
        "skt3" = skewt::rskt(n, df = 3, gamma = 2),
        stop("unknown dist ", dist))
}

## precompute HP2 standardisation constants once per distribution
set.seed(20260524)
STD <- lapply(c("N","CN","t3","skt3"), function(d) {
    hp <- huber_propII(raw_draw(d, 2e5)); list(loc = hp[["location"]],
                                               scale = hp[["scale"]])
})
names(STD) <- c("N","CN","t3","skt3")
sample_std <- function(dist, n) (raw_draw(dist, n) - STD[[dist]]$loc) /
                                 STD[[dist]]$scale

## ---- simulate one replicate of a generator ----------------------------
simulate_rep <- function(gen, seed) {
    parts <- strsplit(gen, "/", fixed = TRUE)[[1]]
    re_dist <- parts[1]; err_dist <- parts[2]
    set.seed(seed)
    d <- sleepstudy
    for (i in seq_len(nS)) {
        u <- sample_std(re_dist, 2L)            # standardised RE coords
        b <- as.numeric(L_b %*% u)
        rows <- d$Subject == subj[i]
        d$Reaction[rows] <- beta_true[1] + beta_true[2] * d$Days[rows] +
            b[1] + b[2] * d$Days[rows] +
            sigma_true * sample_std(err_dist, sum(rows))
    }
    d
}

rho_of <- function(fit) {
    th <- tryCatch(getME(fit, "theta"), error = function(e) NULL)
    if (is.null(th) || length(th) < 3) return(NA_real_)
    L <- matrix(c(th[1], th[2], 0, th[3]), 2)   # lme4 lower-tri theta
    V <- tcrossprod(L)
    V[1, 2] / sqrt(V[1, 1] * V[2, 2])
}
slope_of <- function(fit)
    tryCatch(unname(robustlmm:::.fixef(fit)[2]), error = function(e) NA_real_)

## ---- one replicate: default vs ransac+bisquare ------------------------
one <- function(r, gen) {
    d <- simulate_rep(gen, SEED0 + r)
    fml <- Reaction ~ Days + (Days | Subject)
    ## default monotone RSE (baseline)
    fd <- tryCatch(suppressMessages(suppressWarnings(
        rlmer(fml, d, method = "DASvar"))), error = function(e) NULL)
    ## RANSAC + bisquare (redescending); capture the basin warning
    warned <- FALSE
    fr <- tryCatch(withCallingHandlers(
        suppressMessages(rlmer_ransac(fml, d, K = RANSACK, seed = SEED0 + r,
            method = "DASvar", rho.e = bisquarePsi,
            rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
            rho.b = chgDefaults(smoothPsi, k = 5.14, s = 10),
            rho.sigma.b = chgDefaults(smoothPsi, k = 5.14, s = 10))),
        warning = function(w) {
            if (grepl("basin", conditionMessage(w))) warned <<- TRUE
            invokeRestart("muffleWarning")
        }), error = function(e) NULL)
    rd <- rho_of(fd); rr <- rho_of(fr)
    c(rep = r,
      gen = match(gen, GENS),
      rho_def = rd, rho_ran = rr,
      slope_def = slope_of(fd), slope_ran = slope_of(fr),
      phony_def = as.numeric(is.finite(rd) && abs(rd) > 0.9999),
      phony_ran = as.numeric(is.finite(rr) && abs(rr) > 0.9999),
      warned = as.numeric(warned),
      ok_def = as.numeric(!is.null(fd)), ok_ran = as.numeric(!is.null(fr)))
}

## ---- run --------------------------------------------------------------
t0 <- Sys.time()
grid <- expand.grid(rep = seq_len(NREPS), gen = GENS,
                    stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(grid)), function(i)
    tryCatch(one(grid$rep[i], grid$gen[i]), error = function(e) NULL),
    mc.cores = NCORES)
res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
res <- as.data.frame(res)
res$gen <- GENS[res$gen]
runtime <- as.numeric(Sys.time() - t0, units = "mins")

## ---- aggregate --------------------------------------------------------
agg <- do.call(rbind, lapply(GENS, function(g) {
    s <- res[res$gen == g, ]
    sr <- s[s$ok_ran == 1, ]
    ph <- sr$phony_ran == 1
    data.frame(
        generator   = g,
        n           = nrow(s),
        phony_default = mean(s$phony_def[s$ok_def == 1]),
        phony_ransac  = mean(sr$phony_ran),
        ## basin gate vs the phony outcome on the ransac+bisquare fits
        warn_rate     = mean(sr$warned),
        TPR_warn_if_phony = if (any(ph)) mean(sr$warned[ph]) else NA,
        FPR_warn_if_ok    = if (any(!ph)) mean(sr$warned[!ph]) else NA,
        slope_rmse_def = sqrt(mean((s$slope_def[s$ok_def==1]-beta_true[2])^2)),
        slope_rmse_ran = sqrt(mean((sr$slope_ran - beta_true[2])^2)))
}))

cat(sprintf("\nransacBasin: %d cells, %.1f min, %d cores, K<=%d, %d reps/gen\n",
            nrow(res), runtime, NCORES, RANSACK, NREPS))
print(agg, row.names = FALSE, digits = 3)

## pooled basin-gate confusion (warned x phony) over all generators
sr <- res[res$ok_ran == 1, ]
cat("\nPooled basin-gate confusion (ransac+bisquare):\n")
print(table(warned = sr$warned == 1, phony = sr$phony_ran == 1))
cat(sprintf("\nPooled: P(warn | phony) = %.3f,  P(warn | not phony) = %.3f\n",
            mean(sr$warned[sr$phony_ran==1]), mean(sr$warned[sr$phony_ran==0])))

saveRDS(list(results = res, agg = agg,
             config = list(NREPS = NREPS, K = RANSACK, gens = GENS,
                           beta_true = beta_true, sigma_true = sigma_true,
                           V_b_true = V_b_true, runtime_min = runtime,
                           R = R.version.string,
                           robustlmm = as.character(packageVersion("robustlmm")))),
        file.path("inst", "simulationStudy", "ransacBasin_results.rds"))
cat("\nsaved inst/simulationStudy/ransacBasin_results.rds\n")
