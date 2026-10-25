## WS12-D payoff: does the multi-start RANSAC consensus lower the phony-
## correlation rate of the bisquare redescender vs the single best start?
## Same Koller-Stahel / Phase-C generators as ransacBasin.R; per replicate
## we fit ransac+bisquare with n_starts = 1 (single best start) and
## n_starts = 5 (consensus) and compare the phony rate (|rho-hat| > 0.9999)
## and slope RMSE.
##
## ---- FINDING (200 reps/gen, n_starts = 5, K <= 80) ----------------------
## Re-run 2026-09-29 with version 3.5.1 (RANSAC candidates
## ranked by the centred M-scale; list starts' fixed effects honoured):
## the consensus lowers the phony rate in every generator,
##   N/N 0.145->0.135, N/CN 0.340->0.320, CN/N 0.065->0.060,
##   CN/CN 0.170->0.165, t3/t3 0.210->0.185, skt3/skt3 0.170->0.155;
##   pooled 0.183 -> 0.170, all 16 discordant pairs in its favour, slope
## RMSE unchanged (1.91 vs 1.92). The effect is small; earlier runs gave
## 0.178 -> 0.165 (3.5.0-2) and 0.149 -> 0.027 (June 2026, before the fix
## of init = "ransac" fitting the RANSAC subsample). The vignettes quote
## the shipped results via vignettes/precompute-3.5.0.R.
## ----------------------------------------------------------------------

suppressMessages({ library(robustlmm); library(lme4); library(parallel) })

NREPS <- 200L; NCORES <- 5L; RANSACK <- 80L; NSTART <- 5L
SEED0 <- 20260622L
GENS  <- c("N/N", "N/CN", "CN/N", "CN/CN", "t3/t3", "skt3/skt3")

lm0 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
beta_true <- lme4::fixef(lm0); sigma_true <- sigma(lm0)
L_b <- t(chol(matrix(VarCorr(lm0)$Subject[1:4], 2, 2)))
subj <- unique(sleepstudy$Subject); nS <- length(subj)

huber_propII <- function(x, k = 1.345, max.iter = 200, rel.tol = 1e-9) {
    mu <- median(x); s <- mad(x)
    kappa <- integrate(function(z) pmin(1,k/abs(z))^2*z^2*dnorm(z),-Inf,Inf)$value
    for (it in seq_len(max.iter)) {
        mu_new <- { w <- pmin(1,k/pmax(abs((x-mu)/s),1e-12)); sum(w*x)/sum(w) }
        s_new <- sqrt(mean(pmin(1,k/abs((x-mu_new)/s))^2*(x-mu_new)^2)/kappa)
        if (abs(mu_new-mu)+abs(s_new-s) < rel.tol*(abs(mu)+abs(s))) break
        mu <- mu_new; s <- s_new
    }
    c(location = mu, scale = s)
}
raw_draw <- function(dist, n) switch(dist,
    "N" = rnorm(n), "CN" = ifelse(runif(n)<0.1, rnorm(n,4,1), rnorm(n)),
    "t3" = rt(n,3), "skt3" = skewt::rskt(n,df=3,gamma=2))
set.seed(20260524)
STD <- lapply(c("N","CN","t3","skt3"), function(d) {
    hp <- huber_propII(raw_draw(d,2e5)); list(loc=hp[["location"]],scale=hp[["scale"]]) })
names(STD) <- c("N","CN","t3","skt3")
samp <- function(d,n) (raw_draw(d,n)-STD[[d]]$loc)/STD[[d]]$scale

simulate_rep <- function(gen, seed) {
    p <- strsplit(gen,"/",fixed=TRUE)[[1]]; set.seed(seed); d <- sleepstudy
    for (i in seq_len(nS)) {
        b <- as.numeric(L_b %*% samp(p[1],2L)); rows <- d$Subject==subj[i]
        d$Reaction[rows] <- beta_true[1]+beta_true[2]*d$Days[rows] +
            b[1]+b[2]*d$Days[rows] + sigma_true*samp(p[2],sum(rows))
    }
    d
}
phony_of <- function(fit) {
    r <- tryCatch(robustlmm:::.re_max_abscor(fit), error=function(e) NA_real_)
    if (!is.finite(r)) return(NA); as.numeric(r > 0.9999)
}
slope_of <- function(fit) tryCatch(unname(robustlmm:::.fixef(fit)[2]),
                                   error=function(e) NA_real_)
rho_args <- list(method="DASvar", rho.e=bisquarePsi,
    rho.sigma.e=psi2propII(smoothPsi,k=2.28),
    rho.b=chgDefaults(smoothPsi,k=5.14,s=10),
    rho.sigma.b=chgDefaults(smoothPsi,k=5.14,s=10))
fml <- Reaction ~ Days + (Days | Subject)

one <- function(r, gen) {
    d <- simulate_rep(gen, SEED0 + r)
    f1 <- tryCatch(suppressWarnings(do.call(rlmer_ransac, c(list(fml, d,
        K=RANSACK, n_starts=1L, seed=SEED0+r), rho_args))), error=function(e) NULL)
    fc <- tryCatch(suppressWarnings(do.call(rlmer_ransac, c(list(fml, d,
        K=RANSACK, n_starts=NSTART, seed=SEED0+r), rho_args))), error=function(e) NULL)
    c(gen = match(gen, GENS),
      phony1 = if (is.null(f1)) NA else phony_of(f1),
      phonyc = if (is.null(fc)) NA else phony_of(fc),
      slope1 = if (is.null(f1)) NA else slope_of(f1),
      slopec = if (is.null(fc)) NA else slope_of(fc))
}

t0 <- Sys.time()
grid <- expand.grid(rep=seq_len(NREPS), gen=GENS, stringsAsFactors=FALSE)
res <- mclapply(seq_len(nrow(grid)), function(i)
    tryCatch(one(grid$rep[i], grid$gen[i]), error=function(e) NULL),
    mc.cores=NCORES)
res <- as.data.frame(do.call(rbind, res[!vapply(res, is.null, logical(1))]))
res$gen <- GENS[res$gen]; runtime <- as.numeric(Sys.time()-t0, units="mins")

agg <- do.call(rbind, lapply(GENS, function(g) {
    s <- res[res$gen==g,]
    data.frame(generator=g, n=nrow(s),
        phony_single=mean(s$phony1, na.rm=TRUE),
        phony_consensus=mean(s$phonyc, na.rm=TRUE),
        slope_rmse_single=sqrt(mean((s$slope1-beta_true[2])^2, na.rm=TRUE)),
        slope_rmse_consensus=sqrt(mean((s$slopec-beta_true[2])^2, na.rm=TRUE))) }))
cat(sprintf("\nransacConsensus: %d cells, %.1f min, n_starts=%d, K<=%d, %d reps/gen\n",
            nrow(res), runtime, NSTART, RANSACK, NREPS))
print(agg, row.names=FALSE, digits=3)
cat(sprintf("\nPOOLED phony rate: single = %.3f -> consensus = %.3f\n",
            mean(res$phony1, na.rm=TRUE), mean(res$phonyc, na.rm=TRUE)))
saveRDS(list(results=res, agg=agg, config=list(NREPS=NREPS, NSTART=NSTART,
    K=RANSACK, runtime_min=runtime)),
    file.path("inst","simulationStudy","ransacConsensus_results.rds"))
cat("saved inst/simulationStudy/ransacConsensus_results.rds\n")
