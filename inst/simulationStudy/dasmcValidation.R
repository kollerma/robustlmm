## Validation study for the EXPERIMENTAL Monte-Carlo DAS-tau calibration
## ("DASmc": options robustlmm.dastau.mc / robustlmm.dasmc.nsim /
## robustlmm.dasmc.seed, see ?`robustlmm-options`). The closed-form
## DASvar T_b assumes sphericity of the whitened linearized predictor,
## which leaves a small calibration residual in the fitted correlation
## of non-diagonal blocks; the MC path computes the DAStau fixed point
## self-consistently for any block dimension. Headline cells: balanced
## ar1 repeated measures, NC = 4 levels (block dimension s = 4, i.e.
## beyond the Gauss-Hermite s <= 2 limit), reps = 2, J = 200 subjects,
## rho0 in {-0.5, +0.6}; per seed the estimators lmer (REML reference),
## rlmer DASvar and rlmer DAStau+MC (N_MC = 3e5, robustlmm.dasmc.seed =
## data seed so the residual MC calibration error decorrelates across
## replicates) are fitted on IDENTICAL data; all contrasts are paired.
##
## ---- FINDING (nrep = 200/cell, 2026-07-03; shipped in ----------------
## dasmcValidation_results.rds) -- paired rho-hat offsets vs lmer:
##   rho0 = -0.5: DASvar - lmer = -0.0037 (t = -4.0)
##                DASmc  - lmer = -0.0001 (t = -0.13)
##   rho0 = +0.6: DASvar - lmer = +0.0026 (t = +3.2)
##                DASmc  - lmer = -0.0010 (t = -1.2)
## The MC calibration removes the ~0.004 DASvar sphericity residual
## (|d| <= 0.001 at the pre-registered cell rho0 = -0.5); these are the
## numbers cited in the 3.5.0 feature vignette. Full study (also J =
## 100 and an s = 2 GH-vs-MC comparison): dasmc-validation.R in the
## robustlmm-open-problems archive.
## ----------------------------------------------------------------------

suppressMessages({ library(robustlmm); library(lme4); library(parallel) })

## NREP = 200 for the shipped results; the default here is modest so a
## test run finishes quickly. ~8s per replicate-triple on an M1 Max.
NREP   <- 20L
NCORES <- 5L
J      <- 200L
RHO0   <- c(-0.5, 0.6)
NMC    <- 3e5      # MC draws for the s = 4 DASmc path
SIGMA0 <- 1.5      # true marginal RE sd
RESID  <- 0.7      # residual sd

makeData <- function(seed, nc, reps, J, rho0) {
    set.seed(seed)
    R  <- rho0^abs(outer(1:nc, 1:nc, "-"))
    ch <- chol(SIGMA0^2 * R)
    fl <- factor(1:nc)
    d <- do.call(rbind, lapply(seq_len(J), function(j) {
        b <- as.vector(crossprod(ch, rnorm(nc)))
        do.call(rbind, lapply(seq_len(reps), function(rr) {
            data.frame(subj = j, f = fl, y = 10 + b + rnorm(nc, 0, RESID))
        }))
    }))
    d$subj <- factor(d$subj)
    d
}

sdRho <- function(fit) {
    vc <- VarCorr(fit)$subj
    c(sd  = as.numeric(attr(vc, "stddev")[1]),
      rho = as.numeric(attr(vc, "correlation")[1, 2]))
}

## one replicate: all three estimators on the same data
one <- function(seed, rho0) {
    dat  <- makeData(seed, nc = 4L, reps = 2L, J = J, rho0 = rho0)
    form <- y ~ 1 + ar1(0 + f | subj)
    row_for <- function(label, fitter) tryCatch({
        v <- sdRho(fitter())
        data.frame(seed = seed, rho0 = rho0, est = label,
                   sd = v[["sd"]], rho = v[["rho"]])
    }, error = function(e) NULL)
    rbind(
        row_for("lmer", function()
            lmer(form, data = dat, REML = TRUE,
                 control = lmerControl(calc.derivs = FALSE))),
        row_for("DASvar", function()
            rlmer(form, data = dat, method = "DASvar")),
        row_for("DASmc", function() {
            opts <- options(robustlmm.dastau.mc = TRUE,
                            robustlmm.dasmc.nsim = NMC,
                            robustlmm.dasmc.seed = seed)
            on.exit(options(opts))
            suppressMessages(rlmer(form, data = dat, method = "DAStau"))
        }))
}

t0 <- Sys.time()
## seed bases 2000/4000 reproduce the shipped run's J = 200 cells exactly
SEED_BASE <- c(2000L, 4000L)
results <- do.call(rbind, lapply(seq_along(RHO0), function(k) {
    seeds <- SEED_BASE[k] + seq_len(NREP)
    do.call(rbind, mclapply(seeds, one, rho0 = RHO0[k], mc.cores = NCORES))
}))
runtime <- as.numeric(Sys.time() - t0, units = "mins")

## marginal summaries and paired contrasts (per rho0 cell)
marginal <- function(df) do.call(rbind, by(df, df[c("est", "rho0")],
    function(g) data.frame(est = g$est[1], rho0 = g$rho0[1], n = nrow(g),
        bias_rho = mean(g$rho) - g$rho0[1], se_rho = sd(g$rho)/sqrt(nrow(g)),
        bias_sd = mean(g$sd) - SIGMA0, se_sd = sd(g$sd)/sqrt(nrow(g)))))
paired <- function(df, a, b) do.call(rbind, by(df, df["rho0"], function(g) {
    m <- merge(g[g$est == a, c("seed", "rho", "sd")],
               g[g$est == b, c("seed", "rho", "sd")],
               by = "seed", suffixes = c(".a", ".b"))
    drho <- m$rho.a - m$rho.b
    data.frame(contrast = paste(a, "-", b), rho0 = g$rho0[1],
               n_pairs = nrow(m), d_rho = mean(drho),
               se_d_rho = sd(drho)/sqrt(nrow(m)),
               t_rho = mean(drho)/(sd(drho)/sqrt(nrow(m))),
               d_sd = mean(m$sd.a - m$sd.b))
}))

marg       <- marginal(results)
p_var_lmer <- paired(results, "DASvar", "lmer")
p_mc_lmer  <- paired(results, "DASmc", "lmer")
p_mc_var   <- paired(results, "DASmc", "DASvar")

cat(sprintf("\ndasmcValidation: %d fits, %.1f min, nrep = %d/cell, J = %d\n",
            nrow(results), runtime, NREP, J))
cat("\nmarginal:\n");             print(marg, row.names = FALSE, digits = 3)
cat("\npaired DASvar - lmer:\n"); print(p_var_lmer, row.names = FALSE, digits = 3)
cat("\npaired DASmc - lmer:\n");  print(p_mc_lmer, row.names = FALSE, digits = 3)
cat("\npaired DASmc - DASvar:\n"); print(p_mc_var, row.names = FALSE, digits = 3)

saveRDS(list(results = results, marg = marg, p_var_lmer = p_var_lmer,
             p_mc_lmer = p_mc_lmer, p_mc_var = p_mc_var,
             config = list(NREP = NREP, J = J, NC = 4L, reps = 2L,
                           RHO0 = RHO0, NMC = NMC, SIGMA0 = SIGMA0,
                           RESID = RESID, runtime_min = runtime)),
        file.path("inst", "simulationStudy", "dasmcValidation_results.rds"))
cat("saved inst/simulationStudy/dasmcValidation_results.rds\n")
