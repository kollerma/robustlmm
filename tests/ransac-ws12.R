## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## WS12 regression tests: adaptive K, stratified subsampling, basin gate.

suppressMessages({require(robustlmm); require(lme4)})

## --- 1. adaptive K: early stop on an easy problem; K_used reported ----
set.seed(1)
ra <- ransac_lme4(Reaction ~ Days + (Days | Subject), sleepstudy,
                  K = 400L, K_min = 30L, patience = 30L)
stopifnot(ra$K == 400L,            # cap unchanged
          ra$K_used < 400L,        # stopped early
          length(ra$scales) == 400L,
          is.finite(ra$scale))
## adaptive = FALSE runs exactly K
set.seed(1)
rf <- ransac_lme4(Reaction ~ Days + (Days | Subject), sleepstudy,
                  K = 12L, adaptive = FALSE)
stopifnot(rf$K_used == 12L, sum(!is.na(rf$scales)) == 12L)

## --- 2. cluster subsampling keeps every draw identifiable -------------
## The nonsingular draw adds WHOLE clusters (Koller and Stahel 2017): it is
## a random CLUSTER subsample, not a within-cluster stratified draw. On
## many small clusters every draw still leaves >= 2 grouping levels and a
## full-rank fixed design, and the chosen clusters are a prefix of the RNG
## permutation.
set.seed(7); J <- 40L; m <- 3L; n <- J * m
g <- factor(rep(seq_len(J), each = m)); x <- rnorm(n)
y <- 1 + 0.5 * x + rnorm(J, 0, 1)[as.integer(g)] + rnorm(n)
d <- data.frame(y, x, g)
strata <- robustlmm:::.ransac_strata(y ~ x + (1 | g), d)
stopifnot(!is.null(strata), nlevels(strata) == J)
info <- robustlmm:::.ransac_design_info(y ~ x + (1 | g), d)
set.seed(2)
identifiable <- replicate(50, {
    dr <- robustlmm:::.ransac_nonsingular_subsample(y ~ x + (1 | g), d,
                                                    0.5, strata, info)
    dr$ok && length(unique(g[dr$idx])) >= 2L
})
stopifnot(all(identifiable))                          # every draw usable
## crossed/nested: picks the finest factor; the cluster subsample keeps
## ceiling(sub_frac * n) observations = 9 whole Subjects of sleepstudy.
rs <- ransac_lme4(Reaction ~ Days + (Days | Subject), sleepstudy, K = 5L)
stopifnot(isTRUE(rs$stratified),
          length(rs$subset) == ceiling(0.5 * nrow(sleepstudy)),
          length(unique(sleepstudy$Subject[rs$subset])) == 9L)
## stratify = FALSE -> plain SRS (n_sub = ceiling(sub_frac n))
rsf <- ransac_lme4(Reaction ~ Days + (Days | Subject), sleepstudy,
                   K = 3L, stratify = FALSE)
stopifnot(!isTRUE(rsf$stratified),
          rsf$n_sub == ceiling(0.5 * nrow(sleepstudy)))

## --- 3. basin radius + redescender gate ------------------------------
lm0 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
r <- ransac_basin_radius(lm0)                        # default c = 4.685
stopifnot(is.finite(as.numeric(r)), as.numeric(r) > 0,
          ## scales linearly in c
          abs(as.numeric(ransac_basin_radius(lm0, cc = 2)) -
              as.numeric(r) * 2 / 4.685) < 1e-8,
          ## r* = c sigma / (2 max||x_j||); max row norm = sqrt(82)
          abs(attr(r, "max_xnorm") - sqrt(82)) < 1e-6)
## redescending detector
stopifnot(robustlmm:::.is_redescending(bisquarePsi),
          !robustlmm:::.is_redescending(smoothPsi),
          !robustlmm:::.is_redescending(cPsi))
## monotone rho -> no gate (returns NULL)
rfit <- suppressWarnings(rlmer(Reaction ~ Days + (Days | Subject),
                               sleepstudy, method = "DASvar",
                               rho.e = smoothPsi))
stopifnot(is.null(robustlmm:::.ransac_check_basin(rfit, lm0, smoothPsi)))
## redescending + in-basin (clean fit near start) -> no warning, in_basin
chk <- robustlmm:::.ransac_check_basin(rfit, lm0, bisquarePsi)
stopifnot(isTRUE(chk$in_basin), chk$distance <= chk$radius)
## redescending + far start -> warning, out of basin
ss2 <- sleepstudy; ss2$Reaction <- ss2$Reaction + 30 * ss2$Days
lm_far <- lmer(Reaction ~ Days + (Days | Subject), ss2)
w <- tryCatch(robustlmm:::.ransac_check_basin(rfit, lm_far, bisquarePsi),
              warning = function(w) conditionMessage(w))
stopifnot(is.character(w), grepl("basin", w))

## --- 4. phony random-effects correlation check (WS12 follow-up b) -----
## .re_max_abscor matches the lme4 VarCorr correlation; a diagonal RE has
## none (NA -> no check); the warning fires only above threshold.
rho_fit <- robustlmm:::.re_max_abscor(rfit)
stopifnot(is.finite(rho_fit),
          abs(rho_fit - abs(attr(VarCorr(lm0)$Subject,
                                 "correlation")[1, 2])) < 0.05)
## clean fit (rho ~ 0.07): no phony warning at the 0.99 default
chk_ok <- robustlmm:::.ransac_check_phony(rfit)
stopifnot(!isTRUE(chk_ok$phony))
ok <- withCallingHandlers({ robustlmm:::.ransac_check_phony(rfit); TRUE },
                          warning = function(w) FALSE)
stopifnot(isTRUE(ok))                                  # no warning emitted
## lowering the threshold below the fitted rho fires the warning
wp <- tryCatch(robustlmm:::.ransac_check_phony(rfit, threshold = 0.01),
               warning = function(w) conditionMessage(w))
stopifnot(is.character(wp), grepl("[Pp]hony", wp))
## diagonal random effects -> no correlation parameter -> NA -> no check
fdiag <- suppressWarnings(rlmer(Reaction ~ Days + (1 | Subject),
                                sleepstudy, method = "DASvar"))
stopifnot(is.na(robustlmm:::.re_max_abscor(fdiag)),
          is.null(robustlmm:::.ransac_check_phony(fdiag)))

## --- 5. multi-start consensus (WS12-D) -------------------------------
## ransac_lme4(n_keep) returns up to n_keep distinct candidate starts
rk <- ransac_lme4(Reaction ~ Days + (Days | Subject), sleepstudy,
                  K = 40L, n_keep = 5L, seed = 5)
stopifnot(length(rk$candidates) >= 2L, length(rk$candidates) <= 5L,
          all(vapply(rk$candidates,
                     function(f) methods::is(f, "lmerMod"), logical(1))))
## candidates are distinct (deduped by fixed effects)
bmat <- t(vapply(rk$candidates, lme4::fixef, numeric(2)))
stopifnot(nrow(unique(round(bmat, 6))) == nrow(bmat))
## n_starts = 1 is the previous behaviour: rlmerMod, no consensus attr
f1 <- suppressWarnings(rlmer_ransac(Reaction ~ Days + (Days | Subject),
                                    sleepstudy, K = 30L, method = "DASvar",
                                    seed = 7))
stopifnot(methods::is(f1, "rlmerMod"), is.null(attr(f1, "consensus")))
## n_starts > 1: consensus table; exactly one chosen; chosen interior here
f5 <- suppressWarnings(rlmer_ransac(Reaction ~ Days + (Days | Subject),
                                    sleepstudy, K = 40L, n_starts = 5L,
                                    seed = 7, method = "DASvar"))
co <- attr(f5, "consensus")
stopifnot(methods::is(f5, "rlmerMod"), is.data.frame(co),
          sum(co$chosen) == 1L,
          isTRUE(co$interior[co$chosen]),           # picked an interior fit
          ## the chosen interior fit has the lowest residual scale among
          ## the interior candidates
          co$resid_scale[co$chosen] ==
              min(co$resid_scale[co$interior]))
stopifnot(!isTRUE(robustlmm:::.ransac_check_phony(f5)$phony))

## --- 6. all-starts-phony consensus warning branch ---------------------
## Honesty note (2026-07-03): a MIXED phony/interior consensus (one start
## caught by a phony attractor, another interior, consensus picking the
## interior one) is only reproduced reliably by the full bisquare
## simulation in inst/simulationStudy/ransacConsensus.R; a compact
## deterministic version is knife-edge (which start goes phony depends on
## platform-level numerics), i.e. exactly the flaky test we must not
## ship. What IS deterministic is the all-starts-phony branch: a
## generator whose random slope is EXACTLY the random intercept (true
## rho = 1, near-zero residual noise) drives every start to |rho-hat| = 1,
## so the consensus must (a) warn that all starts are phony, (b) fall
## back to the best residual scale overall, and (c) return a fit the
## phony detector still flags.
set.seed(11L); J6 <- 12L; m6 <- 8L
g6 <- factor(rep(seq_len(J6), each = m6))
x6 <- rep(seq_len(m6) - 1, J6) / 2
u6 <- rnorm(J6, 0, 3)
d6 <- data.frame(y = 1 + 0.5 * x6 + u6[as.integer(g6)] * (1 + x6) +
                     rnorm(J6 * m6, 0, 0.05),
                 x = x6, g = g6)
warns6 <- character(0)
f6 <- withCallingHandlers(
    rlmer_ransac(y ~ x + (x | g), d6, K = 15L, n_starts = 3L,
                 seed = 11L, method = "DASvar"),
    warning = function(w) {
        warns6 <<- c(warns6, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
stopifnot(any(grepl("all 3 RANSAC starts converged to a phony", warns6)))
co6 <- attr(f6, "consensus")
stopifnot(methods::is(f6, "rlmerMod"), is.data.frame(co6),
          !any(co6$interior),                      # every start phony
          sum(co6$chosen) == 1L,
          ## fallback picks the best residual scale overall
          co6$resid_scale[co6$chosen] == min(co6$resid_scale),
          ## and the detector still flags the returned fit
          isTRUE(suppressWarnings(
              robustlmm:::.ransac_check_phony(f6))$phony))

cat("ransac-ws12.R: all checks passed\n")
