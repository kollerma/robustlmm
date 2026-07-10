## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Nonsingular subsampling in ransac_lme4 (Koller and Stahel 2017,
## Algorithm 1). Each subsample is drawn full-rank BY CONSTRUCTION with a
## Gaxpy-variant LU decomposition that skips linearly-dependent
## observations, rather than the old draw-then-check-then-repair. The draw
## reduces to plain (cluster) random subsampling when there is no
## collinearity.
##
## Asserts: (1) the C routine returns a nonsingular p-subset and reduces to
## the first p observations on a full-rank design; (2) the draw ALWAYS
## yields a full-rank fixed design (0 degenerate draws) even with a rare
## factor level, and the fit never drops a coefficient; (3) no collinearity
## -> the draw is a uniform random cluster subsample (matches the
## permutation prefix); (4) determinism under an outer seed; (5) the
## refinement-singularity guard (Remark 2 (i)).

suppressMessages(require(robustlmm))

LU  <- robustlmm:::nonsingularSubsampleLU
ns  <- robustlmm:::.ransac_nonsingular_subsample
di  <- robustlmm:::.ransac_design_info
str <- robustlmm:::.ransac_strata

## ==== (1) C-routine unit ============================================
## Full-rank continuous design: the selected columns are exactly the first
## p of the walk order, and index a full-rank submatrix.
set.seed(1)
X <- cbind(1, rnorm(30), rnorm(30)); p <- ncol(X); Xt <- t(X)
for (s in 1:20) {
    ord <- sample.int(30) - 1L
    r <- LU(Xt, as.integer(ord), 1e-7)
    stopifnot(!r$singular, r$rank == p,
              identical(as.integer(r$selected), as.integer(ord[1:p])),
              r$n_skipped == 0L,
              qr(X[r$selected + 1L, , drop = FALSE])$rank == p)
}
## Collinear structure: a factor whose dummies alias under some orderings.
## Over many permutations the routine still returns a NONSINGULAR p-subset
## (it skips the aliased columns) whenever the full design has full rank.
set.seed(2)
fac <- factor(sample(letters[1:3], 40, replace = TRUE))
Xc  <- model.matrix(~ fac); pc <- ncol(Xc)
skipped_any <- FALSE
for (s in 1:50) {
    ord <- sample.int(40) - 1L
    r <- LU(t(Xc), as.integer(ord), 1e-7)
    stopifnot(!r$singular, r$rank == pc,
              qr(Xc[r$selected + 1L, , drop = FALSE])$rank == pc)
    if (r$n_skipped > 0L) skipped_any <- TRUE
}
stopifnot(skipped_any)   # collinear candidates ARE encountered and skipped
## A genuinely rank-deficient full design (duplicated column) -> singular.
Xd <- cbind(1, X[, 2], X[, 2])
rd <- LU(t(Xd), as.integer(sample.int(30) - 1L), 1e-7)
stopifnot(rd$singular, rd$rank < ncol(Xd))
cat("(1) C routine: nonsingular by construction, reduces to first p: ok\n")

## ==== (2) guarantee: rare factor level, 0 degenerate draws ===========
## Level "d" is present in only 2 of 24 clusters; a uniform cluster
## subsample frequently misses it, but the nonsingular draw must always
## include a cluster carrying it (the LU core needs the d column for full
## rank), so no draw is degenerate and no fit drops a coefficient.
set.seed(11); nb <- 24
g <- factor(rep(1:nb, each = 4))
f <- factor(rep("a", nb * 4), levels = letters[1:4])
f[g %in% 3:13]  <- sample(c("b", "c"), sum(g %in% 3:13),  replace = TRUE)
f[g %in% 14:24] <- sample(c("a", "b", "c"), sum(g %in% 14:24), replace = TRUE)
f[g %in% 1:2]   <- "d"
da <- data.frame(g = g, f = f, x = rnorm(nb * 4))
re <- rnorm(nb, sd = 0.5)
da$y <- 1 + as.integer(da$f) + da$x + rnorm(nb * 4) + re[as.integer(da$g)]
p_exp   <- ncol(model.matrix(~ f + x, da))
info    <- di(y ~ f + x + (1 | g), da)
strata  <- str(y ~ f + x + (1 | g), da)

N <- 60L; degen <- 0L; dropped <- 0L; not_fullrank <- 0L; incl_d <- 0L
for (s in seq_len(N)) {
    set.seed(s)
    dr <- ns(y ~ f + x + (1 | g), da, 0.5, strata, info)
    if (!dr$ok) degen <- degen + 1L
    Xs <- model.matrix(~ f + x, da[dr$idx, , drop = FALSE])
    if (ncol(Xs) != p_exp || qr(Xs)$rank != p_exp) not_fullrank <- not_fullrank + 1L
    if ("d" %in% da$f[dr$idx]) incl_d <- incl_d + 1L
}
stopifnot(degen == 0L, not_fullrank == 0L, incl_d == N)
## and the full ransac_lme4 never returns a dropped-coefficient fit
for (s in 1:20) {
    ra <- suppressWarnings(ransac_lme4(y ~ f + x + (1 | g), da, K = 8L, seed = s))
    if (!is.null(ra$fit) && length(lme4::fixef(ra$fit)) != p_exp)
        dropped <- dropped + 1L
}
stopifnot(dropped == 0L)
cat(sprintf(paste0("(2) rare level: %d/%d draws degenerate, %d not full ",
                   "rank, level-d present in %d/%d, dropped-coef fits %d\n"),
            degen, N, not_fullrank, incl_d, N, dropped))

## ==== (3) no-collinearity equivalence ================================
## On a purely continuous design the draw is a uniform random CLUSTER
## subsample: the selected clusters are exactly the first
## ceil(sub_frac * .) clusters of the same RNG permutation.
set.seed(101)
dc   <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:12, 5)))
infC <- di(y ~ x + (1 | g), dc); strC <- str(y ~ x + (1 | g), dc)
set.seed(202); d1 <- ns(y ~ x + (1 | g), dc, 0.5, strC, infC)
## reconstruct the permutation the draw consumed and the target size
set.seed(202); perm <- sample(levels(as.factor(strC)))
target <- sum(pmax(1L, ceiling(0.5 * as.integer(table(strC)))))   # 36 obs
## clusters of 5 obs -> need ceiling(36/5) = 8 clusters (40 obs)
sel_cl <- sort(unique(as.character(dc$g[d1$idx])))
stopifnot(d1$n_singular == 0L, d1$ok,
          identical(sel_cl, sort(perm[seq_len(8L)])),
          length(d1$idx) == 40L)
## distributional check: each cluster's inclusion frequency ~ 8/12.
inc <- integer(12); M <- 400L
for (s in seq_len(M)) {
    set.seed(1000L + s)
    di2 <- ns(y ~ x + (1 | g), dc, 0.5, strC, infC)
    present <- unique(as.integer(dc$g[di2$idx]))
    inc[present] <- inc[present] + 1L
}
freq <- inc / M
stopifnot(all(abs(freq - 8 / 12) < 0.08))
cat(sprintf("(3) no-collinearity equivalence: cluster inclusion freq in [%.3f, %.3f] ~ %.3f\n",
            min(freq), max(freq), 8 / 12))

## ==== (4) determinism under an outer seed ============================
set.seed(7); a1 <- ns(y ~ x + (1 | g), dc, 0.5, strC, infC)
set.seed(7); a2 <- ns(y ~ x + (1 | g), dc, 0.5, strC, infC)
stopifnot(identical(a1$idx, a2$idx), identical(a1$n_singular, a2$n_singular))
r1 <- ransac_lme4(y ~ f + x + (1 | g), da, K = 8L, seed = 123L)
r2 <- ransac_lme4(y ~ f + x + (1 | g), da, K = 8L, seed = 123L)
stopifnot(identical(r1$subset, r2$subset),
          identical(r1$n_singular, r2$n_singular),
          isTRUE(all.equal(fixef(r1$fit), fixef(r2$fit))))
cat("(4) determinism: ok\n")

## ==== (5) refinement-singularity guard (Remark 2 (i)) ================
## The guard inspects a converged fit's e-side robustness weights and warns
## when the design restricted to the positively-weighted observations is
## rank-deficient. On a clean fit with a positive-weight psi the design
## stays full rank (no warning). Reproducing a genuine redescender-induced
## collapse deterministically is config/seed-sensitive -- the fitted factor
## coefficient tends to absorb level shifts, so the observations that would
## have to be rejected are instead fitted -- so the guard mechanics (weight
## threshold, rank check, warning) are unit-tested directly via `eps`.
set.seed(3)
J <- 12; m <- 6; nn <- J * m
dg <- data.frame(y = rnorm(nn), x = rnorm(nn),
                 fac = factor(ifelse(rep(1:J, each = m) <= 3, "b", "a")),
                 g = factor(rep(1:J, each = m)))
dg$y <- 1 + dg$x + rnorm(J)[as.integer(dg$g)] + rnorm(nn)
fit <- suppressWarnings(rlmer(y ~ fac + x + (1 | g), dg, method = "DASvar"))
## clean fit, default eps: full-rank positive-weight design, no warning
dgn <- robustlmm:::.ransac_refinement_rank(fit)
stopifnot(!dgn$singular, dgn$rank == dgn$p, dgn$n_pos == dgn$n)
ok_nowarn <- withCallingHandlers(
    { robustlmm:::.ransac_check_refinement(fit); TRUE },
    warning = function(w) FALSE)
stopifnot(isTRUE(ok_nowarn))
## force the singular branch: an eps above every weight leaves no
## positive-weight observation -> rank-deficient -> the guard must warn.
dgn2 <- robustlmm:::.ransac_refinement_rank(fit, eps = 2)
stopifnot(dgn2$n_pos == 0L, dgn2$singular)
wmsg <- tryCatch(robustlmm:::.ransac_check_refinement(fit, eps = 2),
                 warning = function(w) conditionMessage(w))
stopifnot(is.character(wmsg), grepl("rank-deficient", wmsg),
          grepl("Remark 2", wmsg))
cat("(5) refinement-singularity guard: ok\n")

## ---- (6) auto-reseed loop: max_tries then warn, accept-on-success -------
## Instrument .ransac_refinement_rank so rlmer_ransac's guard loop sees a
## controlled singular/non-singular verdict, and check the retry behaviour.
ns_env <- asNamespace("robustlmm")
orig_rr <- get(".ransac_refinement_rank", ns_env)
unlockBinding(".ransac_refinement_rank", ns_env)
## (a) always singular: loop runs max_tries times, then warns and returns
calls_a <- 0L
assign(".ransac_refinement_rank",
       function(fit, eps = 1e-4) { calls_a <<- calls_a + 1L
           list(rank = 1L, p = 2L, n_pos = 3L, n = 180L, singular = TRUE) },
       ns_env)
wa <- NULL
fa <- withCallingHandlers(
    rlmer_ransac(Reaction ~ Days + (Days | Subject), sleepstudy,
                 K = 15L, seed = 11, max_tries = 3L),
    warning = function(cnd) { wa <<- c(wa, conditionMessage(cnd))
                              invokeRestart("muffleWarning") })
stopifnot(is(fa, "rlmerMod"),
          calls_a >= 3L,                       # loop re-seeded max_tries times
          any(grepl("re-seeded RANSAC starts", wa)),
          any(grepl("Remark 2", wa)))
## (b) singular once then ok: loop stops early, no re-seed warning
calls_b <- 0L
assign(".ransac_refinement_rank",
       function(fit, eps = 1e-4) { calls_b <<- calls_b + 1L
           list(rank = 2L, p = 2L, n_pos = 180L, n = 180L,
                singular = (calls_b == 1L)) },
       ns_env)
wb <- NULL
fb <- withCallingHandlers(
    rlmer_ransac(Reaction ~ Days + (Days | Subject), sleepstudy,
                 K = 15L, seed = 12, max_tries = 5L),
    warning = function(cnd) { wb <<- c(wb, conditionMessage(cnd))
                              invokeRestart("muffleWarning") })
assign(".ransac_refinement_rank", orig_rr, ns_env)   # restore
lockBinding(".ransac_refinement_rank", ns_env)
stopifnot(is(fb, "rlmerMod"), !any(grepl("re-seeded", wb)))
cat("(6) auto-reseed loop (max_tries + accept-on-success): ok\n")

cat("ransac-nonsingular: all tests passed\n")
