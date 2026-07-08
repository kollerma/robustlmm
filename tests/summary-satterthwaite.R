## WS16 regression tests: summary(fit, df = ...) and the df-default
## ("auto") logic.
##
## Checks the wiring only (the df numbers themselves are validated in
## tests/vcov-varpar.R and the IF-thread1 anchor):
## 0. df = "none" reproduces the historic 3-column table.
## 1. df = "auto" (default) shows the df when the dimension-based cost
##    estimate is within budget, skips it (with a note) when it is not,
##    and ALWAYS shows it when the influence function is already cached
##    (cooks.distance), even with a zero budget.
## 2. df = "satterthwaite" on a single-grouping fit adds "df" and
##    "Pr(>|t|)" columns; p == 2 * pt(-|t|, df), df finite and >= 1.
## 3. Crossed/nested designs ARE covered (WS17): 3a shows a crossed
##    diagonal-RE fit gets its df via the Cameron-Gelbach-Miller
##    multiway covariance, 3b the nested case (top-cluster sandwich).
## 4. The small-J note fires when the smallest grouping factor has
##    < 20 levels.

suppressMessages(require(robustlmm))

makeLong <- function(J, n = 6, seed = 1L, beta = c(1, 0.5)) {
    set.seed(seed)
    id   <- factor(rep(seq_len(J), each = n))
    time <- rep(seq_len(n), J)
    b    <- rnorm(J, 0, 1)[as.integer(id)]
    y    <- beta[1] + beta[2] * time + b + rnorm(J * n, 0, 1)
    data.frame(y, time, id)
}

fit25 <- rlmer(y ~ time + (1 | id), makeLong(25, seed = 11),
               method = "DASvar")

## ---- 0. df = "none": historic 3-column table ----------------------
cd <- summary(fit25, df = "none")$coefficients
stopifnot(identical(colnames(cd),
                    c("Estimate", "Std. Error", "t value")),
          !("Pr(>|t|)" %in% colnames(cd)))

## ---- 1. df = "auto" logic (deterministic: the WORKLOAD is a
## dimension-only formula, independent of machine speed) --------------
## workload is positive
stopifnot(robustlmm:::.satterthwaite_df_workload(fit25) > 0)
## a zero budget forces the auto skip on a fit with no cached IF, with an
## explanatory note (uses a fresh fit so no earlier call has cached it)
local({
    fskip <- rlmer(y ~ time + (1 | id), makeLong(20, seed = 31),
                   method = "DASvar")
    stopifnot(is.null(robustlmm:::.IFfullCached(fskip)))
    op <- options(robustlmm.summary.df.max = 0); on.exit(options(op))
    ss <- summary(fskip)
    stopifnot(!("Pr(>|t|)" %in% colnames(ss$coefficients)),
              !is.null(ss$dfNote),
              grepl("skipped", ss$dfNote))
})
## fit25 (DASvar, n = 150) is within the default 5000 cutoff -> computes
## by default (and this populates the IF cache)
stopifnot(robustlmm:::.satterthwaite_df_workload(fit25) <= 5000,
          is.null(robustlmm:::.IFfullCached(fit25)),
          "Pr(>|t|)" %in% colnames(summary(fit25)$coefficients),
          !is.null(robustlmm:::.IFfullCached(fit25)))
## ... and a cached IF is shown even at budget 0 (cache overrides budget)
local({
    op <- options(robustlmm.summary.df.max = 0); on.exit(options(op))
    stopifnot("Pr(>|t|)" %in% colnames(summary(fit25)$coefficients))
})

## ---- 2. df = "satterthwaite" adds df + p, internally consistent ----
s25 <- summary(fit25, df = "satterthwaite")
cs  <- s25$coefficients
stopifnot(identical(colnames(cs),
                    c("Estimate", "Std. Error", "df",
                      "t value", "Pr(>|t|)")),
          all(is.finite(cs[, "df"])),
          all(cs[, "df"] >= 1),
          is.null(s25$dfNote))                 # J = 25 >= 20: no note
## p-value must equal 2 * pt(-|t|, df) exactly
p_expected <- 2 * pt(-abs(cs[, "t value"]), df = cs[, "df"])
stopifnot(max(abs(cs[, "Pr(>|t|)"] - p_expected)) < 1e-12)
## estimates / SE / t untouched by requesting df
stopifnot(max(abs(cs[, "Estimate"]   - cd[, "Estimate"]))   < 1e-12,
          max(abs(cs[, "Std. Error"] - cd[, "Std. Error"])) < 1e-12,
          max(abs(cs[, "t value"]    - cd[, "t value"]))    < 1e-12)

## ---- 2b. IF cache returns exactly the freshly-computed IF ---------
fitc <- rlmer(y ~ time + (1 | id), makeLong(12, seed = 21),
              method = "DAStau")
IF_fresh  <- robustlmm:::implicitIF_full(fitc, use.cache = FALSE)
IF_cached <- robustlmm:::implicitIF_full(fitc)            # computes + caches
IF_again  <- robustlmm:::implicitIF_full(fitc)            # returns the cache
stopifnot(max(abs(IF_cached$IF_beta  - IF_fresh$IF_beta))  == 0,
          max(abs(IF_cached$IF_sigma - IF_fresh$IF_sigma)) == 0,
          identical(IF_again$Jpar, IF_cached$Jpar),
          ## the converged theta is the cache key and is unchanged
          isTRUE(all.equal(getME(fitc, "theta"),
                           fitc@pp$cache.IFfull.theta,
                           check.attributes = FALSE)))

## ---- 3a. crossed REs (1|g1)+(1|g2): df IS computed via the
## Cameron-Gelbach-Miller multiway covariance (WS17 Phase 2) ----------
set.seed(7)
dd <- expand.grid(g1 = factor(1:14), g2 = factor(1:9), rep = 1:2)
dd$y <- 1 + rnorm(14, 0, 1.2)[dd$g1] + rnorm(9, 0, 1.0)[dd$g2] +
            rnorm(nrow(dd))
fitx <- suppressWarnings(rlmer(y ~ 1 + (1 | g1) + (1 | g2), dd,
                               method = "DASvar"))
stopifnot(is.null(robustlmm:::.top_cluster(fitx)))     # crossed
Ax <- robustlmm:::.vcov_varpar(fitx)
stopifnot(isTRUE(attr(Ax, "multiway")),
          min(eigen(Ax, symmetric = TRUE, only.values = TRUE)$values) >=
              -1e-10)                                  # PSD (after fix)
sx <- summary(fitx, df = "satterthwaite")
stopifnot("Pr(>|t|)" %in% colnames(sx$coefficients),
          all(is.finite(sx$coefficients[, "df"])))

## ---- 3b. nested REs (1 | school/class): df IS computed (WS17) ------
## strong, well-identified variance components so the fit is not at a
## boundary (where the df would correctly be suppressed)
set.seed(8)
dn <- expand.grid(class = factor(1:4), school = factor(1:12), rep = 1:6)
csn <- interaction(dn$school, dn$class, drop = TRUE)
dn$y <- 1 + rnorm(12, 0, 1.5)[dn$school] +
            rnorm(nlevels(csn), 0, 1.2)[csn] + rnorm(nrow(dn))
fitn <- suppressWarnings(rlmer(y ~ 1 + (1 | school/class), dn,
                               method = "DASvar"))
stopifnot(all(getME(fitn, "theta") > 1e-4))   # non-boundary
sn <- summary(fitn, df = "satterthwaite")
stopifnot("Pr(>|t|)" %in% colnames(sn$coefficients),
          all(is.finite(sn$coefficients[, "df"])),
          all(sn$coefficients[, "df"] >= 1))

## ---- 4. small-J note fires for J < 20 -----------------------------
fit10 <- rlmer(y ~ time + (1 | id), makeLong(10, seed = 5),
               method = "DASvar")
s10 <- summary(fit10, df = "satterthwaite")
stopifnot("Pr(>|t|)" %in% colnames(s10$coefficients),
          !is.null(s10$dfNote),
          grepl("< 20", s10$dfNote, fixed = TRUE))

## ---- printing must not error in either mode -----------------------
invisible(capture.output(print(s25)))
invisible(capture.output(print(s10)))
invisible(capture.output(print(sx)))

cat("summary-satterthwaite.R: all checks passed\n")
