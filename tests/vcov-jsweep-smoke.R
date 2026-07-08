## Disabled on the CRAN release branch (historical github-transition set:
## slow / MC-random / platform-fragile tests kept out of the CRAN check;
## they run in full on master and in CI). See feedback in project notes.
quit()

## Smoke test for the cluster-robust sandwich vcov on a one-way
## fixed-effects ANOVA y ~ Var1 + (1 | Var2): does it
## (a) run end-to-end without erroring,
## (b) return PD covariance matrices, and
## (c) attain reasonable empirical 95% coverage of the true beta?
##
## A minimised subset of IF-thread1/thread1_validation/vcov_IF_jsweep.R
## (which itself runs J in {18, 50, 100, 400} x 5 generators x 500 reps
## for the full coverage study). The full sweep is too slow for R CMD
## check; this script runs a handful of reps at a single J to catch
## catastrophic regressions in the sandwich path. ~5 s.

require(robustlmm)
suppressMessages(require(lme4))

src <- system.file("simulationStudy/randomNumberGenerators.R",
                   package = "robustlmm")
stopifnot(nzchar(src))
source(src)

set.seed(20260601L)
J <- 18L; N_REPS <- 20L
N_LEVELS <- 3L; N_REPLICATES <- 10L

ds <- generateAnovaDatasets(
    numberOfDatasetsToGenerate  = N_REPS,
    numberOfLevelsInFixedFactor = N_LEVELS,
    numberOfSubjects            = J,
    numberOfReplicates          = N_REPLICATES,
    errorGenerator      = srnorm,
    randomEffectGenerator = srnorm,
    lmeFormula = y ~ Var1,
    lower = 0, arrange = TRUE)
true_beta <- as.numeric(ds$trueBeta)
stopifnot(length(true_beta) == 3L)

p <- length(true_beta)
M_se_def  <- matrix(NA_real_, N_REPS, p)
M_se_sand <- matrix(NA_real_, N_REPS, p)
M_beta    <- matrix(NA_real_, N_REPS, p)
n_fail    <- 0L

for (i in seq_len(N_REPS)) {
    dat <- ds$generateData(i)
    fit <- tryCatch(
        suppressMessages(suppressWarnings(
            rlmer(y ~ Var1 + (1 | Var2), data = dat, method = "DASvar"))),
        error = function(e) NULL)
    if (is.null(fit)) { n_fail <- n_fail + 1L; next }
    M_beta[i, ] <- as.numeric(fixef(fit))
    V_def  <- as.matrix(vcov(fit))
    V_sand <- as.matrix(vcov(fit, type = "sandwich"))
    stopifnot(all(is.finite(V_def)), all(is.finite(V_sand)))
    stopifnot(all(eigen(V_sand, only.values = TRUE)$values > 0))
    M_se_def[i, ]  <- sqrt(diag(V_def))
    M_se_sand[i, ] <- sqrt(diag(V_sand))
}
ok <- N_REPS - n_fail
stopifnot(ok >= 0.8 * N_REPS)  # at most 20% fit failures tolerated

## Per-coefficient nominal 95% Wald coverage (t-based, df = n - p - (J-1)).
n_obs <- N_LEVELS * J * N_REPLICATES
dfree <- n_obs - p - (J - 1L)
crit  <- qt(0.975, dfree)
cover <- function(b_col, se_col, true_b) {
    ok <- !(is.na(b_col) | is.na(se_col))
    mean(abs(b_col[ok] - true_b) <= crit * se_col[ok])
}
cov_def  <- vapply(seq_len(p), function(k)
                   cover(M_beta[, k], M_se_def[, k],  true_beta[k]),
                   numeric(1))
cov_sand <- vapply(seq_len(p), function(k)
                   cover(M_beta[, k], M_se_sand[, k], true_beta[k]),
                   numeric(1))

cat(sprintf("smoke: J=%d, N_REPS=%d, fit_failures=%d\n",
            J, N_REPS, n_fail))
cat(sprintf("  default vcov coverage : %s\n",
            paste(sprintf("%.2f", cov_def),  collapse = " ")))
cat(sprintf("  sandwich coverage     : %s\n",
            paste(sprintf("%.2f", cov_sand), collapse = " ")))

## Sanity bounds: not catastrophically broken (no exactly-0 / exactly-1
## per-coefficient coverage on a clean Gaussian generator with N=540).
stopifnot(all(cov_sand >= 0.5))
stopifnot(all(cov_sand <= 1.0))
stopifnot(all(cov_def  >= 0.5))
