## WS11 PLAN section 5.4: leverage-contamination breakdown study.
##
## The gating demonstration for the Mallows design-weight feature: under
## (x, y) leverage contamination -- a fraction of observations moved to
## an extreme x and placed on a *different* (contaminating) regression
## line -- lmer and the plain RSE break down in the slope (a good
## leverage point has a small residual, so the psi-functions do not
## catch it), while Mallows-RSE (design.weights = "mcd") downweights the
## extreme-x points and stays near the clean-model estimate.
##
## Output: bias and RMSE of beta1 (slope) and bias of theta for
## lmer / RSE / Mallows-RSE over a grid of contamination fractions.
## "If it does not [stay near the clean bias], the feature does not
## ship." (PLAN-WS11-mallows.md section 5.4)

suppressMessages({
    library(robustlmm); library(lme4); library(parallel)
})

## ---- configuration ------------------------------------------------
J        <- 30L          # subjects
m        <- 6L           # observations per subject
n        <- J * m
beta0    <- 0; beta1 <- 1
sigma_b  <- 1; sigma_e <- 1
slope_c  <- -1           # contaminating slope (opposite sign)
x_lev    <- 7            # leverage x location (~7 robust SD)
eps_grid <- c(0, 0.02, 0.05, 0.10)
nreps    <- 1000L
ncores   <- 10L          # full blast (M1 Max)
seed0    <- 20260618L

makeClean <- function() {
    g <- factor(rep(seq_len(J), each = m))
    x <- rnorm(n)
    b <- rnorm(J, 0, sigma_b)[as.integer(g)]
    y <- beta0 + beta1 * x + b + rnorm(n, 0, sigma_e)
    data.frame(y = y, x = x, g = g)
}

## move a fraction eps of observations to leverage points: extreme x,
## response on the contaminating line beta0 + slope_c * x.
contaminate <- function(d, eps) {
    k <- floor(eps * nrow(d))
    if (k < 1L) return(d)
    idx <- sample.int(nrow(d), k)
    d$x[idx] <- x_lev + rnorm(k, 0, 0.1)
    d$y[idx] <- beta0 + slope_c * d$x[idx] + rnorm(k, 0, sigma_e)
    d
}

grab <- function(f, robust) {
    if (is.null(f)) return(c(b1 = NA_real_, theta = NA_real_))
    b  <- if (robust) robustlmm:::.fixef(f) else lme4::fixef(f)
    c(b1 = unname(b[2]), theta = unname(getME(f, "theta")[1]))
}

fitAll <- function(d) {
    fl <- tryCatch(suppressMessages(suppressWarnings(
              lmer(y ~ x + (1 | g), d, REML = TRUE))), error = function(e) NULL)
    fr <- tryCatch(suppressMessages(suppressWarnings(
              rlmer(y ~ x + (1 | g), d))), error = function(e) NULL)
    fm <- tryCatch(suppressMessages(suppressWarnings(
              rlmer(y ~ x + (1 | g), d, design.weights = "mcd"))),
              error = function(e) NULL)
    rbind(lmer = grab(fl, FALSE), RSE = grab(fr, TRUE),
          Mallows = grab(fm, TRUE))
}

oneRep <- function(r) {
    set.seed(seed0 + r)
    d <- makeClean()
    lapply(eps_grid, function(eps) fitAll(contaminate(d, eps)))
}

t0   <- Sys.time()
reps <- mclapply(seq_len(nreps), oneRep, mc.cores = ncores)
runtime <- as.numeric(Sys.time() - t0, units = "mins")
cat(sprintf("ran %d reps on %d cores in %.1f min\n", nreps, ncores, runtime))

## ---- aggregate ----------------------------------------------------
methods <- c("lmer", "RSE", "Mallows")
agg <- expand.grid(method = methods, eps = eps_grid,
                   stringsAsFactors = FALSE)
agg$bias_b1 <- agg$rmse_b1 <- agg$bias_theta <- agg$n_ok <- NA_real_
for (i in seq_len(nrow(agg))) {
    ei <- match(agg$eps[i], eps_grid)
    b1 <- vapply(reps, function(rr) rr[[ei]][agg$method[i], "b1"],    numeric(1))
    th <- vapply(reps, function(rr) rr[[ei]][agg$method[i], "theta"], numeric(1))
    agg$bias_b1[i]    <- mean(b1, na.rm = TRUE) - beta1
    agg$rmse_b1[i]    <- sqrt(mean((b1 - beta1)^2, na.rm = TRUE))
    agg$bias_theta[i] <- mean(th, na.rm = TRUE) - sigma_b
    agg$n_ok[i]       <- sum(!is.na(b1))
}
agg <- agg[order(agg$eps, agg$method), ]
print(agg, row.names = FALSE, digits = 3)

out <- list(agg = agg,
            config = list(J = J, m = m, beta1 = beta1, slope_c = slope_c,
                          x_lev = x_lev, eps_grid = eps_grid, nreps = nreps,
                          seed0 = seed0, runtime_min = runtime))
saveRDS(out, file.path("inst", "simulationStudy",
                       "leverageBreakdown_results.rds"))
cat("saved inst/simulationStudy/leverageBreakdown_results.rds\n")
