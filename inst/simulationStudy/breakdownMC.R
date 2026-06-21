########################################################
## breakdownMC.R                                       ##
##                                                     ##
## Monte-Carlo breakdown-point study across the three  ##
## default robustlmm designs (Dyestuff, Penicillin,    ##
## sleepstudy) and several contamination strategies   ##
## per design. Complements the single-trajectory       ##
## breakdown.R (Koller 2013, Fig 4.5) with a           ##
## distributional view across multiple strategies and  ##
## fitting methods.                                    ##
########################################################

require(robustlmm)
require(ggplot2)

source(system.file("simulationStudy/contaminationStrategies.R",
                   package = "robustlmm"))

path <- system.file("simulationStudy", package = "robustlmm")
ncores <- max(1L, parallel::detectCores() - 1L)

## --------------------------------------------------------------
## Configuration
## --------------------------------------------------------------
N_REPS <- 30L
EPS_GRID <- c(0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50)
SIGMA_TOL <- 5    # broken if sigma_hat / sigma_true outside [1/5, 5]
BETA_TOL  <- 5    # broken if |beta_hat - beta_true| > 5 * SE(lmer)

## --------------------------------------------------------------
## Designs (canonical fit + truth)
## --------------------------------------------------------------
designs <- list(
    dyestuff = list(
        formula   = Yield ~ 1 + (1 | Batch),
        data      = Dyestuff,
        rho_known = FALSE
    ),
    penicillin = list(
        formula   = diameter ~ 1 + (1 | plate) + (1 | sample),
        data      = Penicillin,
        rho_known = FALSE
    ),
    sleepstudy = list(
        formula   = Reaction ~ Days + (Days | Subject),
        data      = sleepstudy,
        rho_known = TRUE
    )
)

canonical <- lapply(designs, function(d) {
    f <- lme4::lmer(d$formula, data = d$data, REML = TRUE)
    th <- as.numeric(lme4::getME(f, "theta"))
    sigma_e <- sigma(f)
    sigma_b <- if (length(th) == 1L) {
        th * sigma_e
    } else if (length(th) == 2L) {
        th * sigma_e
    } else if (length(th) == 3L) {
        L <- matrix(c(th[1], 0, th[2], th[3]), 2)
        V <- L %*% t(L)
        sigma_e * sqrt(diag(V))
    } else NA_real_
    rho_true <- if (length(th) == 3L) {
        L <- matrix(c(th[1], 0, th[2], th[3]), 2)
        V <- L %*% t(L); V[1, 2] / sqrt(V[1, 1] * V[2, 2])
    } else NA_real_
    list(fit = f,
         beta = unname(lme4::fixef(f)),
         sigma_e = sigma_e,
         sigma_b = sigma_b,
         rho = rho_true,
         se_beta = sqrt(diag(vcov(f))))
})

## --------------------------------------------------------------
## Methods
## --------------------------------------------------------------
##
## All fitters take (formula, data, init_seed) and return a fitted
## model (lmerMod or rlmerMod), or NULL on error.
methods_for <- function(design_name) {
    base <- list(
        lmer            = function(formula, data, init_seed)
            lme4::lmer(formula, data = data, REML = TRUE),
        DAStau          = function(formula, data, init_seed)
            rlmer(formula, data = data, method = "DAStau"),
        DAStau_bisq     = function(formula, data, init_seed) {
            rho.b   <- chgDefaults(smoothPsi, k = if (design_name == "sleepstudy") 5.14 else 1.345)
            rho.s.b <- chgDefaults(smoothPsi, k = if (design_name == "sleepstudy") 5.14 else 2.28)
            rlmer(formula, data = data, method = "DAStau",
                  rho.e = bisquarePsi,
                  rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                  rho.b = rho.b, rho.sigma.b = rho.s.b)
        },
        ransac          = function(formula, data, init_seed)
            rlmer_ransac(formula, data = data, method = "DAStau",
                         K = 50L, seed = init_seed),
        ransac_bisq     = function(formula, data, init_seed) {
            rho.b   <- chgDefaults(smoothPsi, k = if (design_name == "sleepstudy") 5.14 else 1.345)
            rho.s.b <- chgDefaults(smoothPsi, k = if (design_name == "sleepstudy") 5.14 else 2.28)
            rlmer_ransac(formula, data = data, method = "DAStau",
                         K = 50L, seed = init_seed,
                         rho.e = bisquarePsi,
                         rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                         rho.b = rho.b, rho.sigma.b = rho.s.b)
        }
    )
    if (design_name == "sleepstudy") {
        base$DAStau_sizeOBR <- function(formula, data, init_seed)
            rlmer(formula, data = data, method = "DAStau",
                  size_obr = TRUE)
    }
    base
}

## --------------------------------------------------------------
## Parameter extraction and breakdown indicator
## --------------------------------------------------------------
extract_params <- function(fit) {
    if (is.null(fit) || inherits(fit, "try-error")) return(NULL)
    th <- as.numeric(lme4::getME(fit, "theta"))
    sigma_e <- as.numeric(sigma(fit))
    if (length(th) == 3L) {
        L <- matrix(c(th[1], 0, th[2], th[3]), 2)
        V <- L %*% t(L)
        sigma_b <- sigma_e * sqrt(diag(V))
        rho <- V[1, 2] / sqrt(V[1, 1] * V[2, 2])
    } else {
        sigma_b <- th * sigma_e
        rho <- NA_real_
    }
    list(beta = unname(lme4::fixef(fit)),
         sigma_e = sigma_e, sigma_b = sigma_b, rho = rho)
}

broken_indicator <- function(p, truth, rho_known) {
    out <- c(beta = 0L, sigma_e = 0L, sigma_b = 0L, rho = 0L)
    if (is.null(p))
        return(c(any = 1L, beta = 1L, sigma_e = 1L,
                 sigma_b = 1L, rho = 1L))
    if (any(abs(p$beta - truth$beta) > BETA_TOL * truth$se_beta))
        out["beta"] <- 1L
    if (p$sigma_e / truth$sigma_e > SIGMA_TOL ||
        p$sigma_e / truth$sigma_e < 1 / SIGMA_TOL)
        out["sigma_e"] <- 1L
    sb  <- sort(p$sigma_b); sbt <- sort(truth$sigma_b)
    if (any(sb / sbt > SIGMA_TOL | sb / sbt < 1 / SIGMA_TOL))
        out["sigma_b"] <- 1L
    if (rho_known && (is.na(p$rho) || abs(p$rho) > 0.99))
        out["rho"] <- 1L
    c(any = as.integer(max(out)), out)
}

## --------------------------------------------------------------
## Generate one contaminated replicate
## --------------------------------------------------------------
simulate_contaminated <- function(truth_fit, strategy_fn, eps, seed) {
    set.seed(seed)
    sim_y <- as.numeric(simulate(truth_fit, nsim = 1, use.u = FALSE)[, 1])
    data <- truth_fit@frame
    yname <- as.character(formula(truth_fit)[[2]])
    data[[yname]] <- sim_y
    strategy_fn(data, eps, seed = seed + 1L)
}

## --------------------------------------------------------------
## Run one (design, strategy, eps, rep, method) cell
## --------------------------------------------------------------
run_cell <- function(args) {
    design_name <- args$design
    strategy_fn <- args$strategy_fn
    method_fn   <- args$method_fn
    truth       <- args$truth
    eps         <- args$eps
    rep_idx     <- args$rep_idx
    formula     <- args$formula
    seed_base   <- args$seed_base
    data <- simulate_contaminated(truth$fit, strategy_fn, eps,
                                   seed = seed_base)
    fit <- tryCatch(
        suppressMessages(suppressWarnings(
            method_fn(formula, data, init_seed = seed_base + 17L))),
        error = function(e) NULL)
    p <- extract_params(fit)
    list(p = p, success = !is.null(fit))
}

## --------------------------------------------------------------
## Run one design end-to-end
## --------------------------------------------------------------
run_design <- function(design_name, verbose = TRUE) {
    truth <- canonical[[design_name]]
    strategies <- strategiesForDesign(design_name, sigma_e = truth$sigma_e)
    methods <- methods_for(design_name)
    rho_known <- designs[[design_name]]$rho_known
    formula <- designs[[design_name]]$formula

    ## Build the full task grid.
    grid <- expand.grid(
        strategy = names(strategies),
        method   = names(methods),
        eps      = EPS_GRID,
        rep      = seq_len(N_REPS),
        stringsAsFactors = FALSE
    )
    ## Per-cell args.
    cells <- lapply(seq_len(nrow(grid)), function(i) {
        s <- grid$strategy[i]; m <- grid$method[i]
        list(design = design_name,
             strategy_name = s,
             strategy_fn   = strategies[[s]],
             method_name   = m,
             method_fn     = methods[[m]],
             truth = truth, formula = formula,
             eps = grid$eps[i], rep_idx = grid$rep[i],
             seed_base = 1000000L *
                 which(EPS_GRID == grid$eps[i]) +
                 1000L * grid$rep[i] +
                 100L * match(s, names(strategies)) +
                 match(m, names(methods)))
    })

    if (verbose)
        cat(sprintf("\n=== %s: %d cells across %d cores ===\n",
                    design_name, length(cells), ncores))

    t0 <- Sys.time()
    if (ncores > 1L) {
        out <- parallel::mclapply(cells, run_cell, mc.cores = ncores,
                                   mc.preschedule = TRUE)
    } else {
        out <- lapply(cells, run_cell)
    }
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    if (verbose)
        cat(sprintf("  done in %.0fs\n", elapsed))

    ## Compute breakdown indicators per cell.
    grid$broken_any     <- NA_integer_
    grid$broken_beta    <- NA_integer_
    grid$broken_sigma_e <- NA_integer_
    grid$broken_sigma_b <- NA_integer_
    grid$broken_rho     <- NA_integer_
    grid$failed         <- NA_integer_
    for (i in seq_along(out)) {
        bi <- broken_indicator(out[[i]]$p, truth, rho_known)
        grid$broken_any[i]     <- bi["any"]
        grid$broken_beta[i]    <- bi["beta"]
        grid$broken_sigma_e[i] <- bi["sigma_e"]
        grid$broken_sigma_b[i] <- bi["sigma_b"]
        grid$broken_rho[i]     <- bi["rho"]
        grid$failed[i]         <- as.integer(!out[[i]]$success)
    }

    grid$design <- design_name
    grid
}

## --------------------------------------------------------------
## Main: run all designs
## --------------------------------------------------------------
filename <- file.path(path, "datasets_breakdownMC-results.Rdata")
if (file.exists(filename)) {
    load(filename)
    cat(sprintf("Loaded cached results: %d rows from %s\n",
                nrow(breakdownMC_results), basename(filename)))
} else {
    results_per_design <- lapply(names(designs), function(d) run_design(d))
    breakdownMC_results <- do.call(rbind, results_per_design)
    save(breakdownMC_results, file = filename)
    cat(sprintf("Saved results: %d rows to %s\n",
                nrow(breakdownMC_results), basename(filename)))
}

## --------------------------------------------------------------
## Aggregation: breakdown rate per (design, strategy, method, eps)
## --------------------------------------------------------------
breakdownMC_summary <- aggregate(
    cbind(broken_any, broken_beta, broken_sigma_e, broken_sigma_b,
          broken_rho, failed) ~ design + strategy + method + eps,
    data = breakdownMC_results, FUN = mean)
names(breakdownMC_summary)[5:10] <- c("rate_any", "rate_beta",
    "rate_sigma_e", "rate_sigma_b", "rate_rho", "rate_failed")

## --------------------------------------------------------------
## Plots
## --------------------------------------------------------------
plot_breakdownMC_dyestuff <-
    ggplot(subset(breakdownMC_summary, design == "dyestuff"),
           aes(eps, rate_any, color = method)) +
    geom_line() + geom_point() +
    facet_wrap(~ strategy) +
    ylim(0, 1) +
    xlab(expression(epsilon)) + ylab("breakdown rate (any indicator)") +
    ggtitle("Dyestuff (n = 30, 6 batches x 5)")

plot_breakdownMC_penicillin <-
    ggplot(subset(breakdownMC_summary, design == "penicillin"),
           aes(eps, rate_any, color = method)) +
    geom_line() + geom_point() +
    facet_wrap(~ strategy) +
    ylim(0, 1) +
    xlab(expression(epsilon)) + ylab("breakdown rate (any indicator)") +
    ggtitle("Penicillin (n = 144, 24 x 6 crossed)")

plot_breakdownMC_sleepstudy <-
    ggplot(subset(breakdownMC_summary, design == "sleepstudy"),
           aes(eps, rate_any, color = method)) +
    geom_line() + geom_point() +
    facet_wrap(~ strategy) +
    ylim(0, 1) +
    xlab(expression(epsilon)) + ylab("breakdown rate (any indicator)") +
    ggtitle("Sleepstudy (n = 180, 18 subj x 10 days)")

if (interactive()) {
    print(plot_breakdownMC_dyestuff)
    print(plot_breakdownMC_penicillin)
    print(plot_breakdownMC_sleepstudy)
}
