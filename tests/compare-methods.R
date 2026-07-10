## Short results-constancy check for the CRAN build: a representative
## subset of the full battery (which runs on master and in CI) to keep
## the CRAN check time within limits. Keeps a classic random-intercept
## fit (DASvar + DAStau, compared against lmer) and a robust correlated-
## random-effects fit.
require(robustlmm)

fit <- function(formula, data, methods =  c("DASvar", "DAStau"),
                rho.e = cPsi, rho.b = cPsi, ...) {
    fits <- list()
    ## compare with result of lmer if rho arguments are not given
    classic <- ! any(c("rho.e", "rho.b") %in% names(match.call())[-1])
    if (classic) fm <- lmer(formula, data, control=lmerControl(optimizer="bobyqa"))
    for (method in methods) {
        fits[[method]] <- list()
        if (classic) fits[[method]][["lmer"]] <- fm
        cat("\n########", method, "########\n")
        try({cat("Time elapsed:",
                 system.time(m <- rlmer(formula, data, method=method,
                                        rho.e = rho.e, rho.b = rho.b, ...)),
                 "\n")
             fits[[method]][["IRWLS"]] <- m
             ## df = "none": this is an estimation-equality regression
             ## (constant golden output); keep the historic 3-column
             ## table, not the WS16 default Satterthwaite df.
             print(summary(m, df = "none"))
             print(robustlmm:::u.rlmerMod(m), 4)
             if (classic) {
                 ## compare with lmer fit
                 cat("#### Checking equality with lmer... ####\n")
                 cat("Fixed effects: ", all.equal(fixef(fm), fixef(m), tolerance = 1e-4), "\n")
                 ranef.fm <- ranef(fm, condVar=FALSE)# lme4 now has default  condVar=TRUE
                 cat("Random effects:", all.equal(ranef.fm, ranef(m), tolerance = 1e-4,
                                                  check.attributes=FALSE), "\n")
                 cat("Theta:         ", all.equal(theta(fm), theta(m), tolerance = 1e-4), "\n")
                 cat("Sigma:         ", all.equal(sigma(fm), sigma(m), tolerance = 1e-4), "\n")
                 if (packageVersion("lme4") >= "0.99999911.0") {
                     tmp <- all.equal(fm@pp$unsc(), unname(m@pp$unsc()), tolerance = 1e-4)
                     if (!isTRUE(tmp))
                         cat("Unsc:          ", tmp , "\n")
                 }
             }
         })
        fits[[method]][["dnames"]] <- names(fits[[method]])
    }
    cat("\n################################################\n")
    cat("################################################\n")
    cat("################################################\n")
    for (method in methods) {
        cat("\n################ results for",method," ##############\n")
        cmp <- do.call(compare, fits[[method]])
        cmp <- cmp[grep("^rho", rownames(cmp), invert=TRUE),,drop=FALSE]
        print.default(cmp, quote="FALSE")
    }
}

Dyestuff$Yield <- Dyestuff$Yield - 0.5

## classic (REML): simple random intercept, both methods vs lmer
fit(Yield ~ (1 | Batch), Dyestuff)

## robust, correlated random effects
fit(Reaction ~ Days + (Days|Subject), sleepstudy,
    rho.e = smoothPsi, rho.b = smoothPsi,
    methods = c("DASvar"))

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
