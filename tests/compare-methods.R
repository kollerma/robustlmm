## test that should make sure that the results remain constant
require(robustlmm)

fit <- function(formula, data, methods =  c("DASexp", "DAStau"),
                methods.effects = c("IRWLS"), rho.e = cPsi, rho.b = cPsi,
                wExp.b = 0, wExp.e = 0, ...) {
    fits <- list()
    ## compare with result of lmer if rho arguments are not given
    classic <- ! any(c("rho.e", "rho.b") %in% names(match.call())[-1])
    if (classic) fm <- lmer(formula, data, ...)
    for (method in methods) {
        fits[[method]] <- list()
        if (classic) fits[[method]][["lmer"]] <- fm
        for (method.effects in methods.effects) {
            cat("\n########", method, method.effects, "########\n")
            try({cat("Time elapsed:",
                     system.time(m <- rlmer(formula, data, method=method,
                                            method.effects=method.effects,
                                            rho.e = rho.e, rho.b = rho.b,
                                            wExp.b = wExp.b, wExp.e = wExp.e, ...)),
                     "\n")
                 fits[[method]][[method.effects]] <- m
                 print(m)
                 print(robustlmm:::u.rlmerMod(m), 4)
                 if (classic) {
                     ## compare with lmer fit
                     cat("#### Checking equality with lmer... ####\n")
                     cat("Fixed effects: ", all.equal(fixef(fm), fixef(m), tol = 1e-4), "\n")
                     cat("Random effects:", all.equal(ranef(fm), ranef(m), tol = 1e-4), "\n")
                     cat("Theta:         ", all.equal(theta(fm), theta(m), tol = 1e-4), "\n")
                     cat("Sigma:         ", all.equal(sigma(fm), sigma(m), 1e-4), "\n")
                     cat("Deviance:      ", all.equal(deviance(m),
                            unname(fm@devcomp$cmp[ifelse(isREML(m), "REML", "dev")])), "\n")
                     cat("Unsc:          ", all.equal(fm@pp$unsc(),
                                                      unname(m@pp$unsc()), tol = 1e-4), "\n")
                 }
             })
        }
        fits[[method]][["dnames"]] <- names(fits[[method]])
    }
    cat("\n################################################\n")
    cat("################################################\n")
    cat("################################################\n")
    for (method in methods) {
        cat("\n################ results for",method," ##############\n")
        cmp <- do.call(compare, fits[[method]])
        cmp <- cmp[!rownames(cmp) %in% c("rho.e", "wExp.e", "rho.b", "wExp.b"), ]
        print.default(cmp, quote="FALSE")
    }
}

## classic (REML)
fit(Yield ~ (1 | Batch), Dyestuff)
fit(Yield ~ (1 | Batch), Dyestuff2)
fit(diameter ~ (1|plate) + (1|sample), Penicillin)

## classic (ML)
fit(Yield ~ (1 | Batch), Dyestuff, REML=FALSE, methods = c("Opt"))
fit(Yield ~ (1 | Batch), Dyestuff2, REML=FALSE, methods = c("Opt"))
fit(diameter ~ (1|plate) + (1|sample), Penicillin, REML=FALSE, methods = c("Opt"))

## classic (no init)
fit(Yield ~ (1 | Batch), Dyestuff,
    methods.effects = c("IRWLS", "Rcgmin"), init = lmerNoFit)
fit(Yield ~ (1 | Batch), Dyestuff2,
    methods.effects = c("IRWLS", "Rcgmin"), init = lmerNoFit)
fit(diameter ~ (1|plate) + (1|sample), Penicillin,
    methods.effects = c("IRWLS", "Rcgmin"), init = lmerNoFit)

## smoothPsi
fit(Yield ~ (1 | Batch), Dyestuff, rho.e = smoothPsi, rho.b = smoothPsi)
fit(Yield ~ (1 | Batch), Dyestuff2, rho.e = smoothPsi, rho.b = smoothPsi)
fit(diameter ~ (1|plate) + (1|sample), Penicillin, rho.e = smoothPsi, rho.b = smoothPsi)

## smoothPsi Proposal II for estimating sigma
fit(Yield ~ (1 | Batch), Dyestuff, rho.e = smoothPsi, rho.b = smoothPsi,
    wExp.e = 2, wExp.b = 0)
fit(Yield ~ (1 | Batch), Dyestuff2, rho.e = smoothPsi, rho.b = smoothPsi,
    wExp.e = 2, wExp.b = 0)
fit(diameter ~ (1|plate) + (1|sample), Penicillin, rho.e = smoothPsi, rho.b = smoothPsi,
    wExp.e = 2, wExp.b = 0)

## smoothPsi Proposal II
fit(Yield ~ (1 | Batch), Dyestuff, rho.e = smoothPsi, rho.b = smoothPsi,
    wExp.e = 2, wExp.b = 2)
fit(Yield ~ (1 | Batch), Dyestuff2, rho.e = smoothPsi, rho.b = smoothPsi,
    wExp.e = 2, wExp.b = 2)
fit(diameter ~ (1|plate) + (1|sample), Penicillin, rho.e = smoothPsi, rho.b = smoothPsi,
    wExp.e = 2, wExp.b = 2)

## correlated random effects
## omitting tests to speed up R CMD check...
## fit(Reaction ~ Days + (Days|Subject), sleepstudy,
##     wExp.e = 2, wExp.b = 2, methods = c("DASexp"), methods.effects = c("IRWLS"))
## fit(Reaction ~ Days + (Days|Subject), sleepstudy,
##     wExp.e = 2, wExp.b = 2, methods = c("DASexp"), methods.effects = c("IRWLS"),
##     init = lmerNoFit)
## including a 0 variance compontent
## sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])
## fit(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
##     wExp.e = 2, wExp.b = 2, methods = c("DASexp"), methods.effects = c("IRWLS"))
## fit(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
##     wExp.e = 2, wExp.b = 2, methods = c("DASexp"), methods.effects = c("IRWLS"),
##     init = lmerNoFit)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
