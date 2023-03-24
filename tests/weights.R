require(robustlmm)

checkEquality <- function(cW, rW, tolerance = 1e-8) {
    stopifnot(
        all.equal(
            coef(cW),
            coef(rW),
            tolerance = tolerance,
            check.attributes = FALSE
        ),
        all.equal(
            fixef(cW),
            fixef(rW),
            tolerance = tolerance,
            check.attributes = FALSE
        ),
        all.equal(
            ranef(cW) ,
            ranef(rW),
            tolerance = 100 * tolerance,
            check.attributes = FALSE
        ),
        all.equal(
            fitted(cW) ,
            fitted(rW),
            tolerance = tolerance,
            check.attributes = FALSE
        ),
        all.equal(
            predict(cW) ,
            predict(rW),
            tolerance = tolerance,
            check.attributes = FALSE
        ),
        all.equal(
            coef(summary(cW)) ,
            coef(summary(rW)),
            tolerance = 100 * tolerance,
            check.attributes = FALSE
        )
    )
}

testBattery <- function(formula, data, tolerance) {
    nobs <- nrow(data)

    test <- function(weights) {
        cW <- lmer(formula, data, weights = weights)
        rW <-
            rlmer(
                formula,
                data,
                weights = weights,
                rho.e = cPsi,
                rho.b = cPsi,
                init = lmerNoFit
            )
        checkEquality(cW, rW, tolerance)
    }

    test(rep(2, nobs))
    test(rep(0.5, nobs))

    set.seed(133)
    test(runif(nobs))
    test(rexp(nobs))
}

testBattery(Yield ~ 1 | Batch, Dyestuff, 1e-8)

# Skip these to speed up tests
#testBattery(diameter ~ (1 | plate) + (1 | sample),
#            Penicillin, 1e-6)
#testBattery(Reaction ~ Days + (Days | Subject), sleepstudy, 1e-5)

# testBattery(y ~ service * dept + studage + lectage +
#                 (1 | s) + (1 | d), InstEval)
