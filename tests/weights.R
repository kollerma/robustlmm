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
            ## df = "none": compare the 3-column coefficient table against
            ## lmer's; not the WS16 default Satterthwaite df table.
            coef(summary(rW, df = "none")),
            tolerance = 100 * tolerance,
            check.attributes = FALSE
        )
    )
}

testBattery <- function(formula, data, tolerance) {
    nobs <- nrow(data)

    ## both tau methods: DAStau was always pinned here; DASvar joined
    ## 2026-06-12 (WS24) -- its tau_e double-counted the prior weights
    ## (v_e in the tau^2 formula on top of the U_e whitening), biasing
    ## sigma-hat low by up to 15% at cPsi for weights in [0.25, 4],
    ## which nothing detected because this battery only ran the
    ## default method.
    test <- function(weights) {
        ## lme4 (>= upcoming release) looks up `weights` in the formula's
        ## environment rather than the calling frame; make it visible
        ## there (GitHub issue #36, Ben Bolker).
        assign("weights", weights, environment(formula))
        cW <- lmer(formula, data, weights = weights)
        for (method in c("DAStau", "DASvar")) {
            rW <-
                rlmer(
                    formula,
                    data,
                    weights = weights,
                    rho.e = cPsi,
                    rho.b = cPsi,
                    init = lmerNoFit,
                    method = method
                )
            checkEquality(cW, rW, tolerance)
        }
    }

    test(rep(2, nobs))
    test(rep(0.5, nobs))

    set.seed(133)
    test(runif(nobs))
    test(rexp(nobs))
}

testBattery(Yield ~ 1 | Batch, Dyestuff, 1e-8)
testBattery(diameter ~ (1 | plate) + (1 | sample),
            Penicillin, 1e-6)
testBattery(Reaction ~ Days + (Days | Subject), sleepstudy, 1e-5)

# testBattery(y ~ service * dept + studage + lectage +
#                 (1 | s) + (1 | d), InstEval)
