require(robustlmm)

set.seed(1)
DyestuffWithOffset <- within(Dyestuff, {
                           offset <- rnorm(length(Yield))
                           Yield <- Yield + offset
})

testFormula <- function(formula, data) {
    print(summary(fm <- lmer(formula, data, control=lmerControl(optimizer="bobyqa"))))
    print(summary(rm <- rlmer(formula, data, rho.e = cPsi, rho.b = cPsi, init = lmerNoFit)))
    ranef.fm <- ranef(fm, condVar=FALSE)
    stopifnot(all.equal(coef(fm), coef(rm), tolerance = 1e-3, check.attributes = FALSE),
              all.equal(fixef(fm), fixef(rm), tolerance = 1e-3, check.attributes = FALSE),
              all.equal(ranef.fm , ranef(rm), tolerance = 1e-2, check.attributes = FALSE),
              all.equal(fitted(fm) , fitted(rm), tolerance = 1e-3, check.attributes = FALSE),
              all.equal(predict(fm) , predict(rm), tolerance = 1e-3, check.attributes = FALSE))
    invisible(list(fm, rm))
}

testFormula(Yield ~ offset(offset) + (1 | Batch), DyestuffWithOffset)

