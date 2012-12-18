## test the availability of generic functions.
ae <- function(target, current, ...) {
    ret <- all.equal(target, current, ...)
    if (isTRUE(ret)) return(ret)
    print(ret)
    stop("Objects not equal")
}

require(robustlmm)


set.seed(3)
sleepstudy2 <- within(sleepstudy, {
    Group <- letters[1:4]
    Covar <- rnorm(180)
})
rfm <- rlmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
             rho.e = cPsi, rho.b = cPsi, doFit=FALSE)
fm <- lmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2)    
## all three print the same:
print(rfm)
show(rfm)
summary(rfm)

## object information
## ae(df.residual(fm), df.residual(rfm))
ae(formula(fm), formula(rfm))
ae(model.frame(fm), model.frame(rfm))
ae(model.matrix(fm), model.matrix(rfm), check.attr=FALSE)
nobs(rfm)
## ae(getInitial(fm), getInitial(rfm))
ae(terms(fm), terms(rfm))
weights(rfm)

## basic accessors for the results
coef(rfm)
## dummy.coef(rfm)
deviance(rfm)
ae(fitted(fm), fitted(rfm))
ae(lme4::fixef(fm), fixef(rfm))
ranef(rfm)
ae(resid(fm), resid(rfm))
ae(sigma(fm), sigma(rfm))
## weighted.residuals(rfm)

## var-covar methods
check.attr <- is(fm, "merMod")
## VarCorr(rfm)
ae(lme4::VarCorr(fm), VarCorr(rfm), check.attr=check.attr)
ae(vcov(fm), vcov(rfm), check.attr = check.attr, tolerance=1e-4)
## vcov(rfm)

## confidence intervals
## confint(rfm)

## other (deprecated?)
theta(rfm)

## other methods
getInfo(rfm)
compare(fm, rfm)
update(rfm, ~ . + Covar)
