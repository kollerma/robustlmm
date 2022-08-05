require(robustlmm)

asymptoticVariance_generic <- robustlmm:::asymptoticVariance_generic
asymptoticVariance_classic <- robustlmm:::asymptoticVariance_classic
asymptoticVariance_huber_proposal2 <- robustlmm:::asymptoticVariance_huber_proposal2
partialMoment_standardNormal <-
    robustlmm:::partialMoment_standardNormal

## check classic implementation
for (equation in c("location", "scale", "eta", "tau", "mu")) {
    dimension <- if (equation %in% c("location", "scale"))
        1
    else
        2
    stopifnot(all.equal(
        asymptoticVariance_classic(equation, dimension),
        asymptoticVariance_generic(cPsi, equation, dimension)
    ))
}

## check partial moment function
integrate_partialMoment <- function(z, n) {
    int <- integrate(function(x)
        x ^ n * dnorm(x),-Inf, z)
    return(int$value)
}
stopifnot(
    all.equal(
        partialMoment_standardNormal(2, 2),
        integrate_partialMoment(2, 2)
    ),
    all.equal(
        partialMoment_standardNormal(2, 3),
        integrate_partialMoment(2, 3)
    ),
    all.equal(
        partialMoment_standardNormal(2, 4),
        integrate_partialMoment(2, 4)
    ),
    all.equal(
        partialMoment_standardNormal(-2, 2),
        integrate_partialMoment(-2, 2)
    ),
    all.equal(
        partialMoment_standardNormal(-2, 3),
        integrate_partialMoment(-2, 3)
    ),
    all.equal(
        partialMoment_standardNormal(-2, 4),
        integrate_partialMoment(-2, 4)
    )
)

## check Huber Proposal 2 implementation
huberPsiRcppProp2 <- psi2propII(huberPsiRcpp)
for (equation in c("location", "scale", "eta", "tau", "mu")) {
    dimension <- if (equation %in% c("location", "scale"))
        1
    else
        2
    stopifnot(all.equal(
        asymptoticVariance_huber_proposal2(huberPsiRcppProp2, equation, dimension),
        asymptoticVariance_generic(huberPsiRcppProp2, equation, dimension),
        tolerance = 1e-3
    ))
}
