## test kappa functions
require(robustlmm)

calcKappaTau <- robustlmm:::calcKappaTau
calcTau <- robustlmm:::calcTau

E <- function(fun) {
    ## sum(na.omit(fun(rfm@pp$ghz)*rfm@pp$ghw))
    int <- function(x) fun(x)*dnorm(x)
    integrate(int, -Inf, Inf)$value
}

testKappaTau <- function(rho, wExp) {
    rho2 <- if (wExp == 2) psi2propII(rho) else rho
    kappa <- calcKappaTau(rho2, 1)
    fun <- switch(wExp,
                  function(x) rho@psi(x)*x - kappa*rho@wgt(x),
                  function(x) rho@psi(x)^2 - kappa*rho@wgt(x)^2)
    E(fun)
}

testKappa <- function(rho) {
    stopifnot(all.equal(0, testKappaTau(rho, 1), 1e-4),
              all.equal(0, testKappaTau(rho, 2), 1e-4))
}

testKappa(huberPsi)
testKappa(smoothPsi)

E2 <- function(fun, ...) {
    int1 <- function(x, y, ...) fun(x, y, ...)*dnorm(x)*dnorm(y)
    int2 <- Vectorize(function(y, ...) integrate(int1, -Inf, Inf, y=y, ...)$value)
    integrate(int2, -Inf, Inf, ...)$value    
}

rfm <- rlmer(Yield ~ (1|Batch), Dyestuff)

testTau <- function(rho, rho.sigma, wExp, a, s) {
    psi <- rho@psi
    rho.e <- rho.sigma
    rho.sigma2 <- if (wExp == 2) psi2propII(rho.sigma) else rho.sigma
    kappa <- calcKappaTau(rho.sigma2, 1)
    i <- 1
    
    fun <- switch(wExp + 1,
                  { ## wExp.e == 0:
                      wgt <- function(x) ifelse(x == 0, rho.e@Dpsi(0)/2, rho@rho(x)/(x*x))
                      function(w, v, tau) {
                          t <- (v-a[i]*psi(v)+w*s[i])/tau
                          rho.e@rho(t) - kappa*wgt(t)
                      } }, ## wExp.e == 1:
                  function(w, v, tau) {
                      t <- (v-a[i]*psi(v)+w*s[i])/tau
                      rho.e@psi(t)*t - kappa*rho.e@wgt(t)
                  }, ## wExp.e == 2:
                  function(w, v, tau) {
                      t <- (v-a[i]*psi(v)+w*s[i])/tau
                      rho.e@psi(t)^2 - kappa*rho.e@wgt(t)^2
                  })

    tau <- calcTau(a, s, rho, rho.sigma2, rfm@pp, kappa)
    print(tau)
    
    ret <- tau
    for (i in seq_along(a)) ret[i] <- E2(fun, tau = tau[i])

    ret
}

a <- c(0.1, 0.4, 1, 1, 0.8)
s <- c(0.4, 0.5, 0.8, 0.1, 0.1)
## FIXME: increase accuracy of calcTau.
stopifnot(all.equal(rep(0, length(a)), testTau(huberPsi, huberPsi, 1, a, s), 1e-2),
          all.equal(rep(0, length(a)), testTau(huberPsi, huberPsi, 2, a, s), 1e-2))
stopifnot(all.equal(rep(0, length(a)), testTau(smoothPsi, huberPsi, 1, a, s), 1e-2),
          all.equal(rep(0, length(a)), testTau(smoothPsi, huberPsi, 2, a, s), 1e-1))
stopifnot(all.equal(rep(0, length(a)), testTau(smoothPsi, smoothPsi, 1, a, s), 1e-2),
          all.equal(rep(0, length(a)), testTau(smoothPsi, smoothPsi, 2, a, s), 1e-1))
sPsi <- chgDefaults(smoothPsi, k = 1, s = 10)
stopifnot(all.equal(rep(0, length(a)), testTau(smoothPsi, sPsi, 1, a, s), 1e-2),
          all.equal(rep(0, length(a)), testTau(smoothPsi, sPsi, 2, a, s), 1e-2))

