## test .calcE.D.re and .calcEchi
require(robustlmm)

.calcE.D.re <- robustlmm:::.calcE.D.re
.calcEchi <- robustlmm:::.calcEchi
cPsi <- robustlmm:::cPsi

## simple tests for q = 1:
test1 <- function(rho) {
    stopifnot(all.equal(.calcE.D.re(1, rho), rho@EDpsi()),
              all.equal(.calcEchi(1, rho), rho@Epsi2()))
}

test1(cPsi)
test1(smoothPsi)

## check very large q
stopifnot(all.equal(.calcE.D.re(10000, huberPsi), 0),
          all.equal(.calcE.D.re(10000, smoothPsi), 0))

## tests fro q > 1:
test2 <- function(q, rho, npts = 1000000) {
    cat("\nq", q, "\n")
    ## test .calcE.D.re
    integrand <- function(u2, v2, s2) {
        dk <- sqrt(u2 + v2)
        (rho@Dpsi(dk)/dk - rho@psi(dk)/dk^2)*u2/dk + rho@psi(sqrt(s2))/sqrt(s2)
    }
    u2 <- rnorm(npts)^2
    v2 <- rchisq(npts, q-1)
    s2 <- rchisq(npts, q)
    value <- mean(integrand(u2, v2, s2))
    cat("\nImportance Sampling:", value, "\n")
    cat(".calcE.D.re:        ", .calcE.D.re(q, rho), "\n")
    #print(all.equal(value, .calcE.D.re(q, rho)))
    stopifnot(all.equal(value, .calcE.D.re(q, rho), tolerance = 1e-2))

    ## test .calcEchi
    integrand <- function(u2, v2) {
        dk <- sqrt(u2 + v2)
        rho@wgt(dk)^2*u2
    }
    value <- mean(integrand(u2, v2))
    cat("\nImportance Sampling:", value, "\n")
    cat(".calcEchi:          ", .calcEchi(q, rho), "\n")
    stopifnot(all.equal(value, .calcEchi(q, rho), tolerance = 1e-2))
      
}

for (q in 1:10) test2(q, cPsi)
for (q in 1:10) test2(q, smoothPsi)



cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
