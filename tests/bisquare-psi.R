## Test bisquarePsi and makeBisquarePsi.
##
## - constructor returns a psi_func_rcpp
## - default c = 4.685 reproduces classical bisquare values
## - integrand expectations are finite and sensible
## - bisquarePsi can be used as rho.e in rlmer

require(robustlmm)

## --- 1. Construction -----------------------------------------------
pf <- makeBisquarePsi(c = 4.685)
stopifnot(methods::is(pf, "psi_func_rcpp"))
stopifnot(identical(pf@getRcppClass(), character(0)))

## --- 2. Bisquare psi values at canonical points ---------------------
x <- c(0, 1, 2, 4, 4.685, 5)
expected_wgt <- c(1, (1 - (1/4.685)^2)^2,
                  (1 - (2/4.685)^2)^2,
                  (1 - (4/4.685)^2)^2,
                  0, 0)
stopifnot(all.equal(pf@wgt(x), expected_wgt, tolerance = 1e-10))
stopifnot(all.equal(pf@psi(x), x * expected_wgt, tolerance = 1e-10))
stopifnot(pf@rho(0) == 0)
stopifnot(pf@rho(100) == 4.685^2 / 6)

## --- 3. Expectations are finite -------------------------------------
stopifnot(is.finite(pf@Erho()))
stopifnot(is.finite(pf@Epsi2()))
stopifnot(is.finite(pf@EDpsi()))
stopifnot(pf@Erho() > 0)

## --- 4. Default bisquarePsi exists and matches ----------------------
stopifnot(methods::is(bisquarePsi, "psi_func_rcpp"))
stopifnot(all.equal(bisquarePsi@wgt(x), expected_wgt, tolerance = 1e-10))

## --- 5. bisquarePsi can be used as rho.e in rlmer -------------------
##     Requires the psi2propII fix for getRcppClass returning character(0).
fit <- rlmer(Reaction ~ Days + (Days | Subject),
             data = sleepstudy,
             method = "DASvar",
             rho.e = bisquarePsi,
             rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
             rho.b = chgDefaults(smoothPsi, k = 5.14),
             rho.sigma.b = chgDefaults(smoothPsi, k = 5.14))
stopifnot(methods::is(fit, "rlmerMod"))
stopifnot(robustlmm:::.sigma(fit) > 0)
stopifnot(all(is.finite(robustlmm:::.fixef(fit))))

## --- 6. Argument validation ----------------------------------------
err <- tryCatch(makeBisquarePsi(c = -1), error = identity)
stopifnot(inherits(err, "error"))
err <- tryCatch(makeBisquarePsi(c = c(1, 2)), error = identity)
stopifnot(inherits(err, "error"))
err <- tryCatch(makeBisquarePsi(c = NA_real_), error = identity)
stopifnot(inherits(err, "error"))

## --- 7. bit-compat: robustbase-backed bisquarePsi reproduces the old
##     hand-coded values to 1e-10. References captured from the HEAD
##     (pre-robustbase-delegation) build.
xg <- c(-6, -3, -1, -0.5, 0, 0.5, 1, 2, 4, 4.685, 5, 6, 0.123, 3.3)
cc <- 4.685; u <- xg / cc
ref_psi  <- ifelse(abs(u) <= 1, xg * (1 - u^2)^2, 0)
ref_wgt  <- ifelse(abs(u) <= 1, (1 - u^2)^2, 0)
ref_Dpsi <- ifelse(abs(u) <= 1, (1 - u^2) * (1 - 5 * u^2), 0)
ref_rho  <- ifelse(abs(u) <= 1, cc^2 / 6 * (1 - (1 - u^2)^3), cc^2 / 6)
ref_Dwgt <- ifelse(abs(u) <= 1, -4 * xg / cc^2 * (1 - u^2), 0)
stopifnot(all.equal(bisquarePsi@psi(xg),  ref_psi,  tolerance = 1e-10),
          all.equal(bisquarePsi@wgt(xg),  ref_wgt,  tolerance = 1e-10),
          all.equal(bisquarePsi@Dpsi(xg), ref_Dpsi, tolerance = 1e-10),
          all.equal(bisquarePsi@rho(xg),  ref_rho,  tolerance = 1e-10),
          all.equal(bisquarePsi@Dwgt(xg), ref_Dwgt, tolerance = 1e-10))
## E-slots (hardcoded references from the HEAD build)
stopifnot(all.equal(bisquarePsi@Erho(),  0.4368496294409207, tolerance = 1e-10),
          all.equal(bisquarePsi@Epsi2(), 0.6044483627183062, tolerance = 1e-10),
          all.equal(bisquarePsi@EDpsi(), 0.7577759185952667, tolerance = 1e-10))

## --- 8. makeRobustbasePsi general constructor + lqqPsi --------------
stopifnot(methods::is(lqqPsi, "psi_func_rcpp"))
## bisquare via the general constructor equals the hand path
bp2 <- makeRobustbasePsi("bisquare", cc = 4.685)
stopifnot(all.equal(bp2@psi(xg), bisquarePsi@psi(xg), tolerance = 1e-12),
          all.equal(bp2@rho(xg), bisquarePsi@rho(xg), tolerance = 1e-12))
## every family: finite E-slots, psi' matches numeric derivative, rho'==psi
for (fam in c("bisquare", "lqq", "optimal", "hampel", "ggw")) {
    pf <- makeRobustbasePsi(fam)
    stopifnot(methods::is(pf, "psi_func_rcpp"),
              is.finite(pf@Erho()), pf@Erho() > 0,
              is.finite(pf@Epsi2()), is.finite(pf@EDpsi()))
    xx <- c(0.3, 1, 2, 4)
    nd <- (pf@psi(xx + 1e-6) - pf@psi(xx - 1e-6)) / 2e-6
    stopifnot(max(abs(pf@Dpsi(xx) - nd)) < 1e-4)   # ggw has kinks -> looser
    nr <- (pf@rho(xx + 1e-6) - pf@rho(xx - 1e-6)) / 2e-6
    stopifnot(max(abs(pf@psi(xx) - nr)) < 1e-3)
}

## --- 9. lqqPsi runs as rho.e in rlmer (redescender; needs a good init;
##     a default (smoothPsi) fit provides one via doFit staging) --------
fit_lqq <- rlmer(Reaction ~ Days + (Days | Subject),
                 data = sleepstudy, method = "DASvar",
                 rho.e = lqqPsi,
                 rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                 rho.b = chgDefaults(smoothPsi, k = 5.14),
                 rho.sigma.b = chgDefaults(smoothPsi, k = 5.14))
stopifnot(methods::is(fit_lqq, "rlmerMod"),
          robustlmm:::.sigma(fit_lqq) > 0,
          all(is.finite(robustlmm:::.fixef(fit_lqq))))

cat("bisquare-psi: all tests passed\n")
