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
