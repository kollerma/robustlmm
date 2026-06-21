## Test: rlmer accepts user-supplied psi_func objects whose
## getRcppClass returns an empty vector (i.e. constructed via
## psiFunc + a no-op Rcpp slot rather than the built-in Rcpp
## factories).
##
## Background. psi2propII used to dispatch on
##   object@getRcppClass()[1] == "..."
## which raises "missing value where TRUE/FALSE needed" when
## getRcppClass returns character(0) (since [1] is NA). The fix
## guards the length and NA case, and emits a clear error
## directing users to supply rho.sigma.e explicitly when
## providing a custom rho without a backing Rcpp class.

require(robustlmm)

## A Tukey-bisquare psi_func_rcpp constructed via R-level closures.
make_bisquare <- function(c = 4.685) {
    psi_fn <- function(x, c) {
        if (!is.finite(c) || c <= 0) return(x)
        u <- x / c
        ifelse(abs(u) <= 1, x * (1 - u^2)^2, 0)
    }
    wgt_fn <- function(x, c) {
        if (!is.finite(c) || c <= 0) return(rep(1, length(x)))
        u <- x / c
        ifelse(abs(u) <= 1, (1 - u^2)^2, 0)
    }
    Dpsi_fn <- function(x, c) {
        if (!is.finite(c) || c <= 0) return(rep(1, length(x)))
        u <- x / c
        ifelse(abs(u) <= 1, (1 - u^2) * (1 - 5*u^2), 0)
    }
    Dwgt_fn <- function(x, c) {
        if (!is.finite(c) || c <= 0) return(rep(0, length(x)))
        u <- x / c
        ifelse(abs(u) <= 1, -4*x/c^2 * (1 - u^2), 0)
    }
    rho_fn <- function(x, c) {
        if (!is.finite(c) || c <= 0) return(x^2 / 2)
        u <- x / c
        ifelse(abs(u) <= 1, c^2/6 * (1 - (1 - u^2)^3), c^2/6)
    }
    U <- qnorm(1 - 1e-10)
    Erho_fn <- function(c) vapply(c, function(cc)
        2 * integrate(function(z) rho_fn(z, cc)*dnorm(z), 0, U)$value,
        numeric(1))
    Epsi2_fn <- function(c) vapply(c, function(cc)
        2 * integrate(function(z) psi_fn(z, cc)^2*dnorm(z), 0, U)$value,
        numeric(1))
    EDpsi_fn <- function(c) vapply(c, function(cc)
        2 * integrate(function(z) Dpsi_fn(z, cc)*dnorm(z), 0, U)$value,
        numeric(1))
    pf <- robustbase::psiFunc(rho = rho_fn, psi = psi_fn, wgt = wgt_fn,
                              Dpsi = Dpsi_fn, Dwgt = Dwgt_fn,
                              Erho = Erho_fn, Epsi2 = Epsi2_fn,
                              EDpsi = EDpsi_fn,
                              name = sprintf("bisquare (c=%.4g)", c),
                              c = c)
    args <- list("psi_func_rcpp")
    for (sn in slotNames("psi_func")) args[[sn]] <- slot(pf, sn)
    args[["getRcppClass"]] <- function() character(0)
    args[["getInstanceWithOriginalDefaults"]] <- function() NULL
    do.call("new", args)
}

bisq <- make_bisquare(4.685)

## --- 1. psi2propII rejects custom rho cleanly (no NaN crash) ---
err <- tryCatch(psi2propII(bisq), error = identity)
stopifnot(inherits(err, "error"))
## Error message should mention rho.sigma.e (the guidance for users).
stopifnot(grepl("rho.sigma.e", conditionMessage(err)))

## --- 2. rlmer with custom rho.e + explicit rho.sigma.e works ----
fit <- rlmer(Reaction ~ Days + (Days | Subject),
             data = sleepstudy,
             method = "DASvar",
             rho.e = bisq,
             rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
             rho.b       = chgDefaults(smoothPsi, k = 5.14, s = 10),
             rho.sigma.b = chgDefaults(smoothPsi, k = 5.14, s = 10))
stopifnot(robustlmm:::.sigma(fit) > 0)
b <- robustlmm:::.fixef(fit)
stopifnot(length(b) == 2)
stopifnot(is.finite(b[1]), is.finite(b[2]))
th <- getME(fit, "theta")
stopifnot(all(is.finite(th)))
stopifnot(th[1] > 0)  # intercept SD positive

cat("All custom-redescending-psi tests passed.\n")
