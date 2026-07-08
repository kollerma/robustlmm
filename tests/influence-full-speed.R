## Disable test
quit()

## github-transition: SLOW TEST (~45s; a performance regression guard,
## least essential on CRAN) -- when merging to the github-transition
## (CRAN-release) branch, disable it there by adding `quit()` at the top,
## per that branch's convention.
## WS22 regression test: the fast implicitIF_full (analytic d score/dy;
## tau_e and tau_b reused on the non-theta parameter-Jacobian columns,
## since both are theta-only) must equal the original all-finite-
## difference computation. Correctness, not timing, is asserted here
## (the speedup -- ~25x at n=240 -- is documented in PLAN-WS22).

suppressMessages(require(robustlmm))

## reference: the pre-WS22 implementation -- full 2n-column FD for Jy and
## a numDeriv parameter Jacobian with NO tau reuse (every score eval
## re-runs the DAStau (tau | s) alternation).
.oldIF_par <- function(fit, eps = 1e-6) {
    pp <- fit@pp; p <- pp$p; q <- pp$q; n <- pp$n
    L  <- length(getME(fit, "theta"))
    par0 <- c(robustlmm:::.fixef(fit), pp$b.s,
              log(robustlmm:::.sigma(fit)), getME(fit, "theta"))
    y0 <- fit@resp$y
    Jpar <- numDeriv::jacobian(function(par) robustlmm:::.scoreVec(par, fit, y0),
                               par0, method = "Richardson",
                               method.args = list(eps = eps, d = 1e-4, r = 4))
    Jy <- matrix(NA_real_, p + q + 1L + L, n); h <- 1e-4
    for (i in seq_len(n)) {
        yp <- y0; yp[i] <- yp[i] + h; ym <- y0; ym[i] <- ym[i] - h
        Jy[, i] <- (robustlmm:::.scoreVec(par0, fit, yp) -
                    robustlmm:::.scoreVec(par0, fit, ym)) / (2 * h)
    }
    -solve(Jpar, Jy)
}

cmp <- function(fit, label) {
    IFn <- robustlmm:::implicitIF_full(fit)
    sig <- robustlmm:::.sigma(fit)
    new_par <- rbind(IFn$IF_beta, IFn$IF_u, IFn$IF_sigma / sig, IFn$IF_theta)
    ref <- .oldIF_par(fit)
    rel <- max(abs(new_par - ref)) / max(abs(ref))
    cat(sprintf("%-28s fast vs full-FD IF: max rel err %.2e\n", label, rel))
    stopifnot(rel < 1e-5)
}

## a well-conditioned design (J large enough, RE sd away from 0) so the
## reference's plain solve() of Jpar is not at a variance-component
## boundary.
set.seed(8)
J <- 20L; m <- 6L; n <- J * m
g <- factor(rep(seq_len(J), each = m)); x <- rnorm(n)
y <- 1 + 0.5 * x + rnorm(J, 0, 1.5)[as.integer(g)] + rnorm(n)
d <- data.frame(y, x, g)

cmp(suppressWarnings(rlmer(y ~ x + (1 | g), d, method = "DAStau")),  "DAStau")
cmp(suppressWarnings(rlmer(y ~ x + (1 | g), d, method = "DASvar")),  "DASvar")
## weighted (exercises U_e != I in d r_std / d y)
set.seed(9); w <- runif(n, 0.5, 3)
cmp(suppressWarnings(rlmer(y ~ x + (1 | g), d, weights = w, method = "DAStau")),
    "DAStau + prior weights")
## Mallows design weights (exercises eta on the e-side rows)
eta <- pmin(1, 1.3 / abs(x)); eta[!is.finite(eta)] <- 1
cmp(suppressWarnings(rlmer(y ~ x + (1 | g), d, design.weights = eta,
                           method = "DAStau")), "DAStau + design weights")

cat("influence-full-speed.R: all checks passed\n")
