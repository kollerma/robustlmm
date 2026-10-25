## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## makeSScalePsi(): the S-type (bounded-rho) scale psi-function.
require(robustlmm)

ss <- makeSScalePsi()
a <- attr(ss, "sscale")

## 1. Fisher consistency: b = E rho(Z) = 1/2 for the default constant, and
##    rlmer's own consistency constant equals kappa0 = rho^-1(b)^2.
stopifnot(abs(a$b - 0.5) < 1e-6,
          abs(robustlmm:::calcKappaTau(ss, 1) - a$kappa0) < 1e-6,
          abs(robustbase::Mchi(sqrt(a$kappa0), a$c, "bisquare") - a$b) < 1e-10)

## 2. The scale equation's chi(t) = W(t) (t^2 - kappa0) is the bounded
##    rho(t) - b (so the breakdown point is min(b, 1 - b) = 1/2), including
##    at the removable singularity t^2 = kappa0.
tt <- c(0, 0.3, sqrt(a$kappa0) + c(-1e-8, 0, 1e-8), 1, 1.5, 3, 50, 1e6)
chi <- ss@wgt(tt) * (tt^2 - a$kappa0)
stopifnot(max(abs(chi - (robustbase::Mchi(tt, a$c, "bisquare") - a$b))) < 1e-6,
          all(ss@wgt(tt) > 0), all(is.finite(ss@wgt(tt))))

## 3. The analytic weight derivative agrees with finite differences.
x <- c(-2, -0.5, 0.2, sqrt(a$kappa0) + 1e-3, 1.2, 3)
h <- 1e-6
num <- (ss@wgt(x + h) - ss@wgt(x - h)) / (2 * h)
stopifnot(max(abs(ss@Dwgt(x) - num)) < 1e-5)

## 4. rlmer's reweighting iteration for sigma, in the location-scale limit
##    (tau = 1), converges to the M-scale with bounded rho.
set.seed(1)
r <- c(rnorm(160), rnorm(40, 8))
s <- mad(r)
for (i in 1:500) {
    w <- ss@wgt(r / s)
    s1 <- sqrt(sum(w * r^2) / sum(w) / a$kappa0)
    if (abs(s1 - s) < 1e-12) break
    s <- s1
}
mscale <- uniroot(function(s) mean(robustbase::Mchi(r / s, a$c, "bisquare")) - a$b,
                  c(0.1, 20), tol = 1e-12)$root
stopifnot(abs(s1 - mscale) < 1e-8)

## 5. In rlmer: with 20% of the observations shifted by 8 error standard
##    deviations the S-type residual scale stays near 1 where the default
##    scale is inflated, and a bisquare rho.e keeps the outliers rejected.
set.seed(11)
J <- 20; m <- 10
g <- factor(rep(seq_len(J), each = m))
xx <- rnorm(J * m)
y <- 2 + xx + rnorm(J)[g] + rnorm(J * m)
bad <- sample.int(J * m, 40)
y[bad] <- y[bad] + 8
d <- data.frame(g = g, x = xx, y = y)
## (a psi-function built without an Rcpp class needs explicit scale rhos)
f_def <- rlmer(y ~ x + (1 | g), d, method = "DASvar", rho.e = bisquarePsi,
               rho.sigma.e = psi2propII(smoothPsi), rho.b = smoothPsi,
               rho.sigma.b = psi2propII(smoothPsi))
f_S <- rlmer(y ~ x + (1 | g), d, method = "DASvar", rho.e = bisquarePsi,
             rho.sigma.e = sScalePsi, rho.b = smoothPsi,
             rho.sigma.b = psi2propII(smoothPsi))
stopifnot(sigma(f_S) < sigma(f_def),
          abs(fixef(f_S)[1] - 2) < abs(fixef(f_def)[1] - 2))

## 6. A larger constant trades breakdown for efficiency: b = 1/4 at
##    c = 2.937015.
s4 <- makeSScalePsi(2.937015)
stopifnot(abs(attr(s4, "sscale")$b - 0.25) < 1e-5)

## 7. Guards (review M2): the S-type scale is a weight function for the
##    scale equations of the residuals and of one-dimensional blocks only;
##    rlmer refuses it for a correlated block (whose multivariate scale
##    equation has another form) and as rho.e / rho.b, chgDefaults() points
##    to makeSScalePsi(), and NA arguments give NA.
e1 <- tryCatch(rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                     rho.sigma.b = sScalePsi), error = conditionMessage)
e2 <- tryCatch(rlmer(y ~ x + (1 | g), d, rho.e = sScalePsi),
               error = conditionMessage)
e3 <- tryCatch(rlmer(y ~ x + (1 | g), d, rho.b = sScalePsi,
                     rho.sigma.b = sScalePsi), error = conditionMessage)
e4 <- tryCatch(chgDefaults(sScalePsi, 2), error = conditionMessage)
stopifnot(grepl("one-dimensional", e1), grepl("dimension 2", e1),
          grepl("only be used for rho.sigma.e", e2),
          grepl("only be used for rho.sigma.e", e3),
          grepl("makeSScalePsi", e4),
          identical(is.na(ss@wgt(c(1, NA))), c(FALSE, TRUE)),
          identical(is.na(ss@Dwgt(c(1, NA))), c(FALSE, TRUE)))

## 8. With the default method = "DAStau" the S-type residual scale is
##    consistent at the central model (the DAS-tau fixed point uses the same
##    weight function): clean data, 100 groups of 10.
set.seed(3)
J <- 100; m <- 10
g <- factor(rep(seq_len(J), each = m))
xx <- rnorm(J * m)
dd <- data.frame(g = g, x = xx, y = 2 + xx + rnorm(J)[g] + rnorm(J * m))
fT <- rlmer(y ~ x + (1 | g), dd, rho.sigma.e = sScalePsi,
            rho.sigma.b = sScalePsi)
fL <- lme4::lmer(y ~ x + (1 | g), dd)
stopifnot(abs(sigma(fT) / sigma(fL) - 1) < 0.05,
          abs(getME(fT, "theta") / getME(fL, "theta") - 1) < 0.2)

## 9. The setting = "RSEhb" bundle: bisquare rho.e with c = 3.443689,
##    bisquare rho.b and S-type scales for one-dimensional blocks, the
##    multivariate defaults for a correlated block; "RSEh" partially
##    matches; the default start is init = "robust" unless init is given.
fh <- rlmer(y ~ x + (1 | g), d, setting = "RSEhb")
stopifnot(abs(fh@rho.e@tDefs[["cc"]] - 3.443689) < 1e-6,
          grepl("bisquare", tolower(fh@rho.b[[1]]@name)),
          !is.null(attr(fh@rho.sigma.e, "sscale")),
          !is.null(attr(fh@rho.sigma.b[[1]], "sscale")))
fc <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
            setting = "RSEhb")
stopifnot(is.null(attr(fc@rho.sigma.b[[1]], "sscale")),
          !is.null(attr(fc@rho.sigma.e, "sscale")))
fp <- rlmer(y ~ x + (1 | g), d, setting = "RSEh")
stopifnot(isTRUE(all.equal(fixef(fp), fixef(fh))))
fx <- rlmer(y ~ x + (1 | g), d, rho.e = makeBisquarePsi(3.443689),
            rho.b = bisquarePsi, rho.sigma.e = sScalePsi,
            rho.sigma.b = sScalePsi, init = "robust")
fn <- rlmer(y ~ x + (1 | g), d, setting = "RSEhb", init = NULL)
stopifnot(isTRUE(all.equal(fixef(fx), fixef(fh))),
          !isTRUE(all.equal(fixef(fn), fixef(fh))))

cat("sscale-psi.R: all checks passed\n")
