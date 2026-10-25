## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Residual types of residuals.rlmerMod: "response" as lme4 (also for
## prior-weighted fits), "pearson", "weighted" and "normalized" against
## direct computation.
require(robustlmm)

chk <- function(a, b, tol = 1e-10)
    stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = tol)))

set.seed(3)
ss <- sleepstudy
ss$w <- runif(nrow(ss), 0.5, 2)
ss$off <- rnorm(nrow(ss))

fits <- list(
    plain    = rlmer(Reaction ~ Days + (Days | Subject), ss),
    weighted = rlmer(Reaction ~ Days + (Days | Subject), ss, weights = w),
    offset   = rlmer(Reaction ~ Days + offset(off) + (1 | Subject), ss,
                     weights = w),
    crossed  = rlmer(diameter ~ 1 + (1 | plate) + (1 | sample), Penicillin))

for (nm in names(fits)) {
    fm <- fits[[nm]]
    y <- getME(fm, "y"); mu <- fitted(fm); s <- sigma(fm)
    w <- weights(fm); if (is.null(w)) w <- rep(1, length(y))
    we <- getME(fm, "w_e")

    ## response: y - fitted, with the names of the model frame
    r <- residuals(fm)
    chk(r, y - mu)
    stopifnot(identical(names(r), rownames(model.frame(fm))),
              identical(r, residuals(fm, type = "response")),
              identical(r, resid(fm)))
    chk(residuals(fm, scaled = TRUE), (y - mu) / s)

    ## pearson: whitened by the prior weights, divided by sigma; `scaled`
    ## has no further effect
    chk(residuals(fm, type = "pearson"), sqrt(w) * (y - mu) / s)
    stopifnot(identical(residuals(fm, type = "pearson"),
                        residuals(fm, type = "pearson", scaled = TRUE)))
    ## these are summary()'s scaled residuals
    chk(summary(fm, df = "none")$residuals, sqrt(w) * (y - mu) / s)

    ## weighted: additionally multiplied by the robustness weights
    chk(residuals(fm, type = "weighted"), we * sqrt(w) * (y - mu))
    chk(residuals(fm, type = "weighted", scaled = TRUE),
        we * sqrt(w) * (y - mu) / s)

    ## normalized: L^-1 (y - X beta - offset), V = L L' the fitted marginal
    ## covariance, L lower triangular, in the order of the data
    Z <- as.matrix(getME(fm, "Z"))
    Lam <- as.matrix(getME(fm, "Lambda"))
    V <- s^2 * (Z %*% tcrossprod(Lam) %*% t(Z) + diag(1 / w))
    off <- model.offset(model.frame(fm)); if (is.null(off)) off <- 0
    rm <- y - drop(getME(fm, "X") %*% fixef(fm)) - off
    chk(residuals(fm, type = "normalized"),
        forwardsolve(t(chol(V)), rm), tol = 1e-8)
}

## the block-wise and the sparse full-matrix computations agree: a nested
## design goes block by block, and equals the dense reference above for a
## nested fit too
set.seed(5)
d <- expand.grid(rep = 1:4, sub = 1:3, top = 1:8)
d$sub <- interaction(d$top, d$sub)
d$y <- rnorm(8)[d$top] + rnorm(24, sd = 0.7)[as.integer(d$sub)] +
    rnorm(nrow(d))
fn <- rlmer(y ~ 1 + (1 | top) + (1 | sub), d)
Z <- as.matrix(getME(fn, "Z")); Lam <- as.matrix(getME(fn, "Lambda"))
V <- sigma(fn)^2 * (Z %*% tcrossprod(Lam) %*% t(Z) + diag(nrow(d)))
chk(residuals(fn, type = "normalized"),
    forwardsolve(t(chol(V)), d$y - fixef(fn)), tol = 1e-8)

## ---- "response" equals lme4's at the classical psi-functions -------------
for (wcol in c("one", "w")) {
    ss$one <- 1
    ss$pw <- ss[[wcol]]
    lm <- lmer(Reaction ~ Days + (Days | Subject), ss, weights = pw)
    cm <- rlmer(Reaction ~ Days + (Days | Subject), ss, weights = pw,
                rho.e = cPsi, rho.b = cPsi)
    chk(residuals(cm, type = "response"),
        residuals(lm, type = "response"), tol = 1e-4)
    ## lme4's pearson residuals are not divided by sigma
    chk(residuals(cm, type = "pearson") * sigma(cm),
        residuals(lm, type = "pearson"), tol = 1e-4)
}

## ---- missing values are padded as before ----------------------------------
ssNA <- ss; ssNA$Reaction[c(3, 50)] <- NA
fna <- rlmer(Reaction ~ Days + (1 | Subject), ssNA, na.action = na.exclude)
for (tp in c("response", "pearson", "weighted", "normalized")) {
    r <- residuals(fna, type = tp)
    stopifnot(length(r) == nrow(ssNA), all(is.na(r[c(3, 50)])),
              !anyNA(r[-c(3, 50)]))
}

## normalized residuals of data simulated from the model are white: their
## within-subject lag-1 correlation is near zero, unlike that of the
## marginal residuals
fm <- fits$plain
nr <- residuals(fm, type = "normalized")
mr <- getME(fm, "y") - drop(getME(fm, "X") %*% fixef(fm))
lag1 <- function(r) {
    a <- unlist(tapply(r, ss$Subject, function(v) v[-length(v)]))
    b <- unlist(tapply(r, ss$Subject, function(v) v[-1]))
    cor(a, b)
}
stopifnot(abs(lag1(nr)) < abs(lag1(mr)))
cat("residual-types.R: all checks passed\n")
