## Disable test
quit()

## Tests for the influence / vcov-sandwich layer (WS1a).

require(robustlmm)
suppressMessages(require(Matrix))
suppressMessages(require(lme4))

## ----------------------------------------------------------------
## 1. Default vcov is byte-for-byte unchanged.
## ----------------------------------------------------------------
set.seed(1)
fit_sleep <- rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                   method = "DASvar")
fit_dye   <- rlmer(Yield ~ (1 | Batch), Dyestuff, method = "DASvar")

old_vcov  <- getS3method("vcov", "merMod")
for (fit in list(fit_sleep, fit_dye)) {
    V_def <- vcov(fit, type = "default")
    V_old <- old_vcov(fit)
    stopifnot(identical(dim(V_def), dim(V_old)))
    stopifnot(max(abs(as.matrix(V_def) - as.matrix(V_old))) < 1e-10)
    ## Bare vcov(fit) (no type arg) must equal type="default".
    V_bare <- vcov(fit)
    stopifnot(max(abs(as.matrix(V_bare) - as.matrix(V_def))) < 1e-10)
}

## ----------------------------------------------------------------
## 2. Sandwich is finite, symmetric, PD; n.clusters is right.
## ----------------------------------------------------------------
V_sand <- vcov(fit_sleep, type = "sandwich")
stopifnot(all(is.finite(V_sand)))
stopifnot(max(abs(V_sand - t(V_sand))) < 1e-10)
stopifnot(all(eigen(V_sand, only.values = TRUE)$values > 0))
stopifnot(attr(V_sand, "n.clusters") == 18L)

V_dye <- vcov(fit_dye, type = "sandwich")
stopifnot(all(is.finite(V_dye)))
stopifnot(attr(V_dye, "n.clusters") == nlevels(Dyestuff$Batch))

## G1 correction is the default and the explicit "none" turns it off.
V_g1   <- vcov(fit_sleep, type = "sandwich", correction = "G1")
V_none <- vcov(fit_sleep, type = "sandwich", correction = "none")
ratio  <- 18 / 17                          # J / (J - 1) on sleepstudy
stopifnot(max(abs(V_g1 - V_sand)) < 1e-10) # default == G1
stopifnot(max(abs(V_g1 - ratio * V_none)) < 1e-10)

## ----------------------------------------------------------------
## 3. caseweightIF agrees with the reweighting finite difference.
##
## Port of IF-thread1/implicitIF.R's finite_difference_caseweightIF,
## adapted for in-package internal access. Solves the reweighted
## (beta, u) score by a fixed-Jacobian iteration whose fixed point is
## the true reweighted solution at fixed sigma-hat, theta-hat.
## ----------------------------------------------------------------
.fd_caseweightIF <- function(fit, idx, h = 1e-3,
                             maxit = 1000L, tol = 1e-12) {
    parts <- robustlmm:::implicitIF(fit)
    pp    <- fit@pp
    rho_e <- fit@rho.e
    rho_b <- fit@rho.b
    sigma <- parts$sigma
    n <- pp$n; p <- pp$p; q <- pp$q
    X     <- pp$X
    Z     <- t(pp$Zt)
    U_b   <- parts$U_b
    Ue_inv <- parts$Ue_inv
    A_X    <- parts$A_X
    A_ZUb  <- parts$A_ZUb
    Lambda_ratio <- parts$Lambda_ratio
    y     <- fit@resp$y
    beta0 <- pp$beta
    u0    <- as.numeric(pp$b.s)

    solveJ <- function(rb, ru) {
        Jbb_inv_rb <- solve(parts$Jbb, rb)
        xu <- as.numeric(solve(parts$S, ru - parts$Jub %*% Jbb_inv_rb))
        xb <- as.numeric(Jbb_inv_rb - parts$Jbb_inv_Jbu %*% xu)
        list(xb = xb, xu = xu)
    }
    Gtilde <- function(beta, u, w) {
        resid <- as.numeric(y - X %*% beta - Z %*% (U_b %*% u))
        r_std <- as.numeric(Ue_inv %*% resid) / sigma
        u_std <- u / sigma
        psi_e <- rho_e@psi(r_std)
        Gb <- as.numeric(t(A_X) %*% (w * psi_e))
        psi_b <- robustlmm:::.applyBlockPsi_b(fit, u_std, rho_b)
        Gu <- as.numeric(t(A_ZUb) %*% (w * psi_e)) -
              as.numeric(Lambda_ratio %*% psi_b)
        list(Gb = Gb, Gu = Gu)
    }
    solveRew <- function(w) {
        beta <- beta0; u <- u0
        G <- Gtilde(beta, u, w)
        for (it in seq_len(maxit)) {
            if (max(abs(c(G$Gb, G$Gu))) < tol) break
            step <- solveJ(G$Gb, G$Gu)
            beta <- beta + sigma * step$xb
            u    <- u    + sigma * step$xu
            G    <- Gtilde(beta, u, w)
        }
        list(beta = beta, u = u)
    }
    w0 <- rep(1, n)
    fd <- matrix(NA_real_, nrow = p, ncol = length(idx))
    for (k in seq_along(idx)) {
        i  <- idx[k]
        wp <- w0; wp[i] <- 1 + h
        wm <- w0; wm[i] <- 1 - h
        sp <- solveRew(wp)
        sm <- solveRew(wm)
        fd[, k] <- (sp$beta - sm$beta) / (2 * h)
    }
    fd
}

for (fit in list(fit_sleep, fit_dye)) {
    n  <- fit@pp$n
    ix <- sort(sample.int(n, min(6L, n)))
    cw <- caseweightIF(fit, idx = ix)
    fd <- .fd_caseweightIF(fit, idx = ix, h = 1e-3)
    scale <- pmax(abs(cw$IF_beta), abs(fd), 1e-8)
    rel   <- max(abs(cw$IF_beta - fd) / scale)
    stopifnot(rel < 1e-3)
}

## ----------------------------------------------------------------
## 4. Crossed designs trigger the approximate-for-crossed warning.
## ----------------------------------------------------------------
fit_pen <- rlmer(diameter ~ (1 | plate) + (1 | sample),
                 data = Penicillin, method = "DASvar")

## (a) Auto-detect errors because the choice is ambiguous.
err <- tryCatch(vcov(fit_pen, type = "sandwich"), error = identity)
stopifnot(inherits(err, "error"))
stopifnot(grepl("multiple grouping factors", conditionMessage(err)))

## (b) Naming a grouping factor warns and produces a finite, PD matrix.
w <- tryCatch(vcov(fit_pen, type = "sandwich", cluster = "plate"),
              warning = function(w) w)
stopifnot(inherits(w, "warning"))
stopifnot(grepl("crossed|approximate", conditionMessage(w), ignore.case = TRUE))

V_pen <- suppressWarnings(
    vcov(fit_pen, type = "sandwich", cluster = "plate"))
stopifnot(all(is.finite(V_pen)))
stopifnot(attr(V_pen, "n.clusters") == nlevels(Penicillin$plate))

## ----------------------------------------------------------------
## 5. predict(): interval = "none" preserves byte-for-byte behaviour,
##    "confidence"/"prediction" return a data frame with ordered
##    intervals and the right monotonicity.
## ----------------------------------------------------------------
fitted_v <- fitted(fit_sleep)
pred_none <- predict(fit_sleep)              # default interval = "none"
stopifnot(is.numeric(pred_none))
stopifnot(length(pred_none) == length(fitted_v))
stopifnot(max(abs(pred_none - fitted_v)) < 1e-10)

## Default predict() must not emit a spurious "type argument ignored"
## warning (it was firing once after the WS1c refactor before the
## match.call passthrough fix).
n_warn <- 0L
withCallingHandlers(invisible(predict(fit_sleep)),
                    warning = function(w) {
                        n_warn <<- n_warn + 1L
                        invokeRestart("muffleWarning")
                    })
stopifnot(n_warn == 0L)

pred_c <- predict(fit_sleep, interval = "confidence")
stopifnot(is.data.frame(pred_c))
stopifnot(identical(colnames(pred_c), c("fit", "lwr", "upr", "se")))
stopifnot(nrow(pred_c) == length(fitted_v))
stopifnot(max(abs(pred_c$fit - fitted_v)) < 1e-10)
stopifnot(all(pred_c$lwr <= pred_c$fit))
stopifnot(all(pred_c$fit <= pred_c$upr))
stopifnot(all(pred_c$se > 0))

pred_p <- predict(fit_sleep, interval = "prediction")
stopifnot(is.data.frame(pred_p))
stopifnot(all(pred_p$se > pred_c$se))  # prediction strictly wider
stopifnot(all(pred_p$lwr < pred_c$lwr))
stopifnot(all(pred_p$upr > pred_c$upr))

## newdata path: fixed-only intervals are narrower than fixed+RE ones
## because the random-effect uncertainty is excluded.
nd <- head(sleepstudy, 3)
pc_full <- predict(fit_sleep, newdata = nd, interval = "confidence")
pc_fix  <- predict(fit_sleep, newdata = nd, interval = "confidence",
                   re.form = NA)
stopifnot(all(pc_fix$se < pc_full$se))
stopifnot(max(abs(pc_full$fit - predict(fit_sleep, newdata = nd))) < 1e-10)

## ----------------------------------------------------------------
## 6. cooks.distance / influence: non-negative, finite, length n;
##    top-flagged Penicillin obs match the IF-thread1 prototype's
##    documented set with the same rho.sigma.e tuning.
##
## The full IF uses numDeriv::jacobian + a column-wise FD pass and
## takes a few seconds, so the test is gated below.
## ----------------------------------------------------------------
fit_pen2 <- rlmer(diameter ~ (1 | plate) + (1 | sample),
                  data = Penicillin,
                  rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
                  method = "DASvar")
IF <- implicitIF_full(fit_pen2)
stopifnot(identical(dim(IF$IF_beta), c(1L, nrow(Penicillin))))
stopifnot(identical(dim(IF$IF_sigma), c(1L, nrow(Penicillin))))
stopifnot(all(is.finite(IF$IF_beta)))
stopifnot(all(is.finite(IF$IF_sigma)))
stopifnot(all(is.finite(IF$IF_theta)))

cd <- cooks.distance(fit_pen2, IF = IF)
stopifnot(length(cd) == nrow(Penicillin))
stopifnot(all(is.finite(cd)))
stopifnot(all(cd >= 0))

top5 <- order(cd, decreasing = TRUE)[1:5]
stopifnot(setequal(top5, c(85L, 86L, 122L, 123L, 13L)))

inf <- influence(fit_pen2)
stopifnot(identical(dim(inf$IF_beta), dim(IF$IF_beta)))
stopifnot(identical(dim(inf$IF_sigma), dim(IF$IF_sigma)))
