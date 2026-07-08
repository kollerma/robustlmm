## Step B (PLAN-WS7-WS9 §5.4): run the WS7-WS9 test computations against
## the OLD code (master, f2de854) and record the failing magnitudes.
## Printing instead of stopifnot; every section wrapped so all four run.

.libPaths(c("/Users/kollerma/Documents/R/lib-master-old", .libPaths()))
suppressMessages(require(robustlmm))
cat("== Step B: old-code bug confirmations ==\n")
cat(R.version.string, "\n")
cat("robustlmm", as.character(packageVersion("robustlmm")),
    "from", find.package("robustlmm"), "(master f2de854)\n\n")

sec <- function(title, expr) {
    cat("---", title, "---\n")
    tryCatch(expr, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
    cat("\n")
}

## ----------------------------------------------------------------
sec("(i) weighted-GLS sensitivity identity (test 1)", {
    set.seed(42)
    w <- runif(nrow(Dyestuff), 0.5, 3)
    fitW <- rlmer(Yield ~ 1 + (1 | Batch), Dyestuff, weights = w,
                  rho.e = cPsi, rho.b = cPsi, init = lmerNoFit,
                  method = "DASvar")
    parts <- robustlmm:::implicitIF(fitW)
    pp  <- fitW@pp
    X   <- as.matrix(pp$X)
    Z   <- as.matrix(t(pp$Zt))
    U_b <- as.matrix(pp$U_b)
    Omega <- Z %*% U_b %*% t(U_b) %*% t(Z) + diag(1 / w)
    M_gls <- solve(t(X) %*% solve(Omega, X)) %*% t(X) %*% solve(Omega)
    relErr1 <- max(abs(parts$IF_beta - M_gls)) / max(abs(M_gls))
    cat("max rel err IF_beta vs GLS sensitivity:", format(relErr1),
        "(branch: 4.4e-15; threshold 1e-7)\n")

    ## caseweight identity on the same fit, for completeness
    cw <- caseweightIF(fitW)
    resid_raw <- Dyestuff$Yield - fitted(fitW)
    relErr2 <- max(abs(cw$IF_beta - sweep(parts$IF_beta, 2, resid_raw, "*"))) /
        max(abs(cw$IF_beta))
    cat("max rel err caseweightIF identity:", format(relErr2),
        "(branch: 3.8e-15; threshold 1e-7)\n")
})

## ----------------------------------------------------------------
sec("(ii) strict offset comparison (test 4)", {
    set.seed(7)
    DyeOff <- within(Dyestuff, {
        off   <- rnorm(length(Yield), sd = sd(Yield))
        Yield <- Yield + off
    })
    cO <- lmer(Yield ~ 1 + (1 | Batch) + offset(off), DyeOff)
    rO <- rlmer(Yield ~ 1 + (1 | Batch) + offset(off), DyeOff,
                rho.e = cPsi, rho.b = cPsi, init = lmerNoFit)
    cat("lmer  fixef:", format(fixef(cO)), "\n")
    cat("rlmer fixef:", format(fixef(rO)), "\n")
    cat("abs diff fixef:", format(abs(fixef(cO) - fixef(rO))), "\n")
    cat("max abs diff fitted:", format(max(abs(fitted(cO) - fitted(rO)))),
        " (sd(Yield) =", format(sd(DyeOff$Yield)), ")\n")
    cat("max rel diff fitted:",
        format(max(abs(fitted(cO) - fitted(rO))) / sd(DyeOff$Yield)),
        "x sd(Yield)\n")
})

## ----------------------------------------------------------------
sec("(iii) DAStau theta-score root consistency (dastau-fallback sec. 1)", {
    ss <- sleepstudy
    set.seed(3)
    ii <- sample(nrow(ss), 12L)
    ss$Reaction[ii] <- ss$Reaction[ii] + 200
    fitT <- rlmer(Reaction ~ Days + (Days | Subject), ss, method = "DAStau")
    p <- fitT@pp$p; q <- fitT@pp$q
    L <- length(getME(fitT, "theta"))
    par0 <- c(robustlmm:::.fixef(fitT), fitT@pp$b.s,
              log(robustlmm:::.sigma(fitT)), getME(fitT, "theta"))
    s0       <- robustlmm:::.scoreVec(par0, fitT)
    F_theta0 <- s0[p + q + 1L + seq_len(L)]
    parP <- par0
    parP[p + q + 2L] <- parP[p + q + 2L] * 1.05
    sP       <- robustlmm:::.scoreVec(parP, fitT)
    F_thetaP <- sP[p + q + 1L + seq_len(L)]
    ratio <- max(abs(F_theta0)) / max(abs(F_thetaP - F_theta0))
    cat("theta-score at fit:", format(max(abs(F_theta0))),
        "; at 5% theta perturbation:", format(max(abs(F_thetaP - F_theta0))),
        "; ratio:", format(ratio), "(threshold 0.05)\n")
})

## ----------------------------------------------------------------
sec("(iv) dim-3 blocks: rlmer downgrade + direct calcTau.nondiag crash", {
    set.seed(11)
    nG <- 12L; nPer <- 6L; n <- nG * nPer
    g  <- factor(rep(seq_len(nG), each = nPer))
    x1 <- rnorm(n); x2 <- rnorm(n)
    b  <- matrix(rnorm(3L * nG, sd = c(1, 0.8, 0.6)), nrow = nG, byrow = TRUE)
    y  <- 1 + x1 + 0.5 * x2 +
        b[as.integer(g), 1L] + b[as.integer(g), 2L] * x1 +
        b[as.integer(g), 3L] * x2 + rnorm(n, sd = 0.5)
    d3 <- data.frame(y, x1, x2, g)
    warned <- character(0)
    fit3 <- withCallingHandlers(
        rlmer(y ~ x1 + x2 + (x1 + x2 | g), d3, method = "DAStau"),
        warning = function(w) {
            warned <<- c(warned, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    cat("rlmer warnings:", paste(warned, collapse = " | "), "\n")
    cat("fit3@method:", fit3@method, "\n")
    cat("direct calcTau.nondiag call:\n")
    Tb3 <- tryCatch(
        withCallingHandlers(
            robustlmm:::calcTau.nondiag(fit3,
                                        skbs = robustlmm:::.S(fit3),
                                        kappas = robustlmm:::.kappa_b(fit3),
                                        max.iter = 50L),
            warning = function(w) {
                cat("  warning:", conditionMessage(w), "\n")
                invokeRestart("muffleWarning")
            }),
        error = function(e) {
            cat("  ERROR (as predicted):", conditionMessage(e), "\n")
            NULL
        })
    if (!is.null(Tb3)) cat("  unexpectedly returned a value\n")
})

cat("== Step B done ==\n")
