## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Messages before long-running computations (options
## robustlmm.message.if.max and robustlmm.message.boot.max).
require(robustlmm)

msgs <- function(expr) {
    m <- character(0)
    withCallingHandlers(expr, message = function(c) {
        m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") })
    m
}
f0 <- rlmer(Reaction ~ Days + (1 | Subject), sleepstudy)
f1 <- rlmer(Reaction ~ Days + (Days || Subject), sleepstudy)
W <- robustlmm:::.satterthwaite_df_workload(f0)

## ---- defaults: these small computations are silent -------------------------
stopifnot(length(msgs(implicitIF_full(f0, use.cache = FALSE))) == 0,
          length(msgs(suppressWarnings(
              anova(f0, f1, test = "score", nsim = 10, seed = 1)))) == 0)

## ---- influence function ------------------------------------------------------
old <- options(robustlmm.message.if.max = W - 1)
m <- msgs(IF <- implicitIF_full(f0, use.cache = FALSE))
stopifnot(length(m) == 1, grepl("influence function", m),
          grepl("robustlmm.message.if.max", m, fixed = TRUE),
          grepl("summary(df = \"none\")", m, fixed = TRUE))
## the cached result needs no computation and gives no message
suppressMessages(invisible(implicitIF_full(f0)))
stopifnot(length(msgs(implicitIF_full(f0))) == 0)
## at the threshold and with Inf: silent
options(robustlmm.message.if.max = W)
stopifnot(length(msgs(implicitIF_full(f0, use.cache = FALSE))) == 0)
options(robustlmm.message.if.max = Inf)
stopifnot(length(msgs(implicitIF_full(f0, use.cache = FALSE))) == 0)
options(old)

## ---- bootstraps ----------------------------------------------------------------
old <- options(robustlmm.message.boot.max = 1)
mb <- msgs(suppressWarnings(anova(f0, f1, test = "boot", nsim = 10, seed = 1)))
ms <- msgs(suppressWarnings(anova(f0, f1, test = "score", nsim = 10, seed = 1)))
stopifnot(length(mb) == 1, grepl("anova(test = \"boot\") with 10 replicates",
                                 mb, fixed = TRUE),
          grepl("test = \"score\"", mb, fixed = TRUE),
          length(ms) == 1, grepl("anova(test = \"score\")", ms, fixed = TRUE),
          grepl("robustlmm.message.boot.max", ms, fixed = TRUE))
if (requireNamespace("confintROB", quietly = TRUE)) {
    mc <- msgs(suppressWarnings(confint(f0, method = "boot", nsim = 10,
                                        seed = 1)))
    stopifnot(length(mc) == 1, grepl("confint(method = \"boot\")", mc,
                                     fixed = TRUE),
              grepl("method = \"Wald\"", mc, fixed = TRUE))
}
## the result does not depend on the message
options(robustlmm.message.boot.max = Inf)
a1 <- suppressWarnings(anova(f0, f1, test = "score", nsim = 10, seed = 1))
options(robustlmm.message.boot.max = 1)
a2 <- suppressMessages(suppressWarnings(
    anova(f0, f1, test = "score", nsim = 10, seed = 1)))
options(old)
stopifnot(identical(a1[2, "Pr(>=Score)"], a2[2, "Pr(>=Score)"]))
cat("slow-message.R: all checks passed\n")
