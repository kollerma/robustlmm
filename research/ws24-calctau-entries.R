.libPaths(c("/Users/kollerma/Documents/R/ws24-sweep/libs/dev-prefix", .libPaths()))
suppressPackageStartupMessages({library(lme4); library(robustlmm)})
## Why is DAStau broken at cPsi but not at smoothPsi?  Count how many
## times calcTau is entered for the e-side during the fit: once = the
## cold DASvar start is final; more = the warm start (<= 1) re-enters
## and the guard no longer blocks, so the damage is undone.
MASTER_SEED <- 20260611L; J <- 50L; n_per <- 5L; n <- J*n_per
set.seed(MASTER_SEED); W <<- exp(runif(n, log(0.25), log(4))); x <- rnorm(n)
g <- factor(rep(seq_len(J), each=n_per))
set.seed(MASTER_SEED+1L); b <- rnorm(J, sd=0.5); e <- rnorm(n, sd=1/sqrt(W))
dat <- data.frame(y=1+0.5*x+b[as.integer(g)]+e, x=x, g=g)
for (rho in c("cPsi","smoothPsi")) {
  n_e <- 0L; over1 <- integer(0)
  trace(robustlmm:::calcTau, tracer = quote({
      if (length(tau) == 250) { n_e <<- n_e + 1L; over1 <<- c(over1, sum(tau > 1)) }
  }), print = FALSE)
  f <- suppressMessages(suppressWarnings(
      rlmer(y ~ x + (1|g), dat, weights=W, rho.e=get(rho), rho.b=get(rho), method="DAStau")))
  untrace(robustlmm:::calcTau)
  cat(sprintf("%-9s e-side calcTau entries: %d | starts above 1 per entry: %s | final tau==1: %3d | sigma %.6f\n",
      rho, n_e, paste(over1, collapse=","), sum(f@pp$tau_e()==1), sigma(f)))
}
