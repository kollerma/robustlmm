## WS24: direction of the bias.  Prior weights are identified only up to
## a constant -- rescaling w -> c*w must rescale sigma-hat by sqrt(c) and
## leave everything else alone.  The pre-fix DASvar tau^2 (leading term
## v_e instead of 1) breaks that equivariance, so the sign of the bias
## depends on whether the weights sit above or below 1.  "Biased low" is
## therefore only half the story.
##   Rscript ws24-weight-scale.R <libpath> <label> <outdir>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1], .libPaths()))
label <- args[2]; outdir <- args[3]
suppressPackageStartupMessages({library(lme4); library(robustlmm)})
W <- NULL
cs <- c(0.01, 0.1, 0.5, 1, 2, 10, 100)
rows <- lapply(cs, function(cc) {
    W <<- rep(cc, nrow(Dyestuff))
    fl <- lmer(Yield ~ 1 | Batch, Dyestuff, weights = W)
    f <- function(m) suppressMessages(suppressWarnings(
        sigma(rlmer(Yield ~ 1 | Batch, Dyestuff, weights = W, rho.e = cPsi,
                    rho.b = cPsi, init = lmerNoFit, method = m))))
    data.frame(c = cc, lmer = sigma(fl), DASvar = f("DASvar"), DAStau = f("DAStau"))
})
tab <- do.call(rbind, rows)
## everything divided by sqrt(c): must be constant (= the unweighted sigma)
tab$lmer_eq   <- tab$lmer   / sqrt(tab$c)
tab$DASvar_eq <- tab$DASvar / sqrt(tab$c)
tab$DAStau_eq <- tab$DAStau / sqrt(tab$c)
tab$DASvar_relerr <- tab$DASvar / tab$lmer - 1
tab$robustlmm <- as.character(packageVersion("robustlmm"))
saveRDS(tab, file.path(outdir, paste0("weight-scale-", label, ".rds")))
cat("Dyestuff, constant weights w = c, rho = cPsi.  sigma-hat/sqrt(c) must be\n",
    "constant (= 49.5101, the unweighted REML sigma).\n\n", sep = "")
print(tab[, c("c","lmer_eq","DASvar_eq","DAStau_eq","DASvar_relerr")],
      digits = 6, row.names = FALSE)
