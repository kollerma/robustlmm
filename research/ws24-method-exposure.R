## WS24: WHICH method was actually broken, and how often?
##
## DASvar is broken directly (the tau^2 formula it evaluates).  DAStau
## is broken only INDIRECTLY: the DASvar formula is its initial value,
## and calcTau's `while (... && tau[i] <= 1)` guard pins tau at 1 for
## every observation whose (inflated) start exceeds 1.  Whether the fit
## recovers depends on whether the outer loop re-enters calcTau from a
## start that satisfies the guard -- so DAStau's exposure is a
## function of the design, the rho, and the starting values.
##
## This measures it: pre-fix vs post-fix sigma-hat over designs x
## weight regimes x rho x init, for both methods.  No lmer oracle is
## needed -- the question is whether the released answer differs from
## the corrected one.
##   Rscript ws24-method-exposure.R <prefix-lib> <fixed-lib> <outdir>
args <- commandArgs(trailingOnly = TRUE)
LIBS <- c(prefix = args[1], fixed = args[2]); outdir <- args[3]

SEED <- 20260726L
W <- NULL
designs <- list(
    "oneway 50x5"   = list(J = 50, m = 5,  form = y ~ x + (1 | g)),
    "oneway 12x4"   = list(J = 12, m = 4,  form = y ~ x + (1 | g)),
    "slope 20x8"    = list(J = 20, m = 8,  form = y ~ x + (x | g)))
wsets <- list(
    "const 2"            = function(n) rep(2, n),
    "const 0.5"          = function(n) rep(0.5, n),
    "const 100"          = function(n) rep(100, n),
    "logunif[0.25,4]"    = function(n) exp(runif(n, log(0.25), log(4))),
    "logunif[0.01,100]"  = function(n) exp(runif(n, log(0.01), log(100))),
    "rexp"               = function(n) rexp(n))

mkdata <- function(d, wfun, seed) {
    set.seed(seed)
    J <- d$J; m <- d$m; n <- J * m
    g <- factor(rep(seq_len(J), each = m))
    W <<- wfun(n)
    sc <- exp(mean(log(1 / sqrt(W))))
    x <- if (grepl("slope", deparse(d$form)[1])) rep(seq_len(m), J) / m else rnorm(n)
    re <- rnorm(J, sd = 0.7 * sc)[as.integer(g)]
    if (any(grepl("x \\| g", deparse(d$form))))
        re <- re + rnorm(J, sd = 1.2 * sc)[as.integer(g)] * x
    data.frame(y = 1 + 0.5 * x + re + rnorm(n, sd = 1 / sqrt(W)), x = x, g = g)
}

runOne <- function(lib) {
    script <- tempfile(fileext = ".R")
    writeLines(c(
      sprintf('.libPaths(c("%s", .libPaths()))', lib),
      'suppressPackageStartupMessages({library(lme4); library(robustlmm)})',
      sprintf('SEED <- %dL', SEED),
      deparse(bquote(designs <- .(designs))),
      deparse(bquote(wsets <- .(wsets))),
      deparse(bquote(mkdata <- .(mkdata))),
      'W <- NULL',
      'rows <- list()',
      'for (dn in names(designs)) for (wn in names(wsets)) {',
      '  d <- designs[[dn]]; dat <- mkdata(d, wsets[[wn]], SEED + nchar(wn))',
      '  for (rho in c("cPsi","smoothPsi")) for (m in c("DASvar","DAStau"))',
      '    for (ini in c("lmer","lmerNoFit")) {',
      '      f <- tryCatch(suppressMessages(suppressWarnings(',
      '        if (ini == "lmer") rlmer(d$form, dat, weights = W, rho.e = get(rho),',
      '                                 rho.b = get(rho), method = m)',
      '        else rlmer(d$form, dat, weights = W, rho.e = get(rho), rho.b = get(rho),',
      '                   method = m, init = lmerNoFit))), error = function(e) NULL)',
      '      rows[[length(rows)+1]] <- data.frame(design = dn, weights = wn, rho = rho,',
      '        method = m, init = ini,',
      '        sigma = if (is.null(f)) NA_real_ else sigma(f),',
      '        slope = if (is.null(f)) NA_real_ else fixef(f)[2],',
      '        se    = if (is.null(f)) NA_real_ else sqrt(diag(as.matrix(vcov(f))))[2],',
      '        tauEq1 = if (is.null(f)) NA_integer_ else sum(f@pp$tau_e() == 1),',
      '        nWarn = if (is.null(f)) NA_integer_ else length(unlist(f@optinfo$warnings)),',
      '        stringsAsFactors = FALSE)',
      '    }',
      '}',
      'saveRDS(do.call(rbind, rows), commandArgs(TRUE)[1])'), script)
    out <- tempfile(fileext = ".rds")
    system2("Rscript", c(script, out), stdout = FALSE, stderr = FALSE)
    readRDS(out)
}

pre <- runOne(LIBS["prefix"]); post <- runOne(LIBS["fixed"])
key <- c("design","weights","rho","method","init")
m <- merge(pre, post, by = key, suffixes = c("_pre","_post"))
m$sigma_relerr <- m$sigma_pre / m$sigma_post - 1
m$slope_relerr <- m$slope_pre / m$slope_post - 1
m$se_relerr    <- m$se_pre / m$se_post - 1
m$broken <- abs(m$sigma_relerr) > 1e-6
## version stamps of the two builds compared
stampOf <- function(lib)
    as.character(utils::packageVersion("robustlmm", lib.loc = lib))
attr(m, "stamp") <- c(prefix = stampOf(LIBS["prefix"]),
                      fixed  = stampOf(LIBS["fixed"]),
                      lme4   = as.character(utils::packageVersion("lme4")),
                      R      = R.version.string)
saveRDS(m, file.path(outdir, "ws24-method-exposure-results.rds"))

cat("Pre-fix vs post-fix sigma-hat.  'broken' = |rel err| > 1e-6.\n")
cat(sprintf("\n%d arms total (3 designs x 6 weight regimes x 2 rho x 2 methods x 2 inits)\n",
            nrow(m)))
tab <- with(m, table(method = method, broken = broken))
cat("\nBy method:\n"); print(tab)
cat("\nBy method x rho x init (fraction of arms broken):\n")
print(round(with(m, tapply(broken, list(paste(method, rho), init), mean)), 3))
cat("\nWorst |sigma rel err| by method x rho:\n")
print(signif(with(m, tapply(abs(sigma_relerr), list(method, rho), max, na.rm = TRUE)), 3))
cat("\nDAStau arms that ARE broken:\n")
b <- m[m$method == "DAStau" & m$broken, c(key, "sigma_relerr", "se_relerr", "tauEq1_pre")]
print(b, digits = 3, row.names = FALSE)
