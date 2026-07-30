## WS24 / Q1: push the cPsi == lmer oracle past the tests/weights.R grid.
## At rho.e = rho.b = cPsi the RSE reduces to REML, so a weighted fit
## must reproduce lmer exactly, for BOTH tau methods.  Sweeps weight
## ranges out to [0.01, 100], unbalanced designs and a random-slope
## (block) design, and records the achieved tolerance, not just PASS.
## Also records whether the DASvar tau^2 initial value ever exceeds 1
## (which makes calcTau's `while (... && tau[i] <= 1)` guard skip the
## fixed point -- the mechanism behind the DAStau half of the bug).
##   Rscript ws24-identity-grid.R <libpath> <label> <outdir>
args <- commandArgs(trailingOnly = TRUE)
libpath <- args[1]; label <- args[2]; outdir <- args[3]
.libPaths(c(libpath, .libPaths()))
suppressPackageStartupMessages({library(lme4); library(robustlmm)})

SEED <- 20260726L

## ---------------- designs -------------------------------------------
mkBalanced <- function(J = 20, m = 6) {
    g <- factor(rep(seq_len(J), each = m)); n <- J * m
    list(name = "balanced 20x6", form = y ~ x + (1 | g),
         data = data.frame(x = rnorm(n), g = g), n = n)
}
mkUnbalanced <- function(J = 25) {
    m <- sample(1:9, J, replace = TRUE); g <- factor(rep(seq_len(J), m))
    n <- length(g)
    list(name = paste0("unbalanced ", J, " grp (n=", n, ")"),
         form = y ~ x + (1 | g), data = data.frame(x = rnorm(n), g = g), n = n)
}
mkSlope <- function(J = 20, m = 8) {
    g <- factor(rep(seq_len(J), each = m)); n <- J * m
    list(name = "random slope 20x8", form = y ~ x + (x | g),
         data = data.frame(x = rep(seq_len(m), J) / m, g = g), n = n)
}
mkCrossed <- function(a = 8, b = 7) {
    d <- expand.grid(f1 = factor(seq_len(a)), f2 = factor(seq_len(b)))
    n <- nrow(d); d$x <- rnorm(n)
    list(name = "crossed 8x7", form = y ~ x + (1 | f1) + (1 | f2),
         data = d, n = n)
}

## The random-effect SDs are scaled by the typical residual SD
## (1/sqrt(w)) so that the variance-component signal-to-noise -- and
## hence the distance of theta from the singular boundary -- does not
## change with the weight scale.  Without this, `weights = 0.01` alone
## pushes lmer to a boundary fit and the two optimizers stop at
## different points on it, which has nothing to do with tau.
addResponse <- function(des, sb = 0.7) {
    d <- des$data; n <- des$n
    sc <- exp(mean(log(1 / sqrt(WGT))))       # geometric-mean residual SD
    if (!is.null(d$g)) {
        J <- nlevels(d$g); b0 <- rnorm(J, sd = sb * sc)
        re <- b0[as.integer(d$g)]
        if (grepl("slope", des$name))
            re <- re + rnorm(J, sd = 1.2 * sc)[as.integer(d$g)] * d$x
    } else {
        re <- rnorm(nlevels(d$f1), sd = sb * sc)[as.integer(d$f1)] +
              rnorm(nlevels(d$f2), sd = sb * sc)[as.integer(d$f2)]
    }
    d$y <- 1 + 0.5 * d$x + re + rnorm(n, sd = 1 / sqrt(WGT))
    des$data <- d; des
}

WGT <- NULL   # global, so the formula environment can see it

wsets <- list(
    "const 100"        = function(n) rep(100, n),
    "const 0.01"       = function(n) rep(0.01, n),
    "logunif[0.01,100]"= function(n) exp(runif(n, log(0.01), log(100))),
    "logunif[0.1,10]"  = function(n) exp(runif(n, log(0.1), log(10))),
    "two-point{0.01,100}" = function(n) ifelse(runif(n) < .5, 0.01, 100),
    "rexp"             = function(n) rexp(n))

relerr <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-12))

rows <- list()
for (mk in list(mkBalanced, mkUnbalanced, mkSlope, mkCrossed)) {
  set.seed(SEED)
  des0 <- mk()
  for (wn in names(wsets)) {
    set.seed(SEED + nchar(wn))
    WGT <<- wsets[[wn]](des0$n)
    des <- addResponse(des0)
    d <- des$data; form <- des$form
    ## UNWEIGHTED CONTROL: the same design and response fitted with no
    ## prior weights at all.  rlmer and lmer still disagree at the level
    ## of the theta optimiser's own stopping rule, and that baseline is
    ## what the weighted arms have to be judged against -- anything at
    ## or below it is not attributable to the weights or to tau.
    fl <- tryCatch(lmer(form, d, weights = WGT, REML = TRUE,
                        control = lmerControl(check.conv.singular = "ignore")),
                   error = function(e) NULL)
    flU <- tryCatch(lmer(form, d, REML = TRUE,
                         control = lmerControl(check.conv.singular = "ignore")),
                    error = function(e) NULL)
    for (method in c("DASvar", "DAStau")) {
      fr <- tryCatch(suppressMessages(suppressWarnings(
              rlmer(form, d, weights = WGT, rho.e = cPsi, rho.b = cPsi,
                    method = method))), error = function(e) NULL)
      frU <- tryCatch(suppressMessages(suppressWarnings(
              rlmer(form, d, rho.e = cPsi, rho.b = cPsi,
                    method = method))), error = function(e) NULL)
      ctrl <- if (is.null(flU) || is.null(frU)) NA_real_ else
          max(abs(sigma(frU)/sigma(flU) - 1), relerr(fixef(frU), fixef(flU)),
              relerr(sqrt(diag(as.matrix(vcov(frU)))),
                     sqrt(diag(as.matrix(vcov(flU))))),
              max(abs(unlist(ranef(frU)) - unlist(ranef(flU)))) /
                  max(1e-12, max(abs(unlist(ranef(flU))))))
      if (is.null(fl) || is.null(fr)) {
        rows[[length(rows) + 1]] <- data.frame(
          design = des$name, weights = wn, method = method,
          sigma_r = NA, sigma_l = NA, tol_sigma = NA, tol_fixef = NA,
          tol_se = NA, tol_ranef = NA, initTauOver1 = NA, tauEq1 = NA,
          nWarn = NA, singular = NA, unweighted_ctrl = NA,
          stringsAsFactors = FALSE)
        next
      }
      pp <- fr@pp
      ## DASvar-form initial tau^2 at the converged fit
      Btmp <- pp$B()
      tmp  <- Matrix::tcrossprod(pp$Epsi_bpsi_bt, Btmp)
      bt   <- robustlmm:::computeDiagonalOfProduct(as(Btmp, "unpackedMatrix"),
                                                   as(tmp, "unpackedMatrix"))
      tau2init <- 1 - pp$EDpsi_e * 2 * pp$diagA + pp$Epsi2_e * pp$diagAAt + bt
      tt <- pp$tau_e()
      rows[[length(rows) + 1]] <- data.frame(
        design = des$name, weights = wn, method = method,
        sigma_r = sigma(fr), sigma_l = sigma(fl),
        tol_sigma = abs(sigma(fr) / sigma(fl) - 1),
        tol_fixef = relerr(fixef(fr), fixef(fl)),
        tol_se    = relerr(sqrt(diag(as.matrix(vcov(fr)))),
                           sqrt(diag(as.matrix(vcov(fl))))),
        tol_ranef = max(abs(unlist(ranef(fr)) - unlist(ranef(fl)))) /
                    max(1e-12, max(abs(unlist(ranef(fl))))),
        initTauOver1 = sum(tau2init > 1), tauEq1 = sum(tt == 1),
        nWarn = length(unlist(fr@optinfo$warnings)),
        singular = isSingular(fl), unweighted_ctrl = ctrl,
        stringsAsFactors = FALSE)
    }
  }
}
res <- do.call(rbind, rows)
res$stamp_robustlmm <- as.character(packageVersion("robustlmm"))
res$stamp_lme4 <- as.character(packageVersion("lme4"))
res$stamp_R <- R.version.string
res$seed <- SEED
saveRDS(res, file.path(outdir, paste0("identity-grid-", label, ".rds")))
res$worst <- pmax(res$tol_sigma, res$tol_fixef, res$tol_se, res$tol_ranef)
print(res[, c("design","weights","method","worst","unweighted_ctrl",
              "initTauOver1","tauEq1","nWarn","singular")],
      digits = 3, row.names = FALSE)
ok <- !res$singular & !is.na(res$worst)
cat("\nWORST over the", sum(ok), "NON-SINGULAR arms: sigma",
    max(res$tol_sigma[ok]), " fixef", max(res$tol_fixef[ok]),
    " se", max(res$tol_se[ok]), " ranef", max(res$tol_ranef[ok]), "\n")
cat("singular lmer reference in", sum(res$singular, na.rm=TRUE), "arms\n")
cat("worst weighted deviation / unweighted control, by design:\n")
print(round(tapply(res$worst[ok] / res$unweighted_ctrl[ok],
                   res$design[ok], max), 2))
cat("arms whose weighted deviation exceeds the unweighted control:",
    sum(res$worst[ok] > res$unweighted_ctrl[ok]), "of", sum(ok), "\n")
cat("arms where the DASvar-form tau^2 initial value exceeded 1:",
    sum(res$initTauOver1 > 0, na.rm=TRUE), "; tau pinned at 1:",
    sum(res$tauEq1 > 0, na.rm=TRUE), "; fit warnings:", sum(res$nWarn), "\n")
