## Precompute the numbers quoted in the section "Studies for the
## high-breakdown configuration" of vignettes/simulationStudies.Rnw, so that
## the vignette pulls them from data via \Sexpr{} rather than hard-coding
## them. Reads the (build-ignored) result files under research/breakdown and
## writes the small, shipped RDS vignettes/high-breakdown-numbers.rds.
## Needs robustlmm (for the breakdown points of the scale equations).
##
## Re-run from the package root after re-running any of the studies:
##   Rscript vignettes/precompute-high-breakdown.R
##
## Build-ignored (see .Rbuildignore); the studies it reads are not shipped.

suppressPackageStartupMessages(library(robustlmm))
rd <- "research/breakdown"
old <- options(bd.source_only = TRUE)
source(file.path(rd, "bd_summarise.R"), local = TRUE)   # bd23(), dir
options(old)

## --- breakdown point of the scale equation sum w(t)(t^2 - kappa) = 0 ----
## kappa / (kappa + chi(Inf)), chi(t) = w(t)(t^2 - kappa)
scaleBp <- function(rho) {
    kappa <- robustlmm:::calcKappaTau(rho, 1)
    t <- 1e6
    chiInf <- rho@wgt(t) * (t^2 - kappa)
    kappa / (kappa + chiInf)
}
bp <- c(RSEn = scaleBp(psi2propII(smoothPsi)),
        RSEa = scaleBp(psi2propII(smoothPsi, k = 2.28)))

## --- the scale equation alone, location known (BD2) --------------------
th <- readRDS(file.path(rd, "bd2_sigma_theory.rds"))
sig20 <- function(tag) th$sigma[grepl(tag, th$setting) & th$eps == 0.2 &
                                th$Delta == 8]
## --- the breakdown map of the current settings (BD2/BD3) ---------------
s23 <- bd23()
v20 <- s23[s23$type == "V8" & s23$eps == 0.2, ]
current <- c("RSEn", "RSEa", "bisq_RSEn", "bisq_RSEa")
map <- list(
    nreps = max(s23$n_ok),
    sigmaKnownLoc = c(RSEn = sig20("RSEn"), RSEa = sig20("RSEa")),
    ## median residual-scale estimate at 20% vertical outliers, by the
    ## scale function the fit uses (RSEn: k = 1.345, RSEa: k = 2.28)
    medSigmaRSEn = range(v20$med_sigma[v20$est %in% c("RSEn", "bisq_RSEn")]),
    medSigmaRSEa = range(v20$med_sigma[v20$est %in% c("RSEa", "bisq_RSEa")]),
    sdbZero = range(v20$sd_b_zero[v20$est %in% current]),
    sdbZeroLmer = v20$sd_b_zero[v20$est == "lmer"])

## --- starts (BD5) --------------------------------------------------------
fs <- list.files(file.path(rd, "bd5_results"), "^cell_.*[.]rds$",
                 full.names = TRUE)
bd5 <- do.call(rbind, lapply(fs, function(f) {
    o <- readRDS(f)
    if (o$eps == 0 || o$eps > 0.3) return(NULL)
    arr <- simplify2array(o$res)                      # 8 x starts x reps
    do.call(rbind, lapply(dimnames(arr)[[2]], function(st) {
        err <- if (o$type %in% c("V8", "C8")) arr[5, st, ] - 2
               else arr[6, st, ] - 1
        data.frame(type = o$type, eps = o$eps, start = st,
                   broken = mean(abs(err) > 0.5, na.rm = TRUE))
    }))
}))
starts <- tapply(bd5$broken, bd5$start, mean)

## --- the RSEhb configuration: breakdown (BD6b2) ---------------------------
fs <- list.files(file.path(rd, "bd6b2_results"), "^cell_.*[.]rds$",
                 full.names = TRUE)
b6 <- do.call(rbind, lapply(fs, function(f) {
    o <- readRDS(f); r <- o$res
    do.call(rbind, lapply(c("RSEn", "RSEhb"), function(es) {
        int <- r[, paste0(es, "_int")] - 2; sl <- r[, paste0(es, "_slope")] - 1
        err <- if (o$type %in% c("V8", "C8")) int else sl
        z <- r[, paste0(es, "_sd_b")] < 1e-6
        data.frame(type = o$type, eps = o$eps, est = es, n = nrow(r),
                   broken = mean(abs(err) > 0.5, na.rm = TRUE),
                   sdbZero = mean(z, na.rm = TRUE),
                   brokenIfZero = if (any(z, na.rm = TRUE))
                       mean(abs(err[z]) > 0.5, na.rm = TRUE) else NA_real_,
                   nZero = sum(z, na.rm = TRUE),
                   brokenIfNotZero = mean(abs(err[!z]) > 0.5, na.rm = TRUE),
                   nNotZero = sum(!z, na.rm = TRUE))
    }))
}))
b6 <- b6[order(match(b6$type, c("V8", "C8", "L4", "L8")), b6$eps), ]
cont <- b6[b6$eps > 0 & b6$eps <= 0.3, ]
hb <- b6[b6$est == "RSEhb", ]
hbCont <- hb[hb$eps >= 0.2 & hb$eps <= 0.3, ]
v30 <- hb[hb$type == "V8" & hb$eps == 0.3, ]
c30 <- hb[hb$type == "C8" & hb$eps == 0.3, ]
config <- list(
    nreps = max(b6$n),
    table = b6,
    meanBroken = tapply(cont$broken, cont$est, mean),
    cleanBroken = b6$broken[b6$type == "V8" & b6$eps == 0],
    sdbZeroMax = max(hbCont$sdbZero),
    zeroV30 = c(broken = round(v30$brokenIfZero * v30$nZero), n = v30$nZero,
                brokenOther = round(v30$brokenIfNotZero * v30$nNotZero),
                nOther = v30$nNotZero),
    zeroC30 = c(broken = round(c30$brokenIfZero * c30$nZero), n = c30$nZero,
                brokenOther = round(c30$brokenIfNotZero * c30$nNotZero),
                nOther = c30$nNotZero))

## --- the RSEhb configuration: efficiency and inference (BD6c, BD6) -------
cellFile <- function(sub, J, type)
    readRDS(file.path(rd, sub, sprintf("cell_J%d_%s.rds", J, type)))$res
mse <- function(r, nm, which, truth)
    mean((r[, paste0(nm, "_", which)] - truth)^2, na.rm = TRUE)
eff <- function(r, nm) c(int = mse(r, "lmer", "int", 2) / mse(r, nm, "int", 2),
                         slope = mse(r, "lmer", "slope", 1) / mse(r, nm, "slope", 1))
cov <- function(r, nm, kind) c(int = mean(r[, paste0(nm, "_", kind, "_int")], na.rm = TRUE),
                               slope = mean(r[, paste0(nm, "_", kind, "_slope")], na.rm = TRUE))
c20 <- cellFile("bd6c_results", 20, "clean"); c50 <- cellFile("bd6c_results", 50, "clean")
v50 <- cellFile("bd6c_results", 50, "V8"); k50 <- cellFile("bd6c_results", 50, "C8")
dv50 <- cellFile("bd6_results", 50, "clean")                # DASvar (BD6)
inference <- list(
    nrepsClean = nrow(c50), nrepsCont = nrow(v50),
    effHB20 = eff(c20, "HB"), effHB50 = eff(c50, "HB"),
    effRSEn = range(c(eff(c20, "RSEn"), eff(c50, "RSEn"))),
    effRSEnInt = range(c(eff(c20, "RSEn")["int"], eff(c50, "RSEn")["int"])),
    effRSEnSlope = range(c(eff(c20, "RSEn")["slope"], eff(c50, "RSEn")["slope"])),
    covCR2 = range(c(cov(c50, "HB", "covs"), cov(v50, "HB", "covs"),
                     cov(k50, "HB", "covs"))),
    covSattIntC8 = cov(k50, "HB", "covt")[["int"]],
    dasvarSigmaLow = 1 - median(dv50[, "HB_sigma"], na.rm = TRUE) /
        median(c50[, "HB_sigma"], na.rm = TRUE),
    meanSigmaC8 = mean(k50[, "HB_sigma"], na.rm = TRUE),
    meanSdbC8 = mean(k50[, "HB_sd_b"], na.rm = TRUE),
    meanSigmaV8 = mean(v50[, "HB_sigma"], na.rm = TRUE))

numbers <- list(bp = bp, map = map, starts = starts, config = config,
                inference = inference)
saveRDS(numbers, "vignettes/high-breakdown-numbers.rds")
cat("wrote vignettes/high-breakdown-numbers.rds\n")
str(numbers, max.level = 2)
