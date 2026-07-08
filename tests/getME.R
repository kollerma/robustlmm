## testing getME function
require(robustlmm)
formatNum <- function(x, ...)
     format(round(x, 8), trim = TRUE, drop0trailing = TRUE, ...)
options(str = strOptions(formatNum = formatNum))

sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])
rfm <- rlmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
             rho.e = smoothPsi, rho.b=smoothPsi,
             rho.sigma.e = psi2propII(smoothPsi, k=2.28),
             rho.sigma.b = psi2propII(smoothPsi, k=2.28),
             init = lmer(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2,
                         control=lmerControl(optimizer="bobyqa")),
             doFit=FALSE)

## lme4's devcomp$dims vector gained a 13th entry in lme4 2.0-2, which
## would make this output differ between CRAN and development lme4.
## Check a version-stable subset explicitly and drop dims from the
## printed structure so the Rout.save works with both.
checkAndDropDims <- function(value) {
    stopifnot(value$dims[["N"]] == 180, value$dims[["n"]] == 180,
              value$dims[["p"]] == 2, value$dims[["nmp"]] == 178,
              value$dims[["q"]] == 40)
    value$dims <- NULL
    value
}

(nmME <- eval(formals(robustlmm:::getME.rlmerMod)$name))
for (nm in nmME) {
    cat("\nName:", nm, "\n")
    value <- getME(rfm, name=nm)
    if (nm == "theta") {
        value <- value + 1
    } else if (nm == "A") {
        value@x <- value@x + 1
    } else if (nm == "devcomp") {
        value <- checkAndDropDims(value)
    }
    if (substr(nm, 1, 3) == "rho") {
       print(value)
    } else {
       str(value)
    }
}
g.all <- getME(rfm, "ALL")
g.all[grepl("^rho", names(g.all))] <- NULL
g.all[["theta"]] <- g.all[["theta"]] + 1
g.all[["devcomp"]] <- checkAndDropDims(g.all[["devcomp"]])
str(g.all, max.level = 2)
