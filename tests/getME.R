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

(nmME <- eval(formals(robustlmm:::getME.rlmerMod)$name))
for (nm in nmME) {
    cat("\nName:", nm, "\n")
    if (substr(nm, 1, 3) == "rho") {
       print(getME(rfm, name=nm))
    } else {
       str(getME(rfm, name=nm))
    }
}
g.all <- getME(rfm, "ALL")
g.all[grepl("^rho", names(g.all))] <- NULL
str(g.all, max.level = 2)
