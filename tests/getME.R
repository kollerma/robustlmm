## testing getME function
require(robustlmm)

rfm <- rlmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy,
             rho.e = smoothPsi, rho.b=smoothPsi,
             rho.sigma.e = psi2propII(smoothPsi, k=2.28),
             rho.sigma.b = psi2propII(smoothPsi, k=2.28))

(nmME <- eval(formals(getME)$name))
for (nm in nmME) {
    cat("\nName:", nm, "\n")
    str(getME(rfm, name=nm))
}
