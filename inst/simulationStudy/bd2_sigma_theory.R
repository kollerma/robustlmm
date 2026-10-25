## Shipped with robustlmm: a script of the breakdown study behind
## setting = "RSEhb" in vignette("simulationStudies"). Run it with Rscript
## from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## BD2: breakdown point and contamination bias of rlmer's residual-scale
## equation  sum_i eta_i tau_i^2 w(t_i) (t_i^2 - kappa) = 0,  t_i = r_i/(tau_i sigma)
## (DAS-scale.R, updateSigma), in the location-scale limit tau = eta = 1.
## chi(t) = w(t) (t^2 - kappa) is even, chi(0) = -kappa, chi(inf) = lim w(t) t^2.
## For an M-scale equation mean chi(r/s) = 0 with bounded chi the explosion
## breakdown point is |chi(0)| / (|chi(0)| + chi(inf)) = kappa / (kappa + chi(inf))
## (Huber 1981, Sec. 5.4; Maronna, Martin, Yohai and Salibian-Barrera 2019,
## Sec. 2.8.2), and the implosion breakdown point is chi(inf)/(kappa + chi(inf)).
suppressMessages(library(robustlmm))

EZ <- function(f) integrate(function(z) f(z) * dnorm(z), -Inf, Inf, rel.tol = 1e-10)$value
scale_eq <- function(rho) {
    w <- rho@wgt
    kappa <- EZ(function(z) w(z) * z^2) / EZ(function(z) w(z))
    chi <- function(t) w(t) * (t^2 - kappa)
    chi_inf <- chi(1e4)
    list(w = w, kappa = kappa, chi = chi, chi_inf = chi_inf,
         bp = kappa / (kappa + chi_inf), bp_imp = chi_inf / (kappa + chi_inf))
}
## population scale under (1 - eps) N(0, 1) + eps N(Delta, 1)
pop_scale <- function(se, eps, Delta) {
    f <- function(s) (1 - eps) * EZ(function(z) se$chi(z / s)) +
        eps * EZ(function(z) se$chi((z + Delta) / s))
    if (f(1e3) > 0) return(Inf)                    # exploded
    uniroot(f, c(0.05, 1e3), tol = 1e-10)$root
}
settings <- list(
    "RSEn: psi2propII(smoothPsi, k = 1.345)" = psi2propII(smoothPsi, k = 1.345),
    "RSEa: psi2propII(smoothPsi, k = 2.28)"  = psi2propII(smoothPsi, k = 2.28))
res <- list()
for (nm in names(settings)) {
    se <- scale_eq(settings[[nm]])
    cat(sprintf("%s: kappa = %.4f, chi(inf) = %.4f, explosion bp = %.3f, implosion bp = %.3f\n",
                nm, se$kappa, se$chi_inf, se$bp, se$bp_imp))
    grid <- expand.grid(eps = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3), Delta = c(4, 6, 8, 12))
    grid$sigma <- mapply(function(e, D) pop_scale(se, e, D), grid$eps, grid$Delta)
    ## outliers at Delta are rejected by a bisquare(c = 4.685) at this sigma
    ## only if Delta > 4.685 sigma (weight 0 beyond c sigma)
    grid$bisq_rejects <- grid$Delta > 4.685 * grid$sigma
    grid$setting <- nm
    res[[nm]] <- grid
    print(reshape(grid[, c("eps", "Delta", "sigma")], idvar = "eps", timevar = "Delta",
                  direction = "wide"), digits = 3, row.names = FALSE)
}
out <- do.call(rbind, res)
saveRDS(out, "bd2_sigma_theory.rds")
