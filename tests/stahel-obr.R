## Test the size_obr option for block-diagonal V_b.
##
## When size_obr = TRUE the size-controlling weight w_delta in
## rlmer.fit.DAS.nondiag is replaced by the Hampel-OBR form
##   w_tau(d) = min(1, b_tau / |d - s*kappa - a|)
## with b_tau from rho.sigma.b's tuning constant and a solved for
## Fisher consistency under chi^2_s. Default (size_obr = FALSE)
## leaves the existing finite-difference w_delta intact.

require(robustlmm)

## --- 1. Default behaviour is unchanged (regression test) ------------
fit_default <- rlmer(Reaction ~ Days + (Days | Subject),
                     data = sleepstudy,
                     method = "DASvar")
fit_explicit_false <- rlmer(Reaction ~ Days + (Days | Subject),
                            data = sleepstudy,
                            method = "DASvar",
                            size_obr = FALSE)
stopifnot(all.equal(robustlmm:::.fixef(fit_default),
                    robustlmm:::.fixef(fit_explicit_false)))
stopifnot(all.equal(robustlmm:::.sigma(fit_default),
                    robustlmm:::.sigma(fit_explicit_false)))
stopifnot(all.equal(getME(fit_default, "theta"),
                    getME(fit_explicit_false, "theta")))
stopifnot(isFALSE(fit_default@optinfo$size_obr))

## --- 2. size_obr = TRUE fits and produces sensible values -----------
fit_obr <- rlmer(Reaction ~ Days + (Days | Subject),
                 data = sleepstudy,
                 method = "DASvar",
                 size_obr = TRUE)
stopifnot(isTRUE(fit_obr@optinfo$size_obr))
stopifnot(robustlmm:::.sigma(fit_obr) > 0)
theta_obr <- getME(fit_obr, "theta")
stopifnot(all(is.finite(theta_obr)))
stopifnot(theta_obr[1] > 0)   # intercept SD positive
stopifnot(theta_obr[3] > 0)   # slope SD positive

## --- 3. OBR fit close to (but not identical to) default fit ----------
## On the clean sleepstudy at default tuning, OBR vs Proposal-II
## differs only at the 1-2 pp asymptotic-efficiency level; the
## fitted (beta, sigma, theta) should be close but not identical.
beta_diff  <- abs(robustlmm:::.fixef(fit_default) -
                  robustlmm:::.fixef(fit_obr))
sigma_diff <- abs(robustlmm:::.sigma(fit_default) -
                  robustlmm:::.sigma(fit_obr))
theta_diff <- abs(getME(fit_default, "theta") -
                  getME(fit_obr, "theta"))
stopifnot(all(beta_diff < 1))        # very loose; should be << 1
stopifnot(sigma_diff < 5)            # very loose; should be small
stopifnot(all(theta_diff < 0.5))     # very loose

## --- 4. size_obr has no effect on diagonal V_b ---------------------
fit_pen_default <- rlmer(diameter ~ (1 | plate) + (1 | sample),
                         data = Penicillin,
                         method = "DASvar")
fit_pen_obr <- rlmer(diameter ~ (1 | plate) + (1 | sample),
                     data = Penicillin,
                     method = "DASvar",
                     size_obr = TRUE)
stopifnot(all.equal(robustlmm:::.fixef(fit_pen_default),
                    robustlmm:::.fixef(fit_pen_obr)))
stopifnot(all.equal(robustlmm:::.sigma(fit_pen_default),
                    robustlmm:::.sigma(fit_pen_obr)))
stopifnot(all.equal(getME(fit_pen_default, "theta"),
                    getME(fit_pen_obr, "theta")))

cat("All stahel-obr tests passed.\n")
