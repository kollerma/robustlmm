## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Centred residual scale for RANSAC candidate selection
## (.ransac_mscale0: bisquare M-scale about zero), replacing Qn,
## which is invariant under a common shift of the residuals.
require(robustlmm)
ms0 <- robustlmm:::.ransac_mscale0

set.seed(3)
r <- rnorm(400)
## Qn cannot see a shift; the centred scale can.
stopifnot(isTRUE(all.equal(robustbase::Qn(r), robustbase::Qn(r + 3))),
          ms0(r) < ms0(r + 0.5), ms0(r + 0.5) < ms0(r + 3))
## symmetric in the sign of the residuals
stopifnot(isTRUE(all.equal(ms0(-r), ms0(r))))
## Fisher consistent at the standard normal
set.seed(4)
stopifnot(abs(ms0(rnorm(1e5)) - 1) < 0.02)
## bounded when 40% of the residuals are huge (breakdown point 1/2)
rr <- r; rr[1:160] <- 1e6
stopifnot(is.finite(ms0(rr)), ms0(rr) < 5)
## exact fit of more than half of the data gives scale 0
stopifnot(ms0(c(rep(0, 30), rnorm(20))) == 0)
## non-finite residuals are dropped
stopifnot(isTRUE(all.equal(ms0(c(r, NA, Inf)), ms0(r))))

## ransac_lme4 uses the centred scale by default; Qn stays available.
i1 <- ransac_lme4(Reaction ~ Days + (Days | Subject), sleepstudy,
                  K = 30, seed = 11)
i2 <- ransac_lme4(Reaction ~ Days + (Days | Subject), sleepstudy,
                  K = 30, seed = 11, scale_fn = robustbase::Qn)
stopifnot(inherits(i1$fit, "lmerMod"), inherits(i2$fit, "lmerMod"),
          is.finite(i1$scale), is.finite(i2$scale))
## the default's reported score is the centred scale of the chosen fit's
## full-data residuals
pr <- sleepstudy$Reaction - predict(i1$fit, newdata = sleepstudy,
                                    allow.new.levels = TRUE)
stopifnot(isTRUE(all.equal(i1$scale, ms0(pr))))

cat("ransac-centred-scale.R: all checks passed\n")
