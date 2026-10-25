## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## fragility(): first-order contamination fraction driving each bounded
## variance parameter to its bound (fragility.R; paper 2 threshold).
require(robustlmm)

fit <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
             method = "DASvar")
fr <- fragility(fit)
stopifnot(nrow(fr) == 2L,
          identical(fr$bound, c("zero variance", "singular covariance")),
          all(fr$eps > 0), all(is.finite(fr$eps)),
          all(fr$driver %in% levels(sleepstudy$Subject)),
          ## m added copies of the driver are a fraction m / (18 + m)
          isTRUE(all.equal(fr$clusters, fr$eps * 18 / (1 - fr$eps))))

## The first-order prediction matches a refit with one extra copy of the
## driving subject (eps = 1/19).
drv <- fr$driver[1]
cp <- sleepstudy[sleepstudy$Subject == drv, ]
cp$Subject <- "copy"
d2 <- rbind(transform(sleepstudy, Subject = as.character(Subject)), cp)
d2$Subject <- factor(d2$Subject)
f2 <- rlmer(Reaction ~ Days + (Days | Subject), d2, method = "DASvar")
th <- getME(fit, "theta")
dth_actual <- getME(f2, "theta")[1] - th[1]
## eps * IF = -(eps / eps*) * theta for the driving cluster
dth_pred <- -(1 / 19) / fr$eps[1] * th[1]
stopifnot(abs(dth_actual - dth_pred) < 0.1 * abs(dth_pred))

## scalar random intercept
fit2 <- rlmer(Yield ~ (1 | Batch), Dyestuff, method = "DASvar")
fr2 <- fragility(fit2)
stopifnot(nrow(fr2) == 1L, fr2$bound == "zero variance", fr2$eps > 0)

## a component already at zero has eps* = 0
fit3 <- suppressWarnings(rlmer(Yield ~ (1 | Batch), Dyestuff2,
                               method = "DASvar"))
if (getME(fit3, "theta") < 1e-8)
    stopifnot(fragility(fit3)$eps == 0)

## crossed designs are rejected with a clear message
fit4 <- rlmer(diameter ~ (1 | plate) + (1 | sample), Penicillin,
              method = "DASvar")
e <- tryCatch(fragility(fit4), error = function(e) conditionMessage(e))
stopifnot(grepl("crossed", e))

## summary() prints the block when the Satterthwaite df was computed
out <- capture.output(print(summary(fit, df = "satterthwaite"),
                            show.resids = FALSE))
## (both sleepstudy values exceed 1/2 and print as such)
stopifnot(any(grepl("^Fragility", out)),
          any(grepl("Subject.\\(Intercept\\) +eps\\* > 0.5 ", out)))
out0 <- capture.output(print(summary(fit, df = "none"), show.resids = FALSE))
stopifnot(!any(grepl("^Fragility", out0)))

cat("fragility.R: all checks passed\n")
