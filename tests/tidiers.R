## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## tidy(), glance() and augment() for rlmerMod: column names and types,
## agreement with the accessors, single / nested / crossed designs and a
## diag() term.
require(robustlmm)
if (!requireNamespace("generics", quietly = TRUE)) {
    cat("tidiers.R: package 'generics' not available, skipped\n")
    quit()
}
tidy <- generics::tidy; glance <- generics::glance
augment <- generics::augment

chk <- function(a, b, tol = 1e-10)
    stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
haveTibble <- requireNamespace("tibble", quietly = TRUE)
isTbl <- function(x) if (haveTibble) inherits(x, "tbl_df") else
    identical(class(x), "data.frame")

fm <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy)

## ---- tidy: fixed effects ---------------------------------------------------
## the default follows summary(): no df for a fit of this size until the
## influence function is cached
t0 <- tidy(fm, effects = "fixed")
stopifnot(identical(is.na(t0$df),
                    rep(!"df" %in% colnames(coef(summary(fm))), 2)))
td <- tidy(fm, df = "satterthwaite")
stopifnot(isTbl(td),
          identical(names(td), c("effect", "group", "term", "estimate",
                                 "std.error", "statistic", "df", "p.value")),
          is.character(td$effect), is.character(td$group),
          is.character(td$term), is.numeric(td$estimate),
          identical(td$effect, rep(c("fixed", "ran_pars"), c(2, 4))))
fx <- td[td$effect == "fixed", ]
co <- coef(summary(fm, df = "satterthwaite"))
chk(fx$estimate, fixef(fm))
stopifnot(identical(fx$term, names(fixef(fm))), all(is.na(fx$group)))
chk(fx$std.error, co[, "Std. Error"])
chk(fx$statistic, co[, "t value"])
chk(fx$df, co[, "df"])
chk(fx$p.value, co[, "Pr(>|t|)"])

## without the degrees of freedom: NA, as in summary(df = "none")
tn <- tidy(fm, effects = "fixed", df = "none")
stopifnot(all(is.na(tn$df)), all(is.na(tn$p.value)), nrow(tn) == 2)
chk(tn$statistic, co[, "t value"])

## confidence limits come from confint(), with the reference of the
## p-values
tc <- tidy(fm, effects = "fixed", conf.int = TRUE, conf.level = 0.9)
ci <- confint(fm, level = 0.9, df = "satterthwaite")
chk(cbind(tc$conf.low, tc$conf.high), ci)
tcn <- tidy(fm, effects = "fixed", conf.int = TRUE, df = "none")
chk(cbind(tcn$conf.low, tcn$conf.high), confint(fm))

## sandwich: CR2 standard errors and Bell-McCaffrey df, as confint()
ts <- tidy(fm, effects = "fixed", conf.int = TRUE, vcov_type = "sandwich")
cis <- confint(fm, vcov_type = "sandwich")
chk(cbind(ts$conf.low, ts$conf.high), cis)
chk(ts$df, attr(cis, "df"))
chk(ts$std.error,
    sqrt(diag(as.matrix(vcov_sandwich(fm, correction = "CR2")))))

## ---- tidy: random-effect parameters ---------------------------------------
rp <- td[td$effect == "ran_pars", ]
vc <- as.data.frame(VarCorr(fm))
stopifnot(identical(rp$term, c("sd__(Intercept)", "sd__Days",
                               "cor__(Intercept).Days", "sd__Observation")),
          identical(rp$group, c("Subject", "Subject", "Subject", "Residual")),
          all(is.na(rp$std.error)), all(is.na(rp$p.value)))
chk(rp$estimate, vc$sdcor)
chk(rp$estimate[4], sigma(fm))
stopifnot(nrow(tidy(fm, effects = "ran_pars")) == 4,
          identical(names(tidy(fm, effects = "ran_pars", conf.int = TRUE)),
                    names(tidy(fm, conf.int = TRUE))),
          identical(tail(names(tidy(fm, conf.int = TRUE)), 2),
                    c("conf.low", "conf.high")))

## a diag() term has no correlation row
fd <- rlmer(Reaction ~ Days + diag(Days | Subject), sleepstudy)
rd <- tidy(fd, effects = "ran_pars")
stopifnot(identical(rd$term, c("sd__(Intercept)", "sd__Days",
                               "sd__Observation")))
chk(rd$estimate[1:2], attr(VarCorr(fd)[[1]], "stddev"))

## ---- glance -----------------------------------------------------------------
gl <- glance(fm)
stopifnot(isTbl(gl), nrow(gl) == 1,
          identical(names(gl), c("nobs", "sigma", "converged",
                                 "ngrps.Subject", "nobs.downweighted",
                                 "ngrps.downweighted.Subject")),
          gl$nobs == 180, gl$ngrps.Subject == 18, isTRUE(gl$converged),
          gl$nobs.downweighted == sum(getME(fm, "w_e") < 0.5),
          gl$ngrps.downweighted.Subject ==
              sum(getME(fm, "w_b")$Subject[[1]] < 0.5),
          !any(c("logLik", "AIC", "BIC") %in% names(gl)))
chk(gl$sigma, sigma(fm))
stopifnot(glance(fm, threshold = 1)$nobs.downweighted ==
              sum(getME(fm, "w_e") < 1))
## an unconverged fit says so
fu <- suppressWarnings(rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                             max.iter = 2))
stopifnot(!glance(fu)$converged)
## an object that was not fitted has no convergence status
stopifnot(is.na(glance(rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                             doFit = FALSE))$converged))

## ---- augment ----------------------------------------------------------------
au <- augment(fm)
stopifnot(isTbl(au), nrow(au) == 180,
          all(c("Reaction", "Days", "Subject", ".fitted", ".resid", ".wt_e",
                ".wt_b", ".hat") %in% names(au)),
          !any(c(".cooksd", ".se.fit") %in% names(au)))
chk(au$.fitted, fitted(fm))
chk(au$.resid, residuals(fm))
chk(au$.wt_e, getME(fm, "w_e"))
chk(au$.hat, hatvalues(fm))
## the cluster weight of each observation
wb <- getME(fm, "w_b")$Subject
chk(au$.wt_b, pmin(wb[[1]], wb[[2]])[as.integer(sleepstudy$Subject)])

au2 <- augment(fm, data = sleepstudy, se_fit = TRUE, influence = TRUE)
chk(au2$.cooksd, cooks.distance(fm))
chk(au2$.se.fit, suppressWarnings(predict(fm, interval = "confidence"))$se)

## new data: fitted values only
nd <- data.frame(Days = c(0, 5), Subject = c("308", "309"))
an <- suppressWarnings(augment(fm, newdata = nd, se_fit = TRUE))
stopifnot(identical(names(an), c("Days", "Subject", ".fitted", ".se.fit")))
chk(an$.fitted, predict(fm, newdata = nd))

## rows dropped for missing values are padded when the full data is given
ssNA <- sleepstudy; ssNA$Reaction[c(3, 50)] <- NA
fna <- rlmer(Reaction ~ Days + (1 | Subject), ssNA)
ana <- augment(fna, data = ssNA)
stopifnot(nrow(ana) == 180, all(is.na(ana$.fitted[c(3, 50)])),
          !anyNA(ana$.fitted[-c(3, 50)]), nrow(augment(fna)) == 178)
chk(ana$.wt_e[-c(3, 50)], getME(fna, "w_e"))
stopifnot(inherits(tryCatch(augment(fna, data = ssNA[1:20, ]),
                            error = function(e) e), "error"))
## na.exclude: predict() pads once, and the standard errors align
fex <- rlmer(Reaction ~ Days + (1 | Subject), ssNA, na.action = na.exclude)
pe <- predict(fex)
stopifnot(length(pe) == 180, identical(unname(which(is.na(pe))), c(3L, 50L)))
chk(pe[-c(3, 50)], predict(fna))
aex <- suppressWarnings(augment(fex, data = ssNA, se_fit = TRUE))
stopifnot(nrow(aex) == 180, identical(unname(which(is.na(aex$.se.fit))), c(3L, 50L)))

## ---- crossed and nested designs: one .wt_b column per grouping factor ------
fc <- rlmer(diameter ~ 1 + (1 | plate) + (1 | sample), Penicillin)
ac <- augment(fc)
stopifnot(all(c(".wt_b.plate", ".wt_b.sample") %in% names(ac)),
          !".wt_b" %in% names(ac))
wc <- getME(fc, "w_b")
chk(ac$.wt_b.plate, wc$plate[[1]][as.integer(Penicillin$plate)])
chk(ac$.wt_b.sample, wc$sample[[1]][as.integer(Penicillin$sample)])
gc <- glance(fc)
stopifnot(gc$ngrps.plate == 24, gc$ngrps.sample == 6,
          all(c("ngrps.downweighted.plate", "ngrps.downweighted.sample") %in%
              names(gc)))
tcr <- tidy(fc)
stopifnot(identical(tcr$term, c("(Intercept)", "sd__(Intercept)",
                                "sd__(Intercept)", "sd__Observation")),
          identical(tcr$group, c(NA, "plate", "sample", "Residual")))
chk(tcr$estimate[2:4], as.data.frame(VarCorr(fc))$sdcor)

set.seed(5)
d <- expand.grid(rep = 1:4, sub = 1:3, top = 1:8)
d$sub <- interaction(d$top, d$sub)
d$y <- rnorm(8)[d$top] + rnorm(24, sd = 0.7)[as.integer(d$sub)] +
    rnorm(nrow(d))
d$top <- factor(d$top)
fn <- rlmer(y ~ 1 + (1 | top) + (1 | sub), d)
an <- augment(fn)
wn <- getME(fn, "w_b")
chk(an$.wt_b.top, wn$top[[1]][match(d$top, rownames(wn$top))])
chk(an$.wt_b.sub, wn$sub[[1]][match(d$sub, rownames(wn$sub))])
stopifnot(nrow(tidy(fn)) == 4, glance(fn)$ngrps.top == 8)

## two terms for the same grouping factor: the smaller weight per subject
f2 <- rlmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
            sleepstudy)
a2 <- augment(f2)
w2 <- getME(f2, "w_b")
chk(a2$.wt_b, pmin(w2[[1]][[1]], w2[[1]][[2]])[as.integer(sleepstudy$Subject)])
stopifnot(".wt_b" %in% names(a2))
cat("tidiers.R: all checks passed\n")
