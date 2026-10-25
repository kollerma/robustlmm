## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## rlmerControl(): control object, clash with loose arguments, options for
## one fit, retry on non-convergence; early argument checks.
require(robustlmm)

chk <- function(a, b, tol = 1e-10)
    stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
err <- function(expr, pattern = ".") {
    m <- tryCatch({expr; NA_character_}, error = function(e) conditionMessage(e))
    !is.na(m) && grepl(pattern, m)
}
par <- function(f) c(fixef(f), sigma(f), getME(f, "theta"))
warns <- function(expr) {
    w <- character(0)
    val <- withCallingHandlers(expr, warning = function(c) {
        w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
    list(value = val, warnings = w)
}

## ---- the control object ----------------------------------------------------
ctrl <- rlmerControl()
stopifnot(inherits(ctrl, "rlmerControl"), ctrl$rel.tol == 1e-8,
          is.null(ctrl$max.iter), !ctrl$retry,
          identical(ctrl$retry.steps, c("iterations", "DASvar", "robust")),
          length(capture.output(print(ctrl))) > 3)
stopifnot(err(rlmerControl(rel.tol = -1), "rel.tol"),
          err(rlmerControl(max.iter = 1.5), "max.iter"),
          err(rlmerControl(retry = NA), "retry"),
          err(rlmerControl(retry.steps = "magic"), "retry.steps"),
          err(rlmerControl(retry.max.iter.factor = 0.5), "factor"),
          err(rlmerControl(dastau.mc = "yes"), "dastau.mc"),
          err(rlmerControl(lmerControl = list()), "lmerControl"))

## ---- defaults are unchanged ---------------------------------------------------
form <- Reaction ~ Days + (Days | Subject)
fm <- rlmer(form, sleepstudy)
fc <- rlmer(form, sleepstudy, control = rlmerControl())
stopifnot(identical(par(fm), par(fc)))

## ---- loose arguments and control ------------------------------------------------
## the same setting through either route gives the same fit
u1 <- warns(rlmer(form, sleepstudy, max.iter = 3, rel.tol = 1e-6))
u2 <- warns(rlmer(form, sleepstudy,
                  control = rlmerControl(max.iter = 3, rel.tol = 1e-6)))
stopifnot(identical(par(u1$value), par(u2$value)),
          identical(u1$warnings, u2$warnings),
          u2$value@optinfo$conv$opt == 1)
## a loose argument combines with other settings in control
u3 <- warns(rlmer(form, sleepstudy, max.iter = 3,
                  control = rlmerControl(rel.tol = 1e-6)))
stopifnot(identical(par(u1$value), par(u3$value)))
## the same setting in both places is an error
stopifnot(err(rlmer(form, sleepstudy, max.iter = 3,
                    control = rlmerControl(max.iter = 5)), "max.iter"),
          err(rlmer(form, sleepstudy, rel.tol = 1e-6,
                    control = rlmerControl(rel.tol = 1e-6)), "rel.tol"),
          err(rlmer(form, sleepstudy, control = list(max.iter = 3)),
              "rlmerControl"))
## an lmerControl object goes to the classical start, as it did through `...`
fl <- rlmer(form, sleepstudy, control = lmerControl(optimizer = "bobyqa"))
fl2 <- rlmer(form, sleepstudy,
             control = rlmerControl(lmerControl = lmerControl(optimizer = "bobyqa")))
chk(par(fl), par(fm), tol = 1e-4)
stopifnot(identical(par(fl), par(fl2)))
## update() keeps the control argument of the call
stopifnot(!is.null(getCall(u2$value)$control))

## ---- Monte-Carlo DAS-tau settings hold for one fit only -----------------------
set.seed(11)
dd <- data.frame(g = factor(rep(1:25, each = 6)),
                 t = factor(rep(1:3, 50)))
dd$y <- rnorm(25)[dd$g] + rnorm(150)
form3 <- y ~ t + (0 + t | g)                     # a block of dimension 3
old <- options(robustlmm.dastau.mc = NULL)
f.var <- suppressWarnings(rlmer(form3, dd))       # falls back to DASvar
f.mc <- suppressMessages(suppressWarnings(
    rlmer(form3, dd, control = rlmerControl(dastau.mc = TRUE,
                                            dasmc.nsim = 2000))))
stopifnot(f.var@method == "DASvar", f.mc@method == "DAStau",
          is.null(getOption("robustlmm.dastau.mc")))   # option restored
## the same fit through the options
options(robustlmm.dastau.mc = TRUE, robustlmm.dasmc.nsim = 2000)
f.op <- suppressMessages(suppressWarnings(rlmer(form3, dd)))
options(robustlmm.dastau.mc = NULL, robustlmm.dasmc.nsim = NULL)
stopifnot(identical(par(f.mc), par(f.op)))
## control overrides the option
options(robustlmm.dastau.mc = TRUE)
f.off <- suppressWarnings(rlmer(form3, dd,
                                control = rlmerControl(dastau.mc = FALSE)))
options(old)
stopifnot(f.off@method == "DASvar", identical(par(f.off), par(f.var)))

## ---- retry on non-convergence ----------------------------------------------------
## a constructed hard case: too few iterations
hard <- function(...) rlmer(form, sleepstudy, max.iter = 4,
                            control = rlmerControl(...))
h0 <- warns(hard())
stopifnot(h0$value@optinfo$conv$opt == 1, length(h0$warnings) >= 1,
          is.null(h0$value@optinfo$retry))
## more iterations, continued from the unconverged estimate
h1 <- warns(hard(retry = TRUE))
stopifnot(h1$value@optinfo$conv$opt == 0, length(h1$warnings) == 0,
          h1$value@optinfo$retry$attempt == 1,
          h1$value@optinfo$retry$step == "iterations")
chk(par(h1$value), par(fm), tol = 1e-5)
for (out in list(capture.output(print(h1$value)),
                 capture.output(print(summary(h1$value, df = "none")))))
    stopifnot(any(grepl("Convergence code: 0", out, fixed = TRUE)),
              any(grepl("converged on retry attempt 1 of 3 (more iterations)",
                        out, fixed = TRUE)))
## the returned fit carries the original call
stopifnot(identical(getCall(h1$value)$max.iter, 4))
## each fallback on its own; all reach the solution
for (st in c("DASvar", "robust")) {
    h <- warns(rlmer(form, sleepstudy, max.iter = 4,
                     control = rlmerControl(retry = TRUE, retry.steps = st,
                                            retry.max.iter.factor = 20)))
    stopifnot(h$value@optinfo$conv$opt == 0, length(h$warnings) == 0,
              h$value@optinfo$retry$step == st,
              h$value@optinfo$retry$attempt == 1,
              h$value@method == "DAStau")
    chk(par(h$value), par(fm), tol = 1e-4)
}
## the chain moves on when a fallback does not apply or fails: a DASvar
## fit has no DASvar fallback
h2 <- warns(rlmer(form, sleepstudy, max.iter = 4, method = "DASvar",
                  control = rlmerControl(retry = TRUE,
                                         retry.steps = c("DASvar", "iterations"))))
stopifnot(h2$value@optinfo$conv$opt == 0,
          h2$value@optinfo$retry$step == "iterations",
          h2$value@optinfo$retry$attempt == 2)
## max.iter given in control instead of as an argument
h4 <- warns(rlmer(form, sleepstudy,
                  control = rlmerControl(max.iter = 4, retry = TRUE,
                                         retry.steps = "robust",
                                         retry.max.iter.factor = 20)))
stopifnot(h4$value@optinfo$conv$opt == 0,
          h4$value@optinfo$retry$step == "robust")
## nothing converges: the first fit and its warnings are returned
h3 <- warns(rlmer(form, sleepstudy, max.iter = 2,
                  control = rlmerControl(retry = TRUE,
                                         retry.steps = "iterations",
                                         retry.max.iter.factor = 1)))
stopifnot(h3$value@optinfo$conv$opt == 1,
          any(grepl("did not converge", h3$warnings)),
          is.na(h3$value@optinfo$retry$attempt),
          any(grepl("no retry attempt converged",
                    capture.output(print(h3$value)))))
## a converged fit is not retried
stopifnot(is.null(rlmer(form, sleepstudy,
                        control = rlmerControl(retry = TRUE))@optinfo$retry))

## the "tried" list names only the fallbacks that applied
h5 <- warns(rlmer(form, sleepstudy, max.iter = 2, method = "DASvar",
                  control = rlmerControl(retry = TRUE, retry.steps = "DASvar")))
stopifnot(grepl("none applicable", h5$value@optinfo$retry$message))

## ---- arguments that keep working ---------------------------------------------------
## a logical verbose, and update() with a changed control
stopifnot(identical(par(rlmer(form, sleepstudy, verbose = FALSE)), par(fm)))
fu <- update(fc, control = rlmerControl(rel.tol = 1e-6))
stopifnot(fu@optinfo$rel.tol == 1e-6, fc@optinfo$rel.tol == 1e-8)

## ---- early argument checks ---------------------------------------------------------
stopifnot(err(rlmer(form, sleepstudy, rel.tol = 0), "rel.tol"),
          err(rlmer(form, sleepstudy, max.iter = -1), "max.iter"),
          err(rlmer(form, sleepstudy, max.iter = 2.5), "max.iter"),
          err(rlmer(form, sleepstudy, doFit = "no"), "doFit"),
          err(rlmer(form, sleepstudy, verbose = "a"), "verbose"),
          err(vcov(fm, cluster = "Subject"), "sandwich"),
          err(vcov(fm, correction = "CR2"), "sandwich"),
          err(confint(fm, level = 1.2), "level"),
          err(confint(fm, cluster = "Subject"), "sandwich"),
          err(confint(fm, method = "boot", vcov_type = "sandwich"),
              "vcov_type"),
          err(confint(fm, method = "boot", df = "satterthwaite"), "df"),
          err(confint(fm, method = "boot", nsim = 0), "nsim"),
          err(anova(fm, cluster = "Subject"), "sandwich"),
          err(anova(fm, nsim = 0), "nsim"),
          err(anova(fm, fm, test = "score", null = "robust"), "null"),
          err(predict(fm, interval = "confidence", level = 0), "level"))
cat("control.R: all checks passed\n")
