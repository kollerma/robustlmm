## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Convergence problems of the fit are reported: raised as warnings by
## rlmer(), stored in optinfo and shown by print() and summary(). A
## converged fit stays silent.
require(robustlmm)

collect <- function(expr) {
    w <- character(0)
    val <- withCallingHandlers(expr, warning = function(c) {
        w <<- c(w, conditionMessage(c))
        invokeRestart("muffleWarning")
    })
    list(value = val, warnings = w)
}
notConv <- "iterations did not converge, returning unconverged estimate."

## ---- converged fits are silent, in both fitting loops -------------------
r1 <- collect(rlmer(Reaction ~ Days + (1 | Subject), sleepstudy))
r2 <- collect(rlmer(Reaction ~ Days + (Days | Subject), sleepstudy))
for (r in list(r1, r2)) {
    stopifnot(length(r$warnings) == 0,
              r$value@optinfo$conv$opt == 0,
              length(r$value@optinfo$warnings) == 0)
    out <- capture.output(print(r$value), print(summary(r$value)))
    stopifnot(!any(grepl("Convergence code", out)))
}

## ---- an unconverged fit warns once, stores and prints the problem -------
for (form in c(Reaction ~ Days + (1 | Subject),
               Reaction ~ Days + (Days | Subject))) {
    u <- collect(rlmer(form, sleepstudy, max.iter = 2))
    stopifnot(sum(u$warnings == notConv) == 1,          # raised, once
              !anyDuplicated(u$warnings),
              u$value@optinfo$conv$opt == 1,
              notConv %in% unlist(u$value@optinfo$warnings))
    for (out in list(capture.output(print(u$value)),
                     capture.output(print(summary(u$value)))))
        stopifnot(any(grepl("Convergence code: 1 (iterations did not converge)",
                            out, fixed = TRUE)),
                  any(grepl(notConv, out, fixed = TRUE)))
}

## ---- warnings of intermediate iterations are not reported ---------------
## the theta iterations of the first outer steps hit max.iter here, the
## last ones do not, and the fit converges
f <- collect(rlmer(Reaction ~ Days + (1 | Subject), sleepstudy, max.iter = 12))
stopifnot(f$value@optinfo$conv$opt == 0, length(f$warnings) == 0,
          f$value@optinfo$nMuffled > 0)

## ---- a failed inner iteration of the last outer step is reported --------
## the outer loop converges, its last theta iteration does not: code 3
u <- collect(rlmer(Reaction ~ Days + (1 | Subject), sleepstudy, max.iter = 9))
stopifnot(u$value@optinfo$conv$opt == 3,
          any(grepl("theta iterations did not converge", u$warnings)),
          !notConv %in% u$warnings,
          any(grepl("Convergence code: 3", capture.output(print(u$value)),
                    fixed = TRUE)))

## ---- max.iter = 0 evaluates the start: reported as not converged ---------
for (form in c(Reaction ~ Days + (1 | Subject),
               Reaction ~ Days + (Days | Subject))) {
    z <- collect(rlmer(form, sleepstudy, max.iter = 0))
    stopifnot(z$value@optinfo$conv$opt == 1, notConv %in% z$warnings)
}

## ---- a fit stored by an earlier version prints no convergence block ------
## (its optinfo$warnings held the warnings of every intermediate iteration)
old <- f$value
old@optinfo$nMuffled <- NULL
old@optinfo$warnings <- list("theta iterations did not converge.")
stopifnot(!any(grepl("Convergence code", capture.output(print(old)))))

## ---- indexed messages are collapsed --------------------------------------
cw <- robustlmm:::.collapseWarnings(
    c("calculation of tau[1] reached max.it", "other",
      "calculation of tau[7] reached max.it", "other"))
stopifnot(identical(unname(cw),
                    c("calculation of tau reached max.it (2 times)", "other")))

## ---- callers that refit in a loop stay quiet ------------------------------
## the parametric bootstrap of anova() muffles its refits
fm0 <- rlmer(Reaction ~ Days + (1 | Subject), sleepstudy)
fm1 <- rlmer(Reaction ~ Days + (Days || Subject), sleepstudy)
b <- collect(anova(fm0, fm1, test = "boot", nsim = 10, seed = 1))
stopifnot(!any(grepl("did not converge", b$warnings)))
