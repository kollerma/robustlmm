## Precompute the numbers quoted in the section "A robust sensitivity
## analysis for an MMRM" of vignettes/robustlmm-3.5.0.Rnw, so that the
## vignette builds without the mmrm package (whose example data set
## fev_data is used and must not be copied into this package). Writes the
## small, shipped RDS vignettes/mmrm-sensitivity-numbers.rds.
##
## Needs robustlmm, mmrm, emmeans and generics. Re-run from the package
## root:
##   Rscript vignettes/precompute-mmrm.R
##
## Build-ignored (see .Rbuildignore).

suppressPackageStartupMessages({
    library(robustlmm); library(mmrm); library(emmeans)
})
d <- mmrm::fev_data
d$visit <- as.numeric(d$AVISIT) - 2.5          # centred visit number
form <- FEV1 ~ RACE + SEX + ARMCD * AVISIT

## ---- the fits ---------------------------------------------------------------
fitMmrm <- mmrm(update(form, . ~ . + us(AVISIT | USUBJID)), data = d)
fitLmer <- lmer(update(form, . ~ . + (1 | USUBJID)), data = d)
fitRob <- rlmer(update(form, . ~ . + (1 | USUBJID)), data = d)
fitRobSlope <- rlmer(update(form, . ~ . + (visit | USUBJID)), data = d,
                     method = "DASvar")
## the structured alternative is refused: one observation per subject and
## visit does not identify a visit-level covariance next to the residual
csError <- tryCatch({
    rlmer(update(form, . ~ . + cs(0 + AVISIT | USUBJID)), data = d); NA
}, error = function(e) conditionMessage(e))

## ---- treatment difference at the last visit ---------------------------------
last <- function(fit) {
    ct <- as.data.frame(suppressWarnings(
        contrast(emmeans(fit, ~ ARMCD | AVISIT), "revpairwise")))
    ct <- ct[ct$AVISIT == "VIS4", ]
    c(estimate = ct$estimate, se = ct$SE, df = ct$df, p = ct$p.value)
}
## the same contrast in terms of the coefficients, for contrast_test()
L <- c("ARMCDTRT" = 1, "ARMCDTRT:AVISITVIS4" = 1)
ctRow <- function(fit, ...) {
    ct <- suppressWarnings(contrast_test(fit, L, ...))
    c(estimate = ct$Estimate, se = ct$`Std. Error`, df = ct$df,
      p = ct$`Pr(>|t|)`)
}
contrasts <- rbind(
    mmrm = last(fitMmrm),
    lmer = last(fitLmer),
    rlmer = last(fitRob),
    rlmerSandwich = ctRow(fitRob, vcov_type = "sandwich"),
    slope = ctRow(fitRobSlope, ddf = "none"),
    slopeSandwich = ctRow(fitRobSlope, vcov_type = "sandwich"))
## emmeans and contrast_test() agree on the robust fit
stopifnot(isTRUE(all.equal(unname(last(fitRob)), unname(ctRow(fitRob)),
                           tolerance = 1e-5)))

## ---- robustness weights -------------------------------------------------------
wtSummary <- function(fit) {
    aug <- generics::augment(fit)
    low <- aug[order(aug$.wt_e), c("USUBJID", "AVISIT", "ARMCD", "FEV1",
                                   ".fitted", ".wt_e")][1:6, ]
    wb <- tapply(aug$.wt_b, droplevels(aug$USUBJID), min)
    list(nobs = nrow(aug), nsubj = length(wb),
         nLow = sum(aug$.wt_e < 0.5),
         belowByVisit = table(aug$AVISIT[aug$.wt_e < 0.8]),
         obsByVisit = table(aug$AVISIT),
         lowest = as.data.frame(low),
         nSubjBelow1 = sum(wb < 1), minWtB = min(wb),
         minWtBSubject = names(wb)[which.min(wb)])
}
weights <- list(rlmer = wtSummary(fitRob), slope = wtSummary(fitRobSlope))

## ---- influential subjects -------------------------------------------------------
## (the influence function of the random-slope fit takes several minutes)
cooksSummary <- function(fit) {
    cd <- suppressWarnings(suppressMessages(
        cooks.distance(fit, groups = TRUE)))
    list(top = sort(cd, decreasing = TRUE)[1:5], median = stats::median(cd),
         max = max(cd))
}
cooks <- list(rlmer = cooksSummary(fitRob), slope = cooksSummary(fitRobSlope))

## ---- variance by visit ----------------------------------------------------------
sdVisit <- function(V) sqrt(diag(V))
vcS <- VarCorr(fitRobSlope)$USUBJID
tt <- (1:4) - 2.5
variance <- list(
    mmrmSd = sdVisit(mmrm::VarCorr(fitMmrm)),
    lmerSd = c(subject = attr(VarCorr(fitLmer)$USUBJID, "stddev")[[1]],
               residual = sigma(fitLmer)),
    rlmerSd = c(subject = attr(VarCorr(fitRob)$USUBJID, "stddev")[[1]],
                residual = sigma(fitRob)),
    slopeSd = sqrt(vapply(tt, function(t)
        drop(c(1, t) %*% vcS %*% c(1, t)), 0) + sigma(fitRobSlope)^2))

numbers <- list(contrasts = contrasts, weights = weights, cooks = cooks,
                variance = variance, csError = csError,
                converged = c(rlmer = fitRob@optinfo$conv$opt == 0,
                              slope = fitRobSlope@optinfo$conv$opt == 0),
                versions = c(mmrm = as.character(packageVersion("mmrm")),
                             robustlmm = as.character(packageVersion("robustlmm"))))
saveRDS(numbers, "vignettes/mmrm-sensitivity-numbers.rds")
cat("wrote vignettes/mmrm-sensitivity-numbers.rds\n")
str(numbers, max.level = 2)
