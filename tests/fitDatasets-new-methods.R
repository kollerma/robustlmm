## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Smoke tests for the new fitDatasets_rlmer_* wrappers:
##   _DAStau_bisq, _DAStau_sizeOBR, _ransac, _ransac_bisq.

require(robustlmm)

## --- Tiny diagonal one-way dataset ---------------------------------
set.seed(7)
ds_diag <- generateAnovaDatasets(numberOfDatasetsToGenerate = 2,
                                 numberOfLevelsInFixedFactor = 1,
                                 numberOfSubjects = 8,
                                 numberOfReplicates = 5)

## --- 1. bisq on diagonal ------------------------------------------
fits_bisq <- fitDatasets_rlmer_DAStau_bisq(ds_diag)
stopifnot(length(fits_bisq) == 2L)
stopifnot(all(vapply(fits_bisq,
                     function(x) methods::is(x, "rlmerMod"), logical(1))))
stopifnot(attr(fits_bisq[[1]], "label") ==
              "fitDatasets_rlmer_DAStau_bisq")

## --- 2. ransac on diagonal ----------------------------------------
set.seed(11)
fits_ransac <- fitDatasets_rlmer_ransac(ds_diag, K = 10L)
stopifnot(length(fits_ransac) == 2L)
stopifnot(all(vapply(fits_ransac,
                     function(x) methods::is(x, "rlmerMod"), logical(1))))
stopifnot(attr(fits_ransac[[1]], "label") ==
              "fitDatasets_rlmer_ransac")

## --- 3. Predefined tuning parameter table is populated --------------
labels <- c("fitDatasets_rlmer_DAStau_bisq",
            "fitDatasets_rlmer_DAStau_sizeOBR",
            "fitDatasets_rlmer_ransac",
            "fitDatasets_rlmer_ransac_bisq")
for (lbl in labels) {
    tp <- extractPredefinedTuningParameter(lbl, "rho.sigma.e")
    stopifnot(is.numeric(tp) && length(tp) == 1L)
}

## --- 4. Block-diagonal smoke test ---------------------------------
set.seed(101)
ds_block <- generateMixedEffectDatasets(
    numberOfDatasetsToGenerate = 1,
    prepareMixedEffectDataset(
        Reaction ~ Days + (Days | Subject),
        data = sleepstudy,
        groups = cbind(rep(1:10, each = 18), rep(1:18, 10)),
        varcov = list(
            Subject.Intercept. = tcrossprod(rep(1, 10), rep(1, 10)),
            Subject.Days.Intercept. = tcrossprod(rep(1, 10), 0:9) +
                                       tcrossprod(0:9, rep(1, 10)),
            Subject.Days = tcrossprod(0:9, 0:9)),
        lower = c(0, -Inf, 0)))

fits_bisq_b <- fitDatasets_rlmer_DAStau_bisq(ds_block)
stopifnot(methods::is(fits_bisq_b[[1]], "rlmerMod"))

fits_sizeOBR <- fitDatasets_rlmer_DAStau_sizeOBR(ds_block)
stopifnot(methods::is(fits_sizeOBR[[1]], "rlmerMod"))
stopifnot(isTRUE(fits_sizeOBR[[1]]@optinfo$size_obr))

set.seed(33)
fits_rb <- fitDatasets_rlmer_ransac_bisq(ds_block, K = 10L)
stopifnot(methods::is(fits_rb[[1]], "rlmerMod"))

## --- 5. Structured-covariance wrappers (cs, ar1), lme4 >= 2.0-0 ----
## Need lme4 >= 2.0-0; skip otherwise so the section is a no-op.
if (packageVersion("lme4") >= "2.0.0") {

    ## 5a. rewriteRandomEffectsStructure: formula surgery
    rw <- robustlmm:::rewriteRandomEffectsStructure
    stopifnot(
        identical(deparse(rw(y ~ 1 + (0 + f | s), "cs")),
                  "y ~ 1 + cs(0 + f | s)"),
        identical(deparse(rw(y ~ 1 + cs(0 + f | s), "ar1")),
                  "y ~ 1 + ar1(0 + f | s)"),
        ## NA unwraps to the plain (expr | group) form
        identical(deparse(rw(y ~ 1 + ar1(0 + f | s), NA)),
                  "y ~ 1 + (0 + f | s)"),
        ## fixed-effect part is preserved
        identical(deparse(rw(y ~ x + (Days | g), "cs")),
                  "y ~ x + cs(Days | g)"))
    ## more than one random-effects term is an error
    stopifnot(inherits(tryCatch(rw(y ~ (1 | a) + (1 | b), "cs"),
                                error = identity), "error"))

    ## 5b. generate cs-structured data through the framework, but
    ## prepare with an UNSTRUCTURED formula on purpose: the wrappers
    ## must impose the structure themselves.
    set.seed(7)
    nv <- 3L; nsubj <- 40L; rho <- 0.5
    Rcs <- matrix(rho, nv, nv); diag(Rcs) <- 1
    Lt <- t(chol(diag(c(2, 1.5, 1.2)) %*% Rcs %*% diag(c(2, 1.5, 1.2))))
    templ <- expand.grid(f = factor(letters[seq_len(nv)]), rep = 1:5,
                         subj = factor(seq_len(nsubj)))
    bb <- t(Lt %*% matrix(rnorm(nv * nsubj), nv, nsubj))
    templ$y <- 10 + bb[cbind(as.integer(templ$subj), as.integer(templ$f))] +
        rnorm(nrow(templ), sd = 0.8)
    prep_cs <- prepareMixedEffectDataset(y ~ 1 + (0 + f | subj), templ)
    stopifnot(length(prep_cs$trueTheta) == 6L)  # 3x3 Cholesky parametrisation
    set.seed(11)
    ds_str <- generateMixedEffectDatasets(2, prep_cs)

    ## 5c. fitDatasets_rlmer_cs enforces compound symmetry
    fits_cs <- fitDatasets_rlmer_cs(ds_str)
    stopifnot(length(fits_cs) == 2L,
              all(vapply(fits_cs, methods::is, NA, "rlmerMod")),
              attr(fits_cs[[1]], "label") == "fitDatasets_rlmer_cs")
    fcs <- fits_cs[[1]]
    stopifnot(identical(deparse(fcs@call$formula),
                        "y ~ 1 + cs(0 + f | subj)"))
    cm <- attr(VarCorr(fcs)$subj, "correlation")
    offdiag <- cm[upper.tri(cm)]                 # all equal (compound symmetric)
    stopifnot(length(offdiag) == 3L,
              max(abs(offdiag - offdiag[1])) < 1e-6,
              fcs@optinfo$conv$opt == 0)

    ## 5d. fitDatasets_rlmer_ar1 enforces the AR(1) lag structure
    fits_ar1 <- fitDatasets_rlmer_ar1(ds_str)
    stopifnot(attr(fits_ar1[[1]], "label") == "fitDatasets_rlmer_ar1")
    far1 <- fits_ar1[[1]]
    stopifnot(identical(deparse(far1@call$formula),
                        "y ~ 1 + ar1(0 + f | subj)"))
    cma <- unname(attr(VarCorr(far1)$subj, "correlation"))
    sda <- as.numeric(attr(VarCorr(far1)$subj, "stddev"))
    stopifnot(max(abs(sda - sda[1])) < 1e-6,              # homogeneous
              max(abs(cma[1, ] - cma[1, 2]^(0:(nv - 1)))) < 1e-6,  # rho^lag
              far1@optinfo$conv$opt == 0)

    ## 5e. generateRepeatedMeasuresDatasets recovers the requested structure
    set.seed(13)
    g_cs <- generateRepeatedMeasuresDatasets(
        1, numberOfSubjects = 120, numberOfVisits = 3, numberOfReplicates = 8,
        structure = "cs", marginalSd = c(2, 1.5, 1.2), correlation = 0.5,
        trueBeta = 10, trueSigma = 0.8)
    stopifnot(identical(deparse(g_cs$formula), "y ~ 1 + cs(0 + visit | subject)"),
              g_cs$numberOfDatasets == 1L)
    fg <- lmer(g_cs$formula, g_cs$generateData(1), REML = TRUE)
    vcg <- VarCorr(fg)$subject
    ## marginal sds and the common correlation are recovered (MC tolerance)
    stopifnot(max(abs(attr(vcg, "stddev") - c(2, 1.5, 1.2))) < 0.3,
              abs(attr(vcg, "correlation")[2, 1] - 0.5) < 0.1,
              abs(sigma(fg) - 0.8) < 0.1)
    ## diag structure => uncorrelated
    g_d <- generateRepeatedMeasuresDatasets(
        1, numberOfSubjects = 120, numberOfVisits = 3, numberOfReplicates = 8,
        structure = "diag", marginalSd = c(2, 1.5, 1.2))
    stopifnot(identical(deparse(g_d$formula), "y ~ 1 + diag(0 + visit | subject)"))

    cat("structured-covariance wrappers (cs, ar1) + generator passed.\n")
}
