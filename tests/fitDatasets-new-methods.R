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
