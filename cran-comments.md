# robustlmm 3.5.0

## Summary

Major release. The full list of changes is in NEWS.md (which also
lists what already shipped in 3.4-4/3.4-5); highlights:

* New inference tools for `rlmer()` fits: `anova()` methods,
  Satterthwaite degrees of freedom in `summary()` and `emmeans`,
  and cluster-level influence diagnostics -- completing the
  inference layer begun in 3.4-4/3.4-5 (cluster-robust sandwich
  `vcov()`, Wald/bootstrap `confint()`).
* New estimation options (Mallows-type design weights, RANSAC-based
  initial estimator, structured covariance support for
  `diag`/`cs`/`ar1`). Features that are still experimental are clearly
  marked as such in the documentation and are off by default; all
  defaults of existing functionality are unchanged.

## Test environments

* Local: macOS (arm64), R release, lme4 development version.
* The full (slow) test suite runs in CI on GitHub Actions
  (macOS native plus rocker-based Linux containers); the tests shipped
  in the tarball are a trimmed subset so that the CRAN check stays
  within a few minutes. The trimming only disables tests, it does not
  alter any code paths.

## R CMD check results

0 errors | 0 warnings. Any NOTE about `tests/getME.Rout.save`
differences only appears when checking against a development version
of lme4 and does not occur with the CRAN release of lme4.

## Reverse dependencies

CRAN reverse dependencies: confintROB, effects, insight,
marginaleffects, misty. All five were checked against both the current
CRAN robustlmm (3.4-5) and this submission (R CMD check, without
vignettes): no new errors, warnings, or notes appeared with 3.5.0.
The one failure observed (an insight test unrelated to robustlmm,
involving panelr and lme4) is pre-existing and byte-identical under
both versions. In addition, all changes to pre-existing functionality
are additive and default-off, so the interfaces these packages use are
unchanged.
