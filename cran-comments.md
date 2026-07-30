# robustlmm 3.5.0-2

## Summary

Patch release. Two bug fixes to behaviour that shipped in earlier CRAN
versions; no new features, no interface changes, and no change to the
defaults of unaffected functionality.

* **`init = "ransac"` estimated the model on the RANSAC subsample
  instead of on the user's data** (affects 3.4-5 and 3.5.0-1). An
  `lmerMod` passed as `init` supplies the model as well as the starting
  values, and the RANSAC start is by construction a fit on a subsample,
  so the returned object was an estimate on roughly half the
  observations. On the regression test's design, 510 of 1000
  observations and 27 of 50 groups. Four entry points were affected
  (`init = "ransac"`, `rlmer_ransac()`, `fitDatasets_rlmer_ransac()`
  and `fitDatasets_rlmer_ransac_bisq()`); all are fixed and covered by
  a new regression test. Users who ran any of them should refit; NEWS.md
  says so explicitly.

* **`calcTau()` refused to solve its fixed point from a starting value
  above 1.** Defensive: no fit produced by the current code reaches
  that condition (largest starting value observed over a randomised
  search is 0.99948), so this is a no-op on current behaviour. It is
  included because it is the mechanism by which the weighted-`tau` bug
  fixed in 3.5.0 was able to reach `method = "DAStau"`.

Documentation: `?rlmer` and the 3.5.0 vignette have been reworded to
describe the RANSAC start as what it is -- a way to raise the breakdown
point of the starting value, and so to make redescending
psi-functions usable -- rather than by its effect on any particular
failure mode. Simulation figures quoted in the vignette were re-measured
under this version.

## Test environments

* Local: macOS (arm64), R 4.6.0 Patched.
* Linux containers (rocker), R devel / release / oldrel.
* win-builder, R devel and R release.

As in 3.5.0-1, the tests shipped in the tarball are a trimmed subset so
that the CRAN check stays well inside the time limit; the full suite
runs in CI. The trimming disables test files only and does not alter
any code path.

## R CMD check results

0 errors | 0 warnings.

Any NOTE about `tests/getME.Rout.save` differences appears only when
checking against a development version of lme4 and does not occur with
the CRAN release of lme4.

## Reverse dependencies

CRAN reverse dependencies: confintROB, effects, insight,
marginaleffects, misty.

The changed code paths are reached only via `init = "ransac"`,
`rlmer_ransac()` and the two `fitDatasets_rlmer_ransac*` wrappers, and
via `calcTau()`'s starting value. None of the five reverse dependencies
calls any of them: they use `rlmer()` at its defaults together with the
extractor and inference methods, none of which changed. The `calcTau()`
change is additionally a verified no-op on current behaviour.
