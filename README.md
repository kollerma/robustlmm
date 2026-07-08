Robust linear mixed effects models
==================================

[![R-CMD-check](https://github.com/kollerma/robustlmm/actions/workflows/check-standard.yaml/badge.svg?branch=master)](https://github.com/kollerma/robustlmm/actions/workflows/check-standard.yaml)
[![cran version](https://www.r-pkg.org/badges/version/robustlmm)](https://CRAN.R-project.org/package=robustlmm)
[![downloads](https://cranlogs.r-pkg.org/badges/robustlmm)](https://cranlogs.r-pkg.org/badges/robustlmm)
[![total downloads](https://cranlogs.r-pkg.org/badges/grand-total/robustlmm)](https://cranlogs.r-pkg.org/badges/grand-total/robustlmm)
[![Research software impact](http://depsy.org/api/package/cran/robustlmm/badge.svg)](http://depsy.org/package/r/robustlmm)

The R-package `robustlmm` provides functions for estimating linear mixed
effects models in a robust way.

The main workhorse is the function `rlmer`; it is implemented as direct
robust analogue of the popular `lmer` function of the `lme4` package. The
two functions have similar abilities and limitations. A wide range of data
structures can be modeled: mixed effects models with hierarchical as well
as complete or partially crossed random effects structures are
possible. While the `lmer` function is optimized to handle large datasets
efficiently, the computations employed in the `rlmer` function are more
complex and for this reason also more expensive to compute. The two
functions have the same limitations in the support of different random
effect and residual error covariance structures. Both support only diagonal
and unstructured random effect covariance structures.

The `robustlmm` package implements most of the analysis tool chain as is
customary in R. The usual functions such as `summary`, `coef`, `resid`,
etc. are provided as long as they are applicable for this type of models
(see `rlmerMod-class` for a full list). The functions are designed to be as
similar as possible to the ones in the `lme4` package to make switching
between the two packages easy.

Inference is supported via:

- `vcov(fit)` (linearised, the lme4-inherited default) and
  `vcov(fit, type = "sandwich")` (robust cluster sandwich;
  see `?vcov_sandwich`).
- `confint(fit)` returns the closed-form Wald interval, optionally with
  `vcov_type = "sandwich"`. `confint(fit, method = "boot")` and
  `method = "BCa"` delegate to the peer-reviewed `confintROB` package
  (Mason, Cantoni & Ghisletta 2021, 2024), which is in `Suggests`.
- `anova(fit)` for a per-term Wald table; `anova(fit0, fit1)` for nested
  model comparison (Wald restriction for fixed effects, parametric-bootstrap
  quasi-deviance for variance-component tests on the PSD-cone boundary).
- `predict(fit, interval = "confidence" / "prediction")` with confidence
  / prediction intervals.
- `cooks.distance(fit)` for per-observation joint influence on
  $(\hat\beta, \hat\sigma, \hat\theta)$; `caseweightIF(fit)` for the
  case-weight influence function.
- `rlmer(formula, data, init = "ransac")` for a high-breakdown
  RANSAC start, useful when redescending psi-functions risk a phony
  local minimum.

An empirical evaluation of the inference methods (CI coverage, Type-I,
power, under contamination) lives in
`inst/simulationStudy/inferenceStudy_results/` on the development branch.
  
Installation
------------

This R-package is [available on
CRAN](https://CRAN.R-project.org/package=robustlmm). Install it
directly in R with the command

    install.packages("robustlmm")

This package requires `lme4` version at least `2.0-1` and other
packages. Make sure to install them as well.

You can also install the package directly from github:

    install.packages("devtools") ## if not already installed
    require(devtools)
    install_github("kollerma/robustlmm")
    require(robustlmm)
