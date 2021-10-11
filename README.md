Robust linear mixed effects models
==================================

[![Build Status](https://travis-ci.org/kollerma/robustlmm.svg?branch=master)](https://travis-ci.org/kollerma/robustlmm)
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
  
Installation
------------

This R-package is [available on
CRAN](https://CRAN.R-project.org/package=robustlmm). Install it
directly in R with the command

    install.packages("robustlmm")

This package requires `lme4` version at least `1.1` and other
packages. Make sure to install them as well.

You can also install the package directly from github:

    install.packages("devtools") ## if not already installed
    require(devtools)
    install_github("kollerma/robustlmm")
    require(robustlmm)
