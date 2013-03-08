Robust linear mixed models
==========================

This is an R-package for fitting linear mixed effects models in a robust
manner. The method is based on the robustification of the scoring equations
and an application of the Design Adaptive Scale approach. More details
forthcoming.

Usage
-----

The main function is `rlmer`, which can be called just like `lmer` of the
`lme4` package. Note that the `family` argument is missing, only
linear mixed models are supported. See `?rlmer` for examples.

Installation
------------

You can install `robustlmm` directly from github using `install_github` of
the R-package `devtools` (on Windows, make sure to have `Rtools` installed):

    install.packages("devtools") ## if not already installed
    require(devtools)
    install_github("robustlmm", "kollerma")
    require(robustlmm)
