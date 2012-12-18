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

This package builds on the *Eigen and S4* implementation of `lme4`. We
require at least version 0.99999911-0. At the time of writing this version
is only available on R-forge, not on CRAN. To install directly from
R-forge, use:

    install.packages("lme4", repos=c("http://R-Forge.R-project.org",
                                     getOption("repos")))

Once you have installed the newest version of `lme4`, you can install
`robustlmm` directly from github using `install_github` of the R-package
`devtools`:

    install.packages("devtools") ## if not already installed
    require(devtools)
    install_github("robustlmm", "kollerma")
    require(robustlmm)
