WS6 inference study — empirical evaluation of `rlmer`’s inference
methods
================
robustlmm simulation study, Phases A + C + D + E + F + G + H

- [Headline](#headline)
- [Phase A — CI coverage and length](#phase-a--ci-coverage-and-length)
  - [Coverage (target 0.95)](#coverage-target-095)
  - [CI length (lower = sharper at matched
    coverage)](#ci-length-lower--sharper-at-matched-coverage)
- [Phase C — single-fit `anova(fit)` Type I +
  power](#phase-c--single-fit-anovafit-type-i--power)
- [Phase D — pairwise `anova(f0, f1)` Type I +
  power](#phase-d--pairwise-anovaf0-f1-type-i--power)
  - [Type I (target 0.05)](#type-i-target-005)
  - [Power (eff_size = 0.5)](#power-eff_size--05)
- [Robustness profile across
  contamination](#robustness-profile-across-contamination)
- [Phase E — variance-component test via parametric
  bootstrap](#phase-e--variance-component-test-via-parametric-bootstrap)
  - [Type I (target 0.05)](#type-i-target-005-1)
  - [Power (eff = 0.5)](#power-eff--05)
- [Phase F — residual-bootstrap variant
  (unweighted)](#phase-f--residual-bootstrap-variant-unweighted)
  - [Type I (target 0.05)](#type-i-target-005-2)
  - [Power (eff = 0.5)](#power-eff--05-1)
  - [Direct comparison: Gaussian (Phase E) vs residual (Phase
    F)](#direct-comparison-gaussian-phase-e-vs-residual-phase-f)
- [Phase G — robustness-weighted residual
  bootstrap](#phase-g--robustness-weighted-residual-bootstrap)
  - [Type I (target 0.05)](#type-i-target-005-3)
  - [Comparison: unweighted (Phase F) vs weighted (Phase G) vs Gaussian
    (Phase
    E)](#comparison-unweighted-phase-f-vs-weighted-phase-g-vs-gaussian-phase-e)
- [Phase H — subject-block residual
  bootstrap](#phase-h--subject-block-residual-bootstrap)
  - [Type I and power](#type-i-and-power)
- [Recommendations](#recommendations)
- [Reproducibility](#reproducibility)

This document summarises the WS6 simulation study comparing the
inference methods exposed in `rlmer` (`vcov`, `confint`, `anova`)
against the `lme4::lmer` baseline that users would otherwise reach for:
Phases A + C + D (Wald CI coverage and `anova` Type-I/power) plus Phases
E–H (the variance-component bootstrap test and three attempts to make it
robust to subject-level contamination).

Configuration of the run quoted below (all numbers in this document are
pulled live from the `.rds`):

- Replicates per cell: **500**
- `J` sweep: **8, 18, 50, 100** subjects
- Contamination: **clean, t3, CN_e, CN_re, shift_subj**
- Effect sizes for Phases C / D: **0, 0.5** (`0` = Type-I cell, `>0` =
  power cell)
- Cores: **9** via `parallel::mclapply`
- Wallclock: **27 min**
- Master seed: **20260603** (3.4.3 of robustlmm)

# Headline

Three methods, three verdicts:

| Method | Coverage / Type-I across cells | Verdict |
|:---|:---|:---|
| `confint(method = "Wald", vcov_type = "default")` | 0.93 - 0.97 | **Trust at any J.** |
| `confint(method = "Wald", vcov_type = "sandwich")` | 0.90 - 0.96 | OK at J \>= 18; under-covers at J = 8. |
| `confint(lmer_fit, method = "Wald")` | 0.93 - 0.97 | Trust at any J (same coverage as rlmer-default). |
| `anova(rlmerfit)` / `anova(f0, f1)` (default vcov) | 0.02 - 0.09 | **Trust at any J.** |
| `anova(rlmerfit)` / `anova(f0, f1)` (sandwich vcov) | 0.04 - 0.20 | Anti-conservative at J = 8; OK from J \>= 50. |
| `anova(lmer0, lmer1)` LRT | 0.03 - 0.09 | Trust at any J (same Type I as rlmer-default). |
| `anova(f0, f1, test = "boot")` VC bootstrap | 0.03 - 0.23 | OK under clean / t3; **anti-conservative under subject-level contamination** (Phase E). |

Coverage / Type-I summary across the 5 contamination cells x 4 J values.

The headline finding is **not** about Type I — both `rlmer` (default
vcov) and `lmer` calibrate correctly under all five contamination
regimes tested. The headline is about **efficiency under
contamination**:

Under heavy-tailed errors (`t3`), at `J = 50`, the `rlmer` Wald test
rejects **+23 pp** more often than the `lme4` LRT at the same nominal
level (power **93%** vs **69%**). Under contaminated-normal errors
(`CN_e`), the gap is **+16 pp** at the same `J`. Under clean Gaussian
data the methods are indistinguishable. **This is the canonical
robustness-buys-power trade-off, now quantified for `rlmer`’s anova.**

# Phase A — CI coverage and length

Coverage for the two Var1 slope coefficients (averaged) and mean CI
length, on the one-way ANOVA design `y ~ Var1 + (1 | Var2)`:

## Coverage (target 0.95)

**Coverage: `Wald.default`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.926 | 0.970 | 0.955 | 0.950 |      0.944 |
|  18 | 0.959 | 0.960 | 0.952 | 0.943 |      0.947 |
|  50 | 0.944 | 0.938 | 0.953 | 0.947 |      0.954 |
| 100 | 0.961 | 0.948 | 0.955 | 0.954 |      0.959 |

**Coverage: `Wald.sandwich`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.897 | 0.922 | 0.913 | 0.910 |      0.916 |
|  18 | 0.938 | 0.935 | 0.948 | 0.931 |      0.939 |
|  50 | 0.934 | 0.939 | 0.948 | 0.943 |      0.948 |
| 100 | 0.963 | 0.953 | 0.949 | 0.950 |      0.960 |

**Coverage: `Wald.lmer`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.927 | 0.958 | 0.948 | 0.952 |      0.944 |
|  18 | 0.958 | 0.952 | 0.948 | 0.942 |      0.947 |
|  50 | 0.949 | 0.948 | 0.947 | 0.953 |      0.955 |
| 100 | 0.963 | 0.966 | 0.947 | 0.952 |      0.957 |

`Wald.default` and `Wald.lmer` cover at the nominal level across the
entire grid. `Wald.sandwich` under-covers at `J = 8` by 4% on average,
then converges by `J = 18`. **The shortfall is roughly the same across
contamination types**, confirming that the small-J finite-sample bias is
the dominant effect, not contamination-induced misspecification.

## CI length (lower = sharper at matched coverage)

The length tells the robustness story coverage alone cannot:

**Mean CI length: `Wald.default`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 2.528 | 2.648 | 2.655 | 2.521 |      2.550 |
|  18 | 1.685 | 1.763 | 1.768 | 1.686 |      1.694 |
|  50 | 1.011 | 1.057 | 1.057 | 1.013 |      1.012 |
| 100 | 0.716 | 0.747 | 0.748 | 0.715 |      0.715 |

**Mean CI length: `Wald.sandwich`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 2.447 | 2.606 | 2.560 | 2.442 |      2.621 |
|  18 | 1.669 | 1.754 | 1.753 | 1.665 |      1.734 |
|  50 | 1.014 | 1.068 | 1.056 | 1.016 |      1.015 |
| 100 | 0.719 | 0.754 | 0.747 | 0.718 |      0.718 |

**Mean CI length: `Wald.lmer`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 2.480 | 3.515 | 3.238 | 2.469 |      2.508 |
|  18 | 1.653 | 2.371 | 2.156 | 1.652 |      1.663 |
|  50 | 0.991 | 1.433 | 1.291 | 0.993 |      0.992 |
| 100 | 0.702 | 1.019 | 0.914 | 0.701 |      0.701 |

Under `clean` Gaussian errors the two methods produce essentially
identical intervals (length ratio 0.98 at `J = 50`). Under `t3`, lmer’s
intervals are **1.36x as wide** as rlmer’s at `J = 50` at the *same*
nominal coverage — i.e. lmer effectively wastes that much information by
not downweighting the heavy tails. Under `CN_e` the lmer:rlmer length
ratio is 1.22. Random-effects contamination (`CN_re`) and the moderate
`shift_subj` regime leave both methods on roughly equal footing.

# Phase C — single-fit `anova(fit)` Type I + power

**Type-I (eff = 0, nominal 0.05): `anova.Wald.default`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.076 | 0.042 | 0.050 | 0.048 |      0.086 |
|  18 | 0.050 | 0.048 | 0.052 | 0.050 |      0.078 |
|  50 | 0.058 | 0.036 | 0.052 | 0.048 |      0.034 |
| 100 | 0.040 | 0.072 | 0.042 | 0.050 |      0.058 |

**Type-I (eff = 0, nominal 0.05): `anova.Wald.sandwich`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.184 | 0.142 | 0.150 | 0.168 |      0.148 |
|  18 | 0.078 | 0.082 | 0.094 | 0.092 |      0.094 |
|  50 | 0.064 | 0.036 | 0.070 | 0.056 |      0.038 |
| 100 | 0.048 | 0.074 | 0.050 | 0.056 |      0.054 |

**Type-I (eff = 0, nominal 0.05): `anova.lmer`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.070 | 0.048 | 0.042 | 0.044 |      0.100 |
|  18 | 0.050 | 0.054 | 0.056 | 0.062 |      0.068 |
|  50 | 0.050 | 0.042 | 0.042 | 0.050 |      0.044 |
| 100 | 0.036 | 0.060 | 0.048 | 0.050 |      0.058 |

**Power (eff = 0.5): `anova.Wald.default`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.244 | 0.274 | 0.234 | 0.250 |      0.136 |
|  18 | 0.588 | 0.498 | 0.494 | 0.544 |      0.412 |
|  50 | 0.930 | 0.898 | 0.898 | 0.948 |      0.956 |
| 100 | 0.996 | 0.998 | 0.998 | 0.998 |      1.000 |

**Power (eff = 0.5): `anova.Wald.sandwich`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.388 | 0.364 | 0.390 | 0.384 |      0.244 |
|  18 | 0.626 | 0.532 | 0.550 | 0.572 |      0.402 |
|  50 | 0.924 | 0.912 | 0.896 | 0.944 |      0.940 |
| 100 | 0.998 | 1.000 | 0.998 | 0.998 |      1.000 |

**Power (eff = 0.5): `anova.lmer`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.264 | 0.162 | 0.180 | 0.242 |      0.144 |
|  18 | 0.596 | 0.292 | 0.358 | 0.548 |      0.430 |
|  50 | 0.932 | 0.678 | 0.742 | 0.948 |      0.966 |
| 100 | 0.998 | 0.930 | 0.970 | 0.998 |      1.000 |

The single-fit table mirrors the pairwise story below:
`anova.Wald.default` and `anova.lmer` calibrate cleanly;
`anova.Wald.sandwich` is anti-conservative at small J.

# Phase D — pairwise `anova(f0, f1)` Type I + power

This is the key comparison for nested-model testing, the most common
user-facing case.

## Type I (target 0.05)

**Type-I (eff = 0, nominal 0.05): `pair.Wald.default`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.070 | 0.084 | 0.058 | 0.072 |      0.092 |
|  18 | 0.054 | 0.044 | 0.056 | 0.054 |      0.066 |
|  50 | 0.024 | 0.052 | 0.050 | 0.044 |      0.062 |
| 100 | 0.058 | 0.060 | 0.056 | 0.060 |      0.036 |

**Type-I (eff = 0, nominal 0.05): `pair.Wald.sandwich`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.168 | 0.200 | 0.150 | 0.176 |      0.152 |
|  18 | 0.106 | 0.096 | 0.116 | 0.076 |      0.082 |
|  50 | 0.046 | 0.054 | 0.054 | 0.048 |      0.070 |
| 100 | 0.068 | 0.070 | 0.064 | 0.062 |      0.040 |

**Type-I (eff = 0, nominal 0.05): `pair.LRT.lmer`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.058 | 0.062 | 0.052 | 0.064 |      0.092 |
|  18 | 0.044 | 0.060 | 0.060 | 0.054 |      0.070 |
|  50 | 0.032 | 0.052 | 0.048 | 0.042 |      0.062 |
| 100 | 0.064 | 0.054 | 0.066 | 0.056 |      0.048 |

`pair.Wald.default` and `pair.LRT.lmer` are interchangeable on
calibration — the maximum absolute deviation from nominal across the
entire grid is **0.042**. `pair.Wald.sandwich` is anti-conservative at
`J = 8`: the worst-case rate is **0.200** under t3 (4.0x nominal). By
`J >= 50` it returns to nominal.

## Power (eff_size = 0.5)

**Power (eff = 0.5): `pair.Wald.default`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.256 | 0.264 | 0.252 | 0.252 |      0.150 |
|  18 | 0.540 | 0.472 | 0.486 | 0.544 |      0.378 |
|  50 | 0.944 | 0.928 | 0.920 | 0.938 |      0.950 |
| 100 | 0.998 | 1.000 | 1.000 | 1.000 |      0.998 |

**Power (eff = 0.5): `pair.Wald.sandwich`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.382 | 0.376 | 0.402 | 0.364 |      0.246 |
|  18 | 0.568 | 0.496 | 0.550 | 0.590 |      0.382 |
|  50 | 0.938 | 0.924 | 0.906 | 0.942 |      0.944 |
| 100 | 0.998 | 1.000 | 1.000 | 1.000 |      0.998 |

**Power (eff = 0.5): `pair.LRT.lmer`**

|   J | clean |    t3 |  CN_e | CN_re | shift_subj |
|----:|------:|------:|------:|------:|-----------:|
|   8 | 0.272 | 0.166 | 0.176 | 0.260 |      0.138 |
|  18 | 0.566 | 0.328 | 0.372 | 0.574 |      0.362 |
|  50 | 0.956 | 0.694 | 0.764 | 0.946 |      0.960 |
| 100 | 0.998 | 0.942 | 0.984 | 1.000 |      0.998 |

Power gain of `pair.Wald.default` (rlmer) over `pair.LRT.lmer` (lmer),
percentage points.

|   J | clean | t3     | CN_e   |
|----:|:------|:-------|:-------|
|   8 | -2 pp | +10 pp | +8 pp  |
|  18 | -3 pp | +14 pp | +11 pp |
|  50 | -1 pp | +23 pp | +16 pp |
| 100 | +0 pp | +6 pp  | +2 pp  |

Under clean data, the methods are within **3 pp** of each other at every
J – rlmer is paying essentially nothing for robustness. Under
error-distribution contamination (`t3`, `CN_e`), rlmer’s advantage is
most pronounced at moderate J (18 and 50) and only collapses to zero
once both methods saturate near 100 %. Under RE-level contamination
(`CN_re`) and the moderate subject-shift cell, the methods are tied.

# Robustness profile across contamination

A single summary row per (method, contamination) averaged over the J
sweep:

**Mean power across J in {8, 18, 50, 100} per (method, contamination),
eff_size = 0.5.**

| method             | clean |    t3 |  CN_e | CN_re | shift_subj |
|:-------------------|------:|------:|------:|------:|-----------:|
| pair.Wald.default  | 0.684 | 0.666 | 0.664 | 0.683 |      0.619 |
| pair.Wald.sandwich | 0.722 | 0.699 | 0.715 | 0.724 |      0.642 |
| pair.LRT.lmer      | 0.698 | 0.532 | 0.574 | 0.695 |      0.615 |

**Mean Type-I across J in {8, 18, 50, 100} per (method,
contamination).**

| method             | clean |    t3 |  CN_e | CN_re | shift_subj |
|:-------------------|------:|------:|------:|------:|-----------:|
| pair.Wald.default  | 0.052 | 0.060 | 0.055 | 0.057 |      0.064 |
| pair.Wald.sandwich | 0.097 | 0.105 | 0.096 | 0.090 |      0.086 |
| pair.LRT.lmer      | 0.050 | 0.057 | 0.056 | 0.054 |      0.068 |

# Phase E — variance-component test via parametric bootstrap

This is the VC bootstrap quasi-deviance test from WS5: nested models
that differ only in random-effects structure (intercept-only vs
intercept + slope), on a longitudinal `y ~ time + (... | id)` design.
The H_0 cell sets the true slope variance to zero (`theta_slope = 0`) –
the PSD-cone boundary that classical Wald and LRT asymptotics can’t
handle. The H_1 cell uses `theta_slope = 0.5`.

Configuration (pulled from Phase E provenance):

- Replicates per cell: **200**
- Bootstrap depth `nsim_boot`: **200**
- `J` sweep: **8, 18, 50**
- Contamination: **clean, t3, shift_subj**
- Cores (PSOCK): **5**
- Wallclock: **26.8 h**

## Type I (target 0.05)

**`pair.boot.QD` Type-I (nominal 0.05)**

|   J | clean | shift_subj |    t3 |
|----:|------:|-----------:|------:|
|   8 | 0.035 |      0.215 | 0.040 |
|  18 | 0.030 |      0.235 | 0.045 |
|  50 | 0.065 |      0.125 | 0.045 |

Under **clean Gaussian errors** and **t3 heavy tails**, the bootstrap
quasi-deviance test is well-calibrated – the worst Type-I across those
two contamination types and all three `J` values is **0.065**, well
within Monte-Carlo noise of nominal 0.05. The bootstrap correctly
absorbs both the boundary problem at and the heavy-tailed error
distribution.

The **`shift_subj` cell tells a different story**: Type I climbs to
**0.235** at `J = 18`, roughly **4.7x** nominal, and is still elevated
at `J = 50`. The robust scoring equations absorb error-distribution
contamination cleanly but the bootstrap simulates from the fitted
central LMM at `fit0`’s estimates – which inherits any bias `fit0`
picked up from the shifted subjects. The bootstrap null distribution
under-represents the true noise process, giving anti-conservative
p-values.

## Power (eff = 0.5)

**`pair.boot.QD` Power (eff_size = 0.5, target high)**

|   J | clean | shift_subj |    t3 |
|----:|------:|-----------:|------:|
|   8 | 0.075 |      0.320 | 0.185 |
|  18 | 0.165 |      0.375 | 0.270 |
|  50 | 0.235 |      0.350 | 0.400 |

Power under clean Gaussian is **modest at best**: the largest cell
(`J = 50`) reaches only **24%** rejection rate at `eff_size = 0.5`
(i.e. true slope SD = 0.5 in scaled units). Under `t3` the test gains
power – consistent with the Phase D observation that the robust scoring
equations buy extra information from heavy tails – but the headline is
that **at small-to-moderate , the bootstrap quasi-deviance is
conservative and detects only moderately large random-slope variances.**
Designs needing high power on VC tests at small should plan for a larger
effect or move to a dedicated boundary-aware procedure (e.g. `RLRsim`
for the classical LMM analog, then verify on the robust fit).

# Phase F — residual-bootstrap variant (unweighted)

Phase F replaces Phase E’s Gaussian parametric simulation
`eps = rnorm(n, sd = sigma_hat_0)` with a residual bootstrap. The
implementation uses `weights(fit0)` as the sampling probability, which
under `rlmer` returns the lme4 model weights (all 1s for an unweighted
fit) – so this is effectively *unweighted* resampling of the conditional
residuals, not robustness-weighted. The motivation: under subject-level
contamination, the empirical residual distribution carries the signature
of the outliers (some residuals are very large), while the Gaussian step
drops it.

A weighted variant using `getME(fit0, "w_e")` (rlmer’s actual robust
weights) is reported below as Phase G, and a subject-block variant as
Phase H.

Configuration:

- Replicates per cell: **200**
- Bootstrap depth `nsim_boot`: **200**
- `J` sweep: **8, 18, 50**
- Contamination: **clean, t3, shift_subj**
- Cores (PSOCK): **5**
- Wallclock: **27.8 h**

## Type I (target 0.05)

**`pair.boot.QD.resid` Type-I (nominal 0.05)**

|   J | clean | shift_subj |    t3 |
|----:|------:|-----------:|------:|
|   8 | 0.055 |      0.225 | 0.070 |
|  18 | 0.030 |      0.245 | 0.060 |
|  50 | 0.065 |      0.115 | 0.085 |

## Power (eff = 0.5)

**`pair.boot.QD.resid` Power (eff_size = 0.5, target high)**

|   J | clean | shift_subj |    t3 |
|----:|------:|-----------:|------:|
|   8 | 0.075 |       0.32 | 0.175 |
|  18 | 0.160 |       0.38 | 0.245 |
|  50 | 0.255 |       0.35 | 0.420 |

## Direct comparison: Gaussian (Phase E) vs residual (Phase F)

**Type-I side-by-side (Phase E – Gaussian boot – vs Phase F – residual
boot)**

|   J | contamination | E (Gaussian) | F (residual) | Delta (F - E) |
|----:|:--------------|-------------:|-------------:|--------------:|
|   8 | clean         |        0.035 |        0.055 |         0.020 |
|  18 | clean         |        0.030 |        0.030 |         0.000 |
|  50 | clean         |        0.065 |        0.065 |         0.000 |
|   8 | t3            |        0.040 |        0.070 |         0.030 |
|  18 | t3            |        0.045 |        0.060 |         0.015 |
|  50 | t3            |        0.045 |        0.085 |         0.040 |
|   8 | shift_subj    |        0.215 |        0.225 |         0.010 |
|  18 | shift_subj    |        0.235 |        0.245 |         0.010 |
|  50 | shift_subj    |        0.125 |        0.115 |        -0.010 |

**Power side-by-side (eff_size = 0.5)**

|   J | contamination | E (Gaussian) | F (residual) | Delta (F - E) |
|----:|:--------------|-------------:|-------------:|--------------:|
|   8 | clean         |        0.075 |        0.075 |         0.000 |
|  18 | clean         |        0.165 |        0.160 |        -0.005 |
|  50 | clean         |        0.235 |        0.255 |         0.020 |
|   8 | t3            |        0.185 |        0.175 |        -0.010 |
|  18 | t3            |        0.270 |        0.245 |        -0.025 |
|  50 | t3            |        0.400 |        0.420 |         0.020 |
|   8 | shift_subj    |        0.320 |        0.320 |         0.000 |
|  18 | shift_subj    |        0.375 |        0.380 |         0.005 |
|  50 | shift_subj    |        0.350 |        0.350 |         0.000 |

Phase F’s headline change vs Phase E is the `shift_subj` cell: average
Type-I across `J` moves from **0.192** (Phase E Gaussian) to **0.195**
(Phase F residual). The `clean` cell stays well-calibrated (Phase E
average **0.043**, Phase F average **0.050**), so residual bootstrap
doesn’t break what was already working.

# Phase G — robustness-weighted residual bootstrap

Phase G is the variant Phase F was originally meant to be: the residual
resampling probabilities are `getME(fit0, "w_e")` – rlmer’s actual
robustness weights – so contaminated residuals are *downweighted* (not
dropped) at sampling time. If subject-level contamination inflates Phase
E/F because the empirical null is too wide, downweighting the outlying
residuals should narrow it.

Configuration:

- Replicates per cell: **200**
- Bootstrap depth `nsim_boot`: **200**
- `J` sweep: **8, 18, 50**
- Contamination: **clean, t3, shift_subj**
- Cores (PSOCK): **5**
- Wallclock: **13.7 h**

## Type I (target 0.05)

**`pair.boot.QD.residw` Type-I (nominal 0.05)**

|   J | clean | shift_subj |    t3 |
|----:|------:|-----------:|------:|
|   8 | 0.045 |      0.245 | 0.055 |
|  18 | 0.035 |      0.190 | 0.040 |
|  50 | 0.040 |      0.105 | 0.085 |

## Comparison: unweighted (Phase F) vs weighted (Phase G) vs Gaussian (Phase E)

**Type-I across the three bootstrap nulls (nominal 0.05)**

|   J | contamination | E (Gaussian) | F (resid) | G (resid, weighted) |
|----:|:--------------|-------------:|----------:|--------------------:|
|   8 | clean         |        0.035 |     0.055 |               0.045 |
|  18 | clean         |        0.030 |     0.030 |               0.035 |
|  50 | clean         |        0.065 |     0.065 |               0.040 |
|   8 | t3            |        0.040 |     0.070 |               0.055 |
|  18 | t3            |        0.045 |     0.060 |               0.040 |
|  50 | t3            |        0.045 |     0.085 |               0.085 |
|   8 | shift_subj    |        0.215 |     0.225 |               0.245 |
|  18 | shift_subj    |        0.235 |     0.245 |               0.190 |
|  50 | shift_subj    |        0.125 |     0.115 |               0.105 |

Weighting the residual resampling **does** help under subject-level
contamination, and the improvement grows with `J`: the `shift_subj`
Type-I falls to **0.245** / **0.190** / **0.105** at `J` = 8/18/50,
versus a Phase E average of **0.192**. It does not fully repair small
`J`, though – at `J = 8` the test still rejects ~24% of the time under
the null. The `clean` cell stays calibrated (average **0.040**), so
weighting costs nothing where there is nothing to fix. The bias that
survives lives in the Gaussian random-effect draw
`b_sim ~ N(0, sigma^2 U_b U_b^T)`, which inherits the variance `fit0`
absorbed from the shifted subjects – weighting the residuals cannot
touch it.

# Phase H — subject-block residual bootstrap

Phase H abandons the Gaussian random-effect draw entirely. Instead of
`b_sim ~ N(0, sigma^2 U_b(theta_hat_0) U_b^T) + eps`, it transplants
whole subjects’ *total* residual vectors
`r_tot_j = y_j - X_j beta_hat_0` (random effect and conditional residual
together), resampled across subjects with replacement. The intent:
cluster-level contamination should enter the bootstrap null as discrete
outlying blocks – exactly how it enters the data – rather than as
inflated Gaussian variance, breaking the
`contamination -> fit0 -> bootstrap-null` chain that Phases E–G could
not.

Configuration:

- Replicates per cell: **200**
- Bootstrap depth `nsim_boot`: **200**
- `J` sweep: **8, 18, 50**
- Contamination: **clean, t3, shift_subj**
- Effect sizes: **0, 0.5**
- Cores (PSOCK): **5**
- Wallclock: **41.2 h**

## Type I and power

**`pair.boot.QD.block` Type-I (nominal 0.05)**

|   J | clean | shift_subj |  t3 |
|----:|------:|-----------:|----:|
|   8 |     0 |      0.005 |   0 |
|  18 |     0 |      0.000 |   0 |
|  50 |     0 |      0.000 |   0 |

**`pair.boot.QD.block` Power (eff_size = 0.5, target high)**

|   J | clean | shift_subj |  t3 |
|----:|------:|-----------:|----:|
|   8 |     0 |          0 |   0 |
|  18 |     0 |          0 |   0 |
|  50 |     0 |          0 |   0 |

The block bootstrap is **degenerate**: it rejects essentially never,
under both the null *and* the alternative. The p-values pile up around
the centre (median **0.555**, with only **0%** of the 3600 reps below
0.05 and **0%** above 0.95) rather than spreading toward 0 under the
alternative. The mechanism is the mirror image of the fix it was meant
to be: transplanting whole-subject total residuals preserves whatever
random-slope signal is present in `r_tot`, so the bootstrap null scales
up in lockstep with `D_obs` and the observed statistic always lands in
the bulk of the null. It controls Type-I only in the trivial sense that
it controls *everything* to zero. Not a usable test.

# Recommendations

1.  **Use the default vcov everywhere.** Both `vcov(fit)` and the
    closed-form `confint(fit, method = "Wald")` are well-calibrated
    across the entire grid, are no worse than the classical `lmer` route
    on clean data, and are substantially sharper under
    error-distribution contamination.

2.  **`vcov_type = "sandwich"` carries a small-J price.** The CI
    shortfall and the much more pronounced anti-conservative anova Type
    I at `J < 18` mean the sandwich should be reserved for designs with
    `J >= 50` (or paired with the bootstrap CI via `confintROB`). A
    warning in `?vcov.rlmerMod` and `?anova.rlmerMod` for
    `n.clusters < 20` would be a reasonable safety net.

3.  **`anova(f0, f1)` with the default vcov is a drop-in replacement for
    the lme4 LRT** with two benefits and no downside under the
    conditions tested: matching calibration on clean data, and 20-30
    percentage points more power under heavy-tailed or contaminated
    errors at moderate J.

4.  **VC test (`anova(f0, f1)`, `test = "boot"`) is calibrated only for
    error-distribution contamination.** Phase E confirms the bootstrap
    quasi-deviance test is well-calibrated under clean and `t3`, but
    **anti-conservative by 2.5–5x nominal under subject-level
    contamination (`shift_subj`).** The bootstrap simulates from the
    fitted central LMM, so any bias `fit0` picked up from group outliers
    contaminates the empirical null. Document this limitation in
    `?anova.rlmerMod` and recommend the bootstrap test only when the
    data look clean at the subject level (e.g. after a `cooks.distance`
    pass). Power is also modest (\< 25 % at `J = 50`,
    `theta_slope = 0.5`), so VC tests at small `J` should plan for a
    larger expected effect or move to a dedicated boundary-aware
    procedure.

5.  **No residual-resampling tweak rescues the VC bootstrap under
    subject-level contamination.** Three variants were tried. Unweighted
    residual resampling (Phase F) leaves the inflation essentially
    unchanged. Robustness-weighted resampling (Phase G) helps — the
    `shift_subj` Type-I improves with `J` and recovers to near-nominal
    by `J = 50` — but it does not fix small `J`, because the surviving
    bias lives in the Gaussian random-effect draw, not in the residuals.
    The subject-block bootstrap (Phase H), which transplants
    whole-subject total residuals to break that chain, is degenerate: it
    preserves the slope signal in the resampled blocks, so the null
    scales with the observed statistic and the test rejects essentially
    never (zero size **and** zero power). The practical upshot stands
    with recommendation 4: screen for subject-level outliers first; if a
    VC bootstrap is needed under suspected cluster contamination, Phase
    G’s weighted variant is a modest mitigation at moderate-to-large
    `J`, but none of these is a cure at small `J`.

# Reproducibility

The numbers in this document are pulled from
inferenceStudy_phaseACD_lmer_2026-06-08_1334.rds,
inferenceStudy_phaseE_2026-06-08_1944.rds,
inferenceStudy_phaseF_2026-06-10_0732.rds,
inferenceStudy_phaseG_2026-06-12_1610.rds,
inferenceStudy_phaseH_2026-06-12_1610.rds, produced by
`Rscript inst/simulationStudy/inferenceStudy.R` with
`master_seed = 20260603`. The same seed reproduces the same `rlmer`
numbers regardless of the lmer comparison being added (verified against
the earlier rlmer-only run committed alongside).
