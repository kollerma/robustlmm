# robustlmm NEWS

## robustlmm 3.5.0-1 (2026-07-10)

Several foundations of the new inference toolchain -- the cluster
sandwich `vcov`, Wald and bootstrap `confint`, prediction intervals,
the observation-level influence diagnostics and the base RANSAC
initial estimator -- were first released to CRAN in the 3.4-4 / 3.4-5
series (see below). 3.5.0 adds the features listed here and extends
several of the 3.4-5 foundations.

### New user-facing inference methods

- **`anova.rlmerMod`**. Single-fit per-term robust Wald table;
  pairwise Wald restriction test for fixed-effects-only nested
  models; parametric-bootstrap quasi-deviance test (`test = "boot"`)
  for variance-component comparisons on the PSD-cone boundary, with
  a common-scale convention (use `fit0`'s $\hat\sigma$ throughout)
  to avoid the sign-flip artefact when scales differ between
  models. Type-I matches `anova(lmer0, lmer1)` LRT under all tested
  conditions; under error-distribution contamination at moderate $J$
  rlmer-Wald is 20-30 percentage points more powerful.
  The bootstrap path gains an experimental `null = "robust"` option: it
  generates the parametric bootstrap from a contamination-cleaned null
  fit (trimming the clusters the robust fit heavily downweights) instead
  of `fit0`'s raw estimates, which reduces the anti-conservativeness of
  the variance-component bootstrap under group contamination at larger
  $J$ while preserving power, at the cost of mild conservatism on clean
  data. Falls back to the plain parametric null when no (or more than
  half the) clusters are flagged, when trimming would leave the design
  rank-deficient (a dropped contrast level, a collapsed grouping factor,
  or a constant random-slope covariate --- checked with the same
  nonsingular-subsampling test as the RANSAC initial estimator), or when
  the trimmed refit fails. Validated by simulation, not by a
  finite-sample theorem.

- **`anova(fit0, fit1, test = "score")`** (experimental). One-sided
  robust score test for the special case of a *single added independent
  scalar variance component* (e.g. `(1|g)` vs `(1|g) + (0 + x|g)`, or
  diagonal structures adding one component; anything else stops with an
  error). The statistic is a self-normalised sum of psi-bounded
  per-cluster score contributions computed from the robust *null* fit
  only, calibrated by a score-only parametric bootstrap that refits
  just the null model per replicate (~4x cheaper than `test = "boot"`;
  default `nsim = 199`). In simulation (Gaussian balanced designs, one
  scalar tested component) it fixes the deviance bootstrap's
  contamination-driven anti-conservativeness without losing power:
  contaminated-null Type-I 0.035 vs 0.115 (`test = "boot"`) at $J = 50$
  with 10% of clusters shifted, clean-null 0.045, power 0.920 vs 0.900;
  an adversarial sweep (small $J$, uncentered $x$, observation-level
  outliers) stayed within 0.015-0.075 Type-I. No null cleaning is
  needed (`null = "robust"` is ignored) and no downweighted-group
  warning is issued for this path. Contamination aligned with the
  tested direction is indistinguishable from the alternative for any
  test with power; the attached per-cluster contributions
  `attr(., "boot")$s_j` identify the driving clusters (pair with
  `cooks.distance(., groups = )` / `hatvalues`). Simulation-validated,
  not proven; the deviance bootstrap `test = "boot"` remains the
  default variance-component path.

- **`confint(., method = "boot" / "BCa")` default `boot.type` is now
  `"wild"`** (matching `confintROB`'s own recommendation);
  `boot.type = "parametric"` remains available. The bootstrap
  delegation itself first shipped in 3.4-5, with `"parametric"` as
  the default.

- **Cluster-level influence and robust leverage**.
  `cooks.distance(fit, groups = TRUE)` (or `groups = "<factor>"`) gives a
  per-cluster Cook's distance -- the Mahalanobis norm of the per-cluster
  influence function $\Xi = -J_{par}^{-1} S$ on $(\hat\beta, \hat\sigma,
  \hat\theta)$ -- flagging whole-group outliers that no single
  observation reveals, as a general influence diagnostic. (It is not a
  reliable pre-screen for the bootstrap variance-component test of
  `anova`, which is anti-conservative under group contamination: the
  robust fit absorbs a contaminated group into a downweighted random
  effect, so the bounded influence saturates exactly when contamination
  is present; that test's automatic guard therefore screens with the
  random-effect robustness weights instead.)
  `hatvalues(fit)` returns the robust leverage (the self-leverage
  $A_{ii}$ in the random-effect-whitened convolution; reduces to the
  classical `lme4` leverage at `rho = cPsi`), with `groups` to sum it per
  cluster. `cooks.distance(fit)` per observation is unchanged. Cluster
  level requires a single or nested grouping structure (crossed designs
  error).

- **RANSAC refinements** (extends the `ransac_lme4()` /
  `rlmer_ransac()` / `init = "ransac"` estimator first released in
  3.4-5). `ransac_lme4()` gains adaptive `K`
  (early-stop once the best Q_n score plateaus; `K` is now a maximum
  budget, `adaptive = FALSE` restores the fixed loop) and stratified
  subsampling (`stratify = TRUE`, default: each subsample covers every
  level of the finest grouping factor, so the per-subsample `lmer` stays
  fittable on designs with many small clusters). It also gains
  **nonsingular subsampling** (Koller and Stahel 2017, Algorithm 1, the
  method behind the `setting = "KS2014"` default of `robustbase`): with
  categorical predictors a rank-deficient subsample does not always make
  `lmer` error — it silently drops the aliased columns of a dropped factor
  level — so, rather than draw-then-repair, each subsample is now drawn
  full-rank **by construction**. A new C++ routine
  (`nonsingularSubsampleLU`) runs a Gaxpy-variant LU decomposition with
  partial pivoting and column skipping over a random permutation of the
  observations, greedily selecting the first `p` linearly-independent
  observations and skipping any that are collinear with those already
  chosen; the whole clusters carrying that core seed the draw, and further
  clusters are added until the retained data are identifiable (full-rank
  fixed design, ≥ 2 grouping levels, varying random slopes) and at the
  target size. The result is guaranteed identifiable whenever the full
  design is, so no degenerate candidate ever enters the scale competition.
  The new `n_singular` return element now counts the collinear candidate
  observations skipped by the LU across all draws. A purely continuous,
  full-rank design reduces exactly to a uniform random cluster subsample
  (the paper's guarantee). Following Koller and Stahel (2017, Remark 2),
  `rlmer_ransac()` guards the refinement: when a redescending `rho.e`
  zero-weights observations back into a rank-deficient positive-weight
  design it automatically re-seeds and re-fits from a fresh nonsingular
  RANSAC start, up to the new `max_tries` argument (default 5); only if
  every attempt still collapses is the last fit returned with a warning
  (suggesting a positive-weight `rho.e` such as `smoothPsi`). New exported
  `ransac_basin_radius()` returns the support-preservation radius
  `r*(c) = c sigma / (2 max_j ||x_j||)`, where `c` is the rejection point
  of the redescender — the bisquare cutoff, or found numerically for lqq
  and the other redescenders. `rlmer_ransac()` now runs
  three post-fit diagnostics: a fixed-effect basin check (warns when a
  redescending `rho.e` left `r*` of the high-breakdown start), a
  random-effects phony-correlation check (warns, for any `rho`, when the
  fitted RE covariance is near-singular, `|rho-hat| -> 1`), and the
  refinement-singularity check above (Remark 2). A simulation
  study (`inst/simulationStudy/ransacBasin.R`) shows the basin radius is a
  fixed-effect condition and the phony failure is a distinct
  random-effects attractor, so the two checks are reported separately.
  `rlmer_ransac(n_starts = m)` adds a multi-start consensus: it fits from
  the `m` best distinct RANSAC starts and returns the lowest-residual-
  scale fit whose random-effects covariance is interior, recovering the
  good solution when the single best start falls into the phony
  `|rho-hat| -> 1` attractor (per-start summary in
  `attr(fit, "consensus")`); `n_starts = 1` is the previous behaviour.

- **Redescending psi-functions from `robustbase`, and `lqqPsi`.**
  `makeBisquarePsi()` / `bisquarePsi` now delegate their psi, rho, weight
  and derivative evaluations to `robustbase`'s compiled bisquare family
  (`Mpsi`, `Mwgt`, `Mchi`, normalised via `MrhoInf`), so they match
  `lmrob()` exactly; the returned values are identical (to numerical
  tolerance, verified to 1e-10) to the previous hand-coded implementation.
  New exported general constructor `makeRobustbasePsi(family, cc)` builds
  a `psi_func_rcpp` from any of the `robustbase` redescenders —
  `"bisquare"`, `"lqq"`, `"optimal"`, `"hampel"`, `"ggw"` — defaulting the
  tuning to `robustbase::.Mpsi.tuning.default(family)`. The bisquare
  redescends comparatively fast; the **lqq** (linear-quadratic-quadratic)
  psi of Koller and Stahel (2011), the recommended redescender used by
  `robustbase::lmrob.control(setting = "KS2014")`, is now pre-built and
  exported as `lqqPsi` and is recommended over `bisquarePsi` for
  redescending fits (pair it with `init = "ransac"` / `rlmer_ransac()`,
  which a redescender needs for a good start).

- **Satterthwaite degrees of freedom and p-values**. Robust
  influence-function-based Satterthwaite df, consistent across all four
  inference surfaces: `summary` adds `df` and `Pr(>|t|)` columns;
  `confint(fit, df = "satterthwaite")` uses t-quantiles;
  `anova(fit, ddf = "satterthwaite")` reports an F-test; and
  `emmeans`/`emtrends`/`contrast` return finite df instead of `Inf`.
  Unlike `lmerTest`, the df derive from the robust IF-based covariance
  of the variance parameters, so they stay honest under contamination.
  `summary` now shows the df **by default** (`df = "auto"`) when it is
  cheap to compute — either the underlying influence function is already
  cached on the fit, or a deterministic, dimension-only size workload is
  within the cutoff `getOption("robustlmm.summary.df.max", 5000)` — and
  otherwise falls back to the historic t-value table with a one-line note
  on how to request it (`df = "satterthwaite"` forces it, `df = "none"`
  disables it). The cutoff is dimensionless and machine-independent, so
  the same fit behaves identically everywhere; raise
  `robustlmm.summary.df.max` to show the df on larger fits. The expensive
  influence
  function is now cached on the fit, so `summary`, `anova`, `emmeans`,
  `confint`, `cooks.distance` and the sandwich `vcov` share a single
  computation. Requires the default `vcov`. Supports single-factor,
  nested (`(1 | school/class)`) and crossed (`(1 | subject) + (1 |
  item)`) designs with diagonal random effects. Nested designs
  use a one-way cluster sandwich over the coarsest grouping factor (into
  which all nested random effects aggregate); crossed designs use a
  Cameron-Gelbach-Miller multiway cluster-robust covariance
  (`V_a + V_b - V_{a:b}`), projected to the nearest
  positive-semidefinite matrix. In both cases the analytic standard
  errors match the empirical ones in Monte-Carlo (ratios ~1). Designs
  not covered (e.g. crossed random slopes) fall back to the t-value
  table. When a variance component is estimated at 0 (a boundary fit),
  the df is now computed conditional on that component being held at the
  boundary -- it equals the df of the model with that component dropped,
  and Monte-Carlo coverage of the resulting intervals is nominal;
  previously such fits suppressed the df entirely. Only a genuinely
  non-identifiable (singular) fit now suppresses it. The package options
  are documented under `?\`robustlmm-options\``.

### Robustness to high-leverage design points

- **`rlmer(., design.weights = ...)`**. Mallows-type design
  weights $\eta_i \in (0, 1]$ for robustness to high-leverage design
  points (the classical gap of the RSE / GM-estimator: the
  $\psi$-functions bound the influence of $y$-outliers but not of
  extreme covariates). `"mcd"` computes $\eta$ from robust Mahalanobis
  distances of the non-constant columns of `X`; a numeric vector or `NULL`
  (default, exact previous behaviour) are also accepted. The weight
  enters the e-side of every estimating equation ($\beta, u, \sigma,
  \theta$), the model `vcov`, the Satterthwaite degrees of freedom and
  the influence diagnostics; `summary`
  reports how many observations are downweighted. In a leverage-
  contamination study (10% of points moved to a contaminating line at
  $\sim$7 robust SD) the plain RSE breaks down with `lmer` while
  Mallows-RSE stays near the clean estimate ($\sim$4x less slope bias
  at 5% contamination), at a clean-model efficiency cost of about 1% at
  the simulation-backed default tuning ($\gamma = 1$, cutoff $0.975$).
  Active design weights are supported for a single grouping factor;
  `rlmer` stops with an error on models with more than one.

### Structured random-effect covariances (lme4 >= 2.0-0)

- **`diag(f | g)`**, **`cs(f | g)`** and **`ar1(f | g)`** random-effect
  covariance structures are now supported (both the heterogeneous and
  homogeneous / equal-variance forms of `cs` and `ar1`; `diag()` was
  already accepted in 3.4-5, `cs()` and `ar1()` are new). `diag()`
  (uncorrelated random slopes, equivalent to `(f || g)`) is fitted
  directly. `cs()` (compound symmetric: one common correlation) and
  `ar1()` (autoregressive: `Cor(i, j) = rho^|i - j|`) are fitted by
  projecting each block's unstructured scoring-equation update onto the
  structured manifold every iteration -- the marginal variances are kept
  and the structure's single correlation parameter is estimated from the
  block's robust correlation matrix: the average off-diagonal correlation
  for `cs` (the constrained estimator for that linear structure, Anderson
  1973), and for `ar1` the conditional maximiser of the working-Gaussian
  objective over the manifold -- an ECME / conditional-maximisation step
  (Meng & Rubin 1993; Liu & Rubin 1994) that makes the projected fixed
  point a stationary point of the constrained objective and, unlike the
  former log-correlation least-squares fit, recovers a negative `rho`.
  So the algorithm converges to a genuinely structured fit. For a 2x2
  block (a single correlation) `cs()` reproduces the unstructured fit
  to numerical precision (1e-7 in the regression test). `VarCorr()`
  reports the structured `(sd, corr)`
  parametrisation. `DAStau` is limited to blocks of size <= 2 (it falls
  back to `DASvar` for larger blocks, as before); `cs`/`ar1` blocks of
  size >= 3 are therefore fitted with `DASvar`. `cs`/`ar1` remain
  experimental (the projection has no consistency proof) and are
  validated only for a single grouping factor: `rlmer` now **warns** when
  a `cs`/`ar1` term is combined with more than one grouping factor, where
  neither the fit nor its inference is validated. Inference degrades
  gracefully -- `summary(df = "satterthwaite")` reports the robust df
  where the constrained covariance admits it and otherwise falls back to
  the t-table. A projected `cs`/`ar1` block of dimension `nc` requires
  more than `nc` observations per cluster (within-cluster replication);
  on a saturated block (cluster size `<= nc`, e.g. one observation per
  visit for a 4-level `ar1`) the structured fit is not identifiable, and
  `rlmer` now **stops** up front with an informative error instead of
  failing cryptically deep in the projection.

- **Experimental Monte-Carlo DAS-tau calibration**
  (`options(robustlmm.dastau.mc = TRUE)`). `method = "DAStau"` computes
  its consistency factors by Gauss-Hermite quadrature, which is limited
  to random-effect blocks of dimension <= 2; fits containing a larger
  block -- including `cs`/`ar1`/unstructured blocks of size >= 3 -- fall
  back to `DASvar` with a warning. The new option instead computes the
  same self-consistent DAS-tau fixed point by plain Monte-Carlo
  integration for those blocks. The MC sample is drawn once per fit
  (common random numbers, deterministically seeded via
  `robustlmm.dasmc.seed`, moment-matched), so fits are reproducible and
  the caller's `.Random.seed` is untouched; the number of draws is
  `robustlmm.dasmc.nsim` (default `1e5`). Simulation-validated on clean
  Gaussian data (in the companion ar1 study it removes the small
  `DASvar` calibration residual of ~0.004 in the fitted correlation at
  4-level `ar1` blocks, J = 200), but not backed by a finite-sample
  theorem. Default behaviour is unchanged: with the option off (the
  default) the `DASvar` fallback fires exactly as before, and blocks of
  dimension <= 2 keep the Gauss-Hermite path even with the option on
  (the MC path offers no improvement there). See `?"robustlmm-options"`.

### Small-J caveats and warnings

- `vcov_sandwich` emits a warning when `n.clusters < 20`. In a
  simulation study at $J = 8$, the sandwich CI undercovers (~0.89
  vs nominal 0.95) and the Wald `anova` Type-I inflates to ~0.15-0.20
  (3-4x nominal). G1 correction is necessary but not sufficient at
  small $J$; prefer `type = "default"` or pair the sandwich CI with
  the `confintROB` bootstrap. Caveat paragraph added to
  `?vcov.rlmerMod`, `?vcov_sandwich`, `?confint.rlmerMod`, and
  `?anova.rlmerMod`.

- `anova(., test = "boot")` is well-calibrated under clean Gaussian
  and `t3` heavy-tailed errors but anti-conservative under
  subject-level contamination (Type-I 0.13-0.24 vs nominal 0.05 at
  10% subjects shifted +5 sigma). It now **warns
  automatically** when the null fit heavily downweights a
  random-effects group (smallest robustness weight below 0.5) -- a
  direct signal that the fit treated that group as an outlier, so the
  bootstrap null built from it may be poisoned. The robust fit absorbs
  a contaminated group into a downweighted random effect, so
  `cooks.distance` (an influence measure) does not reliably flag it for
  this purpose; the random-effects weights are the effective screen.
  The path remains experimental.

### Bug fixes and improvements

- Random-effect robustness weights (`getME(., "w_b_vector")`,
  `"w_sigma_b_vector"`) are now aligned with the spherical random
  effects (`"b_s"`). Previously they were concatenated per covariance
  block type, which silently permuted them for RE terms expanding into
  several block types -- notably heterogeneous `diag()` terms
  (lme4 >= 2.0-0). Consequences fixed: `fitEffects()` applied the
  permuted weights inside the fitting iterations (estimates for
  `diag()` fits may change slightly), `getME(., "w_b")` /
  `"w_sigma_b"` attributed weights to the wrong levels,
  `plot()` / `qq()` colored points by the wrong weights, and the
  `anova(test = "boot")` group guard and `null = "robust"` trimming
  could flag/trim the wrong clusters.

- `anova(test = "boot", null = "robust")` now rebuilds the bootstrap
  generator's `U_b = Lambda(theta)` at the full cleaned theta vector
  (exact for any number of variance components; previously only
  single-theta fits were rescaled). The trimmed refit is checked for
  parameter compatibility with the full fit, and the anova table
  heading now states which fallback applied (layout unresolved / >50%
  flagged / refit unusable) instead of implying no cluster was
  flagged.

### Regression suite

- `tests/anova.R`: single-fit Wald shape; pairwise Wald matches
  single-fit Wald for the same restriction; VC-test warning steers
  to bootstrap; reproducibility under seed.

### Documentation

- New vignette `robustlmm-3.5.0` ("New Features in robustlmm 3.5.0"):
  a feature tour with worked examples of the sandwich `vcov`, Wald /
  bootstrap `confint`, robust `anova`, prediction intervals,
  Satterthwaite degrees of freedom, observation- and cluster-level
  influence diagnostics, the RANSAC initial estimator, Mallows design
  weights and the structured (`diag`/`cs`/`ar1`) covariances. Every
  example runs on a small built-in dataset.

### Dependencies

- `lemon` dropped from Suggests. Its `geom_pointline` was soft-deprecated
  under ggplot2 4.0; the simulation-study plots now use the equivalent
  `ggh4x::geom_pointpath` (`ggh4x` was already a dependency).

### Internal / infrastructure

- ggplot2 deprecation cleanup for ggplot2 4.0: `plot.rlmerMod` no longer
  uses the deprecated `aes_string()` (now the `.data[[]]` tidy-eval
  idiom), and the simulation-study scripts use `ggh4x::geom_pointpath`
  in place of the soft-deprecated `lemon::geom_pointline`. The vignettes
  build without ggplot2 deprecation warnings.

- `inst/simulationStudy/inferenceStudy.R` and its results writeup:
  the inference simulation study (CI coverage, Type-I, power under
  contamination, and the variance-component bootstrap test). Kept on
  the development branch as the research record and excluded from the
  package build via `.Rbuildignore`; the user-facing conclusions are
  distilled into the `?anova.rlmerMod` / `vcov_sandwich` caveats.

- `generateRepeatedMeasuresDatasets()`: new data generator for balanced
  repeated-measures datasets with a within-subject factor and a
  structured (`cs`/`ar1`/`diag`/unstructured) random-effect covariance,
  matching the interface of `generateMixedEffectDatasets()`.

- `fitDatasets_rlmer_cs()` / `fitDatasets_rlmer_ar1()`: simulation-study
  fitting wrappers that force a `cs` / `ar1` random-effects covariance
  structure (rewriting the single random-effects term regardless of how
  the data were generated), plus `inst/simulationStudy/structuredCovariances.R`,
  a consistency-and-robustness study showing that under random-effects
  contamination the classical (`lmer`) structured correlation collapses
  toward zero while the robust structured fit stays near the truth.
  Study script and results are excluded from the build via
  `.Rbuildignore`.

### Known limitations

- VC test (`anova(f0, f1)` with random-effects difference) requires
  `test = "boot"`; the parametric bootstrap is currently the only
  exposed VC path. The residual-bootstrap variant explored in a
  companion simulation study is sim-only at this stage.
- Multi-fit chains `anova(f0, f1, f2, ...)` are not yet supported;
  only pairwise comparison.

## robustlmm 3.4-4 / 3.4-5 (2026-06-21, released to CRAN)

The features below shipped to CRAN in 3.4-4 (2026-06-21) and its
build fix 3.4-5 (2026-06-21). They are listed here so the changelog
matches what CRAN users of those versions already have; 3.5.0
extends several of them (see above).

### New user-facing inference methods

- **`vcov(., type = "sandwich")`**. Robust cluster-sandwich
  covariance of the fixed effects, $\hat V = \hat A^{-1} \hat B \hat A^{-T}$
  with $\hat A$ the Schur-complement (marginal) $\beta$-Jacobian and
  $\hat B = \sum_j s_j s_j^T$ the per-cluster $\beta$-score
  contributions. Defaults are preserved byte-for-byte (`type =
    "default"` returns the pre-existing lme4-inherited linearised
  vcov). `correction = "G1"` (default) applies the $J/(J-1)$
  small-sample scaling. `cluster` auto-detects the sole grouping
  factor; on crossed designs `cluster = "factor-name"` warns that the
  sandwich is approximate.

- **`confint(., method = "Wald")`**. Closed-form per-coefficient
  Wald CI from `vcov(., type = vcov_type)`. The default `vcov_type =
    "default"` matches `confint(lmer_fit, method = "Wald")` on clean
  data and is sharper under error-distribution contamination (20-40%
  shorter intervals at matched coverage under `t3`
  and `CN_e`). Bare `confint(fit)` previously errored on `dpoMatrix`
  vcov; now returns the Wald interval.

- **`confint(., method = "boot" / "BCa")`**. Delegates to
  `confintROB::confintROB` (Mason, Cantoni & Ghisletta 2021, 2024)
  with `boot.type = "parametric"` (the 3.4-5 default; 3.5.0 changes
  the default to `"wild"`) or `"wild"`. `confintROB`
  is in `Suggests`; `requireNamespace` guard emits a clear error if
  it is not installed.

- **`predict(., interval = "confidence" / "prediction", level =
  0.95)`**. Returns `(fit, lwr, upr, se)` data.frame. Fixed-
  effect SE from `vcov(., type = "sandwich")`; RE contribution from
  the partial influence function of $\hat u$. Default `interval =
  "none"` preserves the pre-existing numeric vector return.

- **`cooks.distance.rlmerMod`**. Per-observation joint
  Mahalanobis influence on $(\hat\beta, \hat\sigma, \hat\theta)$,
  computed from `implicitIF_full`. Backed by `numDeriv::jacobian`
  (added to Imports). Top-flagged Penicillin observations match the
  IF-thread1 prototype's documented set `{85, 86, 122, 123, 13}`.

- **`influence.rlmerMod`**. Per-observation IF matrix wrapper
  around `implicitIF_full`.

- **`caseweightIF`**, **`vcov_sandwich`**, **`implicitIF_full`**
  exported as user-facing helpers (the IF substrate is also useful
  outside the `vcov` / `confint` / `anova` paths).

- **`emmeans` support**: `recover_data` / `emm_basis` methods are
  registered for `rlmerMod`, so `emmeans` / `emtrends` / `contrast`
  work on robust fits (asymptotic, `Inf` degrees of freedom; 3.5.0
  adds finite Satterthwaite df).

- **`init = "ransac"`** for `rlmer()`. String form for the
  RANSAC high-breakdown start: equivalent to `init = ransac_lme4(formula,
  data)$fit` with default `K = 200`, `sub_frac = 0.5`. Errors with
  a clear message if no subsample fit succeeds.

### Regression suite

- `tests/influence.R`: default-vcov-unchanged guard; `caseweightIF`
  vs reweight-FD agrees to `< 1e-3` on Dyestuff/sleepstudy; sandwich
  finite + PD; crossed-design warning; predict intervals;
  `cooks.distance` top-5 match.
- `tests/bootstrapWald.R`: Wald closed-form match; `confintROB`
  delegation including the parm-NULL fix (confintROB returns a 0-row
  matrix on explicit `parm = NULL`; the wrapper omits `parm` in the
  delegated call and subsets rows itself).
- `tests/cache-invalidation.R`: `.scoreVec`'s `.tau_e` / `.Tbk`
  cache-clearing guard for the full IF.
- `tests/vcov-jsweep-smoke.R`: minimal sandwich-coverage smoke at
  $J = 18$.
- `tests/sandwich-vs-bootstrap.R`: sandwich SE vs Monte-Carlo SE
  cross-check on sleepstudy.

### Dependencies

- `numDeriv` promoted from absent to Imports (used by
  `implicitIF_full`).
