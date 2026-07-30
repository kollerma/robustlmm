# WS24 / bug 6 — verification of the weighted DASvar `tau_e` fix

Verification of commit `8ef7e59` (2026-06-12, *"WS24: fix DASvar tau_e
for weighted fits — bug 6"*), which replaced the leading term of the
DASvar `tau^2` formula in `R/AllClass.R` from `v_e` to `1`.

Run 2026-07-26 on R 4.6.0 Patched (2026-05-01 r89994), lme4 2.0.4,
Matrix 1.7.5, macOS arm64. The development build under test stamps
`3.5.0.2` (`master` DESCRIPTION `Version: 03.5.0-2`); the pre-fix build at
`8ef7e59^` stamps `3.4.3`, which is what its own DESCRIPTION reads, so
the `dev-prefix` rows and `ws24-weight-scale-prefix-results.rds` keep
that number deliberately. Scripts and stamped `.rds` results are in
this directory; every script takes the library path of the robustlmm
build under test as its first argument, so the same code runs against
the released CRAN versions and against the development builds.

Summary: **the fix is correct, the new test is a real guard, and the
affected window is wider than the commit message implies — it includes
CRAN 3.4-5, and it reaches `method = "DAStau"`, which the commit
message says is unchanged.** The two methods are not equally exposed,
though: DASvar is wrong on every weighted fit, DAStau only at `cPsi`
(10 of 72 arms; 0 of 36 under the default `smoothPsi`).

---

## 1. Is the fix correct? (`ws24-tau-mc.R`, `ws24-identity-grid.R`)

### The structural argument

`tau_{e,i}` is the standard deviation of the *standardised* residual.
`updateSigma()` (`R/DAS-scale.R`) forms it from `object@resp$wtres`,
and `rlmerResp` defines `wtres <- sqrtrwt * (y - mu)` with
`sqrtrwt = sqrt(weights)`. So the quantity `tau` scales is
`sqrt(w_i) (y_i - mu_i) / sigma`, whose variance is 1 under
`Var(e_i) = sigma^2 / w_i` — not `w_i`. Every other term in the formula
is built from the whitened design: `diagA` and `diagAAt` come from
`calculateA(U_eX, U_eZU_b, ...)`, and `B() = Kt %*% Lambda_b` with
`Kt = -(tcrossprod(U_eX, M_bB) + U_eZU_b %*% M_bb)`. After the fix,
`grep -n v_e R/*.R src/*.cpp` returns only the field declaration and
the construction `U_e <- Diagonal(x = sqrt(v_e))` — the prior weights
enter the numerics in exactly one place, which is the point.

This holds on every code path reaching the line, including the
`design.weights` (eta) branch added later by `fe8cce1`: that branch was
written after the fix and already carries the leading `1`.

### The lmer oracle

At `rho.e = rho.b = cPsi` the RSE reduces to REML, so a weighted fit
must reproduce `lmer` exactly. `tests/weights.R` (3 datasets x 4 weight
vectors x 2 tau methods = 24 arms) passes on the fixed build. Achieved,
not just PASS — worst achieved / allowed over all 24 arms and all six
checks, so the battery clears its own thresholds by a factor of 2.6:

| check | allowed | worst achieved / allowed |
|---|---|---|
| coef | `tolerance` | 0.389 |
| fixef | `tolerance` | 0.0072 |
| ranef | `100 * tolerance` | 0.208 |
| fitted | `tolerance` | 0.389 |
| predict | `tolerance` | 0.389 |
| summary | `100 * tolerance` | 0.139 |

In absolute terms: on Dyestuff (one variance component, `tolerance =
1e-8`) the identity holds to 3.9e-09 on coef/fitted/predict and
7.2e-11 on fixef. The loosest arms are `ranef` on sleepstudy
(2.85e-05, against the 1e-03 the test allows) — the same
multi-variance-component `theta`-optimiser effect quantified below,
not a `tau` effect.

Pushed past that grid (`ws24-identity-grid.R`: 4 designs — balanced
20x6, unbalanced 25 groups, random slope 20x8, crossed 8x7 — x 6 weight
regimes including `const 0.01`, `const 100`, log-uniform on
[0.01, 100] and a two-point {0.01, 100} design x 2 methods, 48 arms,
no singular `lmer` reference):

| quantity | worst deviation from `lmer` over 48 arms |
|---|---|
| sigma | 4.2e-06 |
| fixef | 2.8e-06 |
| SE | 4.6e-05 |
| ranef | 1.2e-04 |

That is *not* 1e-8, and the reason is worth stating precisely: it is the
`theta` optimiser, not `tau`. On the worst arm (random slope, `weights
= 100`), `rlmer` and `lmer` stop 1.6e-05 apart in `theta`. Pinning
`rlmer` to `lmer`'s `theta` and refitting collapses the disagreement:

```
as fitted        theta diff 1.558e-05 | sigma rel 2.163e-06 | ranef rel 1.186e-04 | se rel 4.571e-05
at lmer's theta  theta diff 0.000e+00 | sigma rel 7.166e-07 | ranef rel 1.815e-11 | se rel 7.165e-07
```

Single-variance-component designs stay at or below 1.3e-08 throughout.
The median deviation over all 48 arms (3.7e-08) is *below* the median
of the same designs fitted with no weights at all (7.8e-08).

### The oracle that does not use lmer

The lmer comparison only shows the two agree; it does not show the
constant matches `tau`'s definition. `ws24-tau-mc.R` does that
directly. Design: one-way, J = 30 x 6 (n = 180), prior weights
log-uniform on [0.05, 19.95] (a 400-fold range); `theta` and `sigma`
held at their true values; 16,000 Monte-Carlo replicates; compare the
empirical sd of `sqrt(w_i)(y_i - mu_i)/sigma` against each candidate
formula. MC error on an sd at 16,000 reps is 0.56%.

```
post-fix (leading 1)   max|z|     2.79   mean z    -0.17   |z|>3 in   0 of 180 obs   max rel dev   1.56%
pre-fix  (leading v_e) max|z|  1490.73   mean z   143.74   |z|>3 in 177 of 180 obs   max rel dev 833.37%
```

**Verdict: confirmed.** The leading constant is 1.

---

## 2. Is the new test a real guard? (`weights-q2` battery)

`tests/weights.R` as of `8ef7e59` cannot be run verbatim against
`8ef7e59^` under lme4 2.0.4: it fails in `model.frame` with
`invalid type (closure) for variable '(weights)'` — that is GitHub
issue #36 (lme4 resolves `weights` in the formula's environment), fixed
*after* `8ef7e59` on both `master` and `github-transition`. The battery
was therefore reproduced with the weight vector in the global
environment; nothing else was changed, and the `summary()` comparison
uses `df = "none"` exactly as the current `tests/weights.R` does.

Against the pre-fix build (`8ef7e59^` = `a377685`), **all 12 DASvar
arms fail and all 12 DAStau arms pass**. Against the fixed build all 24
pass. The discriminating power comes entirely from the DASvar arms the
commit added — which is what a regression guard should look like.

| dataset | weights | rlmer sigma | lmer sigma | rel |
|---|---|---|---|---|
| Dyestuff | `rep(2)` | 46.678570 | 70.017854 | −0.33333 |
| Dyestuff | `rep(0.5)` | 57.169319 | 35.008927 | +0.63299 |
| Dyestuff | `runif` | 52.753875 | 36.595303 | +0.44155 |
| Dyestuff | `rexp` | 49.125886 | 40.101169 | +0.22505 |
| Penicillin | `rep(2)` | 0.518222 | 0.777709 | −0.33366 |
| Penicillin | `rep(0.5)` | 0.635919 | 0.388854 | +0.63537 |
| Penicillin | `runif` | 0.587836 | 0.383030 | +0.53470 |
| Penicillin | `rexp` | 0.584307 | 0.600149 | −0.02640 |
| sleepstudy | `rep(2)` | 24.128195 | 36.192115 | −0.33333 |
| sleepstudy | `rep(0.5)` | 29.550716 | 18.096113 | +0.63299 |
| sleepstudy | `runif` | 28.626987 | 18.473944 | +0.54959 |
| sleepstudy | `rexp` | 24.837315 | 24.144872 | +0.02868 |

On the coverage-study dataset of the original diagnosis (J = 50 x 5,
weights log-uniform on [0.255, 3.99]) the pre-fix cPsi DASvar sigma is
**0.753568 against lmer's 0.880631** — the commit's 0.7476 / 0.8806 to
the precision it quoted.

**Verdict: the guard is real.**

### The sign of the bias

The commit message says "biased low", which is half the story.
`ws24-weight-scale.R`: prior weights are identified only up to a common
factor, so `weights = c * w` must rescale sigma-hat by `sqrt(c)` and
change nothing else. Pre-fix it did not — the answer depended on the
scale the weights happened to be supplied on. Dyestuff, constant
weights, cPsi; every column should read 49.5101:

| c | lmer / sqrt(c) | DASvar / sqrt(c) (pre-fix) | rel. err |
|---|---|---|---|
| 0.01 | 49.5101 | 619.6437 | +1151.6% |
| 0.1 | 49.5101 | 239.9870 | +384.7% |
| 0.5 | 49.5101 | 80.8496 | +63.3% |
| 1 | 49.5101 | 49.5101 | −0.0% |
| 2 | 49.5101 | 33.0067 | −33.3% |
| 10 | 49.5101 | 14.1457 | −71.4% |
| 100 | 49.5101 | 4.4328 | −91.0% |

Post-fix all seven agree with `lmer` to 1.6e-07 or better. Weights
below 1 inflate sigma-hat; weights above 1 deflate it. The `w = 0.01`
corner is where the `tau2 < 0.01` floor starts binding, which is why
the error runs to an order of magnitude rather than the ~15% the commit
message records for weights in [0.25, 4].

---

## 3. Which released versions are affected? (`ws24-version-sweep.R`)

Every CRAN release from 3.2-0 to 3.4-5 was installed into its own
library and run through one fixed script: the coverage-study dataset,
`rlmer(..., weights = w, rho.e = rho.b = cPsi, method = ...)` against
`lmer(..., weights = w)`.

| version | released | DASvar sigma | rel.err | DAStau sigma | rel.err |
|---|---|---|---|---|---|
| 3.1-2 | 2022-12-08 | — | — | — | — (errors: *"Argument weights is unsave to use at the moment."*) |
| 3.2-0 | 2023-03-25 | 0.753568 | −14.43% | 0.835585 | −5.12% |
| 3.2-2 | 2023-08-18 | 0.753568 | −14.43% | 0.835585 | −5.12% |
| 3.2-3 | 2023-10-05 | 0.753568 | −14.43% | 0.835585 | −5.12% |
| 3.2-5 | 2023-11-14 | 0.753568 | −14.43% | 0.835585 | −5.12% |
| 3.3-1 | 2023-12-14 | 0.753568 | −14.43% | 0.835585 | −5.12% |
| 3.3-2 | 2025-05-04 | 0.753568 | −14.43% | 0.835585 | −5.12% |
| 3.3-3 | 2025-05-26 | 0.753568 | −14.43% | 0.835585 | −5.12% |
| 3.4-2 | 2026-01-28 | 0.753568 | −14.43% | 0.835585 | −5.12% |
| 3.4-3 | 2026-04-20 | 0.753568 | −14.43% | 0.835585 | −5.12% |
| **3.4-5** | **2026-06-21** | **0.753568** | **−14.43%** | **0.835585** | **−5.12%** |
| 3.5.0-1 | 2026-07-10 | 0.880631 | +2.0e-10 | 0.880631 | +3.0e-10 |

`lmer` reference sigma 0.880631, identical in every version.

Two corrections to the task's archaeology:

1. **3.4-5 is affected.** The task expected 3.4-4 / 3.4-5 to be clean
   because `8ef7e59` is an ancestor of `origin/github-transition`. It
   is — but it became one only later. 3.4-5 was cut on
   `github-transition` at `2ef2ce8` (2026-06-21); the fix reached that
   branch in `8fc9156` (2026-07-04, *"Merge master (3.5.0 feature line)
   into the CRAN release candidate"*). `git merge-base --is-ancestor
   8ef7e59 2ef2ce8` returns false, and the CRAN 3.4-5 tarball still
   contains `v_e - EDpsi_e * 2 * diagA`. There is no 3.4-4 on CRAN.

2. **The window starts cleanly at 3.2-0**, as inferred, and there was
   a guard that came off exactly there. Full lifecycle from the release
   sources:

   | releases | `weights` behaviour |
   |---|---|
   | 1.0 (2013) | accepted, with `warning("Argument weights is untested.")` (`2c1046c`, 2012-07-09) |
   | 1.1 (2013-06) .. 3.1-2 (2022-12) | **refused**: `stop("Argument weights is unsave to use at the moment.")` (`911bb68` escalated the warning to a stop) |
   | 3.2-0 (2023-03-25) .. 3.4-5 | accepted silently, and wrong |
   | 3.5.0 (2026-07) | accepted silently, and correct |

   The guard was removed in `484249e` (2023-03-24), which shipped in
   3.2-0. That commit is titled *"Fix weights argument — W^1/2 is
   multiplied everywhere, not W^-1/2!!"*: Koller found the whitening
   direction inverted, fixed it, and opened the feature up in the same
   change. It renamed `.U_eX = U_e\inv X` to `U_eX = U_e %*% X`
   throughout — but **left the `tau^2` leading term at `v_e`**, byte
   for byte (`git show 484249e^:R/AllClass.R` and
   `git show 484249e:R/AllClass.R` agree on that line). So the
   double-count is a leftover of the `W^-1/2` era that survived the
   very commit that corrected the era, which is the same family as
   bug 1 of the correctness pass (the IF layer whitening with
   `U_e^-1`).

   Verified guard present in the released sources of 2.3, 2.4-5,
   2.5-1, 3.0-4 and 3.1-1, and at runtime in 3.1-2; 3.2-0
   (*"enable offset and weights"*) is where it became live. That is
   **10 CRAN releases over 3 years and 4 months**, 2023-03-25 to
   2026-07-10. First clean release: **3.5.0** (CRAN 3.5.0-1).

   How exposed users actually were, stated fairly: `weights` is
   **undocumented** throughout the window. It is not in `rlmer()`'s
   formals and has no `@param` in `?rlmer`; it reaches `lmer()` via
   `...`. There is no NEWS file in any released tarball before 3.5.0,
   and nothing in the README, DESCRIPTION or vignettes mentions it. So
   nobody was told the argument works — but nobody was told it does
   not, either, which is the change 3.2-0 made, and `tests/weights.R`
   asserted correctness against `lmer` for the whole period. The honest
   reading is "silently enabled, tested, never documented, and wrong",
   not "a documented feature broke".

Build note: 3.1-2, 3.2-0 and 3.2-2 do not compile against Matrix 1.7.5
(`MATRIX_VALID_ge_dense` removed; `M_cholmod_sdmult` signature
changed). They were installed from source with the corresponding
3.2-3 fixes backported into `src/rlmerMatrixUtils.cpp` — a pure
build-compatibility patch that touches no estimator code. Patched
sources are kept outside the repo under
`~/Documents/R/ws24-sweep/patched/`.

---

## 4a. DAStau's starting value — **affected** (`ws24-dastau-init.R`)

The commit message asserts "DAStau unchanged (its calcTau fixed point
never depended on the DASvar init)". On this dataset it did.

`calcTau()` (`R/DAS-scale.R`) iterated

```r
while (!conv && (it <- it + 1) < max.it && tau[i] <= 1) { ... }
if (tau[i] > 1) { tau[i] <- 1; warning("tau[", i, "] > 1, setting it to 1") }
```

`tau <= 1` belongs on the answer, not on the starting point: when the
initial value came in above 1 the loop body never ran, and the clamp
returned exactly 1. The pre-fix DASvar formula is DAStau's initial
value, and the double-counted weights push it above 1 for every
observation with `w_i > 1`:

```
old init tau^2: range [0.2102, 3.7392]; > 1 for 101 of 250 obs; < 0.01 (floor): 0
new init tau^2: range [0.5140, 0.9836]; > 1 for   0 of 250 obs; < 0.01 (floor): 0
```

Holding `(a, s)` fixed and varying only the start, one `calcTau` pass
gives `max |tau(old init) − tau(new init)| = 0.283`, with 101 values
pinned at exactly 1 (worst observation: `w = 3.947`, correct tau
0.7170, returned 1.0000). The `tau2 < 0.01` floor did **not** bind
here — the binding constraint was the `> 1` guard.

There is no recovery, because tracing `calcTau` through a whole fit
shows it is entered **twice**: once for the b-side (length 50) and once
for the e-side (length 250, 101 starts above 1). The fitted
`tau_e` therefore still has 101 entries equal to 1, and cPsi DAStau
sigma comes out 0.835585 against lmer's 0.880631 (−5.1%).

Whether a fit recovers depends on how many outer iterations it takes.
On the same data with `smoothPsi`, DAStau converges to 0.823244 both
pre- and post-fix — the outer loop re-enters `calcTau` from a start
that now satisfies the guard. That is the configuration the commit
spot-checked, and its conclusion was right for that configuration and
wrong in general. The 101 `tau[i] > 1` warnings *were* raised, but
`rlmer()` collects fitting warnings into `fit@optinfo$warnings` and
muffles them, so nothing appeared at the console.

### How much of the bug reaches DAStau (`ws24-method-exposure.R`)

The two methods are not equally broken, and the difference is large
enough to change the refit advice. DASvar evaluates the wrong formula;
DAStau only *starts* from it. Pre-fix vs post-fix `sigma`-hat over
3 designs x 6 weight regimes x 2 rho x 2 starting values, 144 arms:

| method | broken (rel err > 1e-6) | worst \|rel err\| |
|---|---|---|
| DASvar | **70 of 70** | 0.914 (cPsi), 1.09 (smoothPsi) |
| DAStau | **10 of 72** | 0.132 (cPsi), 1.3e-07 (smoothPsi) |

Fraction of DAStau arms broken, by rho x starting value:

| | `init` = lmer (default) | `init = lmerNoFit` |
|---|---|---|
| cPsi | 0.556 | 0 |
| smoothPsi | **0** | **0** |

All ten broken DAStau arms are `cPsi` + default init, on the two
one-way designs, and every one of them has some `w_i > 1` --- `const
0.5` is not among them, because weights below 1 cannot push the
starting value above the guard. The random-slope design does not
appear either.

The mechanism accounts for all of it. Counting e-side `calcTau`
entries during a fit on the check dataset (`ws24-calctau-entries.R`):

```
cPsi      e-side calcTau entries:  1 | starts above 1 per entry: 101              | final tau==1: 101 | sigma 0.835585
smoothPsi e-side calcTau entries: 13 | starts above 1 per entry: 101,0,0,...,0    | final tau==1:   0 | sigma 0.823244
```

At `cPsi` initialised from an `lmer` fit the estimator is already at
its solution, the outer loop converges immediately, and the single
cold `calcTau` call is final — so the 101 clamped values survive into
`sigma`. Under `smoothPsi` the outer loop runs 13 times; only the
first entry sees an inflated start, every later one warm-starts below
1, and the damage undoes itself.

So: **weighted DASvar fits are always wrong; weighted DAStau fits are
wrong essentially only at cPsi**, which is the diagnostic setting used
to check that `rlmer()` reproduces `lmer()`, not a setting anyone fits
robustly with. This is a finite grid, not a proof — the guard could in
principle bite a `smoothPsi` fit that happened to converge in one
outer pass; none was observed.

**Fixed here** (`R/DAS-scale.R`): start the iteration from
`min(tau, 1)` instead of refusing to start. Pinned by a new block in
`tests/tau.R` (the fixed point must not depend on the start). Verified
to be a no-op on current behaviour: the 24-arm weights battery is
bit-identical to the build without it, and a randomised search over 60
fits (varying J, cluster size, weight scale over 4 orders of magnitude,
VC ratio, cPsi/smoothPsi) finds a maximum post-fix initial `tau^2` of
**0.99948** — the guard no longer binds anywhere, so the change is
defensive.

## 4b. The eta (Mallows) analogue — **unaffected** (`ws24-eta-mc.R`)

Structurally eta cannot be double-counted the way `v_e` was: it never
enters the whitening (`U_e` is built from `v_e` alone). It enters the
bread through `D_e = Diagonal(eta * rho_e@EDpsi())`, hence `M`, hence
`A`; and the meat through the explicit `H = diag(eta)` /
`H^2 = diag(eta^2)` factors in the linear and quadratic terms.

Checked numerically by the same MC oracle, in the hardest case — prior
weights **and** design weights active at once — against the three
leading terms that a double count would produce (4000 reps, n = 180,
eta in [0.15, 1], w log-uniform on [0.05, 19.95]):

| arm | leading `1` (shipped) | leading `eta` | leading `v_e` | leading `eta*v_e` |
|---|---|---|---|---|
| eta = 1, w varying | max\|z\| 3.28, 1/180 | 3.28, 1/180 | 754.41, 174/180 | 754.41, 174/180 |
| eta varying, w = 1 | max\|z\| 3.38, 1/180 | 47.28, 169/180 | 3.38, 1/180 | 47.28, 169/180 |
| eta varying, w varying | max\|z\| 3.26, 2/180 | 89.42, 167/180 | 498.66, 174/180 | 406.15, 175/180 |

The shipped leading term is the only one consistent with the empirical
sd in all three arms. Eta enters the whitening or the formula, never
both.

---

## Files

| file | what |
|---|---|
| `ws24-version-sweep.R` / `-results.rds` | §3 table; one run per installed CRAN version |
| `ws24-sweep-table.R` | assembles the per-version `.rds` into the table |
| `ws24-identity-grid.R` / `-results.rds` | §1 extended cPsi == lmer grid, 48 arms + unweighted control |
| `ws24-tau-mc.R` / `-results.rds` | §1 Monte-Carlo validation of the leading constant |
| `ws24-eta-mc.R` / `-results.rds` | §4b eta double-count check |
| `ws24-dastau-init.R` | §4a DAStau starting-value diagnosis |
| `ws24-method-exposure.R` / `-results.rds` | §4a which method was broken, and how often (144 arms) |
| `ws24-calctau-entries.R` | §4a calcTau entry counts, cPsi vs smoothPsi |
| `ws24-weight-scale.R` / `-prefix-results.rds`, `-fixed-results.rds` | §2 weight-scale equivariance |
