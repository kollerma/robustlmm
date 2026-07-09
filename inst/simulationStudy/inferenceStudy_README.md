# inferenceStudy.R — WS6 Phases A + C + D + E + F + G + H

Monte-Carlo evaluation of the inference machinery added in WS1–WS5:

| Phase | Method under test | Question answered |
|------:|-------------------|-------------------|
| A     | `vcov(., type=)` × `confint(., method="Wald", vcov_type=)` | Empirical 95 % coverage of per-coefficient CIs |
| C     | `anova(fit)` per-term Wald | Type-I + power, single fit |
| D     | `anova(fit0, fit1)` pairwise Wald | Type-I + power, fixed-effects-only nested comparison |
| E     | `anova(fit0, fit1, test = "boot")` parametric VC bootstrap | Type-I + power, variance-component test on the PSD-cone boundary |
| F     | Phase E but with residual-bootstrap inner step | Does a residual bootstrap recover Type-I calibration under subject-level contamination? |
| G     | Phase F but resampling with `getME(fit0, "w_e")` probabilities | Does *robustness-weighted* residual resampling (Phase F's original intent) help where unweighted didn't? |
| H     | Phase E but with subject-block residual bootstrap | Does transplanting whole-subject total residuals fix `shift_subj` Type-I, and at what power cost? |

Phase B (`confintROB` bootstrap CIs) is the only remaining out-of-scope phase — it would multiply Phase A's compute by ~5x and needs to be set up as a separate study. Phases E and F have shipped; the results live in `inferenceStudy_results/`.

## Design

For each cell

```
(phase, J ∈ {8, 18, 50, 100}, contamination ∈ {clean, t3, CN_e, CN_re, shift_subj},
 eff_size ∈ {0, 0.5}, n_per_subj = 10, n_levels = 3)
```

we generate `n_reps` datasets with `generateAnovaDatasets` (a one-way ANOVA `y ~ Var1 + (1 | Var2)`; Var1 has `n_levels` levels, Var2 is the subject grouping factor with `J` levels), apply the chosen contamination, fit `rlmer(..., method = "DASvar")`, and record per-replicate the relevant outputs (CIs / chi-square / p-value) for both `vcov_type = "default"` and `vcov_type = "sandwich"`.

`trueBeta = c(1, eff_size, 2 * eff_size)`, so `eff_size = 0` is the Type-I cell for the Var1 levels (intercept is always 1) and `eff_size > 0` is the power cell.

Contamination generators (`randomNumberGenerators.R` + `contaminationStrategies.R`, both already in the package):

- `clean`       : `N / N` (Gaussian errors, Gaussian REs)
- `t3`          : `t_3 / t_3`
- `CN_e`        : contaminated-normal errors, Gaussian REs
- `CN_re`       : Gaussian errors, contaminated-normal REs
- `shift_subj`  : `N / N` plus a fixed 10 % of *subjects* with `+5 σ_e` shifts (matches `breakdownMC.R`'s `strategy_within`)

## Output

`runInferenceStudy()` returns a list

```
list(
  results = list(A = data.frame(...), C = data.frame(...), D = data.frame(...)),
  aggregates = list(coverage_A, type1_C, power_C, type1_D, power_D),
  provenance = list(R_version, robustlmm_version, master_seed, n_cores, runtime_secs, ...))
```

The aggregate tables are one row per (J × contamination × method [× term × eff_size]); columns are the relevant rate (`covered`, `rej05`, `rej01`) and, for Phase A, mean CI `length`.

The same object is saved to `outfile` (an `.rds`) if supplied.

## Running

Interactive / driven by another script:

```r
source(system.file("simulationStudy/inferenceStudy.R", package = "robustlmm"))
out <- runInferenceStudy(
    phases   = c("A", "C", "D"),
    J_values = c(8L, 18L, 50L, 100L),
    n_reps   = 500L,
    contamination = c("clean", "t3", "CN_e", "CN_re", "shift_subj"),
    eff_sizes = c(0, 0.5),
    n_cores  = max(1L, parallel::detectCores() - 1L),
    master_seed = 20260603L,
    outfile = "inferenceStudy_phaseACD_2026-06-03.rds")
```

Or as a stand-alone Rscript, with optional `key=value` overrides:

```bash
Rscript inst/simulationStudy/inferenceStudy.R \
    phases=A,C,D J_values=8,18,50,100 n_reps=500 n_cores=8
```

## Compute budget

Per-fit cost is ~0.3 s for `rlmer(DASvar)` on a sleepstudy-sized cell (`J × n_per_subj × n_levels = 18 × 10 × 3 = 540` rows). Phases A / C / D each cost roughly `n_reps × cells` fits (Phase D doubles it: two fits per rep).

At the recommended config (4 J × 5 contamination × 2 effect sizes × 500 reps × 2 fits) ≈ 80 000 fits × 0.3 s ≈ 6.7 h on one core, ~30 min on a 16-core machine.

## Method codes in the aggregate tables

| String | What it is |
|--------|------------|
| `Wald.default`         | `confint(fit, method = "Wald", vcov_type = "default")` |
| `Wald.sandwich`        | `confint(fit, method = "Wald", vcov_type = "sandwich")` |
| `anova.Wald.default`   | `anova(fit, vcov_type = "default")` per-term |
| `anova.Wald.sandwich`  | `anova(fit, vcov_type = "sandwich")` per-term |
| `pair.Wald.default`    | `anova(fit0, fit1, vcov_type = "default")` |
| `pair.Wald.sandwich`   | `anova(fit0, fit1, vcov_type = "sandwich")` |
| `pair.LRT.lmer`        | `anova(lmer0, lmer1)` LRT (lme4 baseline) |
| `Wald.lmer`            | `confint(lmer_fit, method = "Wald")` (lme4 baseline) |
| `anova.lmer`           | per-term Wald chi-sq from `vcov(lmer_fit)` (lme4 baseline) |
| `pair.boot.QD`         | `anova(fit0, fit1, test = "boot")` quasi-deviance (Phase E) |
| `pair.boot.QD.resid`   | residual-bootstrap variant of Phase E (Phase F, sim-only) |
| `pair.boot.QD.residw`  | robustness-weighted residual bootstrap (Phase G, sim-only) |
| `pair.boot.QD.block`   | subject-block residual bootstrap (Phase H, sim-only) |

## Phases E and F

Phase E (parametric bootstrap quasi-deviance, the WS5 production path) and Phase F (residual-bootstrap variant) target the variance-component test on the PSD-cone boundary. Both use a longitudinal design `y ~ time + (... | id)` via `generateLongitudinalDatasets`, with `trueTheta = c(1, 0, eff_size)` so `eff_size = 0` is the boundary cell (no random slope variance) and `eff_size > 0` is the power cell.

Phase E uses the package's exposed `anova(., test = "boot")` API. Phase F's residual bootstrap is sim-only: `.anova_pair_boot_resid` in this script replaces Phase E's `eps = rnorm(n, sd = sigma_hat_0)` with `eps = sample(residuals(fit0), n, replace = TRUE, prob = weights(fit0))`. The `prob` argument resolves to the lme4 model weights (constant 1 for an unweighted fit), so the implementation is effectively *unweighted* residual resampling. Phase G is the same function with `weighted = TRUE`, which switches `prob` to `getME(fit0, "w_e")` (rlmer's actual robustness weights).

## Phases G and H

Phases G and H respond to the Phase E/F findings: both the parametric and the unweighted-residual bootstrap are anti-conservative under `shift_subj` because the simulated `b_sim ~ N(0, sigma^2 U_b(theta_hat_0) U_b^T)` inherits the random-effect variance that `fit0` absorbed from the shifted subjects — the bootstrap null is biased in the same direction as the data.

- **Phase G** rules out the cheap fix: keep the Gaussian `b_sim` draw, but downweight contaminated residuals at resampling time via `w_e`. The mechanism analysis predicts this is insufficient (the bias lives in `b_sim`, not in `eps`), but it was Phase F's original intent and is worth ruling out.
- **Phase H** (`.anova_pair_boot_block`) breaks the `contamination -> fit0 -> bootstrap null` chain: it drops the Gaussian `b_sim` draw entirely and transplants whole subjects' *total* residual vectors `r_tot_j = y_j - X_j beta_hat_0` (sampled with replacement across subjects). Cluster-level contamination then enters the bootstrap null as discrete outlying blocks — exactly how it enters the data. Anticipated cost: in power cells the unmodeled slope variance is also preserved by block transplantation, so the bootstrap null is no longer slope-free and power should drop towards alpha. Quantifying that Type-I/power trade-off is the point of running both cell types.

Per-cell cost is dominated by the inner bootstrap: at `nsim_boot = 200` and `J = 18` it is roughly 80 minutes on 5 cores. The full 18-cell grid takes ~24 h.

Parallelism note: Phases A/C/D use `parallel::mclapply` (fork-based). Phases E–H use PSOCK clusters because fork-after-rlmer is unsustainable at the inner-bootstrap scale on macOS (the symptom is `sendMaster: ignoring SIGPIPE signal` after the first few outer reps). The dispatch is `phase %in% c("E", "F", "G", "H") ? PSOCK : mclapply` inside `.run_phase_cell`. The PSOCK workers source `inferenceStudy.R` from the *installed* package (`system.file`), so reinstall after editing this script before launching a parallel run.
