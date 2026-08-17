# Diagnose Stage — Detailed Guide

The diagnose stage (`00_diagnose.R`) analyses each feature's distributional
properties and writes a `feature_recommendations.csv` file that can be passed
directly to the fit stage. For studies with multiple features, running diagnose
first is strongly recommended.

---

## What the stage produces

| Output file | Location | Purpose |
|---|---|---|
| `diagnostics_summary.csv` | `output/diagnostics/` | Full per-feature statistics — skewness, kurtosis, batch test results, tier assignment |
| `feature_recommendations.csv` | `output/diagnostics/` | Trimmed recommendations ready to pass to fit with `--feature_families` |

Both files accumulate across runs, keyed by `feature`: re-running diagnose for
one feature (e.g. `--one_feature`) updates only that feature's row in each
file and leaves rows for every other previously diagnosed feature untouched.

---

## Step-by-step workflow

### 1. Run diagnose

```bash
gamlssHarmo diagnose --data data/raw/my_data.csv
```

For count/integer features (lesion counts, cognitive scores), add `--discrete TRUE`:

```bash
gamlssHarmo diagnose --data data/raw/my_data.csv --discrete TRUE
```

### 2. Inspect the recommendations

Open `output/diagnostics/feature_recommendations.csv`. It looks like this:

| feature | family_order | mu_terms | sigma_terms | nu_terms | tau_terms | tier |
|---|---|---|---|---|---|---|
| ThicknessAvg | SHASH;GG;NO | pb(age);sex;random({batch}) | pb(age);sex;random({batch}) | random({batch}) | random({batch}) | 3 |
| VolumeTotal | GG;BCCG;NO | pb(age);sex;random({batch}) | random({batch}) | | | 1 |
| WMH_count | ZINBI;ZANBI;NBI | pb(age);sex;random({batch}) | random({batch}) | random({batch}) | | 3 |

### 3. Edit if needed (optional)

You can open `feature_recommendations.csv` in any spreadsheet editor and manually
override specific rows before passing to fit. Common edits:

- **Change `family_order`** — e.g. put `NO` first if you are confident a feature
  is Gaussian and want faster fitting.
- **Remove nu/tau batch terms** — e.g. clear `nu_terms` for a feature where you
  know the batch does not affect shape, only location/scale.
- **Add a family** — append to the semicolon-separated list, e.g.
  `SHASH;GG;BCT;NO`.

### 4. Pass to fit

```bash
gamlssHarmo fit --data data/raw/my_data.csv \
  --feature_families output/diagnostics/feature_recommendations.csv
```

The fit stage uses the per-feature recommendations from the CSV. Features not
listed in the CSV fall back to the global CLI/config values.

---

## The two modes

### Continuous mode (default)

Used for imaging-derived measures, cognitive scores, and any real-valued feature.

**Steps internally:**
1. Rows with missing feature, age, or sex values are removed.
2. The feature is residualised for age (natural cubic spline, 4 df), sex, and
   batch as fixed effects (`lm(y ~ splines::ns(age, df = 4) + sex + batch)`)
   so the moments reflect the within-batch conditional shape rather than
   between-batch mixing. A spline is used instead of a low-order polynomial
   because a single global polynomial fits poorly across a pooled age range
   spanning childhood growth, adolescent development, adult plateau, and
   late-life decline; a rigid fit can leave residual age curvature that gets
   misattributed to skewness/kurtosis batch effects.
3. A second-stage model of `log(e^2)` on the same design is fitted. This yields
   `sigma_varies_p` (does the scale depend on age or sex?), `sigma_cv2`, and
   studentised residuals. Raw residuals from a homoscedastic fit form a scale
   mixture, which carries excess kurtosis `3 * CV^2(sigma^2)` even when every
   conditional distribution is Gaussian — so `3 * sigma_cv2` is the share of a
   feature's kurtosis attributable to varying scale alone. These three are
   **reported only** and do not drive any decision.
4. Skewness (Type 1, /n) and excess kurtosis are computed on the residuals, and
   also on the studentised residuals (`skewness_stud`, `excess_kurtosis_stud`).
5. A permutation test checks whether the between-batch variance **and** the
   between-batch range of each moment are larger than expected by chance. Labels
   are permuted within age-band x sex strata, at the subject level, and strata
   holding a single batch are dropped — see `--stratify_permutation`.
6. A tier is assigned and a family list is recommended (Cullen-Frey logic).

### Discrete / count mode (`--discrete TRUE`)

Used for integer counts: lesion count, number of episodes, cognitive sub-scores
with a bounded integer range.

**Steps internally:**
1. Feature values are cast to integer via `round()`.
2. A Poisson GLM is fitted (natural cubic spline age, 4 df, + sex) to get
   covariate-adjusted means.
3. **Zero-inflation test** (Hall & Berenhaut 2002): one-sided z-test comparing
   observed zero count to the Poisson-predicted zero count.
4. **Overdispersion test** (Cameron-Trivedi 1990): auxiliary regression of
   `((y - μ)² - y) / μ` on `μ`, no intercept; upper-tail t-test. The division
   by `μ` is what makes the fitted coefficient a constant.
5. **Underdispersion test**: lower-tail variant of the same auxiliary regression.
6. **Pearson dispersion index**: `mean((y - μ)² / μ)` — the covariate-adjusted
   Pearson χ² statistic divided by `n`. This is the `dispersion_index` column and
   is what the "severe OD > 3" gate uses. It is *not* the same as the
   `dispersion_ratio` column, which is the raw marginal `var(y) / mean(y)`.
7. **Tail heaviness**: P95/median ratio; values > 5 suggest Sichel-class families.
   Returns `NA` when the median is 0 (common for sparse counts), and
   `heavy_tail` is FALSE whenever it is `NA`.
8. Permutation tests for between-batch variation in zero rate and dispersion ratio.
9. A count tier and family order are assigned.

---

## Tier classification

### Continuous tiers

The conditions are evaluated **in this order** — the batch-effect test is checked
first, so a feature can be near-Gaussian and still be Tier 3:

| Tier | Conditions | What changes in recommendations |
|---|---|---|
| **3** | Batch effect detected in skewness **or** kurtosis | As Tier 2, **plus** `random({batch})` on nu (if the *skewness* effect fired) and/or tau (if the *kurtosis* effect fired) |
| **1** | No batch shape effect, AND `|skewness| < thresh_skew` AND `|excess kurtosis| < thresh_kurt` | sigma simplified: `random({batch})` only; nu/tau intercept-only |
| **2** | No batch shape effect, but not near-Gaussian | sigma gets full formula `pb(age);sex;random({batch})`; nu/tau intercept-only |

Note that nu and tau are gated **separately**: a Tier 3 feature flagged only on
kurtosis gets `tau_terms = random({batch})` with `nu_terms` empty. Check the
`skew_batch_effect` and `kurt_batch_effect` columns to see which fired.

> **What "batch shape effect" means:** three conditions must all hold for the
> moment in question.
>
> 1. The between-batch **variance** is larger than the permutation null allows
>    (`batch_*_pvalue` < `thresh_pval`).
> 2. The between-batch **range** is larger than the permutation null allows
>    (`batch_*_prange` < `thresh_pval`). This is the effect-size guard, calibrated
>    against the same null rather than a fixed cut-off, because what counts as a
>    large range depends on both the batch sizes and the tail weight of the
>    individual feature.
> 3. The range clears the absolute floor `thresh_skew_range` /
>    `thresh_kurt_range`. This is retained as a backstop but is no longer the
>    discriminating criterion.

### Discrete tiers

| Tier | Conditions | Family class |
|---|---|---|
| **1** | No ZI, no OD, stable dispersion across batches | NBI, ZANBI, ZINBI, PO, ZIP |
| **2** | Global ZI or OD, or dispersion varies by batch | Same families; signals complex count structure |
| **3** | Zero rate varies significantly by batch | ZI/hurdle families; the zero-probability parameter gets `random({batch})` — nu for 3-parameter families (ZINBI, ZANBI, ZIPIG), tau for the 4-parameter Sichel ones |
| **4** | Severe OD (dispersion index > 3) or heavy tail (P95/med > 5), **and already at tier ≥ 2** | Sichel-class: ZASICHEL, ZISICHEL, SICHEL |

Tier 4 is a promotion applied on top of tier 2 or 3, not an independent branch.
A heavy-tailed feature with no zero-inflation, no overdispersion and stable batch
dispersion stays at **tier 1** and is never promoted — and its family order leads
with `NBI` rather than the Sichel class.

---

## `feature_recommendations.csv` column reference

| Column | Format | Meaning |
|---|---|---|
| `feature` | string | Feature column name — must match the data CSV exactly |
| `family_order` | semicolons: `SHASH;GG;NO` | Families tried in order; first convergence wins |
| `mu_terms` | semicolons: `pb(age);sex;random({batch})` | Formula terms for the location parameter |
| `sigma_terms` | semicolons | Formula terms for the scale/dispersion parameter |
| `nu_terms` | semicolons or empty | Formula terms for the skewness/ZI parameter; empty = intercept-only |
| `tau_terms` | semicolons or empty | Formula terms for the kurtosis parameter; empty = intercept-only |
| `tier` | integer 1–4 | Assigned tier (informational; not read by fit) |

**Placeholder substitution:** `{batch}` and `{id}` in formula terms are replaced
at runtime with the values of `--batch_var` and `--id_var`. You can use bare column
names (e.g. `sex`) or smoothers (e.g. `pb(age)`, `random({batch})`).

**Empty vs absent cell:**
- A cell that is present but empty (`""`) → `~ 1` (intercept-only for that parameter).
- A column that is entirely absent from the CSV → the global config/CLI value is used.

---

## Parameter reference

### `--n_perm` (default: 499)

Number of permutations for the between-batch skewness and kurtosis tests.

- **499** is adequate when `thresh_pval = 0.01` (minimum achievable p ≈ 0.002).
- Use **999** if you lower `thresh_pval` to 0.001.
- Fewer permutations speed up diagnose on very large feature sets; the cost is
  noisier p-values near the threshold.

### `--min_batch_n` (default: 50)

Minimum number of observations in a batch for it to be included in the permutation
tests. Batches smaller than this are excluded. The test still runs if at least 2
batches meet the threshold.

- Reduce to **20–30** for studies where all sites are small (rare-disease cohorts).
- If too few batches qualify (< 2), the test returns `NA` and the batch effect is
  not declared; the tier defaults to 1 or 2 based on global moment statistics only.

### `--thresh_skew` (default: 0.5) and `--thresh_kurt` (default: 1.0)

Thresholds for Tier 1 (near-Gaussian) classification. A feature must satisfy
**both** `|skewness| < thresh_skew` and `|excess kurtosis| < thresh_kurt`.

- **Loosen** (increase both) to classify more features as Tier 1, giving simpler
  sigma formulas and faster fitting. Suitable when you trust the data is
  approximately Gaussian.
- **Tighten** (decrease) to be more conservative — all features will get the fuller
  sigma formula. Safer for exploratory analyses.

### `--thresh_pval` (default: 0.01)

p-value threshold for declaring a batch higher-moment effect. Used for **both**
continuous (skewness/kurtosis) and discrete (zero-rate/dispersion) permutation tests.

- Raising to **0.05** detects weaker batch shape effects; more features assigned Tier 3.
- Lowering to **0.001** only flags very strong effects; most features stay Tier 1/2.

### `--thresh_skew_range` (default: 0.3) and `--thresh_kurt_range` (default: 0.5)

Absolute floor on the between-batch range required to declare a skewness (or
kurtosis) batch effect. This is now a backstop rather than the operative
criterion: the effect-size guard is the permutation p-value on the range
(`batch_*_prange`), which adapts to batch sizes and to the tail weight of the
individual feature. A fixed constant cannot do that — on lifespan neuroimaging
data the 95th percentile of the null kurtosis range varies by roughly a factor of
two between features.

- Leave at the defaults unless you have a substantive reason for a hard minimum.
- **Raise** to require a larger practical effect before adding batch terms to nu/tau.

### `--seed` (default: 20260818)

Seed for the permutation tests, recorded in the `seed` column of
`diagnostics_summary.csv`.

Each feature is seeded separately, with the seed derived from the **feature
name** rather than its position in the run. Two consequences:

- Results do not depend on `--n_cores`, or on how work is distributed across
  workers.
- Re-running one feature with `--one_feature` reproduces exactly the p-values it
  had in a full run. Since both output CSVs are merged by feature key, a
  subset re-run therefore updates rows without perturbing them.

### `--stratify_permutation` (default: TRUE)

Controls the null used by the batch shape tests.

With `TRUE`, batch labels are permuted only **within age-band × sex strata**, at
the **subject** level, and strata containing fewer than two batches are dropped.
The share of rows surviving that restriction is written per feature as
`frac_common`.

This matters whenever batch is associated with age, which is the norm for
multi-cohort lifespan data. Under an unconditional permutation every permuted
batch becomes a random slice of the pooled age range, so the null describes a
population that does not exist, and any age-varying residual shape is read as a
batch effect. Permuting subjects rather than rows matters separately, whenever
participants contribute repeated scans: batch is constant within participant, so
shuffling rows treats repeat scans as independent and narrows the null.

Set `FALSE` to recover the unconditional test — useful for comparing against a
previous run. Note that subject-level permutation stays on either way.

Where two batches share no stratum at all, the contrast between them is **not
identified**: batch and age cannot be separated without an untestable assumption
about the age trajectory. Stratification excludes those pairs rather than
answering them, which is why `frac_common` should be reported alongside any
claim about batch effects on shape. Treat features below roughly 0.3 as weakly
identified.

### `--n_age_bands` (default: 10)

Number of age quantile bands used to build the strata, crossed with sex. Only
used when `--stratify_permutation TRUE`.

- **Fewer** bands: coarser age control, more data retained (`frac_common` rises).
- **More** bands: tighter age control, less data retained.

Worth running at two settings to confirm the conclusion is not an artefact of the
granularity.

### `--thresh_zerorate_range` (default: 0.10) — discrete mode only

Minimum between-batch range in zero rate to declare a zero-inflation batch effect
(Tier 3). E.g. `0.10` means the zero rates across batches must span at least 10
percentage points (e.g. 5% in one site vs. 20% in another).

### `--thresh_dispersion_range` (default: 0.50) — discrete mode only

Minimum between-batch range in the variance-to-mean dispersion ratio to declare a
dispersion batch effect.

---

## Interpreting `diagnostics_summary.csv`

Key columns for each feature:

| Column | What to look for |
|---|---|
| `seed` | The per-feature seed. Re-running this feature alone reproduces the row exactly. |
| `skewness`, `excess_kurtosis` | The (s, k) coordinate that chose the family. Large absolute values (> 1 for skewness, > 2 for kurtosis) suggest non-Gaussian; check `tier` |
| `skewness_stud`, `excess_kurtosis_stud` | Same moments on studentised residuals. A large gap from the raw values means they are partly measuring varying scale — or that the variance model is misbehaving for this feature. Compare before trusting either. |
| `sigma_varies_p`, `sigma_cv2` | Does the scale depend on age/sex, and how much kurtosis that accounts for. Compare `3 * sigma_cv2` against `excess_kurtosis`. |
| `frac_common` | **Share of rows the batch test actually used.** Low values mean the batches barely overlap in age for this feature and the test is weakly identified. Treat below ~0.3 as uninterpretable. |
| `batch_skew_pvalue`, `batch_skew_prange` | Both must clear `thresh_pval` for a skewness batch effect |
| `batch_kurt_pvalue`, `batch_kurt_prange` | Same for kurtosis |
| `batch_*_range`, `batch_*_nbatches` | Observed effect size, and how many batches contributed |
| `near_gaussian` | The pooled-moment condition behind Tier 1 vs 2 |
| `tier` | 1 = simplest; 3 = most complex formulas recommended |
| `family_override` | TRUE if you supplied a `--feature_families` CSV that overrode this feature |
| `error_message` | Non-empty if diagnose failed for this feature; check n_obs and column names |

For discrete features, the additional columns `zi_pvalue`, `od_pvalue`, `ud_pvalue`,
`dispersion_index`, and `tail_heaviness` give the raw test statistics, and
`batch_zerorate_prange` / `batch_disp_prange` are the calibrated range guards.

Two cautions when reading these against fitted models:

- `nu_terms` / `tau_terms` record what was **requested**, not what was fitted. The
  specification that converged appears as a suffix in the `distribution` column of
  `model_summary.csv` (e.g. `BCT_no_tau`). Check that before concluding a feature
  carries a batch shape effect.
- `family_order` is a preference ordering, not a support constraint. `SHASH` and
  `NO` sit at the tail of every positive-support branch as fallbacks, and fit
  keeps the first specification that *converges* rather than the best by
  likelihood.

---

## Common scenarios

### I have 100+ features and want to run them in parallel

```bash
gamlssHarmo diagnose --data data/raw/my_data.csv --n_cores 8
```

### Some features are count data, others are continuous

Run diagnose twice with different `--suffix` values:

```bash
gamlssHarmo diagnose --data data/raw/my_data.csv \
  --features_file config/continuous_features.txt --suffix _cont

gamlssHarmo diagnose --data data/raw/my_data.csv \
  --features_file config/count_features.txt --discrete TRUE --suffix _disc
```

Then concatenate the two `feature_recommendations.csv` files (they have the same
columns) and pass the merged file to fit.

### I know the family for some features but not others

Manually create a `feature_families_override.csv` with only the features you want
to override, then pass it to diagnose with `--feature_families`. The stage will
apply your overrides and fill in the rest via the normal diagnostic logic.

### My cohorts are small (< 50 per batch)

Reduce `--min_batch_n` so those batches are included in the permutation tests:

```bash
gamlssHarmo diagnose --data data/raw/my_data.csv --min_batch_n 20
```

Note that with very small batches the moment estimates are noisy; interpret the
batch effect p-values cautiously.
