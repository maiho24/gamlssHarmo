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
2. The feature is residualised for age (quadratic: `lm(y ~ poly(age, 2) + sex)`)
   to remove the covariate effect before distributional assessment.
3. Skewness (Type 1, /n) and excess kurtosis are computed on the residuals.
4. A permutation test checks whether the between-batch variance in skewness
   (and kurtosis) is larger than expected by chance.
5. A tier is assigned and a family list is recommended (Cullen-Frey logic).

### Discrete / count mode (`--discrete TRUE`)

Used for integer counts: lesion count, number of episodes, cognitive sub-scores
with a bounded integer range.

**Steps internally:**
1. Feature values are cast to integer via `round()`.
2. A Poisson GLM is fitted (quadratic age + sex) to get covariate-adjusted means.
3. **Zero-inflation test** (Hall & Berenhaut 2002): one-sided z-test comparing
   observed zero count to the Poisson-predicted zero count.
4. **Overdispersion test** (Cameron-Trivedi 1990): auxiliary regression of
   `(y - μ)² - y` on `μ`; upper-tail t-test.
5. **Underdispersion test**: lower-tail variant of the same auxiliary regression.
6. **Pearson dispersion index**: `Var(y) / mean(y)` adjusted for covariates.
7. **Tail heaviness**: P95/median ratio; values > 5 suggest Sichel-class families.
8. Permutation tests for between-batch variation in zero rate and dispersion ratio.
9. A count tier and family order are assigned.

---

## Tier classification

### Continuous tiers

| Tier | Conditions | What changes in recommendations |
|---|---|---|
| **1** | `|skewness| < thresh_skew` AND `|excess kurtosis| < thresh_kurt` | sigma simplified: `random({batch})` only; nu/tau intercept-only |
| **2** | Not near-Gaussian, but no batch shape effect | sigma gets full formula `pb(age);sex;random({batch})`; nu/tau intercept-only |
| **3** | Batch effect detected in skewness **or** kurtosis | As Tier 2 + nu and/or tau also get `random({batch})` |

> **What "batch shape effect" means:** the between-batch variance in skewness (or
> kurtosis) is significantly larger than expected by chance (p < `thresh_pval`) AND
> the range across batches exceeds `thresh_skew_range` (or `thresh_kurt_range`). Both
> conditions must hold to guard against large-sample false positives.

### Discrete tiers

| Tier | Conditions | Family class |
|---|---|---|
| **1** | No ZI, no OD, stable dispersion across batches | NBI, ZANBI, ZINBI, PO, ZIP |
| **2** | Global ZI or OD, or dispersion varies by batch | Same families; signals complex count structure |
| **3** | Zero rate varies significantly by batch | ZI/hurdle families; nu gets `random({batch})` |
| **4** | Severe OD (dispersion index > 3) or heavy tail (P95/med > 5) | Sichel-class: ZASICHEL, ZISICHEL, SICHEL |

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

Minimum between-batch range required to declare a skewness (or kurtosis) batch
effect, in addition to the permutation p-value. This dual criterion prevents
large-sample artefacts where tiny but statistically significant differences in
distributional shape are flagged.

- **Lower** if you want to detect even subtle between-batch shape differences.
- **Raise** to require a larger practical effect before adding batch terms to nu/tau.

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
| `skewness`, `excess_kurtosis` | Large absolute values (> 1 for skewness, > 2 for kurtosis) suggest non-Gaussian; check `tier` |
| `batch_skew_pvalue`, `batch_skew_range` | Both small/large together indicate a skewness batch effect |
| `batch_kurt_pvalue`, `batch_kurt_range` | Same for kurtosis |
| `tier` | 1 = simplest; 3 = most complex formulas recommended |
| `family_override` | TRUE if you supplied a `--feature_families` CSV that overrode this feature |
| `error_message` | Non-empty if diagnose failed for this feature; check n_obs and column names |

For discrete features, the additional columns `zi_pvalue`, `od_pvalue`, `ud_pvalue`,
`dispersion_index`, and `tail_heaviness` give the raw test statistics.

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
