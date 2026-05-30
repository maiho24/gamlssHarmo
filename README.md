# gamlssHarmo

A reproducible R framework for multi-cohort neuroimaging harmonisation and normative modelling using hierarchical GAMLSS models.

> This software supports the publication:
> **Unified Multi-Cohort Harmonisation and Normative Modelling of Neuroimaging Data via Hierarchical GAMLSS** — Mai P. Ho, et al. (2026). Preprint: [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.08.710422v1)

---

## Quick start

### 1. Clone the repository
```bash
git clone https://github.com/maiho24/gamlssHarmo.git
cd gamlssHarmo
```

### 2. Install dependencies and register the command

**Conda** (recommended)
```bash
conda env create -f environment.yml
conda activate gamlssHarmo
Rscript install_dependencies.R
```

See the [Dependencies](#dependencies) section for the full package list.

### 3. Prepare your data

Your input CSV must contain at minimum:

| Column  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| `id`    | Subject identifier (customisable via `--id_var`)                                   |
| `age`   | Age in any unit — used as-is in the model (customisable via `--age_var`)           |
| `sex`   | Sex covariate, any numeric or factor coding (customisable via `--sex_var`)         |
| `batch` | Batch/site identifier (customisable via `--batch_var`)                             |
| `feat_1`| One or more feature columns to harmonise                                           |

`wave` is optional. If present it is carried through into output files but is not
used as a model predictor.

If your CSV uses different column names, pass them via CLI flags:
```bash
gamlssHarmo fit --data my_data.csv --id_var SubjectID --age_var years --sex_var gender --batch_var site
```

### 4. Run the pipeline

```bash
# Step 1 — Diagnose distributional properties and generate family recommendations
gamlssHarmo diagnose --data data/raw/my_data.csv
# -> Inspect (and optionally edit) output/diagnostics/feature_recommendations.csv

# Step 2 — Fit one GAMLSS model per feature using the recommendations
gamlssHarmo fit --data data/raw/my_data.csv \
                --feature_families output/diagnostics/feature_recommendations.csv

# Step 3 — Harmonise and compute normative scores
gamlssHarmo infer --data data/raw/my_data.csv

# Step 4 — Plot pre/post age trajectories
gamlssHarmo plot --pre  data/raw/my_data.csv \
                 --post output/harmonised/combined_harmonised.csv
```

> **Why use the diagnose step?** For studies with multiple features, running diagnose
> first and passing `feature_recommendations.csv` to fit is strongly recommended.
> It automates family and formula selection per feature, reduces convergence failures,
> and lets you review or manually edit the recommendations before fitting.
> See [`docs/diagnose.md`](docs/diagnose.md) for a detailed guide.

You can also invoke each stage directly via `Rscript scripts/0{0-3}_*.R` — both
forms are equivalent.

For stage-specific help:
```bash
gamlssHarmo diagnose --help
gamlssHarmo fit      --help
gamlssHarmo infer    --help
gamlssHarmo plot     --help
```

---

## Options

### `00_diagnose.R` *(optional but recommended for multiple features)*

Analyses distributional properties of each feature and writes a
`feature_recommendations.csv` that can be passed directly to `--feature_families`
in the fit stage. Skip this stage if you already know which families are appropriate.

For a full explanation of the diagnostic workflow, tier classification, and how to
interpret or edit `feature_recommendations.csv`, see
[`docs/diagnose.md`](docs/diagnose.md).

| Argument              | Type       | Default             | Description                                                                                    |
|-----------------------|------------|---------------------|------------------------------------------------------------------------------------------------|
| `--config`            | path       | `config/params.yml` | Path to params.yml                                                                             |
| `--data`              | path       | *(required)*        | Raw input CSV                                                                                  |
| `--output`            | path       | `output/`           | Base output directory; diagnostics written to `<output>/diagnostics/`                         |
| `--features_file`     | path       | —                   | `.txt` file listing feature column names, one per line; `#` lines ignored                     |
| `--one_feature`       | string     | —                   | Single feature name; overrides `--features_file` and config                                   |
| `--feature_families`  | path       | —                   | CSV with per-feature formula/family overrides applied before diagnosis                        |
| `--suffix`            | string     | —                   | Appended to auto-derived subdirectory names (e.g. `--suffix _disc` → `diagnostics_disc/`)     |
| `--batch_var`         | string     | `batch`             | Column identifying the batch/site variable                                                     |
| `--id_var`            | string     | `id`                | Column name for subject identifier                                                             |
| `--age_var`           | string     | `age`               | Column name for age; renamed to `age` internally                                               |
| `--sex_var`           | string     | `sex`               | Column name for sex covariate; renamed to `sex` internally                                     |
| `--longitudinal`      | TRUE/FALSE | `FALSE`             | If TRUE, adds `random({id})` to the recommended mu formula                                    |
| `--discrete`          | TRUE/FALSE | `FALSE`             | If TRUE, runs zero-inflation/overdispersion tests for count data instead of skewness/kurtosis |
| `--n_perm`            | integer    | `499`               | Permutations for batch-moment tests                                                            |
| `--min_batch_n`       | integer    | `50`                | Minimum observations per batch for permutation tests                                           |
| `--thresh_skew`       | numeric    | `0.5`               | \|skewness\| threshold for Tier 1 (near-Gaussian) classification                              |
| `--thresh_kurt`       | numeric    | `1.0`               | \|excess kurtosis\| threshold for Tier 1 classification                                       |
| `--thresh_pval`       | numeric    | `0.01`              | p-value threshold for batch higher-moment tests                                                |
| `--thresh_skew_range` | numeric    | `0.3`               | Minimum between-batch skewness range to declare a batch skewness effect                        |
| `--thresh_kurt_range` | numeric    | `0.5`               | Minimum between-batch kurtosis range to declare a batch kurtosis effect                        |
| `--thresh_zerorate_range`   | numeric | `0.10`           | [discrete] Minimum between-batch zero-rate range for Tier 3                                   |
| `--thresh_dispersion_range` | numeric | `0.50`           | [discrete] Minimum between-batch dispersion-ratio range                                        |
| `--n_cores`           | integer    | `1`                 | Parallel cores                                                                                 |

> **Threshold guidance:** The defaults suit typical multi-site neuroimaging studies
> (≥ 50 subjects per batch, 5–30 batches). Increase `--n_perm` to `999` for stricter
> p-value thresholds. Reduce `--min_batch_n` only for studies with uniformly small
> cohorts. The `--thresh_skew` / `--thresh_kurt` pair controls how aggressively
> near-Gaussian features are promoted to Tier 1 (simpler families); tighten them if
> you want all features treated as potentially non-Gaussian. Raise `--thresh_*_range`
> to require a larger practical effect before batch terms are recommended in nu/tau.

### `01_fit.R`

| Argument             | Type       | Default             | Description                                                                                                           |
|----------------------|------------|---------------------|-----------------------------------------------------------------------------------------------------------------------|
| `--config`           | path       | `config/params.yml` | Path to params.yml                                                                                                    |
| `--data`             | path       | *(required)*        | Raw input CSV. Must contain age, sex, `<batch_var>`, `<id_var>`, and feature columns.                                |
| `--output`           | path       | `output/`           | Base output directory; `models/` and `logs/` are created beneath it                                                  |
| `--features_file`    | path       | —                   | `.txt` file listing feature column names, one per line; `#` lines ignored                                            |
| `--one_feature`      | string     | —                   | Single feature name; overrides `--features_file` and config                                                          |
| `--feature_families` | path       | —                   | Per-feature recommendations CSV (e.g. from diagnose stage). Columns: `feature`, `family_order`, `mu_terms`, `sigma_terms`, `nu_terms`, `tau_terms` |
| `--suffix`           | string     | —                   | Appended to auto-derived subdirectory names (e.g. `--suffix _disc` → `models_disc/`)                                 |
| `--batch_var`        | string     | `batch`             | Column identifying the batch/site variable                                                                            |
| `--id_var`           | string     | `id`                | Column name for subject identifier                                                                                    |
| `--age_var`          | string     | `age`               | Column name for age; renamed to `age` internally                                                                      |
| `--sex_var`          | string     | `sex`               | Column name for sex covariate; renamed to `sex` internally                                                            |
| `--longitudinal`     | TRUE/FALSE | `FALSE`             | If TRUE, adds a random subject effect to mu                                                                           |
| `--log_transform`    | TRUE/FALSE | `FALSE`             | [continuous only] Apply `log(y+1)` before fitting. Ignored when `--discrete TRUE`.                                   |
| `--discrete`         | TRUE/FALSE | `FALSE`             | If TRUE, treats all features as count data. Uses count families (PO, NBI, ZINBI, …). Skips log transform.            |
| `--family_order`     | string     | `SHASH,GG,NO`       | Comma-separated global family fallback order. Overridden per-feature by `--feature_families`.                         |
| `--n_cores`          | integer    | `1`                 | Parallel cores (`1` = sequential)                                                                                     |

### `02_infer.R`

| Argument          | Type       | Default             | Description                                                                                          |
|-------------------|------------|---------------------|------------------------------------------------------------------------------------------------------|
| `--config`        | path       | `config/params.yml` | Path to params.yml                                                                                   |
| `--data`          | path       | *(required)*        | Input CSV to harmonise                                                                               |
| `--models`        | path       | `<output>/models`   | Fitted models directory; overrides the path derived from `--output`                                 |
| `--output`        | path       | `output/`           | Base output directory; must match the value used during fitting                                      |
| `--features_file` | path       | —                   | `.txt` file listing features to harmonise                                                            |
| `--one_feature`   | string     | —                   | Single feature name; overrides `--features_file` and config                                          |
| `--suffix`        | string     | —                   | Must match the `--suffix` used during fitting. **Mismatch causes infer to look in the wrong models directory and fail immediately.** |
| `--batch_var`     | string     | `batch`             | Column identifying the batch/site variable (must match the fit stage)                                |
| `--id_var`        | string     | `id`                | Column name for subject identifier (must match the fit stage)                                        |
| `--age_var`       | string     | `age`               | Column name for age; renamed to `age` internally                                                     |
| `--sex_var`       | string     | `sex`               | Column name for sex covariate; renamed to `sex` internally                                           |
| `--log_transform` | TRUE/FALSE | `FALSE`             | Overrides the log-transform flag; normally inferred automatically from the per-feature `_meta.csv`   |
| `--normative`     | TRUE/FALSE | `TRUE`              | If TRUE, computes normative z-scores and saves them to `combined_normative.csv`                      |
| `--n_cores`       | integer    | `1`                 | Parallel cores                                                                                       |

The discrete/continuous mode for each feature is determined automatically from the
`_meta.csv` written during fitting — no `--discrete` flag is needed at infer time.

### `03_plot.R`

| Argument          | Type       | Default             | Description                                                                             |
|-------------------|------------|---------------------|-----------------------------------------------------------------------------------------|
| `--config`        | path       | `config/params.yml` | Path to params.yml                                                                      |
| `--pre`           | path       | —                   | Pre-harmonisation CSV (raw input); produces a Pre-Harmonisation panel                  |
| `--post`          | path       | —                   | `combined_harmonised.csv` from infer stage; produces a Post-Harmonisation panel        |
| `--output`        | path       | `output/`           | Base output directory; plots saved to `<output>/plots/`                                |
| `--features_file` | path       | —                   | `.txt` file listing features to plot                                                    |
| `--one_feature`   | string     | —                   | Single feature name; overrides `--features_file` and config                             |
| `--suffix`        | string     | —                   | Appended to `plots/` and `logs/` subdirectory names                                    |
| `--batch_var`     | string     | `batch`             | Column identifying the batch/site variable                                              |
| `--id_var`        | string     | `id`                | Column name for subject identifier                                                      |
| `--age_var`       | string     | `age`               | Column name for age (x-axis); renamed to `age` internally                               |
| `--sex_var`       | string     | `sex`               | Column name for sex (faceting variable); renamed to `sex` internally                    |
| `--group_col`     | string     | `batch`             | Column to colour/group trajectories by (can be `batch`, `sex`, diagnosis, etc.)        |
| `--smooth_method` | string     | `loess`             | `loess` or `gam`                                                                        |
| `--age_bin_width` | number     | `5`                 | Age bin width in years for trajectory summaries                                         |
| `--fix_y_limits`  | TRUE/FALSE | `TRUE`              | Lock y-axis to the 1st–99th percentile range across panels for direct comparison       |

At least one of `--pre` or `--post` must be supplied. Supplying both produces
side-by-side pre/post comparison panels.

---

## Configuration: `config/params.yml`

All CLI arguments have equivalents in `config/params.yml`. CLI always takes priority
over the config file.

### Formula specification

Formulas are specified as term lists per distributional parameter (`mu`, `sigma`,
`nu`, `tau`). Use `{batch}` and `{id}` as placeholders — substituted at runtime
with `batch_var` and `id_var`.

| YAML value                                      | R formula                           | Meaning                         |
|-------------------------------------------------|-------------------------------------|---------------------------------|
| `null`                                          | `~ 1`                               | Intercept only                  |
| `["pb(age)", "sex"]`                            | `~ pb(age) + sex`                   | Fixed effects, no batch         |
| `["random({batch})"]`                           | `~ random(batch)`                   | Batch as random effect          |
| `["pb(age)", "sex", "random({batch})"]`         | `~ pb(age) + sex + random(batch)`   | Full specification              |
| `["{batch}"]`                                   | `~ batch`                           | Batch as fixed effect           |

Term lists are used exactly as written — no implicit terms are added. `null` always
means `~ 1`.

### Supported families

**Continuous** (default mode):

| Params | Families                              |
|--------|---------------------------------------|
| 4      | SHASH, JSU, BCT, BCPE                 |
| 3      | GG, SN1, TF, PE2, BCCG, exGAUS       |
| 2      | NO, LOGNO                             |

**Discrete / count** (`--discrete TRUE`):

| Params | Families                                          |
|--------|---------------------------------------------------|
| 4      | ZISICHEL                                          |
| 3      | ZINBI, ZANBI, ZIPIG, SICHEL, ZASICHEL             |
| 2      | NBI, NBII, ZIP, PIG, DPO                          |
| 1      | PO                                                |

### Family fallback

Families are tried in `family_order` order, stopping at first convergence. If the
primary specification fails, gamlssHarmo progressively simplifies formula complexity
(reducing nu/tau to intercept-only, then sigma, before moving to the next family).

<details>
<summary>Full fallback sequences (click to expand)</summary>

| Group | Fallback sequence |
|---|---|
| 4-param continuous (SHASH, JSU, BCT, BCPE) | `_full` → `_no_tau` → `_no_nu` → `_intercept` |
| 3-param continuous (GG, TF, SN1, PE2, BCCG, exGAUS) | `_full` → `_intercept` |
| 2-param continuous (NO, LOGNO) | single spec |
| 4-param discrete (ZISICHEL) | `_full` → `_no_tau` → `_no_nu` → `_no_nu_no_tau` → `_no_sigma` → `_no_sigma_no_tau` → `_no_sigma_no_nu` → `_intercept` |
| 3-param discrete (ZINBI, ZANBI, ZIPIG, SICHEL, ZASICHEL) | `_full` → `_no_nu` → `_no_sigma` → `_intercept` |
| 2-param discrete (NBI, NBII, ZIP, PIG, DPO) | `_full` → `_no_sigma` |
| 1-param (PO) | single spec |

Spec name glossary: `_full` = all extra params have their full formula; `_no_X` = parameter X reduced to intercept-only (`~ 1`); `_intercept` = all extra params at intercept.

</details>

---

## Examples

### Full recommended workflow
```bash
# 1. Diagnose (optional)
gamlssHarmo diagnose --data data/raw/my_data.csv

# 2. Fit using diagnose recommendations
gamlssHarmo fit --data data/raw/my_data.csv \
  --feature_families output/diagnostics/feature_recommendations.csv

# 3. Harmonise
gamlssHarmo infer --data data/raw/my_data.csv

# 4. Plot
gamlssHarmo plot --pre  data/raw/my_data.csv \
                   --post output/harmonised/combined_harmonised.csv
```

### Skip diagnose — specify families directly
```bash
gamlssHarmo fit --data data/raw/my_data.csv --family_order "GG,SHASH,NO"
```

### Single feature, end-to-end
```bash
Rscript scripts/01_fit.R   --data data/raw/my_data.csv --one_feature ThicknessAvg
Rscript scripts/02_infer.R --data data/raw/my_data.csv --one_feature ThicknessAvg
Rscript scripts/03_plot.R  --pre  data/raw/my_data.csv \
                            --post output/harmonised/combined_harmonised.csv \
                            --one_feature ThicknessAvg
```

### Discrete / count data
```bash
gamlssHarmo diagnose --data data/raw/my_data.csv --discrete TRUE
gamlssHarmo fit      --data data/raw/my_data.csv --discrete TRUE \
  --feature_families output/diagnostics/feature_recommendations.csv
gamlssHarmo infer    --data data/raw/my_data.csv
```

### Custom batch column
```bash
gamlssHarmo fit   --data data/raw/my_data.csv --batch_var site
gamlssHarmo infer --data data/raw/my_data.csv --batch_var site
```

### Namespace outputs with --suffix (useful for comparing runs)
```bash
# --suffix must be identical across fit, infer, and plot for the same run.
# If fit uses --suffix _shash but infer omits it, infer looks for models/
# instead of models_shash/ and fails with "Models directory not found".
gamlssHarmo fit   --data data/raw/my_data.csv --suffix _shash --family_order SHASH
gamlssHarmo infer --data data/raw/my_data.csv --suffix _shash
gamlssHarmo plot  --post output/harmonised_shash/combined_harmonised.csv --suffix _shash
```

### Parallel processing
```bash
gamlssHarmo fit   --data data/raw/my_data.csv --n_cores 8
gamlssHarmo infer --data data/raw/my_data.csv --n_cores 8
```

### Longitudinal data
```bash
gamlssHarmo fit   --data data/raw/my_data.csv --longitudinal TRUE
gamlssHarmo infer --data data/raw/my_data.csv
```

### Post-harmonisation trajectories only
```bash
gamlssHarmo plot --post output/harmonised/combined_harmonised.csv
```

---

## Output files

After running the full pipeline the output directory looks like this:

```
output/
├── diagnostics/
│   ├── diagnostics_summary.csv       # full per-feature distributional statistics
│   └── feature_recommendations.csv  # pass to fit with --feature_families
├── models/
│   ├── feature_<name>/
│   │   ├── <name>_model.rds          # fitted GAMLSS model object
│   │   ├── <name>_summary.txt        # human-readable equations, links, AIC/BIC
│   │   ├── <name>_metrics.csv        # AIC, BIC, MSE, converged specification
│   │   ├── <name>_predictions.csv    # in-sample fitted distributional parameters
│   │   ├── <name>_diagnostics.pdf    # residual diagnostic plots
│   │   ├── <name>_meta.csv           # discrete / log_transform flags
│   │   └── <name>_timing.csv
│   └── model_summary.csv             # cross-feature summary table
├── harmonised/
│   ├── feature_<name>/
│   │   ├── <name>_harmonised.csv     # per-subject: original, harmonised, z-score
│   │   └── <name>_timing.csv
│   ├── combined_harmonised.csv       # wide format — use this for downstream analyses
│   ├── combined_normative.csv        # wide format — normative z-scores
│   └── combined_cdf.csv              # wide format — CDF / mid-p values
├── plots/
│   └── <name>_trajectory.png         # pre/post age trajectory panels
└── logs/
    ├── diagnose_YYYYMMDD_HHMMSS.log
    ├── fit_YYYYMMDD_HHMMSS.log
    ├── infer_YYYYMMDD_HHMMSS.log
    └── plot_YYYYMMDD_HHMMSS.log
```

### Diagnose stage (`output/diagnostics/`)

| File                          | Description                                                       |
|-------------------------------|-------------------------------------------------------------------|
| `diagnostics_summary.csv`     | Full distributional statistics for every feature                  |
| `feature_recommendations.csv` | Trimmed recommendations (family_order, formula terms) for fit stage |

### Fit stage (`output/models/`)

| File                                    | Description                                                          |
|-----------------------------------------|----------------------------------------------------------------------|
| `feature_<n>/<n>_model.rds`             | Fitted GAMLSS model object                                           |
| `feature_<n>/<n>_meta.csv`              | Per-feature flags: `discrete`, `log_transform`                       |
| `feature_<n>/<n>_metrics.csv`           | AIC, BIC, MSE, converged specification                               |
| `feature_<n>/<n>_predictions.csv`       | In-sample fitted parameter values (mu, sigma, nu, tau)               |
| `feature_<n>/<n>_summary.txt`           | Human-readable model equation, links, parameter interpretations      |
| `feature_<n>/<n>_diagnostics.pdf`       | Residual diagnostic plots                                            |
| `feature_<n>/<n>_timing.csv`            | Per-feature timing and status                                        |
| `feature_timings.csv`                   | Running timing table across all features (sequential mode)           |
| `model_summary.csv`                     | Summary across all features                                          |

### Infer stage (`output/harmonised/`)

| File                                    | Description                                                                              |
|-----------------------------------------|------------------------------------------------------------------------------------------|
| `feature_<n>/<n>_harmonised.csv`        | Per-feature output: original value, harmonised value, CDF value, normative z-score, fitted and harmonised parameters per distributional parameter |
| `feature_<n>/<n>_timing.csv`            | Per-feature timing and status                                                            |
| `combined_harmonised.csv`               | Wide format: one row per subject, one column per feature (harmonised values)             |
| `combined_normative.csv`                | Wide format: one row per subject, one normative z-score column per feature               |

Multiple infer runs accumulate into `combined_harmonised.csv` and
`combined_normative.csv` — new feature columns are appended and existing columns
are overwritten by the newer result.

### Plot stage (`output/plots/`)

| File                      | Description                                       |
|---------------------------|---------------------------------------------------|
| `<n>_trajectory.png`      | Age trajectory plot (pre and/or post panels)      |

---

## Dependencies

| Package        | Version   | Purpose                                      |
|----------------|-----------|----------------------------------------------|
| `gamlss`       | >= 5.4    | Model fitting                                |
| `gamlss.dist`  | >= 6.1    | Distribution families (SHASH, GG, NO, etc.)  |
| `dplyr`        | >= 1.1    | Data manipulation                            |
| `ggplot2`      | >= 3.4    | Plotting                                     |
| `scales`       | >= 1.3    | Axis formatting                              |
| `RColorBrewer` | >= 1.1    | Colour palettes                              |
| `yaml`         | >= 2.3    | Config file parsing                          |
| `optparse`     | >= 1.7    | Argument parsing                             |
| `logger`       | >= 0.3    | Structured logging                           |
| `parallel`     | base R    | Parallel processing (no installation needed) |

Exact versions used in development are pinned in `environment.yml`.

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
