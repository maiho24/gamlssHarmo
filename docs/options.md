# Command-line Options

Full argument reference for each pipeline stage. Every CLI argument also has an
equivalent in `config/params.yml`; the CLI always takes priority over the config
file. See the [README](../README.md) for the quick-start workflow and
[`diagnose.md`](diagnose.md) for the diagnostic stage in depth.

> **`--n_cores` and progress reporting:** with `--n_cores 1` (the default),
> `diagnose`/`fit`/`infer` log progress per feature as they go. With
> `--n_cores > 1`, per-feature logs are produced by worker processes and are
> not shown in the main console; instead, a live text progress bar (features
> completed + ETA) is shown while the batch runs. Per-feature detail is still
> written to each feature's own output files (e.g. `<feature>_timing.csv`)
> and can be inspected once the run finishes.

---

## `00_diagnose.R` *(optional but recommended for multiple features)*

Analyses distributional properties of each feature and writes a
`feature_recommendations.csv` that can be passed directly to `--feature_families`
in the fit stage. Skip this stage if you already know which families are appropriate.

For a full explanation of the diagnostic workflow, tier classification, and how to
interpret or edit `feature_recommendations.csv`, see [`diagnose.md`](diagnose.md).

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

## `01_fit.R`

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

## `02_infer.R`

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

## `03_plot.R`

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
