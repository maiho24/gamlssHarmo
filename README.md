# gamlssHarmo

A reproducible R framework for multi-cohort neuroimaging harmonisation and normative
modelling using hierarchical GAMLSS models.

> This software supports the publication:
> **Unified Multi-Cohort Harmonisation and Normative Modelling of Neuroimaging Data via Hierarchical GAMLSS** — Mai P. Ho, et al. (2026)
> Preprint: [bioRxiv link](https://www.biorxiv.org/content/10.64898/2026.03.08.710422v1)

---

## Quick start

### 1. Clone the repository
```bash
git clone https://github.com/maiho24/gamlssHarmo.git
cd gamlssHarmo
```

### 2. Install dependencies

The scripts use whatever R and packages are active in your current shell session.
Install the required packages once using whichever approach suits your setup.

**Plain R** (if you already have a working R installation)
```bash
Rscript install_dependencies.R
```

**Conda**
```bash
conda env create -f environment.yml
conda activate gamlssHarmo
```

See the [Dependencies](#dependencies) section for the full package list.

### 3. Prepare your data

Your input CSV must contain at minimum:

| Column    | Description                                                      |
|-----------|------------------------------------------------------------------|
| `id`      | Subject identifier                                               |
| `age`     | Age (any unit -- used as-is in the model)                        |
| `sex`     | Sex covariate (any numeric or factor coding)                     |
| `cohort`  | Batch/site identifier (or whatever `batch_var` is set to)        |
| `feat_1`  | One or more feature columns to harmonise                         |

`wave` is optional. If present it is carried through into output files but is not
used as a model predictor.

Place your data at `data/raw/my_data.csv` or pass its path via `--data`.

### 4. Run the three stages
```bash
Rscript scripts/01_fit.R   --data data/raw/my_data.csv
Rscript scripts/02_infer.R --data data/raw/my_data.csv
Rscript scripts/03_plot.R  --pre  data/raw/my_data.csv \
                            --post output/harmonised/combined_harmonised.csv
```

Alternatively, use the top-level runner:
```bash
./gamlssHarmo fit   --data data/raw/my_data.csv
./gamlssHarmo infer --data data/raw/my_data.csv
./gamlssHarmo plot  --pre  data/raw/my_data.csv \
                    --post output/harmonised/combined_harmonised.csv
```

For stage-specific help:
```bash
Rscript scripts/01_fit.R   --help
Rscript scripts/02_infer.R --help
Rscript scripts/03_plot.R  --help
```

---

## Options

### `01_fit.R`

| Argument          | Type        | Default              | Description                                                                  |
|-------------------|-------------|----------------------|------------------------------------------------------------------------------|
| `--config`        | path        | `config/params.yml`  | Path to params.yml                                                           |
| `--data`          | path        | *(required)*         | Raw input CSV                                                                |
| `--output`        | path        | `output/`            | Base output directory; `models/` and `logs/` are created beneath it         |
| `--features-file` | path        | —                    | `.txt` file listing feature column names, one per line; `#` lines ignored   |
| `--one-feature`   | string      | —                    | Single feature name; overrides `--features-file` and config                 |
| `--batch_var`     | string      | `cohort`             | Batch column name                                                            |
| `--longitudinal`  | TRUE/FALSE  | `FALSE`              | Use `formulas_longitudinal` from config (adds random subject effect to mu)  |
| `--log_transform` | TRUE/FALSE  | `FALSE`              | Apply `log(y+1)` before z-standardising                                     |
| `--family_order`  | string      | `SHASH,GG,NO`        | Comma-separated family list, tried in order until convergence                |
| `--n_cores`       | integer     | `1`                  | Parallel cores (`1` = sequential)                                            |

### `02_infer.R`

| Argument          | Type        | Default                   | Description                                                                  |
|-------------------|-------------|---------------------------|------------------------------------------------------------------------------|
| `--config`        | path        | `config/params.yml`       | Path to params.yml                                                           |
| `--data`          | path        | *(required)*              | Input CSV to harmonise                                                       |
| `--models`        | path        | `<output>/models`         | Fitted models directory; overrides the path derived from `--output`         |
| `--output`        | path        | `output/`                 | Base output directory; must match the value used during fitting              |
| `--features-file` | path        | —                         | `.txt` file listing features to harmonise                                   |
| `--one-feature`   | string      | —                         | Single feature name; overrides `--features-file` and config                 |
| `--batch_var`     | string      | `cohort`                  | Batch column name (must match fit stage)                                     |
| `--log_transform` | TRUE/FALSE  | from scaling file         | Validated automatically against the `_scaling.csv` saved during fitting     |
| `--normative`     | TRUE/FALSE  | `TRUE`                    | Compute normative z-scores and centiles                                      |
| `--n_cores`       | integer     | `1`                       | Parallel cores                                                               |

### `03_plot.R`

| Argument          | Type        | Default              | Description                                                                        |
|-------------------|-------------|----------------------|------------------------------------------------------------------------------------|
| `--config`        | path        | `config/params.yml`  | Path to params.yml                                                                 |
| `--pre`           | path        | —                    | Pre-harmonisation CSV (raw input); produces a Pre-Harmonisation panel              |
| `--post`          | path        | —                    | `combined_harmonised.csv` from infer stage; produces a Post-Harmonisation panel   |
| `--output`        | path        | `output/`            | Base output directory; plots saved to `<output>/plots/`                           |
| `--features-file` | path        | —                    | `.txt` file listing features to plot                                               |
| `--one-feature`   | string      | —                    | Single feature name; overrides `--features-file` and config                       |
| `--batch_var`     | string      | `cohort`             | Batch column name                                                                  |
| `--group_col`     | string      | `cohort`             | Column to colour/group trajectories by                                             |
| `--smooth_method` | string      | `loess`              | `loess` or `gam`                                                                   |
| `--age_bin_width` | number      | `5`                  | Age bin width in years                                                             |
| `--fix_y_limits`  | TRUE/FALSE  | `TRUE`               | Lock y-axis range across pre/post panels for direct comparison                    |

At least one of `--pre` or `--post` must be supplied.

---

## Configuration: `config/params.yml`

All CLI arguments have equivalents in `config/params.yml`. CLI always takes priority
over the config file.

### Formula specification

Formulas are specified as term lists per distributional parameter (`mu`, `sigma`,
`nu`, `tau`). Use `{batch}` and `{id}` as placeholders -- substituted at runtime
with `batch_var` and `id_var`.

| YAML value                                      | R formula                           | Meaning                         |
|-------------------------------------------------|-------------------------------------|---------------------------------|
| `null`                                          | `~ 1`                               | Intercept only                  |
| `["pb(age)", "sex"]`                            | `~ pb(age) + sex`                   | Fixed effects, no batch         |
| `["random({batch})"]`                           | `~ random(cohort)`                  | Batch as random effect          |
| `["pb(age)", "sex", "random({batch})"]`         | `~ pb(age) + sex + random(cohort)`  | Full specification              |
| `["{batch}"]`                                   | `~ cohort`                          | Batch as fixed effect           |

Term lists are used exactly as written -- no implicit terms are added. `null` always
means `~ 1`.

### Family fallback

Families are tried in `family_order` order, stopping at first convergence. Within
each multi-parameter family, nu/tau specs are tried in decreasing complexity:
```
SHASH/BCT/BCPE/ST4/EGB2 (4-param):  full -> nu_only -> tau_only -> intercept
GG/SN1/TF/ST3/BCCG      (3-param):  nu   -> intercept
All others               (2-param):  fitted directly
```

---

## Examples

### Single feature, end-to-end
```bash
Rscript scripts/01_fit.R   --data data/raw/my_data.csv --one-feature ThicknessAvg
Rscript scripts/02_infer.R --data data/raw/my_data.csv --one-feature ThicknessAvg
Rscript scripts/03_plot.R  --pre  data/raw/my_data.csv \
                            --post output/harmonised/combined_harmonised.csv \
                            --one-feature ThicknessAvg
```

### Custom batch column
```bash
Rscript scripts/01_fit.R   --data data/raw/my_data.csv --batch_var site
Rscript scripts/02_infer.R --data data/raw/my_data.csv --batch_var site
```

### Override family order
```bash
Rscript scripts/01_fit.R --data data/raw/my_data.csv --family_order "NO,GG,SHASH"
```

### Parallel processing
```bash
Rscript scripts/01_fit.R   --data data/raw/my_data.csv --n_cores 8
Rscript scripts/02_infer.R --data data/raw/my_data.csv --n_cores 8
```

### Features from a file
```bash
Rscript scripts/01_fit.R --data data/raw/my_data.csv --features-file config/my_features.txt
```

### Post-harmonisation trajectories only
```bash
Rscript scripts/03_plot.R --post output/harmonised/combined_harmonised.csv
```

### Custom config file
```bash
Rscript scripts/01_fit.R --config config/hpc_params.yml --data data/raw/my_data.csv
```

### Longitudinal data
```bash
Rscript scripts/01_fit.R --data data/raw/my_data.csv --longitudinal TRUE
Rscript scripts/02_infer.R --data data/raw/my_data.csv
```

---

## Output files

### Fit stage (`output/models/`)

| File                                        | Description                              |
|---------------------------------------------|------------------------------------------|
| `feature_<n>/<n>_model.rds`                 | Fitted GAMLSS model object               |
| `feature_<n>/<n>_scaling.csv`               | `y_mean`, `y_sd`, `log_transform` flag   |
| `feature_<n>/<n>_metrics.csv`               | AIC, BIC, MSE, converged distribution   |
| `feature_<n>/<n>_predictions.csv`           | In-sample fitted values                  |
| `feature_<n>/<n>_diagnostics.pdf`           | Residual diagnostic plots                |
| `feature_<n>/<n>_timing.csv`                | Per-feature timing and status            |
| `feature_timings.csv`                       | Running timing table (sequential mode)   |
| `model_summary.csv`                         | Summary across all features              |

### Infer stage (`output/harmonised/`)

| File                                        | Description                                                                |
|---------------------------------------------|----------------------------------------------------------------------------|
| `feature_<n>/<n>_harmonised.csv`            | Per-feature output: harmonised values, CDF diagnostics, normative scores   |
| `feature_<n>/<n>_timing.csv`                | Per-feature timing and status                                              |
| `combined_harmonised.csv`                   | Wide format: one row per subject, one column per feature (harmonised)      |
| `combined_normative.csv`                    | Wide format: one row per subject, one z-score column per feature           |
| `harmonisation_summary.csv`                 | Summary across all features                                                |

Multiple infer runs accumulate into `combined_harmonised.csv` and
`combined_normative.csv` -- new feature columns are appended and existing columns
are overwritten by the newer result.

### Plot stage (`output/plots/`)

| File                          | Description                              |
|-------------------------------|------------------------------------------|
| `<n>_trajectory.png`          | Age trajectory (pre and/or post panels)  |

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