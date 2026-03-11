# gamlssHarmo

A reproducible R framework for multi-cohort neuroimaging harmonisation and normative modelling using hierarchical GAMLSS models.

---

## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/maiho24/gamlssHarmo.git
cd gamlssHarmo
```

### 2. Install dependencies

The scripts use whatever R and packages are active in your current shell session. Install the required packages once using whichever approach suits your setup.

**Plain R** (if you already have a working R installation)

```bash
Rscript install_dependencies.R
```

See the [Dependencies](#dependencies) section for the full package list.

**Conda**

```bash
conda env create -f environment.yml
conda activate gamlssHarmo
```

### 3. Prepare your data

Your input CSV must contain at minimum:

| Column   | Description                                               |
|----------|-----------------------------------------------------------|
| `id`     | Subject identifier                                        |
| `age`    | Age (any unit -- used as-is in the model)                 |
| `sex`    | Sex covariate (any numeric or factor coding)              |
| `cohort` | Batch/site identifier (or whatever `batch_var` is set to) |
| `feat_1` | One or more feature columns to harmonise                  |

`wave` is optional. If present it is carried through into output files but not required.

Place your data at `data/raw/my_data.csv` or pass its path via `--data`.

### 4. Run the three stages

Call the runner scripts directly with `Rscript` from the project root:

```bash
Rscript scripts/01_fit.R   --data data/raw/my_data.csv
Rscript scripts/02_infer.R --data data/raw/my_data.csv
Rscript scripts/03_plot.R  --pre  data/raw/my_data.csv \
                            --post output/harmonised/combined_harmonised.csv
```

For help on any stage:

```bash
Rscript scripts/01_fit.R   --help
Rscript scripts/02_infer.R --help
Rscript scripts/03_plot.R  --help
```

---

## Options

### `01_fit.R`

| Argument | Type | Description |
|---|---|---|
| `--config` | path | Path to params.yml [config/params.yml] |
| `--data` | path | Raw input CSV |
| `--output` | path | Base output directory |
| `--features` | path | Features .txt file |
| `--feature` | string | Single feature name (overrides `--features`) |
| `--batch_var` | string | Batch column name (e.g. `"site"`, `"scanner"`) |
| `--longitudinal` | TRUE/FALSE | Include random subject effect in mu |
| `--log_transform` | TRUE/FALSE | Apply log(y+1) before z-standardising |
| `--family_order` | string | Comma-separated family list e.g. `"SHASH,GG,NO"` |
| `--n_cores` | integer | Parallel cores (1 = sequential) |

### `02_infer.R`

| Argument | Type | Description |
|---|---|---|
| `--config` | path | Path to params.yml |
| `--data` | path | Input CSV |
| `--models` | path | Fitted models directory (from fit stage) |
| `--output` | path | Output directory for harmonised CSVs |
| `--features` | path | Features .txt file |
| `--feature` | string | Single feature name |
| `--batch_var` | string | Batch column name (must match fit stage) |
| `--log_transform` | TRUE/FALSE | Validated against `_scaling.csv` from fit |
| `--normative` | TRUE/FALSE | Compute normative z-scores and centiles |
| `--n_cores` | integer | Parallel cores |

### `03_plot.R`

| Argument | Type | Description |
|---|---|---|
| `--config` | path | Path to params.yml |
| `--pre` | path | Pre-harmonisation CSV (optional) |
| `--post` | path | Combined harmonised CSV from infer stage (optional) |
| `--output` | path | Output directory for PNGs |
| `--features` | path | Features .txt file |
| `--feature` | string | Single feature name |
| `--batch_var` | string | Batch column name |
| `--group_col` | string | Column to colour/group trajectories by |
| `--smooth_method` | string | `loess` or `gam` |
| `--age_bin_width` | number | Age bin width [5] |
| `--fix_y_limits` | TRUE/FALSE | Lock y-axis across pre/post panels [TRUE] |

---

## Configuration: params.yml

All arguments have equivalents in `config/params.yml`. CLI always takes priority over the config file.

### Formula specification

Formulas are specified as term lists per distributional parameter (`mu`, `sigma`, `nu`, `tau`). Use `{batch}` and `{id}` as placeholders -- substituted at runtime with `batch_var` and `id_var`.

| YAML value | R formula | Meaning |
|---|---|---|
| `null` | `~ 1` | Intercept only |
| `["pb(age)", "sex"]` | `~ pb(age) + sex` | Fixed effects, no batch |
| `["random({batch})"]` | `~ random(cohort)` | Batch as random effect only |
| `["pb(age)", "sex", "random({batch})"]` | `~ pb(age) + sex + random(cohort)` | Full specification |
| `["{batch}"]` | `~ cohort` | Batch as fixed effect |

Term lists are used exactly as written -- no implicit terms are added. `null` always means `~ 1`.

### Family fallback

Families are tried in `family_order` order, stopping at first convergence. Within each multi-parameter family, nu/tau specs are tried in decreasing complexity automatically:

```
SHASH:  nu=user+tau=user -> nu=user+tau=~1 -> nu=~1+tau=user -> nu=~1+tau=~1
GG:     nu=user          -> nu=~1
NO:     fits directly
```

---

## Examples

### Single feature, end-to-end

```bash
Rscript scripts/01_fit.R   --data data/raw/my_data.csv --feature ThicknessAvg
Rscript scripts/02_infer.R --data data/raw/my_data.csv --feature ThicknessAvg
Rscript scripts/03_plot.R  --pre  data/raw/my_data.csv \
                            --post output/harmonised/combined_harmonised.csv \
                            --feature ThicknessAvg
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
Rscript scripts/01_fit.R --data data/raw/my_data.csv --features config/my_features.txt
```

### Post-harmonisation trajectories only

```bash
Rscript scripts/03_plot.R --post output/harmonised/combined_harmonised.csv
```

### Custom config file

```bash
Rscript scripts/01_fit.R --config config/hpc_params.yml --data data/raw/my_data.csv
```

---

## Output files

### Fit stage
| File | Description |
|---|---|
| `output/models/feature_<n>/<n>_model.rds` | Fitted GAMLSS model object |
| `output/models/feature_<n>/<n>_scaling.csv` | y_mean, y_sd, log_transform flag |
| `output/models/feature_<n>/<n>_metrics.csv` | AIC, BIC, MSE, converged distribution |
| `output/models/feature_<n>/<n>_predictions.csv` | In-sample fitted values |
| `output/models/feature_<n>/<n>_diagnostics.pdf` | Residual diagnostic plots |
| `output/models/model_summary.csv` | Summary across all features |

### Infer stage
| File | Description |
|---|---|
| `output/harmonised/feature_<n>/<n>_harmonised.csv` | Per-feature full output (harmonised values, normative scores, diagnostics) |
| `output/harmonised/combined_harmonised.csv` | Wide format: one row per subject, one column per feature (harmonised values only) |
| `output/harmonised/combined_normative.csv` | Wide format: one row per subject, z-score and centile columns per feature |
| `output/harmonised/harmonisation_summary.csv` | Summary across all features |

### Plot stage
| File | Description |
|---|---|
| `output/plots/<n>_trajectory.png` | Age trajectory (pre and/or post panels) |

---

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| `gamlss` | >= 5.4 | Model fitting |
| `gamlss.dist` | >= 6.1 | Distribution families (SHASH, GG, NO, etc.) |
| `dplyr` | >= 1.1 | Data manipulation |
| `ggplot2` | >= 3.4 | Plotting |
| `scales` | >= 1.3 | Axis formatting |
| `RColorBrewer` | >= 1.1 | Colour palettes |
| `yaml` | >= 2.3 | Config file parsing |
| `optparse` | >= 1.7 | Argument parsing |
| `logger` | >= 0.3 | Structured logging |
| `parallel` | base R | Parallel processing (no installation needed) |

Exact versions used in development are pinned in `environment.yml`.

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
