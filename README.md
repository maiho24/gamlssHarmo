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

**Option A — Conda (recommended)**
```bash
conda env create -f environment.yml
conda activate gamlssHarmo
Rscript install_dependencies.R
```
The conda environment installs all R packages. `install_dependencies.R` then registers the `gamlssHarmo` CLI command inside the active environment.

**Option B — Plain R (R already installed)**
```bash
Rscript install_dependencies.R
```
Installs any missing R packages from CRAN and registers the `gamlssHarmo` command:
- Linux/macOS: symlink placed in `~/.local/bin/`
- Windows: `gamlssHarmo.cmd` created in the repo root (add the repo root to `PATH` to use `gamlssHarmo` from any directory)

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

---

## Options

Each stage exposes its own command-line arguments. The full per-stage reference —
every flag, type, default, and description for `diagnose`, `fit`, `infer`, and
`plot` — lives in **[`docs/options.md`](docs/options.md)**.

For quick help at the terminal:
```bash
gamlssHarmo diagnose --help
gamlssHarmo fit      --help
gamlssHarmo infer    --help
gamlssHarmo plot     --help
```

---

## Configuring family & formula selection

For each feature, gamlssHarmo needs to know **which distribution families to try** and
**what formula to use for each distributional parameter** (`mu`, `sigma`, `nu`, `tau`).
This can be set at two levels:

| Level | Where | Scope | Use when |
|---|---|---|---|
| **Global rule** | `config/params.yml` (or CLI flags) | One family order + one formula applied to **every** feature | All features share the same model, or as a fallback |
| **Per-feature rules** | a `feature_recommendations.csv` passed via `--feature_families` | One row per feature, overriding the global rule for that feature | Features differ in distributional shape (the usual case for multi-feature studies) |

Precedence: CLI flags override `config/params.yml`, and a per-feature CSV row overrides
both — **column by column**. Any column omitted from the CSV falls back to the global rule.

### Global rule — `config/params.yml`

`family_order` lists the families to try (see [Family fallback](#family-fallback)). The
`formulas` block sets the term list for each parameter, applied to all features. Use
`{batch}` and `{id}` as placeholders — substituted at runtime with `batch_var` and `id_var`.

| YAML value                                      | R formula                           | Meaning                         |
|-------------------------------------------------|-------------------------------------|---------------------------------|
| `null`                                          | `~ 1`                               | Intercept only                  |
| `["pb(age)", "sex"]`                            | `~ pb(age) + sex`                   | Fixed effects, no batch         |
| `["random({batch})"]`                           | `~ random(batch)`                   | Batch as random effect          |
| `["pb(age)", "sex", "random({batch})"]`         | `~ pb(age) + sex + random(batch)`   | Full specification              |
| `["{batch}"]`                                   | `~ batch`                           | Batch as fixed effect           |

Term lists are used exactly as written — no implicit terms are added. `null` means `~ 1`.

### Per-feature rules — `feature_recommendations.csv`

When features differ in shape, a per-feature table is more convenient than a single global
rule. The `diagnose` stage generates this file automatically (see
[`docs/diagnose.md`](docs/diagnose.md)), or you can write it by hand. Pass it to fit with:

```bash
gamlssHarmo fit --data data/my_data.csv \
  --feature_families output/diagnostics/feature_recommendations.csv
```

One row per feature, with these columns:

| Column | Meaning |
|---|---|
| `feature` | Feature column name — must match the data CSV exactly |
| `family_order` | Semicolon-separated families to try, in order |
| `mu_terms`, `sigma_terms`, `nu_terms`, `tau_terms` | Semicolon-separated formula terms per parameter. Empty cell → intercept-only (`~ 1`); omit the whole column → use the global rule for that parameter |
| `tier` | Informational diagnose tier (1–4); not read by fit |

**A ready-to-edit template is included in the repository:
[`examples/feature_recommendations_example.csv`](examples/feature_recommendations_example.csv).**
Copy it, replace the `feature` values with your own column names, and adjust each row. Its
contents:

```csv
feature,family_order,mu_terms,sigma_terms,nu_terms,tau_terms,tier
smri_thick_cdk_locclh,GG;BCCG;SHASH;NO,pb(age);sex;random({batch}),random({batch}),,,1
smri_area_cdk_cdacatelh,SHASH;JSU;SN1;NO,pb(age);sex;random({batch}),pb(age);sex;random({batch}),,,2
smri_area_cdk_fusiformlh,BCT;BCPE;SHASH;GG;BCCG;NO,pb(age);sex;random({batch}),pb(age);sex;random({batch}),random({batch}),random({batch}),3
```

Reading the rows: the first feature is near-Gaussian (Tier 1 — simplified `sigma`,
intercept-only `nu`/`tau`); the second adds age/sex to `sigma` (Tier 2); the third has batch
effects on both skewness and kurtosis (Tier 3 — `nu` and `tau` each get a batch random
effect). A blank cell means intercept-only for that parameter.

### Supported families

**Continuous** (default mode):

| Params | Families                              |
|--------|---------------------------------------|
| 4      | SHASH, JSU, BCT, BCPE, ST3, ST4, EGB2 |
| 3      | GG, SN1, TF, PE2, BCCG, exGAUS       |
| 2      | NO, LOGNO                             |

**Discrete / count** (`--discrete TRUE`):

| Params | Families                                          |
|--------|---------------------------------------------------|
| 4      | ZISICHEL, ZASICHEL                                |
| 3      | ZINBI, ZANBI, ZIPIG, SICHEL                        |
| 2      | NBI, NBII, ZIP, PIG, DPO                          |
| 1      | PO                                                |

### Family fallback

The first model that converges is used. **Continuous** families are searched
**breadth-first by structural richness**: every family is first tried at its richest
specification (nu/tau modelled), and the shape parameters are simplified only once no
family converges at that level. A batch effect on nu/tau is therefore dropped only when
no candidate family can support it — not the moment the top-ranked family struggles.
Within a richness level, families follow `family_order`. `mu` is never simplified, and
`sigma` (which carries the batch scale effect) is never simplified for continuous families.

**Discrete** families are tried per family in `family_order`, exhausting each family's
fallback ladder before moving to the next.

<details>
<summary>Per-family fallback ladders (click to expand)</summary>

The table lists each family's own ladder. For continuous families these ladders are
interleaved across families by richness (all `_full` specs before any `_no_tau`, etc.);
for discrete families each ladder is run to exhaustion before the next family.

| Group | Fallback ladder |
|---|---|
| 4-param continuous (SHASH, JSU, BCT, BCPE, ST3, ST4, EGB2) | `_full` → `_no_tau` → `_no_nu` → `_intercept` |
| 3-param continuous (GG, TF, SN1, PE2, BCCG, exGAUS) | `_full` → `_intercept` |
| 2-param continuous (NO, LOGNO) | single spec |
| 4-param discrete (ZISICHEL, ZASICHEL) | `_full` → `_no_tau` → `_no_nu` → `_no_nu_no_tau` → `_no_sigma` → `_no_sigma_no_tau` → `_no_sigma_no_nu` → `_intercept` |
| 3-param discrete (ZINBI, ZANBI, ZIPIG, SICHEL) | `_full` → `_no_nu` → `_no_sigma` → `_intercept` |
| 2-param discrete (NBI, NBII, ZIP, PIG, DPO) | `_full` → `_no_sigma` |
| 1-param (PO) | single spec |

Spec name glossary: `_full` = all extra params have their full formula; `_no_X` = parameter X reduced to intercept-only (`~ 1`); `_intercept` = all extra params at intercept.

</details>

---

## Examples

### Full recommended workflow
```bash
# 1. Diagnose (optional)
gamlssHarmo diagnose --data data/my_data.csv --output output/

# 2. Fit using diagnose recommendations
gamlssHarmo fit --data data/my_data.csv \
  --output output/ \
  --feature_families output/diagnostics/feature_recommendations.csv

# 3. Harmonise
gamlssHarmo infer --data data/my_data.csv \
                  --output output/ \
                  --model output/models

# 4. Plot
gamlssHarmo plot  --pre  data/my_data.csv \
                  --post output/harmonised/combined_harmonised.csv \
                  --output output/
```

### Skip diagnose — specify families directly

One global family order for all features:
```bash
gamlssHarmo fit --data data/my_data.csv --family_order "GG,SHASH,NO"
```

Or hand-write per-feature rules (copy [`examples/feature_recommendations_example.csv`](examples/feature_recommendations_example.csv) as a template):
```bash
gamlssHarmo fit --data data/my_data.csv \
  --feature_families examples/feature_recommendations_example.csv
```
See [Configuring family & formula selection](#configuring-family--formula-selection) for the column reference.

### Single feature, end-to-end
```bash
gamlssHarmo fit   --data data/my_data.csv --one_feature ThicknessAvg
gamlssHarmo infer --data data/my_data.csv --one_feature ThicknessAvg
gamlssHarmo plot  --pre  data/my_data.csv \
                  --post output/harmonised/combined_harmonised.csv \
                  --one_feature ThicknessAvg
```

### Discrete / count data
```bash
gamlssHarmo diagnose --data data/my_data.csv --discrete TRUE
gamlssHarmo fit      --data data/my_data.csv --discrete TRUE \
  --feature_families output/diagnostics/feature_recommendations.csv
gamlssHarmo infer    --data data/my_data.csv
```

### Custom batch column
```bash
gamlssHarmo fit   --data data/my_data.csv --batch_var site
gamlssHarmo infer --data data/my_data.csv --batch_var site
```

### Namespace outputs with --suffix (useful for comparing runs)
```bash
gamlssHarmo fit   --data data/my_data.csv --suffix _shash --family_order SHASH
gamlssHarmo infer --data data/my_data.csv --suffix _shash
gamlssHarmo plot  --post output/harmonised_shash/combined_harmonised.csv --suffix _shash
```

### Parallel processing
```bash
gamlssHarmo fit   --data data/my_data.csv --n_cores 8
gamlssHarmo infer --data data/my_data.csv --n_cores 8
```
With `--n_cores > 1`, `diagnose`/`fit`/`infer` display a live text progress bar
(feature count + ETA) as workers finish, instead of per-feature log lines —
per-feature logs from worker processes aren't visible in the main console, so
the progress bar is the only live signal during a parallel run. Detailed
per-feature logs (family tried, convergence, AIC, etc.) are still written to
each feature's own output files (e.g. `<feature>_timing.csv`) as before, and
can be inspected once the run finishes, or tailed live from a second terminal.

### Longitudinal data
```bash
gamlssHarmo fit   --data data/my_data.csv --longitudinal TRUE
gamlssHarmo infer --data data/my_data.csv
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

Multiple diagnose runs accumulate into both files by `feature`: re-diagnosing a
feature (e.g. via `--one_feature` or a smaller `--features_file`) updates only
that feature's row; rows for features not included in the current run are
carried forward unchanged, rather than being dropped.

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

`model_summary.csv` is rebuilt each run by rescanning every
`feature_<n>/<n>_timing.csv` under `output/models/` — it always reflects every
feature ever fit into this output directory, not just the features covered by
the current invocation. Re-fitting a feature updates its row in place; fitting
a different subset of features does not remove rows for features fit earlier.

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
| `pbapply`      | >= 1.7    | Live progress bar for `--n_cores > 1` runs   |

Exact versions used in development are pinned in `environment.yml`.

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
