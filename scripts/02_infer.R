#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(gamlss)
  library(gamlss.dist)
  library(dplyr)
  library(logger)
})

get_root <- function() {
  if (exists("SCRIPT_DIR", envir = parent.env(environment()), inherits = TRUE))
    return(normalizePath(get("SCRIPT_DIR", envir = parent.env(environment()),
                             inherits = TRUE)))
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
    return(dirname(dirname(script_path)))
  }
  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)
    if (exists("ofile", envir = f, inherits = FALSE)) {
      script_path <- normalizePath(get("ofile", envir = f), mustWork = FALSE)
      return(dirname(dirname(script_path)))
    }
  }
  normalizePath(".")
}

ROOT <- get_root()
source(file.path(ROOT, "R", "utils.R"))
source(file.path(ROOT, "R", "infer.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

option_list <- list(
  make_option("--config",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to params.yml config file. [default: config/params.yml]"),

  make_option("--data",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to input CSV to harmonise. Can be the same as the training data or a held-out dataset. Must contain the same columns used during fitting. [required]"),

  make_option("--models",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to fitted models directory (the models/ subfolder produced by the fit stage). Overrides the path derived from --output. [default: <output>/models]"),

  make_option("--output",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Base output directory. Subdirectories harmonised/, logs/ are created automatically beneath it. Must match the --output used in the fit stage so models are found automatically. [default: output/]"),

  make_option("--features-file",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to a .txt file listing feature names to harmonise, one per line. Lines starting with '#' are ignored. Overridden by --feature."),

  make_option("--one-feature",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Single feature name to harmonise. Overrides --features and config. Results are merged into combined_harmonised.csv alongside any features from previous runs."),

  make_option("--batch_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name identifying the batch/site/scanner variable. Must match the value used during fitting. [default: cohort]"),

  make_option("--log_transform",
    type    = "character",
    default = NULL,
    metavar = "TRUE/FALSE",
    help    = "Whether log(y+1) was applied during fitting. Validated automatically against the scaling file saved during fit -- only set this if you need to override. [default: FALSE]"),

  make_option("--normative",
    type    = "character",
    default = "TRUE",
    metavar = "TRUE/FALSE",
    help    = "If TRUE, computes normative z-scores for each subject relative to the fitted distribution and saves them to combined_normative.csv (one column per feature, named after the feature). [default: TRUE]"),

  make_option("--n_cores",
    type    = "integer",
    default = NULL,
    metavar = "N",
    help    = "Number of parallel cores to use. 1 runs sequentially. [default: 1]")
)

opt <- parse_args(OptionParser(
  usage       = "gamlssHarmo infer [options]\n       Rscript scripts/02_infer.R [options]",
  option_list = option_list,
  description = paste(
    "Apply fitted GAMLSS models to harmonise features and compute normative scores.",
    "",
    "Harmonisation removes batch/site effects by mapping each observation through",
    "the fitted CDF under the original distribution parameters, then inverting",
    "through the CDF under batch-effect-removed parameters.",
    "",
    "Runs on the same or a new dataset (e.g. a held-out evaluation set). Models",
    "must have been fitted first using the fit stage.",
    "",
    "Output (written to <output>/harmonised/):",
    "  feature_<n>/<n>_harmonised.csv  Per-feature full output with harmonised",
    "                                   values, CDF diagnostics, and normative scores",
    "  combined_harmonised.csv          Wide table: one row per subject, one column",
    "                                   per feature (harmonised values)",
    "  combined_normative.csv           Wide table: one row per subject, one column",
    "                                   per feature (normative z-scores)",
    "",
    "Multiple infer runs accumulate into combined_harmonised.csv and",
    "combined_normative.csv -- new feature columns are appended and existing",
    "columns are overwritten by the newer result.",
    "",
    "Examples:",
    "  gamlssHarmo infer --data data/raw/my_data.csv --output output/",
    "  gamlssHarmo infer --data data/eval.csv --output output/ --feature ThicknessAvg",
    "  gamlssHarmo infer --data data/eval.csv --output output/ --normative FALSE",
    sep = "\n"
  ),
  epilogue = "All options can also be set in config/params.yml. CLI arguments always take priority."
))

cfg_path <- opt$config %||% file.path(ROOT, "config", "params.yml")
cfg      <- load_config(cfg_path)

raw_csv          <- resolve_arg(opt$data, cfg$data$raw_csv)
base_dir         <- normalizePath(resolve_arg(opt$output, cfg$output$base, "output"), mustWork = FALSE)
models_dir       <- resolve_arg(opt$models, file.path(base_dir, "models"))
harm_dir         <- file.path(base_dir, "harmonised")
log_dir          <- file.path(base_dir, "logs")
features_txt     <- resolve_arg(opt$features_file, cfg$features$features_txt)
feature_single   <- resolve_arg(opt$one_feature, cfg$features$feature)
batch_var        <- resolve_arg(opt$batch_var, cfg$model$batch_var, "cohort")
id_var           <-                            cfg$model$id_var %||% "id"
log_transform    <- as.logical(resolve_arg(opt$log_transform, cfg$model$log_transform, FALSE))
normative_scores <- as.logical(opt$normative %||% "TRUE")
n_cores          <- as.integer(resolve_arg(opt$n_cores, cfg$compute$n_cores, 1L))

setup_logging(log_dir, "infer")

log_info("=== gamlssHarmo infer ===")
log_info(paste0("config:          ", cfg_path))
log_info(paste0("data:            ", raw_csv %||% "(not set)"))
log_info(paste0("output base:     ", base_dir))
log_info(paste0("models dir:      ", models_dir))
log_info(paste0("harmonised dir:  ", harm_dir))
log_info(paste0("batch_var:       ", batch_var))
log_info(paste0("normative:       ", normative_scores))
log_info(paste0("n_cores:         ", n_cores))

if (is.null(raw_csv) || !nzchar(raw_csv))
  stop("Provide --data or set data.raw_csv in params.yml")
if (!file.exists(raw_csv))
  stop("Data file not found: ", raw_csv)
if (!dir.exists(models_dir))
  stop("Models directory not found: ", models_dir, "\nRun 'gamlssHarmo fit' first.")

data <- read.csv(raw_csv, stringsAsFactors = FALSE)
log_info(paste0("Loaded: ", nrow(data), " rows x ", ncol(data), " cols"))

# Resolve feature subset: CLI --feature takes priority over --features file,
# which takes priority over config, which means process all available models.
log_info(paste0("opt$feature:   ", opt$one_feature  %||% "(not set)"))
log_info(paste0("opt$features:  ", opt$features_file %||% "(not set)"))

feature_subset <- if (!is.null(opt$one_feature) && nzchar(trimws(opt$one_feature))) {
  trimws(opt$one_feature)
} else if (!is.null(opt$features_file) && nzchar(trimws(opt$features_file))) {
  fpath <- trimws(opt$features_file)
  if (!file.exists(fpath))
    stop("Features file not found: ", fpath)
  read_features_file(fpath)
} else if (!is.null(cfg$features$features_txt) &&
           nzchar(cfg$features$features_txt %||% "") &&
           file.exists(cfg$features$features_txt)) {
  read_features_file(cfg$features$features_txt)
} else {
  NULL
}

if (!is.null(feature_subset))
  log_info(paste0("Feature subset (", length(feature_subset), "): ",
                  paste(feature_subset, collapse = ", ")))

results <- harmonise_all_gamlss_models(
  model_base_dir            = models_dir,
  data                      = data,
  output_dir                = harm_dir,
  batch_var                 = batch_var,
  id_var                    = id_var,
  log_transform             = log_transform,
  generate_normative_scores = normative_scores,
  feature_subset            = feature_subset,
  n_cores                   = n_cores
)

features_processed <- sapply(
  results$results[sapply(results$results, function(x) x$status == "success")],
  function(x) x$feature
)

if (length(features_processed) > 0) {
  combine_harmonised_results(
    harmonised_output_dir     = harm_dir,
    id_var                    = id_var,
    generate_normative_scores = normative_scores,
    feature_subset            = features_processed
  )
} else {
  log_info("No successful features to combine.")
}

log_info(paste0("=== infer complete === success: ", results$successes,
                " | failed: ", results$failures))
log_info(paste0("combined CSV: ", file.path(harm_dir, "combined_harmonised.csv")))