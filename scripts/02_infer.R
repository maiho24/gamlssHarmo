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
    script_path <- normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)
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
    help    = "Path to input CSV to harmonise. Must contain the same columns used during fitting. [required]"),

  make_option("--models",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to fitted models directory (the models/ subfolder from the fit stage). [default: <o>/models]"),

  make_option("--output",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Base output directory. Must match the --output used in the fit stage. [default: output/]"),

  make_option("--features_file",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to a .txt file listing feature names to harmonise. Overridden by --one_feature."),

  make_option("--one_feature",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Single feature name to harmonise. Overrides --features_file and config."),

  make_option("--suffix",
    type    = "character",
    default = NULL,
    metavar = "STRING",
    help    = paste(
      "String appended to every auto-derived subdirectory name.",
      "Must match the --suffix used during the fit stage so models are found.",
      "E.g. --suffix _cont resolves models_cont/ and harmonised_cont/.",
      "Has no effect on paths supplied explicitly via --models.",
      "[default: (empty)]"
    )),

  make_option("--batch_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name identifying the batch/site/scanner variable. Must match the value used during fitting. [default: batch]"),

  make_option("--id_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name for the subject identifier. Must match the value used during fitting. [default: id]"),

  make_option("--age_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name for age. Renamed to 'age' internally. [default: age]"),

  make_option("--sex_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name for sex covariate. Renamed to 'sex' internally. [default: sex]"),

  make_option("--log_transform",
    type    = "character",
    default = NULL,
    metavar = "TRUE/FALSE",
    help    = "Whether log(y+1) was applied during fitting. Validated automatically against the per-feature _meta.csv -- only override if needed. [default: FALSE]"),

  make_option("--normative",
    type    = "character",
    default = "TRUE",
    metavar = "TRUE/FALSE",
    help    = paste(
      "If TRUE, computes normative z-scores and saves them to combined_normative.csv.",
      "For continuous features: z = qnorm(F(y | fitted params)).",
      "For discrete / count features: z = qnorm(mid_p(y)) using the Lancaster (1961)",
      "mid-p correction for discreteness.",
      "[default: TRUE]"
    )),

  make_option("--n_cores",
    type    = "integer",
    default = NULL,
    metavar = "N",
    help    = "Number of parallel cores. 1 runs sequentially. [default: 1]")
)

parser <- OptionParser(
  usage       = "gamlssHarmo infer [options]\n       Rscript scripts/02_infer.R [options]",
  option_list = option_list,
  description = paste(
    "Apply fitted GAMLSS models to harmonise features and compute normative scores.",
    "",
    "The discrete/continuous mode and log_transform flag for each feature are",
    "determined automatically from the per-feature _meta.csv written during",
    "the fit stage. No --discrete flag is required here.",
    "",
    "Output (written to <o>/harmonised/):",
    "  feature_<n>/<n>_harmonised.csv   Per-feature output",
    "  combined_harmonised.csv           Wide table: harmonised values",
    "  combined_normative.csv            Wide table: normative z-scores",
    "",
    "Multiple infer runs accumulate into the combined CSVs; new columns are",
    "appended and existing ones overwritten by the newer result.",
    "",
    "Examples:",
    "  gamlssHarmo infer --data data/raw/my_data.csv --output output/",
    "  gamlssHarmo infer --data data/eval.csv --output output/ --one_feature ThicknessAvg",
    "  gamlssHarmo infer --data data/eval.csv --output output/ --normative FALSE",
    sep = "\n"
  ),
  epilogue = "All options can also be set in config/params.yml. CLI arguments always take priority."
)
opt <- parse_args(parser)

cfg_path <- opt$config %||% file.path(ROOT, "config", "params.yml")
cfg      <- load_config(cfg_path)

suffix           <- resolve_arg(opt$suffix, cfg$output$suffix, "")
raw_csv          <- resolve_arg(opt$data,   cfg$data$raw_csv)
base_dir         <- normalizePath(resolve_arg(opt$output, cfg$output$base, "output"), mustWork = FALSE)
models_dir       <- resolve_arg(opt$models, file.path(base_dir, paste0("models",     suffix)))
harm_dir         <- file.path(base_dir, paste0("harmonised", suffix))
log_dir          <- file.path(base_dir, paste0("logs",       suffix))
batch_var        <- resolve_arg(opt$batch_var, cfg$model$batch_var, "batch")
id_var           <- resolve_arg(opt$id_var,   cfg$model$id_var,   "id")
age_var          <- resolve_arg(opt$age_var,  cfg$model$age_var,  "age")
sex_var          <- resolve_arg(opt$sex_var,  cfg$model$sex_var,  "sex")
log_transform    <- as.logical(resolve_arg(opt$log_transform, cfg$model$log_transform, FALSE))
normative_scores <- as.logical(opt$normative %||% "TRUE")
n_cores          <- as.integer(resolve_arg(opt$n_cores, cfg$compute$n_cores, 1L))

setup_logging(log_dir, "infer")

log_info("=== gamlssHarmo infer ===")
log_info(paste0("config:          ", cfg_path))
log_info(paste0("data:            ", raw_csv %||% "(not set)"))
log_info(paste0("output base:     ", base_dir))
log_info(paste0("suffix:          ", if (nzchar(suffix)) suffix else "(none)"))
log_info(paste0("models dir:      ", models_dir))
log_info(paste0("harmonised dir:  ", harm_dir))
log_info(paste0("batch_var:       ", batch_var))
log_info(paste0("id_var:          ", id_var))
log_info(paste0("age_var:         ", age_var))
log_info(paste0("sex_var:         ", sex_var))
log_info(paste0("normative:       ", normative_scores))
log_info(paste0("n_cores:         ", n_cores))
log_info("discrete mode: determined per-feature from _meta.csv")

if (is.null(raw_csv) || !nzchar(raw_csv)) {
  print_help(parser)
  quit(status = 1)
}
if (!file.exists(raw_csv))
  stop("Data file not found: ", raw_csv)
if (!dir.exists(models_dir))
  stop("Models directory not found: ", models_dir, "\nRun 'gamlssHarmo fit' first.")

data <- read.csv(raw_csv, stringsAsFactors = FALSE)
data <- normalise_column_names(data, age_var, sex_var)
log_info(paste0("Loaded: ", nrow(data), " rows x ", ncol(data), " cols"))

log_info(paste0("opt$one_feature:   ", opt$one_feature  %||% "(not set)"))
log_info(paste0("opt$features_file: ", opt$features_file %||% "(not set)"))

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
  results$results[vapply(results$results,
                         function(x) identical(x$status, "success"),
                         logical(1L))],
  function(x) x$feature
)

if (length(features_processed) > 0L) {
  combine_harmonised_results(
    harmonised_output_dir     = harm_dir,
    id_var                    = id_var,
    batch_var                 = batch_var,
    generate_normative_scores = normative_scores,
    feature_subset            = features_processed
  )
} else {
  log_info("No successful features to combine.")
}

log_info(paste0("=== infer complete === success: ", results$successes,
                " | failed: ", results$failures))
log_info(paste0("combined CSV: ", file.path(harm_dir, "combined_harmonised.csv")))