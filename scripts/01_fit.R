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
source(file.path(ROOT, "R", "fit.R"))

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
    help    = "Path to raw input CSV. Must contain: id, age, sex, <batch_var>, and feature columns. [required]"),

  make_option("--output",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Base output directory. Subdirectories models/, logs/ are created automatically. [default: output/]"),

  make_option("--features_file",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to a .txt file listing feature column names to fit, one per line. Overridden by --one_feature."),

  make_option("--one_feature",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Single feature column name to fit. Overrides --features_file and config."),

  make_option("--feature_families",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = paste(
      "Path to a feature recommendations CSV (e.g. output/diagnostics/feature_recommendations.csv",
      "produced by 'gamlssHarmo diagnose', or a manually authored file).",
      "Required column: feature.",
      "Optional columns: family_order (semicolon-separated, e.g. 'ZINBI;NBI;PO'),",
      "mu_terms, sigma_terms, nu_terms, tau_terms (semicolon-separated formula terms).",
      "Recommendations win over CLI and config for any column that is present.",
      "A column present but empty overrides to intercept-only (~ 1).",
      "A column absent from the CSV falls back to the config value.",
      "discrete is not read from this file; use --discrete to control it globally.",
      "Features not listed in the CSV use CLI/config values entirely."
    )),

  make_option("--suffix",
    type    = "character",
    default = NULL,
    metavar = "STRING",
    help    = paste(
      "String appended to every auto-derived subdirectory name.",
      "E.g. --suffix _disc produces models_disc/ and logs_disc/.",
      "Has no effect on paths supplied explicitly via other flags.",
      "[default: (empty)]"
    )),

  make_option("--batch_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name identifying the batch/site/scanner variable. [default: batch]"),

  make_option("--id_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name for the subject identifier. [default: id]"),

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

  make_option("--longitudinal",
    type    = "character",
    default = NULL,
    metavar = "TRUE/FALSE",
    help    = "If TRUE, adds a random subject effect to mu. [default: FALSE]"),

  make_option("--log_transform",
    type    = "character",
    default = NULL,
    metavar = "TRUE/FALSE",
    help    = "[continuous only] Apply log(y+1) before fitting. Ignored when discrete=TRUE. [default: FALSE]"),

  make_option("--discrete",
    type    = "character",
    default = NULL,
    metavar = "TRUE/FALSE",
    help    = paste(
      "If TRUE, treats all features as count/discrete data.",
      "Skips z-scoring and log transformation. Passes raw integer counts",
      "to gamlss(). Supported families: PO, NBI, NBII, ZIP, ZINBI, ZANBI.",
      "Per-feature 'discrete' flags from --feature_families take priority.",
      "[default: FALSE]"
    )),

  make_option("--family_order",
    type    = "character",
    default = NULL,
    metavar = "LIST",
    help    = paste(
      "Comma-separated global family fallback order. Overrides config and mode defaults.",
      "Continuous default: SHASH,GG,NO. Discrete default: ZINBI,NBI,PO.",
      "Overridden per-feature by --feature_families."
    )),

  make_option("--n_cores",
    type    = "integer",
    default = NULL,
    metavar = "N",
    help    = "Number of parallel cores. 1 runs sequentially. [default: 1]")
)

parser <- OptionParser(
  usage       = "gamlssHarmo fit [options]\n       Rscript scripts/01_fit.R [options]",
  option_list = option_list,
  description = paste(
    "Fit hierarchical GAMLSS harmonisation models, one per feature.",
    "",
    "CONTINUOUS mode (default):",
    "  Features are fitted on the raw scale. Families must operate on the",
    "  full real line (SHASH, JSU, TF, SN1, PE2, NO).",
    "",
    "DISCRETE / COUNT mode (--discrete TRUE):",
    "  Raw integer counts are passed to gamlss().",
    "  Supported families: PO, NBI, NBII, ZIP, ZINBI, ZANBI.",
    "",
    "Output (written to <o>/models/):",
    "  feature_<n>/<n>_model.rds       Fitted model object",
    "  feature_<n>/<n>_meta.csv        discrete and log_transform flags",
    "  feature_<n>/<n>_metrics.csv     AIC, BIC, MSE, distribution used",
    "  feature_<n>/<n>_predictions.csv In-sample fitted values",
    "  feature_<n>/<n>_diagnostics.pdf Residual diagnostic plots",
    "  model_summary.csv               Summary across all features",
    "",
    "Examples:",
    "  gamlssHarmo fit --data data/raw/my_data.csv",
    "  gamlssHarmo fit --data data/raw/my_data.csv --discrete TRUE",
    "  gamlssHarmo fit --data data/raw/my_data.csv \\",
    "    --feature_families output/diagnostics/feature_recommendations.csv",
    sep = "\n"
  ),
  epilogue = "All options can also be set in config/params.yml. CLI arguments always take priority."
)
opt <- parse_args(parser)

cfg_path <- opt$config %||% file.path(ROOT, "config", "params.yml")
cfg      <- load_config(cfg_path)

suffix        <- resolve_arg(opt$suffix, cfg$output$suffix, "")
raw_csv       <- resolve_arg(opt$data,   cfg$data$raw_csv)
base_dir      <- normalizePath(resolve_arg(opt$output, cfg$output$base, "output"), mustWork = FALSE)
models_dir    <- file.path(base_dir, paste0("models", suffix))
log_dir       <- file.path(base_dir, paste0("logs",   suffix))
families_csv  <- resolve_arg(opt$feature_families, cfg$features$feature_families_csv)
batch_var     <- resolve_arg(opt$batch_var,    cfg$model$batch_var,   "batch")
id_var        <- resolve_arg(opt$id_var,       cfg$model$id_var,      "id")
age_var       <- resolve_arg(opt$age_var,      cfg$model$age_var,     "age")
sex_var       <- resolve_arg(opt$sex_var,      cfg$model$sex_var,     "sex")
longitudinal  <- as.logical(resolve_arg(opt$longitudinal,  cfg$model$longitudinal,  FALSE))
log_transform <- as.logical(resolve_arg(opt$log_transform, cfg$model$log_transform, FALSE))
discrete      <- as.logical(resolve_arg(opt$discrete,      cfg$model$discrete,      FALSE))
n_cores       <- as.integer(resolve_arg(opt$n_cores, cfg$compute$n_cores, 1L))

family_order <- if (!is.null(opt$family_order)) {
  trimws(strsplit(opt$family_order, ",")[[1L]])
} else if (discrete) {
  cfg$model$family_order_discrete %||% c("ZINBI", "NBI", "PO")
} else {
  cfg$model$family_order %||% c("SHASH", "GG", "NO")
}

setup_logging(log_dir, "fit")

log_info("=== gamlssHarmo fit ===")
log_info(paste0("config:          ", cfg_path))
log_info(paste0("data:            ", raw_csv %||% "(not set)"))
log_info(paste0("output base:     ", base_dir))
log_info(paste0("suffix:          ", if (nzchar(suffix)) suffix else "(none)"))
log_info(paste0("models dir:      ", models_dir))
log_info(paste0("batch_var:       ", batch_var))
log_info(paste0("id_var:          ", id_var))
log_info(paste0("age_var:         ", age_var))
log_info(paste0("sex_var:         ", sex_var))
log_info(paste0("longitudinal:    ", longitudinal))
log_info(paste0("discrete:        ", discrete))
log_info(paste0("log_transform:   ", log_transform))
log_info(paste0("n_cores:         ", n_cores))
log_info(paste0("feature_families:", families_csv %||% "(not provided)"))

if (is.null(raw_csv) || !nzchar(raw_csv)) {
  print_help(parser)
  quit(status = 1)
}
if (!file.exists(raw_csv))
  stop("Data file not found: ", raw_csv)

data      <- read.csv(raw_csv, stringsAsFactors = FALSE)
data      <- normalise_column_names(data, age_var, sex_var)
meta_cols <- get_meta_cols(batch_var, id_var, data)
log_info(paste0("Loaded: ", nrow(data), " rows x ", ncol(data), " cols"))

features <- if (!is.null(opt$one_feature) && nzchar(opt$one_feature)) {
  opt$one_feature
} else if (!is.null(opt$features_file) && nzchar(opt$features_file)) {
  if (!file.exists(opt$features_file))
    stop("Features file not found: ", opt$features_file)
  read_features_file(opt$features_file)
} else if (!is.null(cfg$features$features_txt) &&
           nzchar(cfg$features$features_txt) &&
           file.exists(cfg$features$features_txt)) {
  read_features_file(cfg$features$features_txt)
} else {
  setdiff(names(data), meta_cols)
}

log_info(paste0("Features: ", length(features)))

feature_recs <- NULL
if (!is.null(families_csv) && nzchar(trimws(families_csv))) {
  if (!file.exists(families_csv))
    stop("Feature families CSV not found: ", families_csv)
  feature_recs <- read_feature_recommendations(families_csv)
  n_match <- sum(features %in% names(feature_recs))
  log_info(paste0("Per-feature recommendations matched to ",
                   n_match, "/", length(features), " features"))
}

formula_terms <- parse_and_validate_formulas(cfg, data, longitudinal,
                                              batch_var = batch_var,
                                              id_var    = id_var)
log_info("Formulas validated")

results <- run_gamlss_harmonisation(
  data          = data,
  features      = features,
  output_dir    = models_dir,
  formula_terms = formula_terms,
  batch_var     = batch_var,
  id_var        = id_var,
  longitudinal  = longitudinal,
  log_transform = log_transform,
  discrete      = discrete,
  family_order  = family_order,
  feature_recs  = feature_recs,
  n_cores       = n_cores
)

log_info(paste0("=== fit complete === success: ", results$successes,
                " | failed: ", results$failures))
log_info(paste0("time: ", results$total_time))