#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
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
source(file.path(ROOT, "R", "diagnose.R"))

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
    help    = "Base output directory. Diagnostics written to <o>/diagnostics/. [default: output/]"),

  make_option("--features_file",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to a .txt file listing feature column names to diagnose, one per line. Lines starting with '#' are ignored."),

  make_option("--one_feature",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Single feature column name to diagnose. Overrides --features_file and config."),

  make_option("--feature_families",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = paste(
      "Path to a CSV with per-feature overrides.",
      "Required column: feature.",
      "Optional columns: family_order (semicolon-separated, e.g. 'GG;SHASH;NO'),",
      "nu_terms, tau_terms (semicolon-separated formula terms,",
      "e.g. 'random({batch})').",
      "mu and sigma formulas are set by the tier logic and cannot be overridden here."
    )),

  make_option("--suffix",
    type    = "character",
    default = NULL,
    metavar = "STRING",
    help    = paste(
      "String appended to every auto-derived subdirectory name.",
      "E.g. --suffix _cont produces diagnostics_cont/ and logs_cont/.",
      "Has no effect on paths supplied explicitly via other flags.",
      "[default: (empty)]"
    )),

  make_option("--batch_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name identifying the batch/site variable. [default: batch]"),

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
    help    = "If TRUE, adds random({id}) to the recommended mu formula. [default: FALSE]"),

  make_option("--discrete",
    type    = "character",
    default = NULL,
    metavar = "TRUE/FALSE",
    help    = paste(
      "If TRUE, treats all features as count/discrete data.",
      "Runs zero-inflation and overdispersion tests in place of skewness/kurtosis.",
      "Recommends from: PO, NBI, NBII, ZIP, ZINBI, ZANBI.",
      "When FALSE (default), runs the continuous Cullen-Frey pipeline."
    )),

  make_option("--n_perm",
    type    = "integer",
    default = NULL,
    metavar = "N",
    help    = "Number of permutations for batch-moment tests. [default: 499]"),

  make_option("--min_batch_n",
    type    = "integer",
    default = NULL,
    metavar = "N",
    help    = "Minimum observations per batch for permutation tests. [default: 50]"),

  make_option("--thresh_skew",
    type    = "numeric",
    default = NULL,
    metavar = "N",
    help    = "[continuous] |skewness| threshold for Tier 1 classification. [default: 0.5]"),

  make_option("--thresh_kurt",
    type    = "numeric",
    default = NULL,
    metavar = "N",
    help    = "[continuous] |excess kurtosis| threshold for Tier 1 classification. [default: 1.0]"),

  make_option("--thresh_pval",
    type    = "numeric",
    default = NULL,
    metavar = "N",
    help    = "P-value threshold for batch higher-moment tests. [default: 0.01]"),

  make_option("--thresh_skew_range",
    type    = "numeric",
    default = NULL,
    metavar = "N",
    help    = "[continuous] Minimum between-batch skewness range to declare a skewness batch effect. [default: 0.3]"),

  make_option("--thresh_kurt_range",
    type    = "numeric",
    default = NULL,
    metavar = "N",
    help    = "[continuous] Minimum between-batch excess kurtosis range. [default: 0.5]"),

  make_option("--thresh_zerorate_range",
    type    = "numeric",
    default = NULL,
    metavar = "N",
    help    = "[discrete] Minimum between-batch zero-rate range for Tier 3. [default: 0.10]"),

  make_option("--thresh_dispersion_range",
    type    = "numeric",
    default = NULL,
    metavar = "N",
    help    = "[discrete] Minimum between-batch dispersion-ratio range. [default: 0.50]"),

  make_option("--n_cores",
    type    = "integer",
    default = NULL,
    metavar = "N",
    help    = "Number of parallel cores. [default: 1]")
)

parser <- OptionParser(
  usage       = "gamlssHarmo diagnose [options]\n       Rscript scripts/00_diagnose.R [options]",
  option_list = option_list,
  description = paste(
    "Diagnose distributional properties of features before fitting.",
    "",
    "CONTINUOUS mode (default):",
    "  Residualises for age/sex, then tests skewness, kurtosis, and batch shape",
    "  effects to assign a complexity tier and recommend a real-line GAMLSS family.",
    "",
    "DISCRETE / COUNT mode (--discrete TRUE):",
    "  Runs zero-inflation and overdispersion tests, then recommends from",
    "  PO, NBI, NBII, ZIP, ZINBI, ZANBI following Hilbe (2014).",
    "",
    "Output (written to <o>/diagnostics/):",
    "  diagnostics_summary.csv       Full statistics for every feature",
    "  feature_recommendations.csv   Trimmed recommendations for 01_fit.R",
    "",
    "Examples:",
    "  gamlssHarmo diagnose --data data/raw/my_data.csv",
    "  gamlssHarmo diagnose --data data/raw/my_data.csv --discrete TRUE",
    "  gamlssHarmo diagnose --data data/raw/my_data.csv --n_perm 999",
    "  gamlssHarmo diagnose --data data/raw/my_data.csv \\",
    "    --feature_families config/my_families.csv",
    sep = "\n"
  ),
  epilogue = "All options can also be set in config/params.yml under a 'diagnostics:' key."
)
opt <- parse_args(parser)

cfg_path <- opt$config %||% file.path(ROOT, "config", "params.yml")
cfg      <- load_config(cfg_path)
dcfg     <- cfg$diagnostics %||% list()

suffix   <- resolve_arg(opt$suffix, cfg$output$suffix, "")

raw_csv      <- resolve_arg(opt$data,  cfg$data$raw_csv)
base_dir     <- normalizePath(resolve_arg(opt$output, cfg$output$base, "output"), mustWork = FALSE)
diag_dir     <- file.path(base_dir, paste0("diagnostics", suffix))
log_dir      <- file.path(base_dir, paste0("logs",        suffix))
families_csv <- resolve_arg(opt$feature_families, cfg$features$feature_families_csv)
batch_var    <- resolve_arg(opt$batch_var,   cfg$model$batch_var,  "batch")
id_var       <- resolve_arg(opt$id_var,      cfg$model$id_var,     "id")
age_var      <- resolve_arg(opt$age_var,     cfg$model$age_var,    "age")
sex_var      <- resolve_arg(opt$sex_var,     cfg$model$sex_var,    "sex")
longitudinal <- as.logical(resolve_arg(opt$longitudinal, cfg$model$longitudinal, FALSE))
discrete     <- as.logical(resolve_arg(opt$discrete,     cfg$model$discrete,     FALSE))
n_perm       <- as.integer(resolve_arg(opt$n_perm,      dcfg$n_perm,      499L))
min_batch_n  <- as.integer(resolve_arg(opt$min_batch_n, dcfg$min_batch_n,   50L))
n_cores      <- as.integer(resolve_arg(opt$n_cores,     cfg$compute$n_cores, 1L))

thresholds <- list(
  skew                 = as.numeric(resolve_arg(opt$thresh_skew,              dcfg$thresh_skew,              0.5)),
  kurt                 = as.numeric(resolve_arg(opt$thresh_kurt,              dcfg$thresh_kurt,              1.0)),
  pval                 = as.numeric(resolve_arg(opt$thresh_pval,              dcfg$thresh_pval,              0.01)),
  min_skew_range       = as.numeric(resolve_arg(opt$thresh_skew_range,        dcfg$thresh_skew_range,        0.3)),
  min_kurt_range       = as.numeric(resolve_arg(opt$thresh_kurt_range,        dcfg$thresh_kurt_range,        0.5)),
  min_zerorate_range   = as.numeric(resolve_arg(opt$thresh_zerorate_range,    dcfg$thresh_zerorate_range,    0.10)),
  min_dispersion_range = as.numeric(resolve_arg(opt$thresh_dispersion_range,  dcfg$thresh_dispersion_range,  0.50))
)

setup_logging(log_dir, "diagnose")

log_info("=== gamlssHarmo diagnose ===")
log_info(paste0("config:          ", cfg_path))
log_info(paste0("data:            ", raw_csv %||% "(not set)"))
log_info(paste0("output base:     ", base_dir))
log_info(paste0("suffix:          ", if (nzchar(suffix)) suffix else "(none)"))
log_info(paste0("diagnostics dir: ", diag_dir))
log_info(paste0("batch_var:       ", batch_var))
log_info(paste0("id_var:          ", id_var))
log_info(paste0("age_var:         ", age_var))
log_info(paste0("sex_var:         ", sex_var))
log_info(paste0("longitudinal:    ", longitudinal))
log_info(paste0("discrete:        ", discrete))
log_info(paste0("n_perm:          ", n_perm))
log_info(paste0("min_batch_n:     ", min_batch_n))
log_info(paste0("n_cores:         ", n_cores))

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

features <- if (!is.null(opt$one_feature) && nzchar(trimws(opt$one_feature))) {
  trimws(opt$one_feature)
} else if (!is.null(opt$features_file) && nzchar(trimws(opt$features_file))) {
  read_features_file(opt$features_file)
} else if (!is.null(cfg$features$features_txt) &&
            nzchar(cfg$features$features_txt %||% "") &&
            file.exists(cfg$features$features_txt)) {
  read_features_file(cfg$features$features_txt)
} else {
  setdiff(names(data), meta_cols)
}

log_info(paste0("Features to diagnose: ", length(features)))

user_overrides <- NULL
if (!is.null(families_csv) && nzchar(trimws(families_csv))) {
  if (!file.exists(families_csv))
    stop("Feature families CSV not found: ", families_csv)
  user_overrides <- read_feature_recommendations_csv(families_csv)
  log_info(paste0("Loaded family overrides for ",
                   length(user_overrides), " features from: ", families_csv))
}

results <- diagnose_all_features(
  data           = data,
  features       = features,
  output_dir     = diag_dir,
  batch_var      = batch_var,
  id_var         = id_var,
  longitudinal   = longitudinal,
  discrete       = discrete,
  n_perm         = n_perm,
  min_batch_n    = min_batch_n,
  thresholds     = thresholds,
  user_overrides = user_overrides,
  n_cores        = n_cores
)

log_info(paste0("=== diagnose complete ===",
                " success: ", results$successes,
                " | failed: ", results$failures))
log_info(paste0("Recommendations: ", results$recommendations_path))
log_info(paste0("Pass to fit with: --feature_families ", results$recommendations_path))