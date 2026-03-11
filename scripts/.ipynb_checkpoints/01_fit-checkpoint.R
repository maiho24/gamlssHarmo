#!/usr/bin/env Rscript
# gamlssHarmo fit -- Fit hierarchical GAMLSS harmonisation models
#
# Options:
#   --config        Path to params.yml  [config/params.yml]
#   --data          Raw input CSV
#   --output        Base output directory (overrides config)
#   --features      Features .txt file
#   --feature       Single feature name (overrides --features)
#   --batch_var     Batch/site column name
#   --longitudinal  TRUE or FALSE
#   --log_transform TRUE or FALSE
#   --family_order  Comma-separated e.g. "SHASH,GG,NO"
#   --n_cores       Parallel cores

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
source(file.path(ROOT, "R", "fit.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

option_list <- list(
  make_option("--config",        type = "character", default = NULL),
  make_option("--data",          type = "character", default = NULL),
  make_option("--output",        type = "character", default = NULL),
  make_option("--features",      type = "character", default = NULL),
  make_option("--feature",       type = "character", default = NULL),
  make_option("--batch_var",     type = "character", default = NULL),
  make_option("--longitudinal",  type = "character", default = NULL),
  make_option("--log_transform", type = "character", default = NULL),
  make_option("--family_order",  type = "character", default = NULL),
  make_option("--n_cores",       type = "integer",   default = NULL)
)

opt <- parse_args(OptionParser(
  usage = "Rscript scripts/01_fit.R [options]",
  option_list = option_list,
  description = "Fit hierarchical GAMLSS harmonisation models."
))

cfg_path <- opt$config %||% file.path(ROOT, "config", "params.yml")
cfg      <- load_config(cfg_path)

raw_csv        <- resolve_arg(opt$data,          cfg$data$raw_csv)
# If --output is supplied on CLI it overrides everything; otherwise use
# cfg$output$models directly so the user controls the exact path in config.
models_dir     <- if (!is.null(opt$output) && nzchar(opt$output)) {
                    file.path(opt$output, "models")
                  } else {
                    cfg$output$models %||% "output/models"
                  }
features_txt   <- resolve_arg(opt$features,      cfg$features$features_txt)
feature_single <- resolve_arg(opt$feature,       cfg$features$feature)
batch_var      <- resolve_arg(opt$batch_var,     cfg$model$batch_var,  "cohort")
id_var         <-                                cfg$model$id_var %||% "id"
longitudinal   <- as.logical(resolve_arg(opt$longitudinal,  cfg$model$longitudinal,  FALSE))
log_transform  <- as.logical(resolve_arg(opt$log_transform, cfg$model$log_transform, FALSE))
n_cores        <- as.integer(resolve_arg(opt$n_cores,       cfg$compute$n_cores,     1L))
family_order   <- if (!is.null(opt$family_order))
  trimws(strsplit(opt$family_order, ",")[[1]]) else
  cfg$model$family_order %||% c("SHASH", "GG", "NO")

log_dir <- cfg$output$logs %||% file.path(dirname(models_dir), "logs")
setup_logging(log_dir, "fit")

log_info("=== gamlssHarmo fit ===")
log_info(paste0("config:        ", cfg_path))
log_info(paste0("data:          ", raw_csv %||% "(not set)"))
log_info(paste0("models dir:    ", models_dir))
log_info(paste0("batch_var:     ", batch_var))
log_info(paste0("longitudinal:  ", longitudinal))
log_info(paste0("log_transform: ", log_transform))
log_info(paste0("family_order:  ", paste(family_order, collapse = " -> ")))
log_info(paste0("n_cores:       ", n_cores))

if (is.null(raw_csv) || !nzchar(raw_csv))
  stop("Provide --data or set data.raw_csv in params.yml")
if (!file.exists(raw_csv))
  stop("Data file not found: ", raw_csv)

data <- read.csv(raw_csv, stringsAsFactors = FALSE)
log_info(paste0("Loaded: ", nrow(data), " rows x ", ncol(data), " cols"))

meta_cols <- get_meta_cols(batch_var, id_var, data)
features  <- resolve_features(
  feature_arg   = feature_single,
  features_file = features_txt,
  data          = data,
  meta_cols     = meta_cols
)
log_info(paste0("Features: ", length(features)))

formula_terms <- parse_and_validate_formulas(cfg, data, longitudinal)
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
  family_order  = family_order,
  n_cores       = n_cores
)

log_info(paste0("=== fit complete === success: ", results$successes,
                " | skipped: ", results$skipped,
                " | failed: ", results$failures))
log_info(paste0("time: ", results$total_time))