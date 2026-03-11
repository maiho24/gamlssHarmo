#!/usr/bin/env Rscript
# gamlssHarmo infer -- Apply harmonisation and compute normative scores
#
# Options:
#   --config        Path to params.yml  [config/params.yml]
#   --data          Input CSV
#   --models        Fitted models directory (from fit stage)
#   --output        Output directory for harmonised CSVs (overrides config)
#   --features      Features .txt file
#   --feature       Single feature name (overrides --features)
#   --batch_var     Batch/site column name (must match fit stage)
#   --log_transform TRUE or FALSE
#   --normative     TRUE or FALSE  [TRUE]
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
source(file.path(ROOT, "R", "infer.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

option_list <- list(
  make_option("--config",        type = "character", default = NULL),
  make_option("--data",          type = "character", default = NULL),
  make_option("--models",        type = "character", default = NULL),
  make_option("--output",        type = "character", default = NULL),
  make_option("--features",      type = "character", default = NULL),
  make_option("--feature",       type = "character", default = NULL),
  make_option("--batch_var",     type = "character", default = NULL),
  make_option("--log_transform", type = "character", default = NULL),
  make_option("--normative",     type = "character", default = "TRUE"),
  make_option("--n_cores",       type = "integer",   default = NULL)
)

opt <- parse_args(OptionParser(
  usage = "Rscript scripts/02_infer.R [options]",
  option_list = option_list,
  description = "Apply GAMLSS harmonisation and compute normative scores."
))

cfg_path <- opt$config %||% file.path(ROOT, "config", "params.yml")
cfg      <- load_config(cfg_path)

raw_csv          <- resolve_arg(opt$data,      cfg$data$raw_csv)
models_dir       <- resolve_arg(opt$models,    cfg$output$models %||% "output/models")
output_dir       <- if (!is.null(opt$output) && nzchar(opt$output)) {
                      file.path(opt$output, "harmonised")
                    } else {
                      cfg$output$harmonised %||% "output/harmonised"
                    }
features_txt     <- resolve_arg(opt$features,  cfg$features$features_txt)
feature_single   <- resolve_arg(opt$feature,   cfg$features$feature)
batch_var        <- resolve_arg(opt$batch_var, cfg$model$batch_var,  "cohort")
id_var           <-                            cfg$model$id_var %||% "id"
log_transform    <- as.logical(resolve_arg(opt$log_transform, cfg$model$log_transform, FALSE))
normative_scores <- as.logical(opt$normative %||% "TRUE")
n_cores          <- as.integer(resolve_arg(opt$n_cores, cfg$compute$n_cores, 1L))

log_dir <- cfg$output$logs %||% file.path(dirname(output_dir), "logs")
setup_logging(log_dir, "infer")

log_info("=== gamlssHarmo infer ===")
log_info(paste0("config:          ", cfg_path))
log_info(paste0("data:            ", raw_csv %||% "(not set)"))
log_info(paste0("models dir:      ", models_dir))
log_info(paste0("harmonised dir:  ", output_dir))
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

feature_subset <- if (!is.null(feature_single) && nzchar(feature_single)) {
  feature_single
} else if (!is.null(features_txt) && nzchar(features_txt) && file.exists(features_txt)) {
  read_features_file(features_txt)
} else {
  NULL
}

results <- harmonise_all_gamlss_models(
  model_base_dir            = models_dir,
  data                      = data,
  output_dir                = output_dir,
  batch_var                 = batch_var,
  id_var                    = id_var,
  log_transform             = log_transform,
  generate_normative_scores = normative_scores,
  feature_subset            = feature_subset,
  n_cores                   = n_cores
)

combine_harmonised_results(
  harmonised_output_dir     = output_dir,
  id_var                    = id_var,
  generate_normative_scores = normative_scores
)

log_info(paste0("=== infer complete === success: ", results$successes,
                " | failed: ", results$failures))
log_info(paste0("combined CSV: ", file.path(output_dir, "combined_harmonised.csv")))