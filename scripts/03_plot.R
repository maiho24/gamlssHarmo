#!/usr/bin/env Rscript
# gamlssHarmo plot -- Plot pre/post age trajectories
#

#
# At least one of --pre or --post must be supplied.
# Supplying both produces side-by-side pre/post comparison panels.
#
# Options:
#   --config        Path to params.yml  [config/params.yml]
#   --pre           Pre-harmonisation CSV
#   --post          Combined harmonised CSV (from infer stage)
#   --output        Output directory for plots
#   --features      Features .txt file
#   --feature       Single feature name (overrides --features)
#   --batch_var     Batch/site column name
#   --group_col     Column to colour trajectories by
#   --smooth_method loess or gam  [loess]
#   --age_bin_width Age bin width in years  [5]
#   --fix_y_limits  TRUE or FALSE  [TRUE]

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  library(logger)
})

get_root <- function() {
  if (exists("SCRIPT_DIR", envir = parent.env(environment()), inherits = TRUE))
    return(normalizePath(get("SCRIPT_DIR", envir = parent.env(environment()),
                             inherits = TRUE)))
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0)
    return(normalizePath(dirname(dirname(sub("^--file=", "", file_arg[1])))))
  normalizePath(".")
}

ROOT <- get_root()
source(file.path(ROOT, "R", "utils.R"))
source(file.path(ROOT, "R", "plot.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

option_list <- list(
  make_option("--config",        type = "character", default = NULL),
  make_option("--pre",           type = "character", default = NULL),
  make_option("--post",          type = "character", default = NULL),
  make_option("--output",        type = "character", default = NULL),
  make_option("--features",      type = "character", default = NULL),
  make_option("--feature",       type = "character", default = NULL),
  make_option("--batch_var",     type = "character", default = NULL),
  make_option("--group_col",     type = "character", default = NULL),
  make_option("--smooth_method", type = "character", default = NULL),
  make_option("--age_bin_width", type = "numeric",   default = NULL),
  make_option("--fix_y_limits",  type = "character", default = NULL)
)

opt <- parse_args(OptionParser(
  usage = "Rscript scripts/03_plot.R [options]",
  option_list = option_list,
  description = "Plot pre/post age trajectories."
))

cfg_path <- opt$config %||% file.path(ROOT, "config", "params.yml")
cfg      <- load_config(cfg_path)

pre_csv        <- resolve_arg(opt$pre,           cfg$data$raw_csv)
post_csv       <- resolve_arg(opt$post,          cfg$data$harmonised_csv)
output_dir     <- resolve_arg(opt$output,        cfg$output$plots,       "output/plots")
features_txt   <- resolve_arg(opt$features,      cfg$features$features_txt)
feature_single <- resolve_arg(opt$feature,       cfg$features$feature)
batch_var      <- resolve_arg(opt$batch_var,     cfg$model$batch_var,    "cohort")
id_var         <-                                cfg$model$id_var %||%   "id"
group_col      <- resolve_arg(opt$group_col,     cfg$plot$group_col,     "cohort")
smooth_method  <- resolve_arg(opt$smooth_method, cfg$plot$smooth_method, "loess")
age_bin_width  <- as.numeric(resolve_arg(opt$age_bin_width, cfg$plot$age_bin_width, 5))
fix_y_limits   <- as.logical(resolve_arg(opt$fix_y_limits,  cfg$plot$fix_y_limits,  TRUE))

log_dir <- cfg$output$logs %||% file.path(dirname(output_dir), "logs")
setup_logging(log_dir, "plot")

log_info("=== gamlssHarmo plot ===")
log_info(paste0("pre:           ", pre_csv  %||% "(not provided)"))
log_info(paste0("post:          ", post_csv %||% "(not provided)"))
log_info(paste0("group_col:     ", group_col))
log_info(paste0("smooth_method: ", smooth_method))
log_info(paste0("age_bin_width: ", age_bin_width))

has_pre  <- !is.null(pre_csv)  && nzchar(pre_csv)  && file.exists(pre_csv)
has_post <- !is.null(post_csv) && nzchar(post_csv) && file.exists(post_csv)

if (!has_pre && !has_post)
  stop("At least one of --pre or --post must point to an existing file.")

pre_data  <- if (has_pre)  read.csv(pre_csv,  stringsAsFactors = FALSE) else NULL
post_data <- if (has_post) read.csv(post_csv, stringsAsFactors = FALSE) else NULL

meta_cols <- get_meta_cols(batch_var, id_var, pre_data %||% post_data)

all_features <- if (!is.null(pre_data)) {
  setdiff(names(pre_data), meta_cols)
} else if (!is.null(post_data)) {
  harm_cols <- grep("\\.harmonised_value$", names(post_data), value = TRUE)
  sub("\\.harmonised_value$", "", harm_cols)
} else {
  character(0)
}

features <- if (!is.null(feature_single) && nzchar(feature_single)) {
  feature_single
} else if (!is.null(features_txt) && nzchar(features_txt) && file.exists(features_txt)) {
  read_features_file(features_txt)
} else {
  all_features
}

log_info(paste0("Features to plot: ", length(features)))

plot_trajectories(
  features      = features,
  pre_data      = pre_data,
  post_data     = post_data,
  output_dir    = output_dir,
  batch_var     = batch_var,
  id_var        = id_var,
  group_col     = group_col,
  smooth_method = smooth_method,
  age_bin_width = age_bin_width,
  fix_y_limits  = fix_y_limits
)

log_info(paste0("=== plot complete === saved to: ", output_dir))