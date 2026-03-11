#!/usr/bin/env Rscript
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
source(file.path(ROOT, "R", "plot.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

option_list <- list(
  make_option("--config",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to params.yml config file. [default: config/params.yml]"),

  make_option("--pre",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to pre-harmonisation CSV (raw input data). When supplied, a Pre-Harmonisation panel is shown. At least one of --pre or --post must be provided."),

  make_option("--post",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to combined_harmonised.csv produced by the infer stage. When supplied, a Post-Harmonisation panel is shown. Supplying both --pre and --post produces side-by-side comparison panels."),

  make_option("--output",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Base output directory. Plots are saved to <o>/plots/. [default: output/]"),

  make_option("--features-file",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to a .txt file listing feature names to plot, one per line. Lines starting with '#' are ignored. Overridden by --feature."),

  make_option("--one-feature",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Single feature name to plot. Overrides --features and config."),

  make_option("--batch_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name identifying the batch/site variable. Used to locate the batch column in pre-harmonisation data. [default: cohort]"),

  make_option("--group_col",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name to colour and group age trajectories by (e.g. cohort, sex, diagnosis). [default: cohort]"),

  make_option("--smooth_method",
    type    = "character",
    default = NULL,
    metavar = "METHOD",
    help    = "Smoothing method for trajectory lines. Options: loess, gam. [default: loess]"),

  make_option("--age_bin_width",
    type    = "numeric",
    default = NULL,
    metavar = "N",
    help    = "Width of age bins in years used to compute mean and SE for trajectory points. Smaller values give finer resolution but noisier estimates. [default: 5]"),

  make_option("--fix_y_limits",
    type    = "character",
    default = NULL,
    metavar = "TRUE/FALSE",
    help    = "If TRUE, locks the y-axis range to the 1st-99th percentile of values across all panels, making pre/post panels directly comparable. Set to FALSE to allow free y-axes. [default: TRUE]")
)

opt <- parse_args(OptionParser(
  usage       = "gamlssHarmo plot [options]\n       Rscript scripts/03_plot.R [options]",
  option_list = option_list,
  description = paste(
    "Plot age trajectories for features before and/or after harmonisation.",
    "",
    "Data are binned by age (--age_bin_width), then mean +/- 1 SE is plotted",
    "per group per bin, with a smooth curve overlaid. Facets are split by sex",
    "and by pre/post harmonisation panel.",
    "",
    "At least one of --pre or --post must be supplied:",
    "  --pre only   Shows raw data trajectories only",
    "  --post only  Shows harmonised data trajectories only",
    "  Both         Shows side-by-side pre/post comparison panels",
    "",
    "Output (written to <o>/plots/):",
    "  <feature>_trajectory.png   One PNG per feature (16x10 inches, 300 dpi)",
    "",
    "Examples:",
    "  gamlssHarmo plot --pre data/raw/my_data.csv --output output/",
    "  gamlssHarmo plot --post output/harmonised/combined_harmonised.csv --output output/",
    "  gamlssHarmo plot --pre data/raw/my_data.csv --post output/harmonised/combined_harmonised.csv --output output/",
    "  gamlssHarmo plot --pre data/raw/my_data.csv --output output/ --feature ThicknessAvg --group_col sex",
    sep = "\n"
  ),
  epilogue = "All options can also be set in config/params.yml. CLI arguments always take priority."
))

cfg_path <- opt$config %||% file.path(ROOT, "config", "params.yml")
cfg      <- load_config(cfg_path)

pre_csv        <- resolve_arg(opt$pre,  cfg$data$raw_csv)
post_csv       <- resolve_arg(opt$post, cfg$data$harmonised_csv)
base_dir       <- normalizePath(resolve_arg(opt$output, cfg$output$base, "output"), mustWork = FALSE)
plots_dir      <- file.path(base_dir, "plots")
log_dir        <- file.path(base_dir, "logs")
features_txt   <- resolve_arg(opt$features_file, cfg$features$features_txt)
feature_single <- resolve_arg(opt$one_feature,   cfg$features$feature)
batch_var      <- resolve_arg(opt$batch_var,     cfg$model$batch_var,    "cohort")
id_var         <-                                cfg$model$id_var %||%   "id"
group_col      <- resolve_arg(opt$group_col,     cfg$plot$group_col,     "cohort")
smooth_method  <- resolve_arg(opt$smooth_method, cfg$plot$smooth_method, "loess")
age_bin_width  <- as.numeric(resolve_arg(opt$age_bin_width, cfg$plot$age_bin_width, 5))
fix_y_limits   <- as.logical(resolve_arg(opt$fix_y_limits,  cfg$plot$fix_y_limits,  TRUE))

setup_logging(log_dir, "plot")

log_info("=== gamlssHarmo plot ===")
log_info(paste0("pre:           ", pre_csv  %||% "(not provided)"))
log_info(paste0("post:          ", post_csv %||% "(not provided)"))
log_info(paste0("output base:   ", base_dir))
log_info(paste0("plots dir:     ", plots_dir))
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

post_meta  <- c("id", "batch", "age", "sex", "wave", id_var, batch_var)
pre_feats  <- if (!is.null(pre_data))  setdiff(names(pre_data),  meta_cols)  else NULL
post_feats <- if (!is.null(post_data)) setdiff(names(post_data), post_meta)  else NULL

all_features <- if (!is.null(pre_feats) && !is.null(post_feats)) {
  only_pre  <- setdiff(pre_feats,  post_feats)
  only_post <- setdiff(post_feats, pre_feats)
  if (length(only_pre) > 0)
    log_info(paste0("Skipping (not harmonised, pre only): ",
                    paste(only_pre, collapse = ", ")))
  if (length(only_post) > 0)
    log_info(paste0("Skipping (post only, not in raw data): ",
                    paste(only_post, collapse = ", ")))
  intersect(pre_feats, post_feats)
} else {
  c(pre_feats, post_feats)
}

# Resolve features: CLI --feature takes priority over --features file,
# which takes priority over config, which falls back to all data columns.
features <- if (!is.null(opt$one_feature) && nzchar(opt$one_feature)) {
  opt$one_feature
} else if (!is.null(opt$features_file) && nzchar(opt$features_file)) {
  if (!file.exists(opt$features_file))
    stop("Features file not found: ", opt$features_file)
  read_features_file(opt$features_file)
} else if (!is.null(cfg$features$features_txt) && nzchar(cfg$features$features_txt) &&
           file.exists(cfg$features$features_txt)) {
  read_features_file(cfg$features$features_txt)
} else {
  all_features
}

log_info(paste0("Features to plot: ", length(features)))

plot_trajectories(
  features      = features,
  pre_data      = pre_data,
  post_data     = post_data,
  output_dir    = plots_dir,
  batch_var     = batch_var,
  id_var        = id_var,
  group_col     = group_col,
  smooth_method = smooth_method,
  age_bin_width = age_bin_width,
  fix_y_limits  = fix_y_limits
)

log_info(paste0("=== plot complete === saved to: ", plots_dir))