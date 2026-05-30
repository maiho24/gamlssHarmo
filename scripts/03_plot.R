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
    help    = "Path to pre-harmonisation CSV (raw input data). At least one of --pre or --post must be provided."),

  make_option("--post",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to combined_harmonised.csv produced by the infer stage. Supplying both --pre and --post produces side-by-side comparison panels."),

  make_option("--output",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Base output directory. Plots are saved to <o>/plots<suffix>/. [default: output/]"),

  make_option("--suffix",
    type    = "character",
    default = NULL,
    metavar = "STRING",
    help    = "String appended to plots/ and logs/ subdirectory names. E.g. --suffix _cont writes to plots_cont/. [default: (empty)]"),

  make_option("--features_file",
    type    = "character",
    default = NULL,
    metavar = "PATH",
    help    = "Path to a .txt file listing feature names to plot, one per line. Lines starting with '#' are ignored. Overridden by --one_feature."),

  make_option("--one_feature",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Single feature name to plot. Overrides --features_file and config."),

  make_option("--batch_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name identifying the batch/site variable in pre-harmonisation data. [default: batch]"),

  make_option("--id_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name for the subject identifier. [default: id]"),

  make_option("--age_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name for age (x-axis of trajectory plots). Renamed to 'age' internally. [default: age]"),

  make_option("--sex_var",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name for sex (faceting variable in trajectory plots). Renamed to 'sex' internally. [default: sex]"),

  make_option("--group_col",
    type    = "character",
    default = NULL,
    metavar = "NAME",
    help    = "Column name to colour and group age trajectories by (e.g. batch, sex, diagnosis). [default: batch]"),

  make_option("--smooth_method",
    type    = "character",
    default = NULL,
    metavar = "METHOD",
    help    = "Smoothing method for trajectory lines. Options: loess, gam. [default: loess]"),

  make_option("--age_bin_width",
    type    = "numeric",
    default = NULL,
    metavar = "N",
    help    = "Width of age bins in years for trajectory summaries. [default: 5]"),

  make_option("--fix_y_limits",
    type    = "character",
    default = NULL,
    metavar = "TRUE/FALSE",
    help    = "If TRUE, locks y-axis to 1st-99th percentile across all panels for comparability. [default: TRUE]")
)

parser <- OptionParser(
  usage       = "gamlssHarmo plot [options]\n       Rscript scripts/03_plot.R [options]",
  option_list = option_list,
  description = paste(
    "Plot age trajectories for features before and/or after harmonisation.",
    "",
    "At least one of --pre or --post must be supplied.",
    "Supplying both produces side-by-side pre/post comparison panels,",
    "with Pre-Harmonisation always on the left and Post-Harmonisation on the right.",
    "",
    "Output (written to <o>/plots<suffix>/):",
    "  <feature>_trajectory.png   One PNG per feature (16x10 inches, 300 dpi)",
    "",
    "Examples:",
    "  gamlssHarmo plot --pre data/raw/my_data.csv --output output/",
    "  gamlssHarmo plot --post output/harmonised/combined_harmonised.csv --output output/",
    "  gamlssHarmo plot --pre data/raw/my_data.csv \\",
    "    --post output/harmonised/combined_harmonised.csv --output output/",
    sep = "\n"
  ),
  epilogue = "All options can also be set in config/params.yml. CLI arguments always take priority."
)
opt <- parse_args(parser)

cfg_path <- opt$config %||% file.path(ROOT, "config", "params.yml")
cfg      <- load_config(cfg_path)

suffix         <- resolve_arg(opt$suffix, cfg$output$suffix, "")
base_dir       <- normalizePath(resolve_arg(opt$output, cfg$output$base, "output"), mustWork = FALSE)

pre_csv        <- resolve_arg(opt$pre,  cfg$data$raw_csv)
post_csv       <- opt$post

plots_dir      <- file.path(base_dir, paste0("plots", suffix))
log_dir        <- file.path(base_dir, paste0("logs",  suffix))
features_txt   <- resolve_arg(opt$features_file,  cfg$features$features_txt)
feature_single <- resolve_arg(opt$one_feature,    cfg$features$feature)
batch_var      <- resolve_arg(opt$batch_var,      cfg$model$batch_var,    "batch")
id_var         <- resolve_arg(opt$id_var,         cfg$model$id_var,       "id")
age_var        <- resolve_arg(opt$age_var,        cfg$model$age_var,      "age")
sex_var        <- resolve_arg(opt$sex_var,        cfg$model$sex_var,      "sex")
group_col      <- resolve_arg(opt$group_col,      cfg$plot$group_col,     "batch")
smooth_method  <- resolve_arg(opt$smooth_method,  cfg$plot$smooth_method, "loess")
age_bin_width  <- as.numeric(resolve_arg(opt$age_bin_width, cfg$plot$age_bin_width, 5))
fix_y_limits   <- as.logical(resolve_arg(opt$fix_y_limits,  cfg$plot$fix_y_limits,  TRUE))

setup_logging(log_dir, "plot")

log_info("=== gamlssHarmo plot ===")
log_info(paste0("pre:           ", pre_csv  %||% "(not provided)"))
log_info(paste0("post:          ", post_csv %||% "(not provided)"))
log_info(paste0("output base:   ", base_dir))
log_info(paste0("suffix:        ", if (nzchar(suffix)) suffix else "(none)"))
log_info(paste0("plots dir:     ", plots_dir))
log_info(paste0("batch_var:     ", batch_var))
log_info(paste0("id_var:        ", id_var))
log_info(paste0("age_var:       ", age_var))
log_info(paste0("sex_var:       ", sex_var))
log_info(paste0("group_col:     ", group_col))
log_info(paste0("smooth_method: ", smooth_method))
log_info(paste0("age_bin_width: ", age_bin_width))

has_pre  <- !is.null(pre_csv)  && nzchar(pre_csv)  && file.exists(pre_csv)
has_post <- !is.null(post_csv) && nzchar(post_csv) && file.exists(post_csv)

if (!has_pre && !has_post) {
  print_help(parser)
  quit(status = 1)
}

pre_data  <- if (has_pre)  read.csv(pre_csv,  stringsAsFactors = FALSE) else NULL
post_data <- if (has_post) read.csv(post_csv, stringsAsFactors = FALSE) else NULL
if (!is.null(pre_data))
  pre_data <- normalise_column_names(pre_data, age_var, sex_var)

ref_data  <- if (!is.null(pre_data)) pre_data else post_data
meta_cols <- get_meta_cols(batch_var, id_var, ref_data)

features <- if (!is.null(feature_single) && nzchar(trimws(feature_single))) {
  trimws(feature_single)
} else if (!is.null(features_txt) && nzchar(trimws(features_txt))) {
  if (!file.exists(features_txt))
    stop("Features file not found: ", features_txt)
  read_features_file(features_txt)
} else if (!is.null(cfg$features$features_txt) &&
           nzchar(cfg$features$features_txt %||% "") &&
           file.exists(cfg$features$features_txt)) {
  read_features_file(cfg$features$features_txt)
} else {
  setdiff(names(ref_data), meta_cols)
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

log_info("=== plot complete ===")