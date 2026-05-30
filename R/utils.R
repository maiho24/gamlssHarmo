# R/utils.R

`%||%` <- function(a, b) if (!is.null(a)) a else b

format_time <- function(seconds) {
  if (seconds < 60)        paste0(round(seconds, 1), " seconds")
  else if (seconds < 3600) paste0(round(seconds / 60, 1), " minutes")
  else                     paste0(round(seconds / 3600, 2), " hours")
}

load_config <- function(config_path) {
  if (!file.exists(config_path))
    stop("Config file not found: ", config_path)
  yaml::read_yaml(config_path)
}

resolve_arg <- function(cli_val, config_val, default = NULL) {
  if (!is.null(cli_val) && !identical(cli_val, "")) return(cli_val)
  if (!is.null(config_val))                          return(config_val)
  default
}

make_group_palette <- function(n) {
  if (n <= 8) RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)]
  else        grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
}

read_features_file <- function(path) {
  if (!file.exists(path))
    stop("Features file not found: ", path)
  lines <- readLines(path)
  lines <- lines[!grepl("^\\s*#", lines)]
  lines[nzchar(trimws(lines))]
}

resolve_features <- function(feature_arg   = NULL,
                             features_file = NULL,
                             data          = NULL,
                             meta_cols     = NULL) {
  if (!is.null(feature_arg) && nzchar(feature_arg)) {
    message("Using single feature: ", feature_arg)
    return(feature_arg)
  }
  if (!is.null(features_file) && nzchar(features_file)) {
    feats <- read_features_file(features_file)
    message("Loaded ", length(feats), " features from: ", features_file)
    return(feats)
  }
  if (!is.null(data)) {
    feats <- setdiff(names(data), meta_cols)
    message("Using all ", length(feats), " non-metadata columns")
    return(feats)
  }
  stop("Cannot resolve feature list: provide --feature, --features, or --data")
}

substitute_placeholders <- function(term, batch_var, id_var = NULL) {
  term <- gsub("{batch}", batch_var, term, fixed = TRUE)
  if (!is.null(id_var))
    term <- gsub("{id}", id_var, term, fixed = TRUE)
  term
}

terms_to_formula <- function(terms, batch_var, id_var = NULL) {
  if (is.null(terms) || length(terms) == 0)
    return(as.formula("~ 1"))
  resolved <- vapply(terms, substitute_placeholders,
                     FUN.VALUE = character(1),
                     batch_var = batch_var, id_var = id_var)
  as.formula(paste("~", paste(resolved, collapse = " + ")))
}

# Smoother/spline wrappers recognised by gamlss and R spline packages.
# Shared by extract_covariate_cols() and validate_formula_terms().
.SMOOTHER_FNS <- c("pb", "random", "cs", "fp", "re", "cy", "lo", "tr", "nn",
                   "poly", "ns", "bs", "s", "te", "ti")

extract_covariate_cols <- function(terms, batch_var = NULL, id_var = NULL,
                                   data = NULL) {
  if (is.null(terms) || length(terms) == 0L) return(character(0L))

  skip_tokens <- c("1", batch_var, id_var)
  cols        <- character(0L)

  for (raw in unlist(terms)) {
    term  <- substitute_placeholders(raw, batch_var %||% "{batch}",
                                     id_var  %||% "{id}")
    inner <- sub(
      paste0("^(", paste(.SMOOTHER_FNS, collapse = "|"), ")\\((.+)\\)$"),
      "\\2", term
    )
    if (grepl("^[A-Za-z_\\.][A-Za-z0-9_\\.]*$", inner) &&
        !inner %in% skip_tokens)
      cols <- c(cols, inner)
  }

  if (!is.null(data)) cols <- intersect(cols, names(data))
  unique(cols)
}

validate_formula_terms <- function(formula_spec, data, param_name,
                                   batch_var, id_var) {
  problems <- character(0)

  for (raw_term in formula_spec) {
    term  <- substitute_placeholders(raw_term, batch_var, id_var)
    inner <- gsub(
      paste0("^(", paste(.SMOOTHER_FNS, collapse = "|"), ")\\((.+)\\)$"),
      "\\2", term
    )
    if (grepl("^[A-Za-z_\\.][A-Za-z0-9_\\.]*$", inner) &&
        !inner %in% names(data) && inner != "1")
      problems <- c(problems,
                    paste0("  [", param_name, "] '", raw_term,
                           "' -> column '", inner, "' not found in data"))
  }

  if (length(problems) > 0)
    stop("Formula validation failed:\n", paste(problems, collapse = "\n"))
  invisible(TRUE)
}

parse_and_validate_formulas <- function(config, data, longitudinal,
                                        batch_var = NULL, id_var = NULL) {
  model_cfg <- config$model
  batch_var <- batch_var %||% model_cfg$batch_var
  id_var    <- id_var    %||% model_cfg$id_var
  spec_key  <- if (isTRUE(longitudinal) && !is.null(model_cfg$formulas_longitudinal))
    "formulas_longitudinal" else "formulas"

  specs <- model_cfg[[spec_key]]
  if (is.null(specs))
    stop("No formula spec found under '", spec_key, "' in config")

  result <- list()
  for (p in c("mu", "sigma", "nu", "tau")) {
    terms <- specs[[p]]
    if (!is.null(terms) && length(terms) > 0)
      validate_formula_terms(terms, data, p, batch_var, id_var)
    result[[p]] <- terms
  }
  result
}

get_meta_cols <- function(batch_var, id_var, data = NULL) {
  base <- c(id_var, "age", "sex", batch_var)
  if (!is.null(data) && "wave" %in% names(data))
    base <- c(base, "wave")
  unique(base)
}

normalise_column_names <- function(data, age_var, sex_var) {
  rename_col <- function(df, from, to) {
    if (from == to) return(df)
    if (!from %in% names(df))
      stop("Column '", from, "' not found in data. ",
           "Check --age_var / --sex_var. Available columns: ",
           paste(names(df), collapse = ", "))
    names(df)[names(df) == from] <- to
    logger::log_info(paste0("Renamed column '", from, "' -> '", to, "'"))
    df
  }
  data <- rename_col(data, age_var, "age")
  data <- rename_col(data, sex_var, "sex")
  data
}

get_n_params <- function(family_name) {
  base_name <- strsplit(family_name, "_")[[1L]][1L]
  # 4-parameter families
  if (base_name %in% c("SHASH", "JSU", "BCT", "BCPE", "ST4", "EGB2",
                        "ZISICHEL"))
    return(4L)
  # 3-parameter families
  if (base_name %in% c("GG", "SN1", "TF", "PE2", "ST3", "BCCG",
                        "ZINBI", "ZANBI",
                        "SICHEL", "ZASICHEL", "ZIPIG", "exGAUS"))
    return(3L)
  # 1-parameter families
  if (base_name %in% c("PO", "EXP"))
    return(1L)
  # Default: 2-parameter (NBI, NBII, ZIP, PIG, DPO, LOGNO, GA, WEI, IG,
  #                       LO, NO, ...)
  2L
}

setup_logging <- function(log_dir, stage_name) {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_dir,
                        paste0(stage_name, "_",
                               format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  logger::log_formatter(logger::formatter_paste)
  logger::log_appender(logger::appender_tee(log_file))
  invisible(log_file)
}