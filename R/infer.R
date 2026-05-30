# R/infer.R
#
# Harmonisation and normative scoring functions for gamlssHarmo.
#
# CONTINUOUS mode (discrete = FALSE):
#   Harmonisation uses the exact probability integral transform (PIT):
#     map y -> F(y | fitted params) -> Q(p | harmonised params).
#   Normative score: qnorm(cdf_value).
#
# DISCRETE / COUNT mode (discrete = TRUE):
#   Raw integer counts; no scaling. Mid-p correction (Lancaster 1961):
#     mid_p(y) = P(Y < y) + 0.5 * P(Y = y)
#   Normative score: qnorm(mid_p_value).
#   Harmonised values remain discrete integer counts.
#
# discrete and log_transform are read from <feature>_meta.csv written by fit.R.

get_distribution_functions <- function(family_name) {
  base_name <- strsplit(family_name, "_")[[1L]][1L]
  dist_map  <- list(
    # Continuous ---------------------------------------------------------------
    SHASH  = list(p_func = gamlss.dist::pSHASH,  q_func = gamlss.dist::qSHASH,  d_func = NULL),
    JSU    = list(p_func = gamlss.dist::pJSU,    q_func = gamlss.dist::qJSU,    d_func = NULL),
    GG     = list(p_func = gamlss.dist::pGG,     q_func = gamlss.dist::qGG,     d_func = NULL),
    NO     = list(p_func = gamlss.dist::pNO,     q_func = gamlss.dist::qNO,     d_func = NULL),
    LOGNO  = list(p_func = gamlss.dist::pLOGNO,  q_func = gamlss.dist::qLOGNO,  d_func = NULL),
    GA     = list(p_func = gamlss.dist::pGA,     q_func = gamlss.dist::qGA,     d_func = NULL),
    TF     = list(p_func = gamlss.dist::pTF,     q_func = gamlss.dist::qTF,     d_func = NULL),
    PE2    = list(p_func = gamlss.dist::pPE2,    q_func = gamlss.dist::qPE2,    d_func = NULL),
    SN1    = list(p_func = gamlss.dist::pSN1,    q_func = gamlss.dist::qSN1,    d_func = NULL),
    ST3    = list(p_func = gamlss.dist::pST3,    q_func = gamlss.dist::qST3,    d_func = NULL),
    ST4    = list(p_func = gamlss.dist::pST4,    q_func = gamlss.dist::qST4,    d_func = NULL),
    BCT    = list(p_func = gamlss.dist::pBCT,    q_func = gamlss.dist::qBCT,    d_func = NULL),
    BCPE   = list(p_func = gamlss.dist::pBCPE,   q_func = gamlss.dist::qBCPE,   d_func = NULL),
    BCCG   = list(p_func = gamlss.dist::pBCCG,   q_func = gamlss.dist::qBCCG,   d_func = NULL),
    WEI    = list(p_func = gamlss.dist::pWEI,    q_func = gamlss.dist::qWEI,    d_func = NULL),
    IG     = list(p_func = gamlss.dist::pIG,     q_func = gamlss.dist::qIG,     d_func = NULL),
    LO     = list(p_func = gamlss.dist::pLO,     q_func = gamlss.dist::qLO,     d_func = NULL),
    EGB2   = list(p_func = gamlss.dist::pEGB2,   q_func = gamlss.dist::qEGB2,   d_func = NULL),
    EXP    = list(p_func = gamlss.dist::pEXP,    q_func = gamlss.dist::qEXP,    d_func = NULL),
    # qexGAUS uses numerical inversion (uniroot) -- slower than closed-form families.
    exGAUS = list(p_func = gamlss.dist::pexGAUS, q_func = gamlss.dist::qexGAUS, d_func = NULL),

    # Discrete -----------------------------------------------------------------
    PO     = list(p_func = gamlss.dist::pPO,     q_func = gamlss.dist::qPO,     d_func = gamlss.dist::dPO),
    NBI    = list(p_func = gamlss.dist::pNBI,    q_func = gamlss.dist::qNBI,    d_func = gamlss.dist::dNBI),
    NBII   = list(p_func = gamlss.dist::pNBII,   q_func = gamlss.dist::qNBII,   d_func = gamlss.dist::dNBII),
    ZIP    = list(p_func = gamlss.dist::pZIP,    q_func = gamlss.dist::qZIP,    d_func = gamlss.dist::dZIP),
    ZINBI  = list(p_func = gamlss.dist::pZINBI,  q_func = gamlss.dist::qZINBI,  d_func = gamlss.dist::dZINBI),
    ZANBI  = list(p_func = gamlss.dist::pZANBI,  q_func = gamlss.dist::qZANBI,  d_func = gamlss.dist::dZANBI),
    PIG    = list(p_func = gamlss.dist::pPIG,    q_func = gamlss.dist::qPIG,    d_func = gamlss.dist::dPIG),
    ZIPIG  = list(p_func = gamlss.dist::pZIPIG,  q_func = gamlss.dist::qZIPIG,  d_func = gamlss.dist::dZIPIG),
    SICHEL   = list(p_func = gamlss.dist::pSICHEL,   q_func = gamlss.dist::qSICHEL,   d_func = gamlss.dist::dSICHEL),
    ZISICHEL = list(p_func = gamlss.dist::pZISICHEL, q_func = gamlss.dist::qZISICHEL, d_func = gamlss.dist::dZISICHEL),
    ZASICHEL = list(p_func = gamlss.dist::pZASICHEL, q_func = gamlss.dist::qZASICHEL, d_func = gamlss.dist::dZASICHEL),
    DPO    = list(p_func = gamlss.dist::pDPO,    q_func = gamlss.dist::qDPO,    d_func = gamlss.dist::dDPO)
  )
  if (!base_name %in% names(dist_map))
    stop("Distribution '", base_name, "' not recognised. ",
         "Supported: ", paste(sort(names(dist_map)), collapse = ", "))
  dist_map[[base_name]]
}

# gamlss >=5.0-6 emits "New way of prediction in random()" via cat(); only
# capture.output() can suppress direct stdout writes.
predict_silent <- function(...) {
  val <- NULL
  capture.output(val <- suppressWarnings(predict(...)), type = "output")
  val
}

# mid_p(y) = P(Y < y) + 0.5 * P(Y = y) = F(y-1) + 0.5*f(y)
# Clamped to (eps, 1-eps) so qnorm() stays finite.
compute_midp_values <- function(y, params, dist_funcs) {
  if (is.null(dist_funcs$d_func))
    stop("d_func is NULL -- mid-p requires a PMF function (discrete families only)")

  y_int              <- as.integer(round(y))
  p_lower            <- do.call(dist_funcs$p_func, c(list(q = y_int - 1L), params))
  p_lower[y_int == 0L] <- 0
  p_mass             <- do.call(dist_funcs$d_func, c(list(x = y_int), params))
  mid_p              <- p_lower + 0.5 * p_mass
  pmin(pmax(mid_p, .Machine$double.eps), 1 - .Machine$double.eps)
}

apply_inverse_link <- function(model, param_name, linear_pred) {
  link_name <- model[[paste0(param_name, ".link")]]
  switch(link_name,
    "log"      = exp(linear_pred),
    "identity" = linear_pred,
    "logit"    = plogis(linear_pred),
    "inverse"  = 1 / linear_pred,
    "sqrt"     = linear_pred^2,
    "1/mu^2"   = 1 / sqrt(linear_pred),
    "probit"   = pnorm(linear_pred),
    {
      inv_fn <- paste0(link_name, ".linkinv")
      if (exists(inv_fn)) get(inv_fn)(linear_pred)
      else stop("Unknown link '", link_name, "' for ", param_name)
    }
  )
}

has_random_in_fitted <- function(model, param_name, batch_var) {
  coef_names <- names(coef(model, what = param_name))
  any(grepl(paste0("random\\(", batch_var, "\\)"), coef_names))
}

has_fixed_batch_in_fitted <- function(model, param_name, batch_var) {
  coef_names <- tryCatch(names(coef(model, what = param_name)),
                         error = function(e) character(0L))
  any(grepl(paste0("^", batch_var), coef_names))
}

compute_harmonised_param <- function(model, model_data, param_name, batch_var) {
  logger::log_info(paste0("  Harmonising: ", param_name))

  has_random <- has_random_in_fitted(model, param_name, batch_var)
  has_fixed  <- has_fixed_batch_in_fitted(model, param_name, batch_var)

  if (!has_random && !has_fixed) {
    logger::log_info(paste0("  No batch term in ", param_name, " -- using original"))
    return(predict_silent(model, what = param_name, type = "response",
                          newdata = model_data))
  }

  pred_link  <- predict_silent(model, what = param_name, type = "link",
                               newdata = model_data)
  pred_terms <- tryCatch(
    predict_silent(model, what = param_name, type = "terms",
                   newdata = model_data),
    error = function(e) NULL
  )

  if (!is.null(pred_terms) && is.matrix(pred_terms)) {
    col_idx <- if (has_random) {
      which(grepl(paste0("random(", batch_var, ")"), colnames(pred_terms), fixed = TRUE))
    } else {
      which(grepl(batch_var, colnames(pred_terms), fixed = TRUE))
    }
    if (length(col_idx) > 0L) {
      batch_effect   <- rowSums(pred_terms[, col_idx, drop = FALSE])
      pred_link_zero <- pred_link - batch_effect
      return(apply_inverse_link(model, param_name, pred_link_zero))
    }
    logger::log_warn(paste0("  Batch column not found in terms for ", param_name,
                            "; available: ", paste(colnames(pred_terms), collapse = ", "),
                            " -- batch effect NOT removed"))
  } else {
    logger::log_warn(paste0("  type='terms' prediction failed for ", param_name,
                            " -- batch effect NOT removed"))
  }

  apply_inverse_link(model, param_name, pred_link)
}

load_meta <- function(model_file, feature_name, log_transform_arg) {
  meta_file <- file.path(dirname(model_file),
                         paste0(feature_name, "_meta.csv"))
  if (!file.exists(meta_file)) {
    warning("Meta file not found for '", feature_name,
            "' -- assuming log_transform=", log_transform_arg, ", discrete=FALSE")
    return(list(log_transform = log_transform_arg, discrete = FALSE))
  }

  mc       <- read.csv(meta_file, stringsAsFactors = FALSE)
  saved_lt <- as.logical(mc$log_transform[1L])
  saved_d  <- if ("discrete" %in% names(mc)) as.logical(mc$discrete[1L]) else FALSE

  if (!saved_d && !identical(saved_lt, log_transform_arg))
    warning("log_transform mismatch for '", feature_name,
            "': fitted with ", saved_lt, " -- using fitted value.")

  list(log_transform = saved_lt, discrete = saved_d)
}

# Output columns in <feature>_harmonised.csv:
#   <id_var>, <batch_var>, age, sex, [wave]
#   original_value, harmonised_value, harmonisation_effect
#   cdf_value    F(y | fitted params) or mid_p(y) for discrete
#   cdf_valid    FALSE where cdf/mid-p was degenerate
#   distribution_family, discrete, log_transform
#   <p>_fitted, <p>_harmonised  (per active parameter)
#   normative_z_score

harmonise_gamlss_feature <- function(model_file, data, output_dir,
                                     batch_var                 = "batch",
                                     id_var                    = "id",
                                     log_transform             = FALSE,
                                     generate_normative_scores = TRUE) {
  feature_name <- gsub("_model\\.rds$", "", basename(model_file))
  start_time   <- Sys.time()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  out_csv    <- file.path(output_dir, paste0(feature_name, "_harmonised.csv"))
  out_timing <- file.path(output_dir, paste0(feature_name, "_timing.csv"))

  tryCatch({
    logger::log_info(paste0("Harmonising: ", feature_name))

    meta          <- load_meta(model_file, feature_name, log_transform)
    log_transform <- meta$log_transform
    discrete      <- meta$discrete

    model       <- readRDS(model_file)
    family_name <- model$family[1L]
    n_params    <- get_n_params(family_name)
    dist_funcs  <- get_distribution_functions(family_name)

    logger::log_info(paste0("  Family: ", family_name,
                            " | discrete: ", discrete,
                            " | log_transform: ", log_transform))

    required_cols <- c(feature_name, "age", "sex", id_var, batch_var)
    missing_cols  <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0L)
      stop("Missing columns: ", paste(missing_cols, collapse = ", "))

    structural         <- c(feature_name, "y", "age", "sex", id_var, batch_var)
    model_formula_vars <- tryCatch(
      unique(unlist(lapply(
        list(model$mu.formula, model$sigma.formula,
             model$nu.formula, model$tau.formula),
        function(f) if (!is.null(f)) all.vars(f) else character(0L)
      ))),
      error = function(e) character(0L)
    )
    extra_cols <- intersect(setdiff(model_formula_vars, structural), names(data))

    has_wave  <- "wave" %in% names(data)
    keep_cols <- unique(c(feature_name, "age", "sex", batch_var, id_var,
                          if (has_wave) "wave", extra_cols))

    if (length(extra_cols) > 0L)
      logger::log_info(paste0("  Extra covariates retained: ",
                              paste(extra_cols, collapse = ", ")))

    model_data <- data[, keep_cols, drop = FALSE]
    model_data <- model_data[
      Reduce("&", lapply(keep_cols, function(col) !is.na(model_data[[col]]))), ]
    model_data <- model_data[model_data[[feature_name]] >= 0, ]

    model_data[[batch_var]] <- factor(model_data[[batch_var]])
    model_data[[id_var]]    <- factor(model_data[[id_var]])

    if (nrow(model_data) == 0L) stop("No valid observations after filtering")

    original_y <- model_data[[feature_name]]

    if (discrete) {
      model_data$y <- as.integer(round(original_y))
    } else {
      model_data$y <- if (log_transform) log(original_y + 1) else original_y
    }

    param_names <- c("mu",
                     if (n_params >= 2L) "sigma",
                     if (n_params >= 3L) "nu",
                     if (n_params >= 4L) "tau")

    random_status <- setNames(
      lapply(param_names, function(p) has_random_in_fitted(model, p, batch_var)),
      param_names
    )
    fixed_status <- setNames(
      lapply(param_names, function(p) has_fixed_batch_in_fitted(model, p, batch_var)),
      param_names
    )
    if (!any(unlist(random_status)) && !any(unlist(fixed_status)))
      stop("No batch (", batch_var, ") terms found in any parameter -- ",
           "ensure the formula contains random({batch}) or {batch}")

    logger::log_info(paste0("  Batch effect mode: ",
                            if (any(unlist(random_status))) "random" else "fixed"))

    params_fitted <- lapply(setNames(param_names, param_names), function(p)
      predict_silent(model, what = p, type = "response", newdata = model_data))

    if (discrete) {
      pit_values <- compute_midp_values(model_data$y, params_fitted, dist_funcs)
      cdf_valid  <- is.finite(pit_values) & pit_values > 0 & pit_values < 1
    } else {
      pit_raw    <- do.call(dist_funcs$p_func, c(list(q = model_data$y), params_fitted))
      cdf_valid  <- is.finite(pit_raw) & pit_raw > 0 & pit_raw < 1
      pit_values <- pmin(pmax(pit_raw, .Machine$double.eps), 1 - .Machine$double.eps)
    }

    n_invalid <- sum(!cdf_valid)
    if (n_invalid > 0L)
      logger::log_info(paste0("  ", n_invalid, " PIT values outside (0,1) -> NA"))

    params_harmonised <- lapply(setNames(param_names, param_names), function(p)
      compute_harmonised_param(model, model_data, p, batch_var))

    harmonised_std             <- do.call(dist_funcs$q_func,
                                          c(list(p = pit_values), params_harmonised))
    harmonised_std[!cdf_valid] <- NA_real_

    if (discrete) {
      harmonised_orig <- harmonised_std
    } else {
      harmonised_orig <- if (log_transform) exp(harmonised_std) - 1 else harmonised_std
    }

    out_df <- data.frame(
      age                  = model_data$age,
      sex                  = model_data$sex,
      original_value       = original_y,
      harmonised_value     = harmonised_orig,
      harmonisation_effect = harmonised_orig - original_y,
      cdf_value            = pit_values,
      cdf_valid            = cdf_valid,
      distribution_family  = family_name,
      discrete             = discrete,
      log_transform        = log_transform
    )
    out_df[[id_var]]    <- model_data[[id_var]]
    out_df[[batch_var]] <- model_data[[batch_var]]

    if (has_wave) out_df$wave <- model_data$wave

    for (p in param_names) {
      out_df[[paste0(p, "_fitted")]]     <- params_fitted[[p]]
      out_df[[paste0(p, "_harmonised")]] <- params_harmonised[[p]]
    }

    if (generate_normative_scores) {
      z_scores             <- qnorm(pit_values)
      z_scores[!cdf_valid] <- NA_real_
      out_df$normative_z_score <- z_scores
    }

    write.csv(out_df, out_csv, row.names = FALSE)

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    logger::log_info(paste0("  Done in ", format_time(elapsed)))

    write.csv(data.frame(feature        = feature_name,
                         status         = "success",
                         distribution   = family_name,
                         discrete       = discrete,
                         start_time     = as.character(start_time),
                         end_time       = as.character(Sys.time()),
                         time_secs      = elapsed,
                         n_observations = nrow(out_df),
                         n_invalid_cdf  = n_invalid),
              out_timing, row.names = FALSE)

    list(status = "success", feature = feature_name, distribution = family_name,
         discrete = discrete, processing_time = elapsed,
         n_observations = nrow(out_df), output_file = out_csv)

  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    logger::log_error(paste0("Error harmonising ", feature_name, ": ", e$message))
    write.csv(data.frame(feature       = feature_name,
                         status        = "error",
                         error_message = e$message,
                         time_secs     = elapsed),
              out_timing, row.names = FALSE)
    list(status = "error", feature = feature_name,
         error = e$message, processing_time = elapsed)
  })
}

harmonise_all_gamlss_models <- function(model_base_dir, data,
                                        output_dir                = "output/harmonised",
                                        batch_var                 = "batch",
                                        id_var                    = "id",
                                        log_transform             = FALSE,
                                        generate_normative_scores = TRUE,
                                        feature_subset            = NULL,
                                        n_cores                   = 1L) {
  overall_start <- Sys.time()

  model_files <- list.files(model_base_dir, pattern = "_model\\.rds$",
                            recursive = TRUE, full.names = TRUE)
  if (length(model_files) == 0L)
    stop("No model files found in: ", model_base_dir)

  if (!is.null(feature_subset)) {
    feature_subset <- trimws(feature_subset)
    model_feat     <- gsub("_model\\.rds$", "", basename(model_files))
    model_files    <- model_files[model_feat %in% feature_subset]
    missing        <- setdiff(feature_subset, model_feat)
    if (length(missing) > 0L)
      logger::log_info(paste0("  Warning: no model found for: ",
                              paste(missing, collapse = ", ")))
  }

  if (length(model_files) == 0L)
    stop("No model files matched the requested features. ",
         "Check that --features names match the feature columns used during fitting, ",
         "and that --models points to the correct directory.")

  logger::log_info(paste0("Harmonising ", length(model_files), " features"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  process_one <- function(mf) {
    feat <- gsub("_model\\.rds$", "", basename(mf))
    harmonise_gamlss_feature(
      model_file                = mf,
      data                      = data,
      output_dir                = file.path(output_dir, paste0("feature_", feat)),
      batch_var                 = batch_var,
      id_var                    = id_var,
      log_transform             = log_transform,
      generate_normative_scores = generate_normative_scores
    )
  }

  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, {
      library(gamlss)
      library(gamlss.dist)
      library(dplyr)
      library(logger)
      logger::log_formatter(logger::formatter_paste)
    })
    parallel::clusterExport(cl,
      c("data", "output_dir", "batch_var", "id_var", "log_transform",
        "generate_normative_scores", "process_one"),
      envir = environment()
    )
    parallel::clusterExport(cl,
      c("harmonise_gamlss_feature", "load_meta",
        "get_distribution_functions", "get_n_params",
        "has_random_in_fitted", "has_fixed_batch_in_fitted",
        "predict_silent", "compute_harmonised_param",
        "apply_inverse_link", "compute_midp_values",
        "format_time", "%||%"),
      envir = globalenv()
    )
    results <- parallel::parLapply(cl, model_files, process_one)
    parallel::stopCluster(cl)
  } else {
    results <- vector("list", length(model_files))
    for (i in seq_along(model_files)) {
      logger::log_info(paste0("Feature ", i, "/", length(model_files)))
      results[[i]] <- process_one(model_files[i])
      elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
      logger::log_info(paste0("  ", round(100 * i / length(model_files), 1L),
                              "% | est. remaining: ",
                              format_time(elapsed / i * (length(model_files) - i))))
    }
  }

  statuses  <- vapply(results, `[[`, character(1L), "status")
  successes <- sum(statuses == "success")
  failures  <- sum(statuses == "error")

  list(results   = results,
       successes = successes,
       failures  = failures)
}

# ---------------------------------------------------------------------------
# Combine per-feature harmonised CSVs into wide combined output files.
# ---------------------------------------------------------------------------

combine_harmonised_results <- function(harmonised_output_dir,
                                       id_var                    = "id",
                                       batch_var                 = "batch",
                                       generate_normative_scores = TRUE,
                                       feature_subset            = NULL) {
  feat_dirs <- list.dirs(harmonised_output_dir, full.names = TRUE, recursive = FALSE)
  feat_dirs <- feat_dirs[grepl("^feature_", basename(feat_dirs))]

  if (!is.null(feature_subset)) {
    keep      <- paste0("feature_", feature_subset)
    feat_dirs <- feat_dirs[basename(feat_dirs) %in% keep]
  }

  if (length(feat_dirs) == 0L) {
    logger::log_info("combine_harmonised_results: no feature directories found")
    return(invisible(NULL))
  }

  meta_cols <- c(id_var, batch_var, "age", "sex", "wave")

  harm_wide <- NULL
  norm_wide <- NULL
  cdf_wide  <- NULL

  for (fd in feat_dirs) {
    feat_name <- sub("^feature_", "", basename(fd))
    csv_path  <- file.path(fd, paste0(feat_name, "_harmonised.csv"))

    if (!file.exists(csv_path)) {
      logger::log_info(paste0("  combine: missing ", csv_path, " -- skipping"))
      next
    }

    df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE),
                   error = function(e) {
                     logger::log_info(paste0("  combine: error reading ", csv_path,
                                             " -- ", e$message))
                     NULL
                   })
    if (is.null(df)) next

    key_cols <- intersect(c(id_var, batch_var, "age", "sex", "wave"), names(df))

    harm_col <- df[["harmonised_value"]]
    if (is.null(harm_col)) {
      logger::log_info(paste0("  combine: 'harmonised_value' missing in ", csv_path))
      next
    }

    slice_h <- df[, key_cols, drop = FALSE]
    slice_h[[feat_name]] <- harm_col

    if (is.null(harm_wide)) {
      harm_wide <- slice_h
      harm_wide$.row_index_h <- seq_len(nrow(harm_wide))
    } else {
      harm_wide <- merge(harm_wide, slice_h, by = key_cols, all = TRUE)
    }

    if (generate_normative_scores && "normative_z_score" %in% names(df)) {
      slice_n <- df[, key_cols, drop = FALSE]
      slice_n[[feat_name]] <- df[["normative_z_score"]]

      if (is.null(norm_wide)) {
        norm_wide <- slice_n
        norm_wide$.row_index_n <- seq_len(nrow(norm_wide))
      } else {
        norm_wide <- merge(norm_wide, slice_n, by = key_cols, all = TRUE)
      }
    }

    if ("cdf_value" %in% names(df)) {
      slice_c <- df[, key_cols, drop = FALSE]
      slice_c[[feat_name]] <- df[["cdf_value"]]

      if (is.null(cdf_wide)) {
        cdf_wide <- slice_c
        cdf_wide$.row_index_c <- seq_len(nrow(cdf_wide))
      } else {
        cdf_wide <- merge(cdf_wide, slice_c, by = key_cols, all = TRUE)
      }
    }
  }

  if (!is.null(harm_wide)) {
    if (".row_index_h" %in% names(harm_wide)) {
      harm_wide <- harm_wide[order(harm_wide$.row_index_h, na.last = TRUE), , drop = FALSE]
      harm_wide$.row_index_h <- NULL
    }
    out_harm <- file.path(harmonised_output_dir, "combined_harmonised.csv")
    write.csv(harm_wide, out_harm, row.names = FALSE)
    logger::log_info(paste0("Combined harmonised: ", out_harm))
  }

  if (!is.null(norm_wide)) {
    if (".row_index_n" %in% names(norm_wide)) {
      norm_wide <- norm_wide[order(norm_wide$.row_index_n, na.last = TRUE), , drop = FALSE]
      norm_wide$.row_index_n <- NULL
    }
    out_norm <- file.path(harmonised_output_dir, "combined_normative.csv")
    write.csv(norm_wide, out_norm, row.names = FALSE)
    logger::log_info(paste0("Combined normative:  ", out_norm))
  }

  if (!is.null(cdf_wide)) {
    if (".row_index_c" %in% names(cdf_wide)) {
      cdf_wide <- cdf_wide[order(cdf_wide$.row_index_c, na.last = TRUE), , drop = FALSE]
      cdf_wide$.row_index_c <- NULL
    }
    out_cdf <- file.path(harmonised_output_dir, "combined_cdf.csv")
    write.csv(cdf_wide, out_cdf, row.names = FALSE)
    logger::log_info(paste0("Combined CDF:        ", out_cdf))
  }

  invisible(list(harmonised = harm_wide, normative = norm_wide, cdf = cdf_wide))
}