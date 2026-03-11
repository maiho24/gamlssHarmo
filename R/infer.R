get_distribution_functions <- function(family_name) {
  base_name <- strsplit(family_name, "_")[[1]][1]
  dist_map  <- list(
    SHASH = list(p_func = gamlss.dist::pSHASH, q_func = gamlss.dist::qSHASH),
    GG    = list(p_func = gamlss.dist::pGG,    q_func = gamlss.dist::qGG),
    NO    = list(p_func = gamlss.dist::pNO,    q_func = gamlss.dist::qNO),
    LOGNO = list(p_func = gamlss.dist::pLOGNO, q_func = gamlss.dist::qLOGNO),
    GA    = list(p_func = gamlss.dist::pGA,    q_func = gamlss.dist::qGA),
    TF    = list(p_func = gamlss.dist::pTF,    q_func = gamlss.dist::qTF),
    BCT   = list(p_func = gamlss.dist::pBCT,   q_func = gamlss.dist::qBCT),
    BCPE  = list(p_func = gamlss.dist::pBCPE,  q_func = gamlss.dist::qBCPE),
    BCCG  = list(p_func = gamlss.dist::pBCCG,  q_func = gamlss.dist::qBCCG),
    SN1   = list(p_func = gamlss.dist::pSN1,   q_func = gamlss.dist::qSN1),
    ST3   = list(p_func = gamlss.dist::pST3,   q_func = gamlss.dist::qST3),
    ST4   = list(p_func = gamlss.dist::pST4,   q_func = gamlss.dist::qST4),
    WEI   = list(p_func = gamlss.dist::pWEI,   q_func = gamlss.dist::qWEI),
    IG    = list(p_func = gamlss.dist::pIG,    q_func = gamlss.dist::qIG),
    LO    = list(p_func = gamlss.dist::pLO,    q_func = gamlss.dist::qLO),
    EGB2  = list(p_func = gamlss.dist::pEGB2,  q_func = gamlss.dist::qEGB2),
    EXP   = list(p_func = gamlss.dist::pEXP,   q_func = gamlss.dist::qEXP)
  )
  if (!base_name %in% names(dist_map))
    stop("Distribution '", base_name, "' not recognised")
  dist_map[[base_name]]
}

get_n_params <- function(family_name) {
  base_name <- strsplit(family_name, "_")[[1]][1]
  if (base_name %in% c("SHASH", "BCT", "BCPE", "ST4", "EGB2")) return(4L)
  if (base_name %in% c("GG", "SN1", "TF", "ST3", "BCCG"))      return(3L)
  2L
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

compute_harmonised_param <- function(model, model_data, param_name, batch_var) {
  logger::log_info(paste0("  Harmonising: ", param_name))

  if (!has_random_in_fitted(model, param_name, batch_var)) {
    logger::log_info(paste0("  No random(", batch_var, ") in ", param_name,
                            " -- using original"))
    return(predict(model, what = param_name, type = "response", newdata = model_data))
  }

  pred_link  <- predict(model, what = param_name, type = "link", newdata = model_data)
  pred_terms <- tryCatch(
    predict(model, what = param_name, type = "terms", newdata = model_data),
    warning = function(w) {
      if (grepl("safe.*prediction", w$message, ignore.case = TRUE)) NULL else stop(w)
    }
  )

  if (!is.null(pred_terms) && is.matrix(pred_terms)) {
    col_idx <- which(grepl(paste0("random.*", batch_var),
                           colnames(pred_terms), ignore.case = TRUE))
    if (length(col_idx) == 1 && nrow(pred_terms) == length(pred_link)) {
      logger::log_info(paste0("  Removed random(", batch_var, ") via terms matrix"))
      return(apply_inverse_link(model, param_name, pred_link - pred_terms[, col_idx]))
    }
  }

  logger::log_info(paste0("  Falling back to manual coefficient subtraction"))
  all_coefs  <- coef(model, what = param_name)
  rand_names <- names(all_coefs)[grepl(paste0("random\\(", batch_var, "\\)"),
                                       names(all_coefs))]
  harmonised_link <- pred_link
  for (batch_level in levels(model_data[[batch_var]])) {
    match_idx <- which(grepl(paste0("random.*", batch_level, "\\)"), rand_names))
    if (length(match_idx) > 0) {
      effect  <- all_coefs[rand_names[match_idx[1]]]
      obs_idx <- which(model_data[[batch_var]] == batch_level)
      harmonised_link[obs_idx] <- harmonised_link[obs_idx] - effect
      logger::log_info(paste0("  '", batch_level, "' effect: ", round(effect, 4)))
    }
  }
  apply_inverse_link(model, param_name, harmonised_link)
}

load_scaling <- function(model_file, feature_name, log_transform_arg) {
  scaling_file <- file.path(dirname(model_file),
                            paste0(feature_name, "_scaling.csv"))
  if (!file.exists(scaling_file))
    stop("Scaling file not found: ", scaling_file, "\nRe-run fitting.")

  sc       <- read.csv(scaling_file, stringsAsFactors = FALSE)
  saved_lt <- as.logical(sc$log_transform[1])
  if (!identical(saved_lt, log_transform_arg))
    warning("log_transform mismatch for '", feature_name,
            "': fitted with ", saved_lt, " -- using fitted value.")

  list(y_mean = sc$y_mean[1], y_sd = sc$y_sd[1], log_transform = saved_lt)
}

harmonise_gamlss_feature <- function(model_file, data, output_dir,
                                     batch_var                 = "cohort",
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

    sc            <- load_scaling(model_file, feature_name, log_transform)
    y_mean        <- sc$y_mean
    y_sd          <- sc$y_sd
    log_transform <- sc$log_transform

    model       <- readRDS(model_file)
    family_name <- model$family[1]
    n_params    <- get_n_params(family_name)
    dist_funcs  <- get_distribution_functions(family_name)

    logger::log_info(paste0("  Family: ", family_name,
                            " | log_transform: ", log_transform))

    required_cols <- c(feature_name, "age", "sex", id_var, batch_var)
    missing_cols  <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0)
      stop("Missing columns: ", paste(missing_cols, collapse = ", "))

    has_wave  <- "wave" %in% names(data)
    keep_cols <- unique(c(feature_name, "age", "sex", batch_var, id_var,
                          if (has_wave) "wave"))

    model_data <- data[, keep_cols, drop = FALSE]
    model_data <- model_data[
      Reduce("&", lapply(keep_cols, function(col) !is.na(model_data[[col]]))), ]
    model_data <- model_data[model_data[[feature_name]] >= 0, ]

    model_data[[batch_var]] <- factor(model_data[[batch_var]])
    model_data[[id_var]]    <- factor(model_data[[id_var]])

    if (nrow(model_data) == 0) stop("No valid observations after filtering")

    original_y <- model_data[[feature_name]]
    y_scaled   <- if (log_transform) log(original_y + 1) else original_y
    model_data$y <- (y_scaled - y_mean) / y_sd

    param_names <- c("mu",
                     if (n_params >= 2) "sigma",
                     if (n_params >= 3) "nu",
                     if (n_params >= 4) "tau")

    random_status <- setNames(
      lapply(param_names, function(p) has_random_in_fitted(model, p, batch_var)),
      param_names
    )

    if (!any(unlist(random_status)))
      stop("No random(", batch_var, ") terms found in any parameter")

    params_fitted <- lapply(setNames(param_names, param_names), function(p)
      predict(model, what = p, type = "response", newdata = model_data))

    cdf_values <- do.call(dist_funcs$p_func,
                          c(list(q = model_data$y), params_fitted))
    cdf_valid  <- cdf_values >= 0 & cdf_values <= 1
    n_invalid  <- sum(!cdf_valid)
    if (n_invalid > 0)
      logger::log_info(paste0("  ", n_invalid, " CDF values outside [0,1] -> NA"))

    params_harmonised <- lapply(setNames(param_names, param_names), function(p)
      compute_harmonised_param(model, model_data, p, batch_var))

    harmonised_std             <- do.call(dist_funcs$q_func,
                                          c(list(p = cdf_values), params_harmonised))
    harmonised_std[!cdf_valid] <- NA_real_

    harmonised_raw  <- harmonised_std * y_sd + y_mean
    harmonised_orig <- if (log_transform) exp(harmonised_raw) - 1 else harmonised_raw

    out_df <- data.frame(
      id                   = model_data[[id_var]],
      batch                = model_data[[batch_var]],
      age                  = model_data$age,
      sex                  = model_data$sex,
      original_value       = original_y,
      harmonised_value     = harmonised_orig,
      harmonisation_effect = harmonised_orig - original_y,
      cdf_value            = cdf_values,
      cdf_valid            = cdf_valid,
      distribution_family  = family_name,
      log_transform        = log_transform
    )

    if (has_wave) out_df$wave <- model_data$wave

    for (p in param_names) {
      out_df[[paste0(p, "_fitted")]]     <- params_fitted[[p]]
      out_df[[paste0(p, "_harmonised")]] <- params_harmonised[[p]]
    }

    if (generate_normative_scores) {
      z_scores              <- qnorm(cdf_values)
      z_scores[!cdf_valid] <- NA_real_
      out_df$normative_z_score <- z_scores
      out_df$normative_centile <- pnorm(z_scores) * 100
    }

    write.csv(out_df, out_csv, row.names = FALSE)

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    logger::log_info(paste0("  Done in ", format_time(elapsed)))

    write.csv(data.frame(feature        = feature_name,
                         status         = "success",
                         distribution   = family_name,
                         start_time     = as.character(start_time),
                         end_time       = as.character(Sys.time()),
                         time_secs      = elapsed,
                         n_observations = nrow(out_df),
                         n_invalid_cdf  = n_invalid,
                         log_transform  = log_transform),
              out_timing, row.names = FALSE)

    list(status = "success", feature = feature_name, distribution = family_name,
         processing_time = elapsed, n_observations = nrow(out_df),
         output_file = out_csv)

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
                                        batch_var                 = "cohort",
                                        id_var                    = "id",
                                        log_transform             = FALSE,
                                        generate_normative_scores = TRUE,
                                        feature_subset            = NULL,
                                        n_cores                   = 1) {
  overall_start <- Sys.time()

  model_files <- list.files(model_base_dir, pattern = "_model\\.rds$",
                            recursive = TRUE, full.names = TRUE)
  if (length(model_files) == 0)
    stop("No model files found in: ", model_base_dir)

  if (!is.null(feature_subset)) {
    pattern     <- paste(paste0("feature_", feature_subset, "/"), collapse = "|")
    model_files <- model_files[grepl(pattern, model_files)]
  }

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

  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, {
      library(gamlss); library(gamlss.dist); library(dplyr); library(logger)
    })
    parallel::clusterExport(cl,
      c("data", "output_dir", "batch_var", "id_var", "log_transform",
        "generate_normative_scores", "harmonise_gamlss_feature", "load_scaling",
        "get_distribution_functions", "get_n_params", "has_random_in_fitted",
        "compute_harmonised_param", "apply_inverse_link", "format_time"),
      envir = environment()
    )
    results <- parallel::parLapply(cl, model_files, process_one)
    parallel::stopCluster(cl)
  } else {
    results <- vector("list", length(model_files))
    for (i in seq_along(model_files)) {
      logger::log_info(paste0("Feature ", i, "/", length(model_files)))
      results[[i]] <- process_one(model_files[i])
      elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
      logger::log_info(paste0("  ", round(100 * i / length(model_files), 1),
                              "% | est. remaining: ",
                              format_time(elapsed / i * (length(model_files) - i))))
    }
  }

  successes <- sum(sapply(results, function(x) x$status == "success"))
  failures  <- length(results) - successes
  elapsed   <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))

  summary_df <- data.frame(
    feature        = sapply(results, "[[", "feature"),
    status         = sapply(results, "[[", "status"),
    distribution   = sapply(results, function(x)
      if (x$status == "success") x$distribution else NA_character_),
    time_seconds   = sapply(results, "[[", "processing_time"),
    n_observations = sapply(results, function(x)
      if (x$status == "success") x$n_observations else NA_integer_)
  )
  write.csv(summary_df,
            file.path(output_dir, "harmonisation_summary.csv"), row.names = FALSE)

  invisible(list(results    = results,
                 successes  = successes,
                 failures   = failures,
                 total_time = format_time(elapsed),
                 output_dir = output_dir))
}

combine_harmonised_results <- function(harmonised_output_dir,
                                       id_var                    = "id",
                                       generate_normative_scores = TRUE) {
  logger::log_info(paste0("Combining harmonised results (wide format)"))

  files <- list.files(harmonised_output_dir, pattern = "_harmonised\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0)
    stop("No harmonised CSV files found in: ", harmonised_output_dir)

  meta_cols  <- c(id_var, "batch", "age", "sex", "wave")
  norm_cols  <- c("normative_z_score", "normative_centile")

  all_data   <- lapply(files, read.csv, stringsAsFactors = FALSE)
  found_meta <- Reduce(intersect, lapply(all_data, function(d)
    intersect(meta_cols, names(d))))

  n_rows <- nrow(all_data[[1]])
  if (!all(sapply(all_data, nrow) == n_rows))
    stop("Per-feature CSVs have different row counts -- ensure the same input ",
         "data and filtering were used across all features.")

  wide_harm <- all_data[[1]][, found_meta, drop = FALSE]
  wide_norm <- all_data[[1]][, found_meta, drop = FALSE]

  for (i in seq_along(files)) {
    dat       <- all_data[[i]]
    feat_name <- gsub("_harmonised\\.csv$", "", basename(files[i]))

    if ("harmonised_value" %in% names(dat)) {
      h <- dat[, "harmonised_value", drop = FALSE]
      names(h) <- feat_name
      wide_harm <- cbind(wide_harm, h)
    }

    if (generate_normative_scores) {
      avail <- norm_cols[norm_cols %in% names(dat)]
      if (length(avail) > 0) {
        n <- dat[, avail, drop = FALSE]
        names(n) <- paste0(feat_name, ".", avail)
        wide_norm <- cbind(wide_norm, n)
      }
    }
  }

  harm_path <- file.path(harmonised_output_dir, "combined_harmonised.csv")
  write.csv(wide_harm, harm_path, row.names = FALSE)
  logger::log_info(paste0("Saved: ", harm_path,
                          " (", nrow(wide_harm), " rows, ",
                          length(files), " features)"))

  if (generate_normative_scores) {
    norm_path <- file.path(harmonised_output_dir, "combined_normative.csv")
    write.csv(wide_norm, norm_path, row.names = FALSE)
    logger::log_info(paste0("Saved: ", norm_path,
                            " (", nrow(wide_norm), " rows, ",
                            length(files), " features)"))
  }

  invisible(list(harmonised = wide_harm,
                 normative  = if (generate_normative_scores) wide_norm else NULL))
}