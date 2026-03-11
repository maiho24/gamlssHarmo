try_gamlss <- function(mu.formula, sigma.formula, nu.formula, tau.formula,
                       data, family, n_cyc = 200) {
  model <- try(gamlss(
    formula       = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula    = nu.formula,
    tau.formula   = tau.formula,
    data          = data,
    family        = family,
    control       = gamlss.control(n.cyc = n_cyc, trace = FALSE,
                                   mu.trace = FALSE, sigma.trace = FALSE)
  ), silent = TRUE)
  if (inherits(model, "try-error") || !isTRUE(model$converged)) return(NULL)
  model
}

build_family_specs <- function(family_name, nu_f, tau_f) {
  four_param  <- c("SHASH", "BCT", "BCPE", "ST4", "EGB2")
  three_param <- c("GG", "SN1", "TF", "ST3", "BCCG")
  intercept   <- as.formula("~ 1")

  if (family_name %in% four_param) {
    return(list(
      list(name = paste0(family_name, "_full"),
           nu.formula = nu_f,      tau.formula = tau_f),
      list(name = paste0(family_name, "_nu_only"),
           nu.formula = nu_f,      tau.formula = intercept),
      list(name = paste0(family_name, "_tau_only"),
           nu.formula = intercept, tau.formula = tau_f),
      list(name = paste0(family_name, "_intercept"),
           nu.formula = intercept, tau.formula = intercept)
    ))
  }
  if (family_name %in% three_param) {
    return(list(
      list(name = paste0(family_name, "_nu"),
           nu.formula = nu_f,      tau.formula = NULL),
      list(name = paste0(family_name, "_intercept"),
           nu.formula = intercept, tau.formula = NULL)
    ))
  }
  list(list(name = family_name, nu.formula = NULL, tau.formula = NULL))
}

get_family_fn <- function(family_name) {
  fn <- tryCatch(get(family_name, envir = asNamespace("gamlss.dist")),
                 error = function(e) NULL)
  if (is.null(fn))
    stop("Unknown GAMLSS family: '", family_name, "'")
  fn
}

fit_gamlss_for_feature <- function(data, feature_name, model_dir,
                                   formula_terms, batch_var, id_var,
                                   longitudinal  = FALSE,
                                   log_transform = FALSE,
                                   family_order  = c("SHASH", "GG", "NO")) {
  feature_start <- Sys.time()
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    model       = file.path(model_dir, paste0(feature_name, "_model.rds")),
    summary     = file.path(model_dir, paste0(feature_name, "_summary.txt")),
    predictions = file.path(model_dir, paste0(feature_name, "_predictions.csv")),
    metrics     = file.path(model_dir, paste0(feature_name, "_metrics.csv")),
    timing      = file.path(model_dir, paste0(feature_name, "_timing.csv")),
    diagnostics = file.path(model_dir, paste0(feature_name, "_diagnostics.pdf")),
    scaling     = file.path(model_dir, paste0(feature_name, "_scaling.csv"))
  )

  tryCatch({
    logger::log_info(paste0("Processing: ", feature_name))

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

    n_obs     <- nrow(model_data)
    n_removed <- nrow(data) - n_obs
    logger::log_info(paste0("  n = ", n_obs, " (", n_removed, " removed)"))

    if (n_obs < 100)
      stop("Insufficient data (", n_obs, " < 100)")

    n_batches <- length(unique(model_data[[batch_var]]))
    if (n_batches < 2) {
      logger::log_info(paste0("  Skipping -- only 1 batch level"))
      elapsed <- as.numeric(difftime(Sys.time(), feature_start, units = "secs"))
      write.csv(data.frame(feature = feature_name, reason = "single batch",
                           timestamp = as.character(Sys.time())),
                paths$timing, row.names = FALSE)
      return(list(status = "skipped_single_batch", feature = feature_name,
                  processing_time = elapsed))
    }

    if (log_transform) model_data[[feature_name]] <- log(model_data[[feature_name]] + 1)
    y_mean <- mean(model_data[[feature_name]])
    y_sd   <- sd(model_data[[feature_name]])
    model_data$y <- (model_data[[feature_name]] - y_mean) / y_sd

    write.csv(data.frame(feature = feature_name, log_transform = log_transform,
                         y_mean = y_mean, y_sd = y_sd,
                         timestamp = as.character(Sys.time())),
              paths$scaling, row.names = FALSE)

    mu_f    <- update(terms_to_formula(formula_terms$mu,    batch_var, id_var), y ~ .)
    sigma_f <- terms_to_formula(formula_terms$sigma, batch_var, id_var)
    nu_f    <- terms_to_formula(formula_terms$nu,    batch_var, id_var)
    tau_f   <- terms_to_formula(formula_terms$tau,   batch_var, id_var)

    logger::log_info(paste0("  mu:    ", deparse(mu_f)))
    logger::log_info(paste0("  sigma: ", deparse(sigma_f)))

    fit_start  <- Sys.time()
    model      <- NULL
    final_spec <- NULL

    for (fam_name in family_order) {
      fam_fn <- get_family_fn(fam_name)
      for (spec in build_family_specs(fam_name, nu_f, tau_f)) {
        logger::log_info(paste0("  Trying: ", spec$name))
        model <- try_gamlss(mu.formula    = mu_f,
                            sigma.formula = sigma_f,
                            nu.formula    = spec$nu.formula,
                            tau.formula   = spec$tau.formula,
                            data          = model_data,
                            family        = fam_fn())
        if (!is.null(model)) {
          logger::log_info(paste0("  Converged: ", spec$name,
                                  " | AIC: ", round(AIC(model), 2)))
          final_spec <- spec$name
          break
        }
      }
      if (!is.null(model)) break
    }

    if (is.null(model))
      stop("All model specifications failed to converge")

    fit_time        <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))
    model$call$data <- model_data

    saveRDS(model, paths$model)
    tryCatch(
      capture.output(summary(model), file = paths$summary),
      error   = function(e) logger::log_info(paste0("  Summary skipped: ", e$message)),
      warning = function(w) logger::log_info(paste0("  Summary warning: ", w$message))
    )

    mu_pred    <- predict(model, what = "mu",    type = "response", newdata = model_data)
    sigma_pred <- predict(model, what = "sigma", type = "response", newdata = model_data)
    nu_pred    <- tryCatch(
      predict(model, what = "nu",  type = "response", newdata = model_data),
      error = function(e) rep(NA_real_, nrow(model_data))
    )
    tau_pred   <- tryCatch(
      predict(model, what = "tau", type = "response", newdata = model_data),
      error = function(e) rep(NA_real_, nrow(model_data))
    )

    mu_original <- mu_pred * y_sd + y_mean
    if (log_transform) mu_original <- exp(mu_original) - 1

    pred_df <- data.frame(
      id             = model_data[[id_var]],
      batch          = model_data[[batch_var]],
      age            = model_data$age,
      sex            = model_data$sex,
      observed_value = model_data$y,
      mu             = mu_original,
      sigma          = sigma_pred * y_sd,
      nu             = nu_pred,
      tau            = tau_pred,
      scaling_mean   = y_mean,
      scaling_sd     = y_sd,
      log_transform  = log_transform
    )
    if (has_wave) pred_df$wave <- model_data$wave
    write.csv(pred_df, paths$predictions, row.names = FALSE)

    resid_vals <- residuals(model)
    write.csv(data.frame(
      feature           = feature_name,
      n_observations    = n_obs,
      n_batches         = n_batches,
      distribution      = final_spec,
      AIC               = AIC(model),
      BIC               = BIC(model),
      MSE               = mean(resid_vals^2),
      fitting_time_secs = fit_time,
      converged         = TRUE,
      longitudinal_mode = longitudinal,
      log_transform     = log_transform
    ), paths$metrics, row.names = FALSE)

    tryCatch({
      grDevices::pdf(paths$diagnostics, width = 12, height = 8)
      plot(model)
      grDevices::dev.off()
    }, error = function(e) logger::log_info(paste0("  Diagnostics skipped: ", e$message)))

    feature_time <- as.numeric(difftime(Sys.time(), feature_start, units = "secs"))
    logger::log_info(paste0("  Done in ", format_time(feature_time)))

    write.csv(data.frame(
      feature           = feature_name,
      status            = "success",
      distribution      = final_spec,
      start_time        = as.character(feature_start),
      end_time          = as.character(Sys.time()),
      total_time_secs   = feature_time,
      fit_time_secs     = fit_time,
      longitudinal_mode = longitudinal,
      log_transform     = log_transform
    ), paths$timing, row.names = FALSE)

    list(status = "success", feature = feature_name, distribution = final_spec,
         processing_time = feature_time, n_observations = n_obs, n_batches = n_batches)

  }, error = function(e) {
    logger::log_error(paste0("Error in ", feature_name, ": ", e$message))
    feature_time <- as.numeric(difftime(Sys.time(), feature_start, units = "secs"))
    write.csv(data.frame(feature       = feature_name,
                         status        = "error",
                         error_message = e$message,
                         start_time    = as.character(feature_start),
                         end_time      = as.character(Sys.time()),
                         total_time_secs = feature_time),
              paths$timing, row.names = FALSE)
    list(status = "error", feature = feature_name,
         error = e$message, processing_time = feature_time)
  })
}

run_gamlss_harmonisation <- function(data, features, output_dir, formula_terms,
                                     batch_var     = "cohort",
                                     id_var        = "id",
                                     longitudinal  = FALSE,
                                     log_transform = FALSE,
                                     family_order  = c("SHASH", "GG", "NO"),
                                     n_cores       = 1) {
  overall_start <- Sys.time()
  logger::log_info(paste0("Fitting ", length(features), " features"))

  model_dir <- file.path(output_dir, "models")
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  fit_one <- function(feat) {
    fit_gamlss_for_feature(
      data          = data,
      feature_name  = feat,
      model_dir     = file.path(model_dir, paste0("feature_", feat)),
      formula_terms = formula_terms,
      batch_var     = batch_var,
      id_var        = id_var,
      longitudinal  = longitudinal,
      log_transform = log_transform,
      family_order  = family_order
    )
  }

  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, { library(gamlss); library(dplyr); library(logger) })
    parallel::clusterExport(cl,
      c("data", "formula_terms", "batch_var", "id_var", "longitudinal",
        "log_transform", "family_order", "model_dir", "fit_gamlss_for_feature",
        "try_gamlss", "build_family_specs", "get_family_fn",
        "terms_to_formula", "substitute_placeholders", "format_time"),
      envir = environment()
    )
    results <- parallel::parLapply(cl, features, fit_one)
    parallel::stopCluster(cl)
  } else {
    results <- vector("list", length(features))
    timings <- data.frame()
    for (i in seq_along(features)) {
      logger::log_info(paste0("Feature ", i, "/", length(features), ": ", features[i]))
      results[[i]] <- fit_one(features[i])
      timings <- rbind(timings, data.frame(
        feature      = features[i],
        status       = results[[i]]$status,
        distribution = if (results[[i]]$status == "success")
          results[[i]]$distribution else NA_character_,
        time_seconds = results[[i]]$processing_time
      ))
      write.csv(timings, file.path(output_dir, "feature_timings.csv"), row.names = FALSE)
      elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
      logger::log_info(paste0("  ", round(100 * i / length(features), 1),
                              "% | est. remaining: ",
                              format_time(elapsed / i * (length(features) - i))))
    }
  }

  names(results) <- features
  successes <- sum(sapply(results, function(x) x$status == "success"))
  skipped   <- sum(sapply(results, function(x) x$status == "skipped_single_batch"))
  failures  <- length(features) - successes - skipped
  total_t   <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))

  summary_df <- data.frame(
    feature        = sapply(results, "[[", "feature"),
    status         = sapply(results, "[[", "status"),
    distribution   = sapply(results, function(x)
      if (x$status == "success") x$distribution else NA_character_),
    time_seconds   = sapply(results, "[[", "processing_time"),
    n_observations = sapply(results, function(x)
      if (x$status == "success") x$n_observations else NA_integer_)
  )
  write.csv(summary_df, file.path(output_dir, "model_summary.csv"), row.names = FALSE)

  invisible(list(results    = results,
                 successes  = successes,
                 skipped    = skipped,
                 failures   = failures,
                 total_time = format_time(total_t),
                 output_dir = output_dir))
}