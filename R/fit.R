# R/fit.R
#
# GAMLSS fitting functions for gamlssHarmo.
#
# Supports two data modes controlled by the `discrete` argument:
#
# CONTINUOUS (discrete = FALSE):
#   Features are fitted on the raw scale (no z-standardisation). An optional
#   log(y+1) transform is applied first if log_transform = TRUE.
#   Families: SHASH, JSU, GG, BCCG, BCT, BCPE, LOGNO, TF, SN1, PE2, NO,
#             exGAUS (continuous, 3-param, for reaction-time / cognitive data).
#
# DISCRETE / COUNT (discrete = TRUE):
#   Raw integer counts are passed directly to gamlss().
#   Families: PO, NBI, NBII, ZIP, ZINBI, ZANBI,
#             PIG, ZIPIG (Poisson-Inverse Gaussian family),
#             SICHEL, ZISICHEL, ZASICHEL (Sichel family; heavy-tailed counts),
#             DPO (Double Poisson; handles under-dispersion).

# ---------------------------------------------------------------------------
# Low-level GAMLSS wrapper
# ---------------------------------------------------------------------------

try_gamlss <- function(mu.formula, sigma.formula, nu.formula, tau.formula,
                       data, family, n_cyc = 200, c_crit = 0.001) {
  model <- try(gamlss(
    formula       = mu.formula,
    sigma.formula = sigma.formula,
    nu.formula    = nu.formula,
    tau.formula   = tau.formula,
    data          = data,
    family        = family,
    control       = gamlss.control(n.cyc = n_cyc, c.crit = c_crit, trace = FALSE,
                                   mu.trace = FALSE, sigma.trace = FALSE)
  ), silent = TRUE)
  if (inherits(model, "try-error") || !isTRUE(model$converged)) return(NULL)
  model
}

# ---------------------------------------------------------------------------
# Helper: does a formula contain a random() term?
# ---------------------------------------------------------------------------

formula_has_random <- function(f) {
  if (is.null(f)) return(FALSE)
  grepl("random", deparse1(f))
}

# ---------------------------------------------------------------------------
# Build a family's fallback ladder, richest spec first.
#
# Each spec is a named list: $name plus $sigma.formula/$nu.formula/$tau.formula.
# Later specs reduce nu/tau (and, for discrete, sigma) to intercept; mu is never
# simplified. Per-family sequences are tabulated in README.md ("Family fallback").
#
# Naming: _full = all extra params modelled; _no_X = param X at intercept;
# _intercept = all extra params at intercept.
# ---------------------------------------------------------------------------

build_family_specs <- function(family_name, nu_f, tau_f,
                               sigma_f  = NULL,
                               discrete = FALSE) {

  continuous_four   <- c("SHASH", "JSU", "BCT", "BCPE", "ST3", "ST4", "EGB2")
  continuous_three  <- c("GG", "SN1", "TF", "PE2", "BCCG", "exGAUS")

  discrete_four     <- c("ZISICHEL", "ZASICHEL")
  discrete_three_zi <- c("ZINBI", "ZANBI", "ZIPIG")
  discrete_three_sh <- c("SICHEL")
  discrete_two      <- c("NBI", "NBII", "ZIP", "PIG", "DPO")
  one_param         <- c("PO", "EXP")

  intercept <- as.formula("~ 1")

  sigma_fallback <- if (formula_has_random(sigma_f)) intercept else sigma_f

  sf    <- sigma_f
  sf_fb <- sigma_fallback

  if (family_name %in% continuous_four) {
    return(list(
      list(name          = paste0(family_name, "_full"),
           sigma.formula = sf,    nu.formula = nu_f,      tau.formula = tau_f),
      list(name          = paste0(family_name, "_no_tau"),
           sigma.formula = sf,    nu.formula = nu_f,      tau.formula = intercept),
      list(name          = paste0(family_name, "_no_nu"),
           sigma.formula = sf,    nu.formula = intercept, tau.formula = tau_f),
      list(name          = paste0(family_name, "_intercept"),
           sigma.formula = sf,    nu.formula = intercept, tau.formula = intercept)
    ))
  }

  if (family_name %in% continuous_three) {
    return(list(
      list(name          = paste0(family_name, "_full"),
           sigma.formula = sf, nu.formula = nu_f,      tau.formula = NULL),
      list(name          = paste0(family_name, "_intercept"),
           sigma.formula = sf, nu.formula = intercept, tau.formula = NULL)
    ))
  }

  if (family_name %in% discrete_four) {
    return(list(
      list(name          = paste0(family_name, "_full"),
           sigma.formula = sf,    nu.formula = nu_f,      tau.formula = tau_f),
      list(name          = paste0(family_name, "_no_tau"),
           sigma.formula = sf,    nu.formula = nu_f,      tau.formula = intercept),
      list(name          = paste0(family_name, "_no_nu"),
           sigma.formula = sf,    nu.formula = intercept, tau.formula = tau_f),
      list(name          = paste0(family_name, "_no_nu_no_tau"),
           sigma.formula = sf,    nu.formula = intercept, tau.formula = intercept),
      list(name          = paste0(family_name, "_no_sigma"),
           sigma.formula = sf_fb, nu.formula = nu_f,      tau.formula = tau_f),
      list(name          = paste0(family_name, "_no_sigma_no_tau"),
           sigma.formula = sf_fb, nu.formula = nu_f,      tau.formula = intercept),
      list(name          = paste0(family_name, "_no_sigma_no_nu"),
           sigma.formula = sf_fb, nu.formula = intercept, tau.formula = tau_f),
      list(name          = paste0(family_name, "_intercept"),
           sigma.formula = sf_fb, nu.formula = intercept, tau.formula = intercept)
    ))
  }

  if (family_name %in% discrete_three_zi) {
    return(list(
      list(name          = paste0(family_name, "_full"),
           sigma.formula = sf,    nu.formula = nu_f,      tau.formula = NULL),
      list(name          = paste0(family_name, "_no_nu"),
           sigma.formula = sf,    nu.formula = intercept, tau.formula = NULL),
      list(name          = paste0(family_name, "_no_sigma"),
           sigma.formula = sf_fb, nu.formula = nu_f,      tau.formula = NULL),
      list(name          = paste0(family_name, "_intercept"),
           sigma.formula = sf_fb, nu.formula = intercept, tau.formula = NULL)
    ))
  }

  if (family_name %in% discrete_three_sh) {
    return(list(
      list(name          = paste0(family_name, "_full"),
           sigma.formula = sf,    nu.formula = nu_f,      tau.formula = NULL),
      list(name          = paste0(family_name, "_no_nu"),
           sigma.formula = sf,    nu.formula = intercept, tau.formula = NULL),
      list(name          = paste0(family_name, "_no_sigma"),
           sigma.formula = sf_fb, nu.formula = nu_f,      tau.formula = NULL),
      list(name          = paste0(family_name, "_intercept"),
           sigma.formula = sf_fb, nu.formula = intercept, tau.formula = NULL)
    ))
  }

  if (family_name %in% discrete_two) {
    return(list(
      list(name          = paste0(family_name, "_full"),
           sigma.formula = sf,    nu.formula = NULL, tau.formula = NULL),
      list(name          = paste0(family_name, "_no_sigma"),
           sigma.formula = sf_fb, nu.formula = NULL, tau.formula = NULL)
    ))
  }

  if (family_name %in% one_param) {
    return(list(
      list(name          = family_name,
           sigma.formula = intercept, nu.formula = NULL, tau.formula = NULL)
    ))
  }

  list(list(name          = family_name,
            sigma.formula = sf, nu.formula = NULL, tau.formula = NULL))
}

# ---------------------------------------------------------------------------
# Retrieve a GAMLSS family function by name from gamlss.dist.
# ---------------------------------------------------------------------------

get_family_fn <- function(family_name) {
  fn <- tryCatch(get(family_name, envir = asNamespace("gamlss.dist")),
                 error = function(e) NULL)
  if (is.null(fn))
    stop("Unknown GAMLSS family: '", family_name, "'")
  fn
}

read_feature_recommendations <- function(path) {
  if (!file.exists(path))
    stop("Feature recommendations CSV not found: ", path)

  df <- read.csv(path, stringsAsFactors = FALSE)

  if (!"feature" %in% names(df))
    stop("feature_recommendations CSV must contain a 'feature' column")

  MISSING <- structure(list(), class = "missing_column")

  parse_semi <- function(x) {
    if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(NULL)
    trimws(strsplit(x, ";")[[1L]])
  }

  read_col <- function(r, col) {
    if (col %in% names(r)) parse_semi(r[[col]]) else MISSING
  }

  setNames(lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    list(
      family_order = if ("family_order" %in% names(r)) parse_semi(r$family_order) else NULL,
      mu_terms     = read_col(r, "mu_terms"),
      sigma_terms  = read_col(r, "sigma_terms"),
      nu_terms     = read_col(r, "nu_terms"),
      tau_terms    = read_col(r, "tau_terms")
    )
  }), df$feature)
}

resolve_feature_formula_terms <- function(feature_name, formula_terms,
                                           feature_rec) {
  if (is.null(feature_rec))
    return(list(formula_terms = formula_terms, family_order = NULL))

  MISSING <- structure(list(), class = "missing_column")
  is_missing <- function(x) inherits(x, "missing_column")

  ft <- formula_terms

  if (!is_missing(feature_rec$mu_terms))    ft$mu    <- as.list(feature_rec$mu_terms    %||% list())
  if (!is_missing(feature_rec$sigma_terms)) ft$sigma <- as.list(feature_rec$sigma_terms %||% list())
  if (!is_missing(feature_rec$nu_terms))    ft$nu    <- as.list(feature_rec$nu_terms    %||% list())
  if (!is_missing(feature_rec$tau_terms))   ft$tau   <- as.list(feature_rec$tau_terms   %||% list())

  list(formula_terms = ft,
       family_order  = feature_rec$family_order)
}

# ---------------------------------------------------------------------------
# Write a human-readable mathematical description of the fitted model.
# ---------------------------------------------------------------------------

family_description <- function(base_name) {
  desc <- list(
    NO       = "Normal: Y_i ~ N(mu_i, sigma_i^2)",
    TF       = "t-family: Y_i ~ t(mu_i, sigma_i, nu_i)  [nu_i = degrees of freedom]",
    SHASH    = "Sinh-arcsinh: Y_i ~ SHASH(mu_i, sigma_i, nu_i, tau_i)  [Jones & Pewsey 2009]",
    JSU      = "Johnson SU: Y_i ~ JSU(mu_i, sigma_i, nu_i, tau_i)  [Johnson 1949]",
    SN1      = "Skew-Normal type 1: Y_i ~ SN1(mu_i, sigma_i, nu_i)",
    PE2      = "Power Exponential type 2: Y_i ~ PE2(mu_i, sigma_i, nu_i)  [nu_i = tail index]",
    GG       = "Generalised Gamma: Y_i ~ GG(mu_i, sigma_i, nu_i)  [positive support]",
    LOGNO    = "Log-Normal: log(Y_i) ~ N(mu_i, sigma_i^2)  [positive support]",
    BCCG     = "Box-Cox Cole-Green: Y_i ~ BCCG(mu_i, sigma_i, nu_i)  [positive support]",
    BCT      = "Box-Cox t: Y_i ~ BCT(mu_i, sigma_i, nu_i, tau_i)  [positive support]",
    BCPE     = "Box-Cox Power Exponential: Y_i ~ BCPE(mu_i, sigma_i, nu_i, tau_i)",
    exGAUS   = "Ex-Gaussian: Y_i ~ N(mu_i, sigma_i^2) + Exp(nu_i)\n  mu_i = normal mean, sigma_i = normal SD, nu_i = exponential rate\n  Common in reaction-time and cognitive test score distributions",
    PO       = "Poisson: Y_i ~ Poisson(mu_i)\n  E[Y] = Var[Y] = mu_i",
    NBI      = "Negative Binomial I: Y_i ~ NBI(mu_i, sigma_i)\n  E[Y] = mu_i,  Var[Y] = mu_i + sigma_i * mu_i^2",
    NBII     = "Negative Binomial II: Y_i ~ NBII(mu_i, sigma_i)\n  E[Y] = mu_i,  Var[Y] = mu_i * (1 + sigma_i * mu_i)",
    ZIP      = "Zero-Inflated Poisson: Y_i ~ (1-sigma_i)*Poisson(mu_i) + sigma_i*I(0)\n  sigma_i = P(structural zero)  [logit link in gamlss.dist]",
    ZINBI    = "Zero-Inflated NBI: Y_i ~ nu_i*I(0) + (1-nu_i)*NBI(mu_i, sigma_i)\n  nu_i = P(structural zero)\n  sigma_i = NBI overdispersion",
    ZANBI    = "Zero-Adjusted NBI (hurdle): P(Y=0) = nu_i,  P(Y=y|y>0) ~ NBI(mu_i,sigma_i) truncated\n  nu_i = P(Y = 0)\n  sigma_i = NBI overdispersion",
    PIG      = "Poisson-Inverse Gaussian: Y_i ~ PIG(mu_i, sigma_i)\n  E[Y] = mu_i,  heavier right tail than NBI",
    ZIPIG    = "Zero-Inflated PIG: Y_i ~ nu_i*I(0) + (1-nu_i)*PIG(mu_i, sigma_i)\n  nu_i = P(structural zero)",
    SICHEL   = "Sichel: Y_i ~ Sichel(mu_i, sigma_i, nu_i)\n  Poisson mixture with inverse-Gaussian mixing; nu_i = IG shape\n  E[Y] = mu_i,  heavier tail than NBI/PIG",
    ZISICHEL = "Zero-Inflated Sichel: Y_i ~ tau_i*I(0) + (1-tau_i)*Sichel(mu_i, sigma_i, nu_i)\n  tau_i = P(structural zero)\n  nu_i = IG shape,  sigma_i = Sichel dispersion",
    ZASICHEL = "Zero-Adjusted Sichel (hurdle): P(Y=0) = tau_i,  P(Y=y|y>0) ~ Sichel(mu_i,sigma_i,nu_i) truncated\n  nu_i = IG shape,  tau_i = P(Y = 0)\n  Preferred when true zeros are clinically distinct",
    DPO      = "Double Poisson: Y_i ~ DPO(mu_i, sigma_i)\n  sigma_i < 1 = overdispersion,  sigma_i > 1 = underdispersion\n  Only family that handles underdispersion natively"
  )
  desc[[base_name]] %||% paste0(base_name, ": no description available")
}

param_description <- function(base_name) {
  list(
    NO       = list(mu = "mean",              sigma = "standard deviation"),
    TF       = list(mu = "location",          sigma = "scale",               nu = "degrees of freedom"),
    SHASH    = list(mu = "location",          sigma = "scale",               nu = "skewness",    tau = "kurtosis"),
    JSU      = list(mu = "location",          sigma = "scale",               nu = "skewness",    tau = "kurtosis"),
    SN1      = list(mu = "location",          sigma = "scale",               nu = "skewness"),
    PE2      = list(mu = "location",          sigma = "scale",               nu = "tail index"),
    GG       = list(mu = "mean",              sigma = "dispersion",          nu = "shape"),
    LOGNO    = list(mu = "log-mean",          sigma = "log-scale"),
    BCCG     = list(mu = "median",            sigma = "approx CV",           nu = "Box-Cox power"),
    BCT      = list(mu = "median",            sigma = "approx CV",           nu = "Box-Cox power", tau = "degrees of freedom"),
    BCPE     = list(mu = "median",            sigma = "approx CV",           nu = "Box-Cox power", tau = "tail index"),
    exGAUS   = list(mu = "normal mean",       sigma = "normal SD",           nu = "exponential rate"),
    PO       = list(mu = "mean count"),
    NBI      = list(mu = "mean count",        sigma = "overdispersion"),
    NBII     = list(mu = "mean count",        sigma = "overdispersion"),
    ZIP      = list(mu = "Poisson mean",      sigma = "P(structural zero)"),
    ZINBI    = list(mu = "NBI mean",          sigma = "NBI overdispersion",  nu = "P(structural zero)"),
    ZANBI    = list(mu = "NBI mean",          sigma = "NBI overdispersion",  nu = "P(Y = 0)"),
    PIG      = list(mu = "mean count",        sigma = "dispersion"),
    ZIPIG    = list(mu = "PIG mean",          sigma = "PIG dispersion",      nu = "P(structural zero)"),
    SICHEL   = list(mu = "mean count",        sigma = "dispersion",          nu = "IG shape"),
    ZISICHEL = list(mu = "Sichel mean",       sigma = "Sichel dispersion",   nu = "IG shape",    tau = "P(structural zero)"),
    ZASICHEL = list(mu = "Sichel mean",       sigma = "Sichel dispersion",   nu = "IG shape",    tau = "P(Y = 0)"),
    DPO      = list(mu = "mean count",        sigma = "dispersion index")
  )[[base_name]] %||% list(mu = "mu", sigma = "sigma", nu = "nu", tau = "tau")
}

write_model_equation <- function(path, feature_name, final_spec,
                                  model, mu_f, sigma_f, nu_f, tau_f,
                                  n_obs, n_batches, fit_time,
                                  log_transform, discrete) {
  base_name   <- strsplit(final_spec, "_")[[1L]][1L]
  n_params    <- get_n_params(final_spec)
  param_names <- c("mu",
                   if (n_params >= 2L) "sigma",
                   if (n_params >= 3L) "nu",
                   if (n_params >= 4L) "tau")

  link_of <- function(p) {
    tryCatch(model[[paste0(p, ".link")]], error = function(e) "?")
  }

  fmt_formula <- function(p, f) {
    rhs  <- deparse1(f)
    rhs  <- sub("^~\\s*", "", rhs)
    link <- link_of(p)
    sprintf("  %s(%s_i)  =  %s", link, p, rhs)
  }

  param_lines <- vapply(param_names, function(p) {
    f <- switch(p, mu = mu_f, sigma = sigma_f, nu = nu_f, tau = tau_f)
    fmt_formula(p, f)
  }, character(1L))

  param_interp <- vapply(param_names, function(p) {
    link  <- link_of(p)
    pdesc <- param_description(base_name)[[p]] %||% p
    sprintf("  %-6s  link: %-10s  interpretation: %s", p, link, pdesc)
  }, character(1L))

  rand_params <- param_names[vapply(param_names, function(p) {
    f <- switch(p, mu = mu_f, sigma = sigma_f, nu = nu_f, tau = tau_f)
    grepl("random", deparse1(f))
  }, logical(1L))]

  rand_lines <- if (length(rand_params) > 0L) {
    c("Random batch effects (independently estimated per parameter):",
      vapply(rand_params, function(p)
        sprintf("  u_%s[batch_i] ~ N(0, tau_%s^2)   [in %s equation]", p, p, p),
        character(1L)))
  } else {
    "No random batch effects detected."
  }

  transform_note <- if (discrete)
    "Response scale: raw integer counts (no transformation)"
  else if (log_transform)
    "Response scale: log(y + 1)  [log_transform = TRUE]"
  else
    "Response scale: raw values (no transformation, no z-standardisation)"

  lines <- c(
    paste0("Feature:         ", feature_name),
    paste0("Fitted:          ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Specification:   ", final_spec),
    paste0("Base family:     ", base_name),
    "",
    "Distribution",
    "------------",
    family_description(base_name),
    "",
    "Parameters",
    "----------",
    param_interp,
    "",
    "Fitted equations",
    "----------------",
    param_lines,
    "",
    "  Smoothers:",
    "    pb(age)  = P-spline smoother on age",
    "    random() = random intercept per batch level (BLUPs)",
    "",
    rand_lines,
    "",
    "Response",
    "--------",
    transform_note,
    "",
    "Fit statistics",
    "--------------",
    sprintf("  n observations:  %d", n_obs),
    sprintf("  n batches:       %d", n_batches),
    sprintf("  AIC:             %.3f", AIC(model)),
    sprintf("  BIC:             %.3f", BIC(model)),
    sprintf("  Log-likelihood:  %.3f", logLik(model)),
    sprintf("  Fit time:        %s", format_time(fit_time))
  )

  writeLines(lines, con = path)
  invisible(path)
}

# ---------------------------------------------------------------------------
# Fit a single feature.
# ---------------------------------------------------------------------------

fit_gamlss_for_feature <- function(data, feature_name, model_dir,
                                   formula_terms, batch_var, id_var,
                                   longitudinal       = FALSE,
                                   log_transform      = FALSE,
                                   discrete           = FALSE,
                                   family_order       = c("SHASH", "GG", "NO"),
                                   feature_rec        = NULL,
                                   c_crit             = 0.001) {
  feature_start <- Sys.time()
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    model       = file.path(model_dir, paste0(feature_name, "_model.rds")),
    summary     = file.path(model_dir, paste0(feature_name, "_summary.txt")),
    predictions = file.path(model_dir, paste0(feature_name, "_predictions.csv")),
    metrics     = file.path(model_dir, paste0(feature_name, "_metrics.csv")),
    timing      = file.path(model_dir, paste0(feature_name, "_timing.csv")),
    diagnostics = file.path(model_dir, paste0(feature_name, "_diagnostics.pdf")),
    meta        = file.path(model_dir, paste0(feature_name, "_meta.csv"))
  )

  tryCatch({
    logger::log_info(paste0("Processing: ", feature_name,
                            if (discrete) " [discrete]" else ""))

    resolved     <- resolve_feature_formula_terms(feature_name, formula_terms,
                                                   feature_rec)
    eff_formula  <- resolved$formula_terms
    eff_families <- resolved$family_order %||% family_order
    eff_discrete <- discrete

    if (!is.null(resolved$family_order))
      logger::log_info(paste0("  family_order from recommendations: ",
                              paste(eff_families, collapse = " -> ")))

    required_cols <- c(feature_name, "age", "sex", batch_var, id_var)
    missing_cols  <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0L)
      stop("Missing columns: ", paste(missing_cols, collapse = ", "))

    extra_cols <- extract_covariate_cols(eff_formula, batch_var, id_var, data)
    has_wave   <- "wave" %in% names(data)
    keep_cols  <- unique(c(feature_name, "age", "sex", batch_var, id_var,
                           if (has_wave) "wave", extra_cols))

    model_data <- data[, keep_cols, drop = FALSE]
    model_data <- model_data[
      Reduce("&", lapply(keep_cols, function(col) !is.na(model_data[[col]]))), ]
    # Counts and log(y+1)-transformed responses must be non-negative; raw
    # continuous features may legitimately take real (incl. negative) values.
    if (eff_discrete || log_transform)
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
      # Covers rebuild_model_summary()'s columns so the skip still lands in
      # model_summary.csv.
      write.csv(data.frame(feature        = feature_name,
                           status         = "skipped_single_batch",
                           distribution   = NA_character_,
                           discrete       = NA,
                           n_observations = n_obs,
                           n_batches      = n_batches,
                           error_message  = "single batch",
                           time_secs      = elapsed),
                paths$timing, row.names = FALSE)
      return(list(status = "skipped_single_batch", feature = feature_name,
                  n_observations = n_obs, n_batches = n_batches,
                  processing_time = elapsed))
    }

    if (eff_discrete) {
      model_data$y       <- as.integer(round(model_data[[feature_name]]))
      log_transform_used <- FALSE
    } else {
      if (log_transform)
        model_data[[feature_name]] <- log(model_data[[feature_name]] + 1)
      model_data$y       <- model_data[[feature_name]]
      log_transform_used <- log_transform
    }

    write.csv(data.frame(feature       = feature_name,
                         discrete      = eff_discrete,
                         log_transform = log_transform_used),
              paths$meta, row.names = FALSE)

    mu_f    <- update(terms_to_formula(eff_formula$mu,    batch_var, id_var), y ~ .)
    sigma_f <- terms_to_formula(eff_formula$sigma, batch_var, id_var)
    nu_f    <- terms_to_formula(eff_formula$nu,    batch_var, id_var)
    tau_f   <- terms_to_formula(eff_formula$tau,   batch_var, id_var)

    logger::log_info(paste0("  mu:    ", deparse1(mu_f)))
    logger::log_info(paste0("  sigma: ", deparse1(sigma_f)))
    logger::log_info(paste0("  nu:    ", deparse1(nu_f)))
    logger::log_info(paste0("  tau:   ", deparse1(tau_f)))

    fit_start   <- Sys.time()
    model       <- NULL
    final_spec  <- NULL
    eff_sigma_f <- sigma_f

    candidates <- list()
    for (fam_name in eff_families) {
      fam_fn <- get_family_fn(fam_name)
      for (spec in build_family_specs(fam_name, nu_f, tau_f,
                                      sigma_f  = sigma_f,
                                      discrete = eff_discrete)) {
        spec$family_fn <- fam_fn
        candidates[[length(candidates) + 1L]] <- spec
      }
    }

    # Continuous: order breadth-first by shape richness so the fullest nu/tau
    # structure is attempted across all families before any is simplified; a
    # batch effect on nu/tau is thus dropped only when no family supports it.
    # Discrete keeps its per-family order.
    if (!eff_discrete) {
      active     <- function(f) !is.null(f) && length(all.vars(f)) > 0L
      richness   <- vapply(candidates, function(s)
                           active(s$nu.formula) + active(s$tau.formula),
                           numeric(1L))
      candidates <- candidates[order(-richness, seq_along(richness))]
    }

    for (spec in candidates) {
      logger::log_info(paste0("  Trying: ", spec$name))
      eff_sigma_f <- spec$sigma.formula %||% sigma_f
      model <- try_gamlss(mu.formula    = mu_f,
                          sigma.formula = eff_sigma_f,
                          nu.formula    = spec$nu.formula,
                          tau.formula   = spec$tau.formula,
                          data          = model_data,
                          family        = spec$family_fn(),
                          c_crit        = c_crit)
      if (!is.null(model)) {
        logger::log_info(paste0("  Converged: ", spec$name,
                                " | AIC: ", round(AIC(model), 2)))
        final_spec <- spec$name
        break
      }
    }

    if (is.null(model)) {
      hint <- ""
      if (!eff_discrete &&
          all(model_data$y == floor(model_data$y), na.rm = TRUE) &&
          all(model_data$y >= 0, na.rm = TRUE) &&
          mean(model_data$y == 0, na.rm = TRUE) > 0.05)
        hint <- "\nHint: the feature looks like non-negative integer count data. Try re-running with --discrete TRUE."
      else if (!eff_discrete && length(eff_families) <= 3)
        hint <- "\nHint: try passing --feature_families from the diagnose stage, or expand --family_order."
      stop("All model specifications failed to converge for '", feature_name, "'.", hint)
    }

    fit_time        <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))
    model$call$data <- model_data

    saveRDS(model, paths$model)
    write_model_equation(
      path          = paths$summary,
      feature_name  = feature_name,
      final_spec    = final_spec,
      model         = model,
      mu_f          = mu_f,
      sigma_f       = eff_sigma_f,
      nu_f          = model$nu.formula  %||% as.formula("~ 1"),
      tau_f         = model$tau.formula %||% as.formula("~ 1"),
      n_obs         = n_obs,
      n_batches     = n_batches,
      fit_time      = fit_time,
      log_transform = log_transform_used,
      discrete      = eff_discrete
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

    mu_original <- if (!eff_discrete && log_transform_used)
      exp(mu_pred) - 1
    else
      mu_pred

    pred_df <- data.frame(
      age            = model_data$age,
      sex            = model_data$sex,
      observed_value = if (log_transform_used) exp(model_data$y) - 1 else model_data$y,
      mu             = mu_original,
      sigma          = sigma_pred,
      nu             = nu_pred,
      tau            = tau_pred,
      discrete       = eff_discrete,
      log_transform  = log_transform_used
    )
    pred_df[[id_var]]    <- model_data[[id_var]]
    pred_df[[batch_var]] <- model_data[[batch_var]]
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
      discrete          = eff_discrete,
      longitudinal_mode = longitudinal,
      log_transform     = log_transform_used
    ), paths$metrics, row.names = FALSE)

    tryCatch({
      pdf(paths$diagnostics, width = 10, height = 8)
      plot(model)
      dev.off()
    }, error = function(e) {
      logger::log_info(paste0("  Diagnostics plot failed: ", e$message))
      if (dev.cur() > 1) dev.off()
    })

    elapsed <- as.numeric(difftime(Sys.time(), feature_start, units = "secs"))
    write.csv(data.frame(feature        = feature_name,
                         status         = "success",
                         distribution   = final_spec,
                         discrete       = eff_discrete,
                         start_time     = as.character(feature_start),
                         end_time       = as.character(Sys.time()),
                         time_secs      = elapsed,
                         n_observations = n_obs,
                         n_batches      = n_batches),
              paths$timing, row.names = FALSE)

    list(status         = "success",
         feature        = feature_name,
         distribution   = final_spec,
         discrete       = eff_discrete,
         n_observations = n_obs,
         n_batches      = n_batches,
         processing_time = elapsed)

  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), feature_start, units = "secs"))
    logger::log_error(paste0("Error fitting ", feature_name, ": ", e$message))
    write.csv(data.frame(feature       = feature_name,
                         status        = "error",
                         error_message = e$message,
                         time_secs     = elapsed),
              paths$timing, row.names = FALSE)
    list(status = "error", feature = feature_name,
         error  = e$message, processing_time = elapsed)
  })
}

# ---------------------------------------------------------------------------
# Rebuild model_summary.csv from every feature_<name>/<name>_timing.csv so it
# reflects all features ever fit here, not just the current run's subset.
# Timing files are read by column name, so the differing per-status schemas
# (success / error / skipped) are tolerated; absent columns become NA.
# ---------------------------------------------------------------------------

rebuild_model_summary <- function(output_dir) {
  cols <- c("feature", "status", "distribution", "discrete",
           "n_observations", "n_batches", "error_message", "time_secs")

  feat_dirs <- list.dirs(output_dir, full.names = TRUE, recursive = FALSE)
  feat_dirs <- feat_dirs[grepl("^feature_", basename(feat_dirs))]

  rows <- lapply(feat_dirs, function(fd) {
    feat_name   <- sub("^feature_", "", basename(fd))
    timing_path <- file.path(fd, paste0(feat_name, "_timing.csv"))
    if (!file.exists(timing_path)) return(NULL)

    t <- tryCatch(read.csv(timing_path, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (is.null(t) || nrow(t) == 0L) return(NULL)
    t <- t[1L, , drop = FALSE]

    row <- data.frame(feature = feat_name, stringsAsFactors = FALSE)
    for (col in setdiff(cols, "feature"))
      row[[col]] <- if (col %in% names(t)) t[[col]] else NA
    row[, cols, drop = FALSE]
  })
  rows <- Filter(Negate(is.null), rows)

  if (length(rows) == 0L) {
    empty <- as.data.frame(setNames(lapply(cols, function(x) character(0)), cols),
                           stringsAsFactors = FALSE)
    return(empty)
  }

  summary_df <- do.call(rbind, rows)
  summary_df[order(summary_df$feature), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Fit all features.
# ---------------------------------------------------------------------------

run_gamlss_harmonisation <- function(data, features, output_dir, formula_terms,
                                     batch_var          = "batch",
                                     id_var             = "id",
                                     longitudinal       = FALSE,
                                     log_transform      = FALSE,
                                     discrete           = FALSE,
                                     family_order       = c("SHASH", "GG", "NO"),
                                     feature_recs       = NULL,
                                     c_crit             = 0.001,
                                     n_cores            = 1) {
  overall_start <- Sys.time()
  logger::log_info(paste0("Fitting ", length(features), " features",
                          if (discrete) " [discrete mode]" else ""))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(feature_recs))
    logger::log_info(paste0("Per-feature recommendations loaded for ",
                            length(feature_recs), " features"))

  fit_one <- function(feat) {
    fit_gamlss_for_feature(
      data          = data,
      feature_name  = feat,
      model_dir     = file.path(output_dir, paste0("feature_", feat)),
      formula_terms = formula_terms,
      batch_var     = batch_var,
      id_var        = id_var,
      longitudinal  = longitudinal,
      log_transform = log_transform,
      discrete      = discrete,
      family_order  = family_order,
      feature_rec   = if (!is.null(feature_recs)) feature_recs[[feat]] else NULL,
      c_crit        = c_crit
    )
  }

  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, {
      library(gamlss)
      library(gamlss.dist)
      library(dplyr)
      library(logger)
      logger::log_formatter(logger::formatter_paste)
    })
    parallel::clusterExport(cl,
      c("data", "formula_terms", "batch_var", "id_var", "longitudinal",
        "log_transform", "discrete", "family_order", "feature_recs",
        "output_dir", "fit_one", "c_crit"),
      envir = environment()
    )
    parallel::clusterExport(cl,
      c("fit_gamlss_for_feature", "try_gamlss", "build_family_specs",
        "formula_has_random", "get_family_fn",
        "resolve_feature_formula_terms",
        "write_model_equation", "family_description", "param_description",
        "get_n_params", "extract_covariate_cols",
        "terms_to_formula", "substitute_placeholders", "format_time",
        ".SMOOTHER_FNS", "%||%"),
      envir = globalenv()
    )
    results <- pbapply::pblapply(features, function(feat) {
      tryCatch(fit_one(feat),
               error = function(e)
                 list(status = "error", feature = feat,
                      error = e$message, processing_time = 0))
    }, cl = cl)
    parallel::stopCluster(cl)
  } else {
    results <- vector("list", length(features))
    timings <- data.frame()
    for (i in seq_along(features)) {
      logger::log_info(paste0("Feature ", i, "/", length(features),
                              ": ", features[i]))
      results[[i]] <- fit_one(features[i])
      timings <- rbind(timings, data.frame(
        feature       = features[i],
        status        = results[[i]]$status,
        distribution  = if (results[[i]]$status == "success")
          results[[i]]$distribution else NA_character_,
        discrete      = if (results[[i]]$status == "success")
          results[[i]]$discrete else NA,
        error_message = if (results[[i]]$status == "error")
          results[[i]]$error else NA_character_,
        time_seconds  = results[[i]]$processing_time
      ))
      write.csv(timings, file.path(output_dir, "feature_timings.csv"),
                row.names = FALSE)
      elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
      logger::log_info(paste0("  ", round(100 * i / length(features), 1),
                              "% | est. remaining: ",
                              format_time(elapsed / i * (length(features) - i))))
    }
  }

  statuses  <- vapply(results, `[[`, character(1L), "status")
  successes <- sum(statuses == "success")
  failures  <- sum(statuses == "error")
  skipped   <- sum(statuses != "success" & statuses != "error")

  # This run's results only, for the "Failed features" log; the on-disk
  # model_summary.csv is rebuilt from all feature dirs (rebuild_model_summary()).
  run_summary_df <- do.call(rbind, lapply(results, function(r) {
    data.frame(
      feature        = r$feature          %||% NA_character_,
      status         = r$status           %||% NA_character_,
      error_message  = r$error            %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))

  if (failures > 0L) {
    failed_df <- run_summary_df[!is.na(run_summary_df$status) &
                                  run_summary_df$status == "error", , drop = FALSE]
    logger::log_info("Failed features:")
    for (i in seq_len(nrow(failed_df))) {
      logger::log_info(paste0("  ", failed_df$feature[i], ": ",
                              failed_df$error_message[i]))
    }
  }

  summary_df   <- rebuild_model_summary(output_dir)
  summary_path <- file.path(output_dir, "model_summary.csv")
  write.csv(summary_df, summary_path, row.names = FALSE)
  logger::log_info(paste0("Model summary: ", summary_path,
                          " (", nrow(summary_df), " feature(s) on disk)"))

  total_time <- format_time(as.numeric(difftime(Sys.time(), overall_start,
                                                units = "secs")))
  list(results    = results,
       successes  = successes,
       failures   = failures,
       skipped    = skipped,
       total_time = total_time)
}