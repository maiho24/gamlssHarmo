# R/diagnose.R
#
# Pre-fitting diagnostic functions for gamlssHarmo.
#
# CONTINUOUS mode (discrete = FALSE, default):
#   1. Residualise for age (natural cubic spline), sex, and batch (fixed
#      effects) via OLS so the moments reflect the within-batch conditional
#      shape
#   2. Compute skewness and excess kurtosis on residuals
#   3. Permutation tests for between-batch variance in skewness and kurtosis
#   4. Tier assignment and Cullen-Frey family recommendation
#
# DISCRETE / COUNT mode (discrete = TRUE):
#   1. Cast feature to integer before all computation
#   2. Fit a Poisson GLM (natural cubic spline age + sex) for covariate-adjusted means
#   3. Zero-inflation test: one-sided z-test (Hall & Berenhaut 2002)
#   4. Overdispersion test: Cameron-Trivedi (1990) auxiliary regression
#   5. Underdispersion test: Cameron-Trivedi lower-tail variant
#   6. Pearson dispersion index: NBI-class (D~1-3) vs Sichel-class (D>3)
#   7. Tail heaviness: P95/median ratio; values > 5 suggest Sichel over NBI
#   8. Permutation tests for between-batch variation in zero rate and
#      dispersion ratio
#   9. Count tier and family recommendation:
#        No ZI,  No OD          -> NBI, ZANBI, ZINBI, PO, ZIP
#        No ZI, mild/mod OD     -> ZANBI, ZINBI, NBI, PIG, NBII, PO
#        No ZI, severe/heavy OD -> ZASICHEL, ZISICHEL, SICHEL, ZANBI, ZINBI, PIG, NBI, NBII, PO
#       Yes ZI,  No OD          -> ZANBI, ZINBI, ZIP, NBI, PO
#       Yes ZI, mild/mod OD     -> ZANBI, ZINBI, NBI, PIG, ZIP, PO
#       Yes ZI, severe/heavy OD -> ZASICHEL, ZISICHEL, ZANBI, ZINBI, SICHEL, PIG, NBI, ZIP, PO
#       Underdispersion         -> DPO, PO, NBI
#
#   Hurdle (ZANBI, ZASICHEL) listed before mixture (ZINBI, ZISICHEL).
#
# Output:
#   diagnostics_summary.csv      full statistics
#   feature_recommendations.csv  trimmed recommendations for 01_fit.R

# Natural-cubic-spline df for age adjustment (residualisation and count GLMs).
# A spline rather than poly(age, 2): a global quadratic misfits the pooled
# lifespan and leaks residual age curvature into the skew/kurtosis batch tests.
AGE_SPLINE_DF <- 4L

# ---------------------------------------------------------------------------
# Moment helpers (Type 1, /n, Cullen-Frey convention)
# ---------------------------------------------------------------------------

compute_skewness <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (s < .Machine$double.eps) return(NA_real_)
  (sum((x - m)^3L) / n) / s^3L
}

compute_excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 4L) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (s < .Machine$double.eps) return(NA_real_)
  (sum((x - m)^4L) / n) / s^4L - 3.0
}

# ---------------------------------------------------------------------------
# Per-feature RNG seed.
#
# Derived from the feature name rather than its position so that re-running a
# subset (--one_feature) reproduces the p-values that feature had in a full run.
# Collisions between features are harmless: the analyses are independent.
# ---------------------------------------------------------------------------

feature_seed <- function(base_seed, feature_name) {
  chars <- utf8ToInt(feature_name)
  as.integer((base_seed + sum(chars * seq_along(chars))) %% 2147483647)
}

# ---------------------------------------------------------------------------
# Age x sex strata for the permutation null.
# ---------------------------------------------------------------------------

build_strata <- function(age, sex, n_bands = 10L) {
  br <- unique(quantile(age, seq(0, 1, length.out = n_bands + 1L), na.rm = TRUE))
  if (length(br) < 2L) return(rep("all", length(age)))
  paste0("b", cut(age, br, include.lowest = TRUE, labels = FALSE),
         "_", as.character(sex))
}

# ---------------------------------------------------------------------------
# Permutation test for between-batch variation in a distributional moment.
#
# Batch is typically confounded with age in a multi-cohort design, so labels are
# shuffled only within a stratum; strata holding a single batch carry no batch
# contrast and are dropped, which restricts the test to shared age support.
# Batch is constant within subject, so the exchangeable unit is the subject.
# Passing strata = NULL and subjects = NULL recovers the unconditional row-level
# test. The null range is accumulated alongside the null variance so the
# effect-size guard can be calibrated rather than fixed.
# ---------------------------------------------------------------------------

permutation_test_batch_moment <- function(x, batches, moment_fn,
                                           n_perm      = 499L,
                                           min_batch_n = 50L,
                                           subjects    = NULL,
                                           strata      = NULL) {
  ok <- is.finite(x) & !is.na(batches)
  if (!is.null(strata))   ok <- ok & !is.na(strata)
  if (!is.null(subjects)) ok <- ok & !is.na(subjects)

  x        <- x[ok]
  batches  <- as.character(batches[ok])
  strata   <- if (is.null(strata))   rep("all", length(x))      else as.character(strata[ok])
  subjects <- if (is.null(subjects)) as.character(seq_along(x)) else as.character(subjects[ok])

  na_out <- function(n_used, frac, vals = NULL)
    list(observed_variance = NA_real_,
         p_value           = NA_real_,
         p_range           = NA_real_,
         n_batches_used    = n_used,
         per_batch_values  = vals,
         per_batch_range   = NA_real_,
         frac_common       = frac)

  su <- data.frame(subject = subjects, batch = batches, stratum = strata,
                   stringsAsFactors = FALSE)
  su <- su[!duplicated(su$subject), ]

  n_bat  <- tapply(su$batch, su$stratum, function(z) length(unique(z)))
  su     <- su[su$stratum %in% names(n_bat)[n_bat >= 2L], ]
  row_ok <- subjects %in% su$subject
  frac   <- mean(row_ok)

  x <- x[row_ok]; batches <- batches[row_ok]; subjects <- subjects[row_ok]
  if (length(x) < 2L) return(na_out(0L, frac))

  tab  <- table(batches)
  keep <- names(tab)[tab >= min_batch_n]
  if (length(keep) < 2L) return(na_out(length(keep), frac))

  mask <- batches %in% keep
  x_s  <- x[mask]
  b_s  <- batches[mask]
  su   <- su[su$batch %in% keep, ]
  idx  <- match(subjects[mask], su$subject)

  per_batch <- vapply(keep, function(b)
    moment_fn(x_s[b_s == b]), numeric(1L))
  per_batch <- per_batch[is.finite(per_batch)]

  if (length(per_batch) < 2L) return(na_out(length(per_batch), frac, per_batch))

  obs_var   <- var(per_batch)
  obs_range <- diff(range(per_batch))

  st_idx <- split(seq_len(nrow(su)), su$stratum)
  perm   <- replicate(n_perm, {
    b_new <- su$batch
    for (ii in st_idx) if (length(ii) > 1L) b_new[ii] <- sample(su$batch[ii])
    b_perm <- b_new[idx]
    vals   <- vapply(keep, function(b)
      moment_fn(x_s[b_perm == b]), numeric(1L))
    vals   <- vals[is.finite(vals)]
    if (length(vals) < 2L) c(NA_real_, NA_real_)
    else                   c(var(vals), diff(range(vals)))
  })

  pv <- perm[1L, ][is.finite(perm[1L, ])]
  pr <- perm[2L, ][is.finite(perm[2L, ])]

  list(observed_variance = obs_var,
       p_value           = if (length(pv) > 0L) (1 + sum(pv >= obs_var))   / (1 + length(pv)) else NA_real_,
       p_range           = if (length(pr) > 0L) (1 + sum(pr >= obs_range)) / (1 + length(pr)) else NA_real_,
       n_batches_used    = length(per_batch),
       per_batch_values  = per_batch,
       per_batch_range   = obs_range,
       frac_common       = frac)
}

# ---------------------------------------------------------------------------
# Residualise a continuous feature for age (natural cubic spline) and sex.
#
# A second-stage model of log(e^2) on the same design supplies studentised
# residuals and a direct test of whether sigma depends on age or sex. Raw
# residuals from a homoscedastic fit are a scale mixture, which carries excess
# kurtosis 3 * CV^2(sigma^2) even when every conditional distribution is
# Gaussian; studentising removes that from the shape statistics.
# ---------------------------------------------------------------------------

residualise_feature <- function(data, feature_name, batch_var = NULL) {
  df <- data.frame(
    y   = data[[feature_name]],
    age = data[["age"]],
    sex = data[["sex"]]
  )
  has_batch <- !is.null(batch_var) && batch_var %in% names(data)
  if (has_batch) df$batch <- as.factor(data[[batch_var]])

  df <- df[is.finite(df$y) & is.finite(df$age) & !is.na(df$sex), ]

  plain <- function(e)
    list(residuals = e, studentised = e, sigma_p = NA_real_, sigma_cv2 = NA_real_)

  if (nrow(df) < 10L) {
    warning("Too few observations to residualise '", feature_name, "'")
    return(plain(df$y - mean(df$y, na.rm = TRUE)))
  }

  use_batch <- has_batch && !anyNA(df$batch) &&
               nlevels(droplevels(df$batch)) > 1L
  form <- if (use_batch) y ~ splines::ns(age, df = AGE_SPLINE_DF) + sex + batch
          else           y ~ splines::ns(age, df = AGE_SPLINE_DF) + sex

  fit <- tryCatch(lm(form, data = df), error = function(e) NULL)
  if (is.null(fit)) return(plain(df$y - mean(df$y)))
  e <- residuals(fit)

  df$lv  <- log(pmax(e^2, .Machine$double.eps))
  v_full <- tryCatch(lm(update(form, lv ~ .), data = df), error = function(e) NULL)
  v_null <- tryCatch(if (use_batch) lm(lv ~ batch, data = df) else lm(lv ~ 1, data = df),
                     error = function(e) NULL)
  if (is.null(v_full) || is.null(v_null)) return(plain(e))

  sigma2 <- exp(fitted(v_full))
  list(residuals   = e,
       studentised = e / sqrt(sigma2),
       sigma_p     = tryCatch(anova(v_null, v_full)[2L, "Pr(>F)"],
                              error = function(e) NA_real_),
       sigma_cv2   = var(sigma2) / mean(sigma2)^2)
}

# ---------------------------------------------------------------------------
# Zero-inflation test (Hall & Berenhaut 2002).
# ---------------------------------------------------------------------------

test_zero_inflation <- function(x, age, sex) {
  df <- data.frame(y = as.integer(x), age = age, sex = sex)
  df <- df[df$y >= 0L & is.finite(df$age) & !is.na(df$sex), ]
  n  <- nrow(df)

  if (n < 20L)
    return(list(statistic     = NA_real_, p_value       = NA_real_,
                obs_zero_rate = NA_real_, exp_zero_rate = NA_real_))

  pois_fit <- tryCatch(
    glm(y ~ splines::ns(age, df = AGE_SPLINE_DF) + sex, data = df, family = poisson(link = "log")),
    error = function(e) NULL
  )
  mu_hat <- if (!is.null(pois_fit)) fitted(pois_fit) else rep(mean(df$y), n)

  obs_zeros  <- sum(df$y == 0L)
  p_zero_vec <- exp(-mu_hat)
  exp_zeros  <- sum(p_zero_vec)
  var_zeros  <- sum(p_zero_vec * (1 - p_zero_vec))

  if (var_zeros < .Machine$double.eps)
    return(list(statistic     = NA_real_, p_value       = NA_real_,
                obs_zero_rate = obs_zeros / n, exp_zero_rate = exp_zeros / n))

  z       <- (obs_zeros - exp_zeros) / sqrt(var_zeros)
  p_value <- pnorm(z, lower.tail = FALSE)

  list(statistic     = z,
       p_value       = p_value,
       obs_zero_rate = obs_zeros / n,
       exp_zero_rate = exp_zeros / n)
}

# ---------------------------------------------------------------------------
# Overdispersion test (Cameron-Trivedi 1990, upper-tail).
# ---------------------------------------------------------------------------

test_overdispersion <- function(x, age, sex) {
  df <- data.frame(y = as.integer(x), age = age, sex = sex)
  df <- df[df$y >= 0L & is.finite(df$age) & !is.na(df$sex), ]
  n  <- nrow(df)

  ybar       <- mean(df$y)
  disp_ratio <- if (ybar > 0) var(df$y) / ybar else NA_real_

  if (n < 20L)
    return(list(statistic = NA_real_, p_value = NA_real_,
                dispersion_ratio = disp_ratio))

  pois_fit <- tryCatch(
    glm(y ~ splines::ns(age, df = AGE_SPLINE_DF) + sex, data = df, family = poisson(link = "log")),
    error = function(e) NULL
  )
  mu_hat <- if (!is.null(pois_fit)) fitted(pois_fit) else rep(ybar, n)

  z       <- ((df$y - mu_hat)^2 - df$y) / mu_hat
  aux_fit <- tryCatch(lm(z ~ mu_hat - 1), error = function(e) NULL)

  if (is.null(aux_fit))
    return(list(statistic = NA_real_, p_value = NA_real_,
                dispersion_ratio = disp_ratio))

  s       <- summary(aux_fit)$coefficients
  t_stat  <- s[1L, "t value"]
  p_value <- pt(t_stat, df = n - 1L, lower.tail = FALSE)

  list(statistic        = t_stat,
       p_value          = p_value,
       dispersion_ratio = disp_ratio)
}

# ---------------------------------------------------------------------------
# Underdispersion test (Cameron-Trivedi lower-tail variant).
# ---------------------------------------------------------------------------

test_underdispersion <- function(x, age, sex) {
  df <- data.frame(y = as.integer(x), age = age, sex = sex)
  df <- df[df$y >= 0L & is.finite(df$age) & !is.na(df$sex), ]
  n  <- nrow(df)

  if (n < 20L)
    return(list(statistic = NA_real_, p_value = NA_real_))

  pois_fit <- tryCatch(
    glm(y ~ splines::ns(age, df = AGE_SPLINE_DF) + sex, data = df, family = poisson(link = "log")),
    error = function(e) NULL
  )
  if (is.null(pois_fit))
    return(list(statistic = NA_real_, p_value = NA_real_))

  mu_hat  <- fitted(pois_fit)
  z_vals  <- ((df$y - mu_hat)^2 - df$y) / mu_hat
  aux_fit <- tryCatch(lm(z_vals ~ mu_hat - 1), error = function(e) NULL)

  if (is.null(aux_fit))
    return(list(statistic = NA_real_, p_value = NA_real_))

  alpha   <- coef(aux_fit)["mu_hat"]
  se      <- sqrt(diag(vcov(aux_fit)))["mu_hat"]

  if (is.na(se) || se < .Machine$double.eps)
    return(list(statistic = NA_real_, p_value = NA_real_))

  t_stat  <- alpha / se
  p_value <- pt(t_stat, df = n - 1L, lower.tail = TRUE)
  list(statistic = t_stat, p_value = p_value)
}

# ---------------------------------------------------------------------------
# Pearson dispersion index (covariate-adjusted).
# ---------------------------------------------------------------------------

compute_dispersion_index <- function(x, age, sex) {
  df <- data.frame(y = as.integer(x), age = age, sex = sex)
  df <- df[df$y >= 0L & is.finite(df$age) & !is.na(df$sex), ]
  if (nrow(df) < 20L) return(NA_real_)

  pois_fit <- tryCatch(
    glm(y ~ splines::ns(age, df = AGE_SPLINE_DF) + sex, data = df, family = poisson(link = "log")),
    error = function(e) NULL
  )
  mu_hat <- if (!is.null(pois_fit)) fitted(pois_fit) else rep(mean(df$y), nrow(df))
  mean(((df$y - mu_hat)^2) / mu_hat, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Tail heaviness: P95/median ratio on integer counts.
# ---------------------------------------------------------------------------

compute_tail_heaviness <- function(x) {
  x   <- as.integer(x[is.finite(x) & x >= 0L])
  if (length(x) < 20L) return(NA_real_)
  med <- median(x)
  if (med == 0L) return(NA_real_)
  quantile(x, 0.95) / med
}

# ---------------------------------------------------------------------------
# Continuous family suggestion (Cullen-Frey).
# ---------------------------------------------------------------------------

suggest_families <- function(skewness, excess_kurtosis,
                              is_positive_only, has_near_zero,
                              skew_batch_effect = FALSE,
                              kurt_batch_effect = FALSE) {
  s <- skewness
  k <- excess_kurtosis

  if (is.na(s) || is.na(k)) {
    base <- c("SHASH", "JSU", "GG", "NO")

  } else if (is_positive_only) {
    if (abs(s) < 0.5 && abs(k) < 1.0) {
      base <- c("GG", "BCCG", "SHASH", "NO")
    } else if (s > 0.5) {
      base <- if (k >= 1.5 * s^2) c("GG", "LOGNO", "BCCG", "SHASH")
              else                 c("GG", "BCCG", "LOGNO", "SHASH")
    } else if (s < -0.5) {
      base <- c("BCT", "BCCG", "GG", "SHASH")
    } else {
      base <- c("GG", "SHASH", "BCCG", "LOGNO")
    }

  } else if (has_near_zero) {
    if (abs(s) < 0.5 && abs(k) < 1.0)  base <- c("SHASH", "NO", "TF")
    else if (abs(s) < 0.5 && k > 1.0)  base <- c("TF", "SHASH", "NO")
    else if (abs(s) < 0.5 && k <= -1.0) base <- c("NO", "PE2")
    else if (abs(s) > 0.5)             base <- c("SHASH", "JSU", "SN1", "NO")
    else                                base <- c("SHASH", "JSU", "TF", "NO")

  } else {
    if (abs(s) < 0.5 && abs(k) < 1.0)  base <- c("NO", "TF", "SHASH")
    else if (abs(s) < 0.5 && k > 1.0)  base <- c("TF", "SHASH", "NO")
    else if (abs(s) < 0.5 && k <= -1.0) base <- c("NO", "PE2")
    else if (abs(s) > 0.5)             base <- c("SHASH", "JSU", "SN1", "NO")
    else                                base <- c("SHASH", "JSU", "TF", "NO")
  }

  if (isTRUE(kurt_batch_effect)) {
    lead <- if (is_positive_only) c("BCT", "BCPE", "SHASH") else c("SHASH", "JSU")
    base <- unique(c(lead, base))
  } else if (isTRUE(skew_batch_effect)) {
    lead <- if (is_positive_only) c("BCT", "BCCG", "GG", "SHASH")
            else                  c("SHASH", "JSU", "SN1", "TF")
    base <- unique(c(lead, base))
  }

  base
}

# ---------------------------------------------------------------------------
# Count family suggestion.
# ---------------------------------------------------------------------------

suggest_count_families <- function(zero_inflated, overdispersed,
                                    underdispersed   = FALSE,
                                    dispersion_index = NA_real_,
                                    tail_heaviness   = NA_real_) {

  if (isTRUE(underdispersed))
    return(c("DPO", "PO", "NBI"))

  severe_od  <- isTRUE(overdispersed) &&
                !is.na(dispersion_index) && dispersion_index > 3.0
  heavy_tail <- !is.na(tail_heaviness) && tail_heaviness > 5.0

  if (!zero_inflated && !overdispersed)
    return(c("NBI", "ZANBI", "ZINBI", "PO", "ZIP"))

  if (!zero_inflated && overdispersed) {
    if (severe_od || heavy_tail)
      return(c("ZASICHEL", "ZISICHEL", "SICHEL", "ZANBI", "ZINBI",
               "PIG", "NBI", "NBII", "PO"))
    return(c("ZANBI", "ZINBI", "NBI", "PIG", "NBII", "PO"))
  }

  if (zero_inflated && !overdispersed)
    return(c("ZANBI", "ZINBI", "ZIP", "NBI", "PO"))

  if (severe_od || heavy_tail)
    return(c("ZASICHEL", "ZISICHEL", "ZANBI", "ZINBI", "SICHEL",
             "PIG", "NBI", "ZIP", "PO"))

  c("ZANBI", "ZINBI", "NBI", "PIG", "ZIP", "PO")
}

# ---------------------------------------------------------------------------
# Continuous tier assignment.
# ---------------------------------------------------------------------------

assign_tier <- function(skewness, excess_kurtosis,
                         batch_skew_pval,  batch_skew_range,
                         batch_kurt_pval,  batch_kurt_range,
                         thresholds,
                         batch_skew_prange = NA_real_,
                         batch_kurt_prange = NA_real_) {

  near_gaussian <- !is.na(skewness) &&
                   abs(skewness) < thresholds$skew &&
                   !is.na(excess_kurtosis) &&
                   abs(excess_kurtosis) < thresholds$kurt

  # The permuted range gives a null the fixed min_*_range cannot: it adapts to
  # batch sizes and to the tail weight of the feature. The constant is retained
  # as an absolute floor.
  gate <- function(pval, prange, rng, min_rng)
    !is.na(pval)   && pval   < thresholds$pval &&
    !is.na(prange) && prange < thresholds$pval &&
    !is.na(rng)    && rng   >= min_rng

  skew_batch <- gate(batch_skew_pval, batch_skew_prange,
                     batch_skew_range, thresholds$min_skew_range)

  kurt_batch <- gate(batch_kurt_pval, batch_kurt_prange,
                     batch_kurt_range, thresholds$min_kurt_range)

  tier <- if (skew_batch || kurt_batch) 3L else if (!near_gaussian) 2L else 1L

  list(tier              = tier,
       near_gaussian     = near_gaussian,
       skew_batch_effect = skew_batch,
       kurt_batch_effect = kurt_batch)
}

# ---------------------------------------------------------------------------
# Count tier assignment.
#
# Tier 1 -- No ZI, no OD, stable batch dispersion.
#           mu + sigma (non-PO) get random(batch); nu/tau intercept-only.
# Tier 2 -- Global OD/ZI, or dispersion varies by batch.
#           Same formula as Tier 1; signals more complex count structure.
# Tier 3 -- Zero rate varies by batch. random(batch) on the ZI/hurdle zero
#           parameter (nu for 3-param families, tau for 4-param Sichel).
# Tier 4 -- Severe OD or heavy tail.   As Tier 2/3; Sichel-class families.
# ---------------------------------------------------------------------------

assign_count_tier <- function(zero_inflated, overdispersed,
                               batch_zerorate_pval,   batch_zerorate_range,
                               batch_dispersion_pval, batch_dispersion_range,
                               thresholds,
                               dispersion_index = NA_real_,
                               tail_heaviness   = NA_real_,
                               batch_zerorate_prange   = NA_real_,
                               batch_dispersion_prange = NA_real_) {

  gate <- function(pval, prange, rng, min_rng)
    !is.na(pval)   && pval   < thresholds$pval &&
    !is.na(prange) && prange < thresholds$pval &&
    !is.na(rng)    && rng   >= min_rng

  zerorate_batch <- gate(batch_zerorate_pval, batch_zerorate_prange,
                         batch_zerorate_range, thresholds$min_zerorate_range)

  dispersion_batch <- gate(batch_dispersion_pval, batch_dispersion_prange,
                           batch_dispersion_range, thresholds$min_dispersion_range)

  severe_od  <- isTRUE(overdispersed) &&
                !is.na(dispersion_index) && dispersion_index > 3.0
  heavy_tail <- !is.na(tail_heaviness) && tail_heaviness > 5.0

  tier <- if (zerorate_batch)                3L
          else if (zero_inflated  ||
                   overdispersed  ||
                   dispersion_batch)         2L
          else                               1L

  if (tier >= 2L && (severe_od || heavy_tail)) tier <- 4L

  list(tier             = tier,
       zero_inflated    = zero_inflated,
       overdispersed    = overdispersed,
       zerorate_batch   = zerorate_batch,
       dispersion_batch = dispersion_batch,
       severe_od        = severe_od,
       heavy_tail       = heavy_tail)
}

# ---------------------------------------------------------------------------
# Continuous formula recommendations.
# ---------------------------------------------------------------------------

build_recommended_formulas <- function(tier_result, longitudinal) {

  mu_terms <- if (longitudinal)
    c("pb(age)", "sex", "random({batch})", "random({id})")
  else
    c("pb(age)", "sex", "random({batch})")

  sigma_terms <- if (tier_result$tier >= 2L)
    c("pb(age)", "sex", "random({batch})")
  else
    c("random({batch})")

  nu_terms  <- if (tier_result$tier >= 3L && tier_result$skew_batch_effect)
    "random({batch})" else NULL

  tau_terms <- if (tier_result$tier >= 3L && tier_result$kurt_batch_effect)
    "random({batch})" else NULL

  list(mu    = mu_terms,
       sigma = sigma_terms,
       nu    = nu_terms,
       tau   = tau_terms)
}

# ---------------------------------------------------------------------------
# Count formula recommendations.
#
# sigma gets random(batch) for all non-PO families (age/sex on overdispersion
# rarely converges). The structural-zero probability gets random(batch) at
# Tier 3 (zero rate varies by batch): on nu for 3-param ZI/hurdle families
# (ZINBI, ZANBI, ZIPIG), on tau for 4-param ones (ZISICHEL, ZASICHEL, whose nu
# is the IG shape). SICHEL nu is intercept-only always.
# ---------------------------------------------------------------------------

build_count_recommended_formulas <- function(tier_result, longitudinal,
                                              primary_family) {
  mu_terms <- if (longitudinal)
    c("pb(age)", "sex", "random({batch})", "random({id})")
  else
    c("pb(age)", "sex", "random({batch})")

  sigma_terms <- if (primary_family != "PO")
    c("random({batch})") else NULL

  zi_nu_families  <- c("ZINBI", "ZANBI", "ZIPIG")     # 3-param: nu = zero prob
  zi_tau_families <- c("ZISICHEL", "ZASICHEL")        # 4-param: tau = zero prob
  zerorate_batch  <- isTRUE(tier_result$tier >= 3L)

  nu_terms <- if (primary_family %in% zi_nu_families) {
    if (zerorate_batch) c("random({batch})") else character(0L)
  } else if (primary_family %in% c(zi_tau_families, "SICHEL")) {
    character(0L)
  } else {
    NULL
  }

  tau_terms <- if (primary_family %in% zi_tau_families && zerorate_batch)
    c("random({batch})") else NULL

  list(mu    = mu_terms,
       sigma = sigma_terms,
       nu    = nu_terms,
       tau   = tau_terms)
}

# ---------------------------------------------------------------------------
# Read a user-supplied feature recommendations CSV.
# ---------------------------------------------------------------------------

read_feature_recommendations_csv <- function(path) {
  if (!file.exists(path))
    stop("Feature recommendations CSV not found: ", path)

  df <- read.csv(path, stringsAsFactors = FALSE)

  if (!"feature" %in% names(df))
    stop("Feature recommendations CSV must contain a 'feature' column")

  parse_semi <- function(x) {
    if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(NULL)
    trimws(strsplit(x, ";")[[1L]])
  }

  setNames(lapply(seq_len(nrow(df)), function(i) {
    r <- df[i, ]
    list(
      family_order = if ("family_order" %in% names(r)) parse_semi(r$family_order) else NULL,
      nu_terms     = if ("nu_terms"     %in% names(r)) parse_semi(r$nu_terms)     else NULL,
      tau_terms    = if ("tau_terms"    %in% names(r)) parse_semi(r$tau_terms)    else NULL
    )
  }), df$feature)
}

# ---------------------------------------------------------------------------
# Diagnose a single feature.
# ---------------------------------------------------------------------------

diagnose_feature <- function(data, feature_name, batch_var, id_var,
                              longitudinal, discrete, n_perm, min_batch_n,
                              thresholds, user_override = NULL,
                              seed = NULL, stratify = TRUE, n_age_bands = 10L) {
  start <- Sys.time()
  # kind is pinned so a worker's RNG type cannot change the stream.
  if (!is.null(seed)) set.seed(seed, kind = "Mersenne-Twister")
  logger::log_info(paste0("Diagnosing: ", feature_name,
                          if (discrete) " [discrete]" else ""))

  tryCatch({
    required <- c(feature_name, "age", "sex", batch_var)
    missing  <- setdiff(required, names(data))
    if (length(missing) > 0L)
      stop("Missing columns: ", paste(missing, collapse = ", "))

    # is.finite() rather than !is.na(): residualise_feature() drops non-finite
    # values too, and residuals are paired positionally with batch, subject and
    # stratum below.
    df <- data[is.finite(data[[feature_name]]) &
                 is.finite(data[["age"]]) &
                 !is.na(data[["sex"]]), ]

    # Cast to integer before any computation so all downstream checks
    # and tests operate on actual integer values.
    if (discrete)
      df[[feature_name]] <- as.integer(round(df[[feature_name]]))

    rng <- range(df[[feature_name]], na.rm = TRUE)

    has_zeros        <- any(df[[feature_name]] == 0, na.rm = TRUE)
    is_positive_only <- rng[1L] > 0 && !has_zeros
    near_zero        <- !is_positive_only && rng[1L] == 0

    if (has_zeros && !discrete)
      logger::log_info(paste0("  ", feature_name,
                              ": zeros present -- positive-support families excluded"))

    # Counts must be non-negative; continuous features may be real-valued.
    if (discrete) df <- df[df[[feature_name]] >= 0, ]

    n_obs     <- nrow(df)
    n_batches <- length(unique(df[[batch_var]]))

    if (n_batches < 5L)
      logger::log_info(paste0("  Note: only ", n_batches, " batches -- ",
                              "consider fixed batch effects ({batch} instead of random({batch}))"))

    if (n_obs < 30L)
      stop("Insufficient observations after filtering (", n_obs, " < 30)")

    strata_vec <- if (isTRUE(stratify))
      build_strata(df[["age"]], df[["sex"]], n_age_bands) else NULL
    subj_vec   <- if (!is.null(id_var) && id_var %in% names(df))
      df[[id_var]] else NULL

    fmt_terms <- function(x) {
      if (is.null(x))      return(NA_character_)
      if (length(x) == 0L) return("")
      paste(trimws(unlist(x)), collapse = ";")
    }

    # =========================================================================
    # DISCRETE / COUNT PATH
    # =========================================================================
    if (discrete) {

      y   <- df[[feature_name]]
      age <- df[["age"]]
      sex <- df[["sex"]]

      skew <- compute_skewness(y)
      kurt <- compute_excess_kurtosis(y)
      logger::log_info(paste0("  skew=", round(skew, 3L),
                              "  kurt=", round(kurt, 3L)))

      zi_test <- test_zero_inflation(y, age, sex)
      logger::log_info(paste0(
        "  ZI z=",  round(zi_test$statistic    %||% NA, 3L),
        "  p=",     round(zi_test$p_value       %||% NA, 4L),
        "  obs0=",  round(zi_test$obs_zero_rate %||% NA, 3L),
        "  exp0=",  round(zi_test$exp_zero_rate %||% NA, 3L)
      ))

      od_test <- test_overdispersion(y, age, sex)
      logger::log_info(paste0(
        "  OD t=",  round(od_test$statistic        %||% NA, 3L),
        "  p=",     round(od_test$p_value           %||% NA, 4L),
        "  D=",     round(od_test$dispersion_ratio  %||% NA, 3L)
      ))

      ud_test <- test_underdispersion(y, age, sex)
      logger::log_info(paste0(
        "  UD t=",  round(ud_test$statistic %||% NA, 3L),
        "  p=",     round(ud_test$p_value   %||% NA, 4L)
      ))

      disp_index <- compute_dispersion_index(y, age, sex)
      tail_h     <- compute_tail_heaviness(y)
      logger::log_info(paste0(
        "  dispersion_index=", round(disp_index %||% NA, 3L),
        "  tail_heaviness=",   round(tail_h     %||% NA, 3L)
      ))

      zi_sig <- !is.na(zi_test$p_value) && zi_test$p_value < thresholds$pval
      od_sig <- !is.na(od_test$p_value) && od_test$p_value < thresholds$pval
      ud_sig <- !is.na(ud_test$p_value) && ud_test$p_value < thresholds$pval

      zerorate_fn   <- function(x) mean(as.integer(x) == 0L)
      zerorate_perm <- permutation_test_batch_moment(
        y, df[[batch_var]], zerorate_fn, n_perm, min_batch_n,
        subjects = subj_vec, strata = strata_vec
      )
      logger::log_info(paste0(
        "  zero-rate batch perm p=",
        round(zerorate_perm$p_value %||% NA, 4L),
        "  range=",
        round(zerorate_perm$per_batch_range %||% NA, 3L)
      ))

      dispersion_fn <- function(x) {
        xi <- as.integer(x); mu <- mean(xi)
        if (mu > 0) var(xi) / mu else NA_real_
      }
      dispersion_perm <- permutation_test_batch_moment(
        y, df[[batch_var]], dispersion_fn, n_perm, min_batch_n,
        subjects = subj_vec, strata = strata_vec
      )
      logger::log_info(paste0(
        "  dispersion batch perm p=",
        round(dispersion_perm$p_value %||% NA, 4L),
        "  range=",
        round(dispersion_perm$per_batch_range %||% NA, 3L)
      ))

      tier_res <- assign_count_tier(
        zero_inflated          = zi_sig,
        overdispersed          = od_sig,
        batch_zerorate_pval    = zerorate_perm$p_value,
        batch_zerorate_range   = zerorate_perm$per_batch_range,
        batch_dispersion_pval  = dispersion_perm$p_value,
        batch_dispersion_range = dispersion_perm$per_batch_range,
        thresholds             = thresholds,
        dispersion_index       = disp_index,
        tail_heaviness         = tail_h,
        batch_zerorate_prange   = zerorate_perm$p_range,
        batch_dispersion_prange = dispersion_perm$p_range
      )
      logger::log_info(paste0("  tier=",        tier_res$tier,
                              "  ZI=",          tier_res$zero_inflated,
                              "  OD=",          tier_res$overdispersed,
                              "  severe_od=",   tier_res$severe_od,
                              "  heavy_tail=",  tier_res$heavy_tail))

      if (!is.null(user_override) && !is.null(user_override$family_order)) {
        rec_families <- unlist(user_override$family_order)
        logger::log_info(paste0("  family_order override: ",
                                paste(rec_families, collapse = " -> ")))
      } else {
        rec_families <- suggest_count_families(
          zero_inflated    = tier_res$zero_inflated,
          overdispersed    = tier_res$overdispersed,
          underdispersed   = ud_sig,
          dispersion_index = disp_index,
          tail_heaviness   = tail_h
        )
      }

      primary_family <- rec_families[1L]
      rec_formulas   <- build_count_recommended_formulas(tier_res, longitudinal,
                                                          primary_family)
      if (!is.null(user_override)) {
        if (!is.null(user_override$nu_terms))
          rec_formulas$nu  <- user_override$nu_terms
        if (!is.null(user_override$tau_terms))
          rec_formulas$tau <- user_override$tau_terms
      }

      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      logger::log_info(paste0("  done in ", format_time(elapsed)))

      return(list(
        status                  = "success",
        feature                 = feature_name,
        discrete                = TRUE,
        n_obs                   = n_obs,
        n_batches               = n_batches,
        is_positive_only        = is_positive_only,
        skewness                = round(skew, 4L),
        excess_kurtosis         = round(kurt, 4L),
        zi_statistic            = round(zi_test$statistic        %||% NA, 4L),
        zi_pvalue               = round(zi_test$p_value          %||% NA, 4L),
        obs_zero_rate           = round(zi_test$obs_zero_rate    %||% NA, 4L),
        exp_zero_rate           = round(zi_test$exp_zero_rate    %||% NA, 4L),
        od_statistic            = round(od_test$statistic        %||% NA, 4L),
        od_pvalue               = round(od_test$p_value          %||% NA, 4L),
        dispersion_ratio        = round(od_test$dispersion_ratio %||% NA, 4L),
        ud_statistic            = round(ud_test$statistic        %||% NA, 4L),
        ud_pvalue               = round(ud_test$p_value          %||% NA, 4L),
        dispersion_index        = round(disp_index               %||% NA, 4L),
        tail_heaviness          = round(tail_h                   %||% NA, 4L),
        batch_zerorate_pvalue   = round(zerorate_perm$p_value           %||% NA, 4L),
        batch_zerorate_prange   = round(zerorate_perm$p_range           %||% NA, 4L),
        batch_zerorate_range    = round(zerorate_perm$per_batch_range   %||% NA, 4L),
        batch_zerorate_var      = round(zerorate_perm$observed_variance %||% NA, 4L),
        batch_zerorate_nbatches = zerorate_perm$n_batches_used   %||% NA_integer_,
        batch_disp_pvalue       = round(dispersion_perm$p_value           %||% NA, 4L),
        batch_disp_prange       = round(dispersion_perm$p_range           %||% NA, 4L),
        batch_disp_range        = round(dispersion_perm$per_batch_range   %||% NA, 4L),
        batch_disp_var          = round(dispersion_perm$observed_variance %||% NA, 4L),
        batch_disp_nbatches     = dispersion_perm$n_batches_used %||% NA_integer_,
        frac_common             = round(zerorate_perm$frac_common %||% NA, 4L),
        seed                    = seed %||% NA_integer_,
        tier                    = tier_res$tier,
        zero_inflated           = tier_res$zero_inflated,
        overdispersed           = tier_res$overdispersed,
        underdispersed          = ud_sig,
        severe_od               = tier_res$severe_od,
        heavy_tail              = tier_res$heavy_tail,
        zerorate_batch_effect   = tier_res$zerorate_batch,
        dispersion_batch_effect = tier_res$dispersion_batch,
        family_override         = !is.null(user_override) &&
                                   !is.null(user_override$family_order),
        family_order            = paste(rec_families, collapse = ";"),
        mu_terms                = fmt_terms(rec_formulas$mu),
        sigma_terms             = fmt_terms(rec_formulas$sigma),
        nu_terms                = fmt_terms(rec_formulas$nu),
        tau_terms               = fmt_terms(rec_formulas$tau),
        processing_time         = elapsed
      ))
    }

    # =========================================================================
    # CONTINUOUS PATH
    # =========================================================================

    res_fit <- residualise_feature(df, feature_name, batch_var)
    resid   <- res_fit$residuals
    skew    <- compute_skewness(resid)
    kurt    <- compute_excess_kurtosis(resid)
    skew_st <- compute_skewness(res_fit$studentised)
    kurt_st <- compute_excess_kurtosis(res_fit$studentised)

    logger::log_info(paste0("  skew=", round(skew, 3L),
                            "  kurt=", round(kurt, 3L),
                            "  studentised skew=", round(skew_st, 3L),
                            "  kurt=", round(kurt_st, 3L)))
    logger::log_info(paste0("  sigma varies p=", signif(res_fit$sigma_p %||% NA, 3L),
                            "  CV2(sigma^2)=",   round(res_fit$sigma_cv2 %||% NA, 4L)))

    skew_perm <- permutation_test_batch_moment(
      resid, df[[batch_var]], compute_skewness, n_perm, min_batch_n,
      subjects = subj_vec, strata = strata_vec
    )
    logger::log_info(paste0(
      "  skew batch perm p=", round(skew_perm$p_value %||% NA, 4L),
      "  p_range=",           round(skew_perm$p_range %||% NA, 4L),
      "  range=",             round(skew_perm$per_batch_range %||% NA, 3L),
      "  frac_common=",       round(skew_perm$frac_common %||% NA, 3L)
    ))

    kurt_perm <- permutation_test_batch_moment(
      resid, df[[batch_var]], compute_excess_kurtosis, n_perm, min_batch_n,
      subjects = subj_vec, strata = strata_vec
    )
    logger::log_info(paste0(
      "  kurt batch perm p=", round(kurt_perm$p_value %||% NA, 4L),
      "  p_range=",           round(kurt_perm$p_range %||% NA, 4L),
      "  range=",             round(kurt_perm$per_batch_range %||% NA, 3L)
    ))

    tier_res <- assign_tier(
      skewness          = skew,
      excess_kurtosis   = kurt,
      batch_skew_pval   = skew_perm$p_value,
      batch_skew_range  = skew_perm$per_batch_range,
      batch_kurt_pval   = kurt_perm$p_value,
      batch_kurt_range  = kurt_perm$per_batch_range,
      thresholds        = thresholds,
      batch_skew_prange = skew_perm$p_range,
      batch_kurt_prange = kurt_perm$p_range
    )
    logger::log_info(paste0("  tier=", tier_res$tier))

    if (!is.null(user_override) && !is.null(user_override$family_order)) {
      rec_families <- unlist(user_override$family_order)
      logger::log_info(paste0("  family_order override: ",
                              paste(rec_families, collapse = " -> ")))
    } else {
      rec_families <- suggest_families(skew, kurt, is_positive_only, near_zero,
                                       skew_batch_effect = tier_res$skew_batch_effect,
                                       kurt_batch_effect = tier_res$kurt_batch_effect)
    }

    rec_formulas <- build_recommended_formulas(tier_res, longitudinal)

    if (!is.null(user_override)) {
      if (!is.null(user_override$nu_terms))
        rec_formulas$nu  <- user_override$nu_terms
      if (!is.null(user_override$tau_terms))
        rec_formulas$tau <- user_override$tau_terms
    }

    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    logger::log_info(paste0("  done in ", format_time(elapsed)))

    list(
      status              = "success",
      feature             = feature_name,
      discrete            = FALSE,
      n_obs               = n_obs,
      n_batches           = n_batches,
      is_positive_only    = is_positive_only,
      skewness            = round(skew, 4L),
      excess_kurtosis     = round(kurt, 4L),
      skewness_stud       = round(skew_st %||% NA, 4L),
      excess_kurtosis_stud = round(kurt_st %||% NA, 4L),
      sigma_varies_p      = signif(res_fit$sigma_p   %||% NA, 4L),
      sigma_cv2           = round(res_fit$sigma_cv2 %||% NA, 4L),
      batch_skew_pvalue   = round(skew_perm$p_value           %||% NA, 4L),
      batch_skew_prange   = round(skew_perm$p_range           %||% NA, 4L),
      batch_skew_range    = round(skew_perm$per_batch_range    %||% NA, 4L),
      batch_skew_var      = round(skew_perm$observed_variance  %||% NA, 4L),
      batch_skew_nbatches = skew_perm$n_batches_used %||% NA_integer_,
      batch_kurt_pvalue   = round(kurt_perm$p_value           %||% NA, 4L),
      batch_kurt_prange   = round(kurt_perm$p_range           %||% NA, 4L),
      batch_kurt_range    = round(kurt_perm$per_batch_range    %||% NA, 4L),
      batch_kurt_var      = round(kurt_perm$observed_variance  %||% NA, 4L),
      batch_kurt_nbatches = kurt_perm$n_batches_used %||% NA_integer_,
      frac_common         = round(skew_perm$frac_common %||% NA, 4L),
      seed                = seed %||% NA_integer_,
      tier                = tier_res$tier,
      near_gaussian       = tier_res$near_gaussian,
      skew_batch_effect   = tier_res$skew_batch_effect,
      kurt_batch_effect   = tier_res$kurt_batch_effect,
      family_override     = !is.null(user_override) &&
                             !is.null(user_override$family_order),
      family_order        = paste(rec_families, collapse = ";"),
      mu_terms            = fmt_terms(rec_formulas$mu),
      sigma_terms         = fmt_terms(rec_formulas$sigma),
      nu_terms            = fmt_terms(rec_formulas$nu),
      tau_terms           = fmt_terms(rec_formulas$tau),
      processing_time     = elapsed
    )

  }, error = function(e) {
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    logger::log_error(paste0("Error diagnosing ", feature_name, ": ", e$message))
    list(status          = "error",
         feature         = feature_name,
         error_message   = e$message,
         processing_time = elapsed)
  })
}

# ---------------------------------------------------------------------------
# Diagnose all features.
# ---------------------------------------------------------------------------

diagnose_all_features <- function(data, features, output_dir,
                                   batch_var        = "batch",
                                   id_var           = "id",
                                   longitudinal     = FALSE,
                                   discrete         = FALSE,
                                   n_perm           = 499L,
                                   min_batch_n      = 50L,
                                   thresholds       = list(skew                 = 0.5,
                                                           kurt                 = 1.0,
                                                           pval                 = 0.01,
                                                           min_skew_range       = 0.3,
                                                           min_kurt_range       = 0.5,
                                                           min_zerorate_range   = 0.10,
                                                           min_dispersion_range = 0.50),
                                   user_overrides   = NULL,
                                   n_cores          = 1L,
                                   seed             = 20260818L,
                                   stratify         = TRUE,
                                   n_age_bands      = 10L) {

  overall_start <- Sys.time()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  logger::log_info(paste0("Diagnosing ", length(features), " features",
                          if (discrete) " [discrete mode]" else ""))

  # Seeding per feature rather than once up front keeps results independent of
  # n_cores and of how work is distributed across workers.
  diagnose_one <- function(i) {
    feat <- features[i]
    feat_seed <- if (is.null(seed)) NULL else feature_seed(seed, feat)
    diagnose_feature(
      data          = data,
      feature_name  = feat,
      batch_var     = batch_var,
      id_var        = id_var,
      longitudinal  = longitudinal,
      discrete      = discrete,
      n_perm        = n_perm,
      min_batch_n   = min_batch_n,
      thresholds    = thresholds,
      user_override = if (!is.null(user_overrides)) user_overrides[[feat]] else NULL,
      seed          = feat_seed,
      stratify      = stratify,
      n_age_bands   = n_age_bands
    )
  }

  run_one <- function(i) {
    tryCatch(diagnose_one(i),
             error = function(e)
               list(status = "error", feature = features[i],
                    error_message = e$message, processing_time = 0))
  }

  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    # Covers the unseeded path only; when seed is set, diagnose_feature() pins
    # both the seed and the RNG kind, which overrides this.
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    parallel::clusterEvalQ(cl, {
      library(logger)
      logger::log_formatter(logger::formatter_paste)
    })
    parallel::clusterExport(cl,
      c("data", "features", "batch_var", "id_var", "longitudinal", "discrete",
        "n_perm", "min_batch_n", "thresholds", "user_overrides",
        "seed", "stratify", "n_age_bands", "diagnose_one"),
      envir = environment()
    )
    parallel::clusterExport(cl,
      c("diagnose_feature", "residualise_feature", "AGE_SPLINE_DF",
        "compute_skewness", "compute_excess_kurtosis",
        "permutation_test_batch_moment", "build_strata", "feature_seed",
        "test_zero_inflation", "test_overdispersion",
        "test_underdispersion", "compute_dispersion_index",
        "compute_tail_heaviness",
        "assign_tier", "assign_count_tier",
        "suggest_families", "suggest_count_families",
        "build_recommended_formulas", "build_count_recommended_formulas",
        "format_time", "%||%"),
      envir = globalenv()
    )
    results <- pbapply::pblapply(seq_along(features), run_one, cl = cl)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(seq_along(features), run_one)
  }

  successes <- sum(vapply(results, function(x) x$status == "success", logical(1L)))
  failures  <- length(results) - successes

  fields_discrete   <- c("zi_statistic", "zi_pvalue", "obs_zero_rate", "exp_zero_rate",
                          "od_statistic", "od_pvalue", "dispersion_ratio",
                          "ud_statistic", "ud_pvalue",
                          "dispersion_index", "tail_heaviness",
                          "batch_zerorate_pvalue", "batch_zerorate_prange",
                          "batch_zerorate_range",
                          "batch_zerorate_var", "batch_zerorate_nbatches",
                          "batch_disp_pvalue", "batch_disp_prange",
                          "batch_disp_range",
                          "batch_disp_var", "batch_disp_nbatches",
                          "zero_inflated", "overdispersed", "underdispersed",
                          "severe_od", "heavy_tail",
                          "zerorate_batch_effect", "dispersion_batch_effect")
  fields_continuous <- c("skewness_stud", "excess_kurtosis_stud",
                          "sigma_varies_p", "sigma_cv2",
                          "batch_skew_pvalue", "batch_skew_prange",
                          "batch_skew_range",
                          "batch_skew_var", "batch_skew_nbatches",
                          "batch_kurt_pvalue", "batch_kurt_prange",
                          "batch_kurt_range",
                          "batch_kurt_var", "batch_kurt_nbatches",
                          "near_gaussian", "skew_batch_effect", "kurt_batch_effect")
  fields_shared     <- c("frac_common", "seed")

  summary_rows <- lapply(results, function(r) {
    base <- data.frame(
      status           = r$status           %||% NA_character_,
      feature          = r$feature          %||% NA_character_,
      discrete         = r$discrete         %||% NA,
      n_obs            = r$n_obs            %||% NA_integer_,
      n_batches        = r$n_batches        %||% NA_integer_,
      is_positive_only = r$is_positive_only %||% NA,
      skewness         = r$skewness         %||% NA_real_,
      excess_kurtosis  = r$excess_kurtosis  %||% NA_real_,
      tier             = r$tier             %||% NA_integer_,
      stringsAsFactors = FALSE
    )
    all_fields <- c(fields_discrete, fields_continuous, fields_shared,
                    "family_override", "family_order",
                    "mu_terms", "sigma_terms", "nu_terms", "tau_terms")
    for (f in all_fields) base[[f]] <- r[[f]] %||% NA
    base
  })

  all_cols   <- unique(unlist(lapply(summary_rows, names)))
  summary_df <- do.call(rbind, lapply(summary_rows, function(r) {
    missing_cols <- setdiff(all_cols, names(r))
    for (col in missing_cols) r[[col]] <- NA
    r[, all_cols, drop = FALSE]
  }))
  # Merge into the existing CSV so re-diagnosing a subset updates only those
  # rows instead of dropping the rest (diagnose has no per-feature files to rescan).
  summary_path <- file.path(output_dir, "diagnostics_summary.csv")
  summary_df   <- upsert_rows_by_key(summary_df, summary_path, key_col = "feature")
  write.csv(summary_df, summary_path, row.names = FALSE)
  logger::log_info(paste0("Saved: ", summary_path,
                          " (", nrow(summary_df), " feature(s) on disk)"))

  success_results <- Filter(function(x) x$status == "success", results)

  recs <- if (length(success_results) == 0L) {
    logger::log_info("No successful features -- writing empty recommendations CSV")
    data.frame(feature      = character(0L),
               family_order = character(0L),
               mu_terms     = character(0L),
               sigma_terms  = character(0L),
               nu_terms     = character(0L),
               tau_terms    = character(0L),
               tier         = integer(0L),
               stringsAsFactors = FALSE)
  } else {
    do.call(rbind, lapply(success_results, function(r)
      data.frame(
        feature      = r$feature,
        family_order = r$family_order %||% "",
        mu_terms     = r$mu_terms     %||% "",
        sigma_terms  = r$sigma_terms  %||% "",
        nu_terms     = r$nu_terms,
        tau_terms    = r$tau_terms,
        tier         = r$tier         %||% NA_integer_,
        stringsAsFactors = FALSE
      )
    ))
  }
  recs_path <- file.path(output_dir, "feature_recommendations.csv")
  recs      <- upsert_rows_by_key(recs, recs_path, key_col = "feature")
  write.csv(recs, recs_path, row.names = FALSE)
  logger::log_info(paste0("Saved: ", recs_path,
                          " (", nrow(recs), " feature(s) on disk)"))

  elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
  logger::log_info(paste0("Diagnostics complete | success: ", successes,
                           " | failed: ", failures,
                           " | time: ", format_time(elapsed)))

  invisible(list(
    results              = results,
    summary              = summary_df,
    recommendations      = recs,
    successes            = successes,
    failures             = failures,
    total_time           = format_time(elapsed),
    summary_path         = summary_path,
    recommendations_path = recs_path
  ))
}