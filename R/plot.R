compute_value_limits <- function(values, percentiles = c(0.01, 0.99), pad_frac = 0.05) {
  vals <- values[is.finite(values)]
  if (length(vals) == 0) return(list(lo = NA_real_, hi = NA_real_))
  q_lo <- quantile(vals, percentiles[1])
  q_hi <- quantile(vals, percentiles[2])
  pad  <- (q_hi - q_lo) * pad_frac
  list(lo = q_lo - pad, hi = q_hi + pad)
}

prepare_trajectory_data <- function(pre_data = NULL, post_data = NULL,
                                    feature_name, batch_var = "cohort",
                                    id_var = "id") {
  if (is.null(pre_data) && is.null(post_data))
    stop("At least one of pre_data or post_data must be supplied")

  out <- NULL

  if (!is.null(pre_data)) {
    required <- c(id_var, batch_var, "age", "sex", feature_name)
    missing  <- setdiff(required, names(pre_data))
    if (length(missing) > 0)
      stop("pre_data missing columns: ", paste(missing, collapse = ", "))

    out <- rbind(out, data.frame(
      id        = pre_data[[id_var]],
      group     = as.character(pre_data[[batch_var]]),
      age       = pre_data$age,
      sex       = pre_data$sex,
      value     = pre_data[[feature_name]],
      data_type = "Pre-Harmonisation",
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(post_data)) {
    # Wide format: feature columns are named <feature>.<stat>
    harm_col <- paste0(feature_name, ".harmonised_value")
    if (!harm_col %in% names(post_data))
      stop("Column '", harm_col, "' not found in post_data")
    required <- c(id_var, "batch", "age", "sex", harm_col)
    missing  <- setdiff(required, names(post_data))
    if (length(missing) > 0)
      stop("post_data missing columns: ", paste(missing, collapse = ", "))

    out <- rbind(out, data.frame(
      id        = post_data[[id_var]],
      group     = as.character(post_data$batch),
      age       = post_data$age,
      sex       = post_data$sex,
      value     = post_data[[harm_col]],
      data_type = "Post-Harmonisation",
      stringsAsFactors = FALSE
    ))
  }

  out <- out[is.finite(out$value) & is.finite(out$age), ]
  if (nrow(out) == 0)
    stop("No finite observations for '", feature_name, "'")
  out
}

create_age_bin_summary <- function(data, age_bin_width = 5) {
  data$age_bin        <- floor(data$age / age_bin_width) * age_bin_width
  data$age_bin_center <- data$age_bin + age_bin_width / 2

  data %>%
    dplyr::group_by(group, sex, data_type, age_bin, age_bin_center) %>%
    dplyr::summarise(
      n_obs      = dplyr::n(),
      mean_value = mean(value, na.rm = TRUE),
      sd_value   = sd(value,   na.rm = TRUE),
      se_value   = sd(value,   na.rm = TRUE) / sqrt(dplyr::n()),
      q25        = quantile(value, 0.25, na.rm = TRUE),
      q75        = quantile(value, 0.75, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::filter(n_obs >= 5)
}

create_trajectory_plot <- function(bin_summary, feature_name, output_dir,
                                   group_col = "cohort", smooth_method = "loess",
                                   value_limits = NULL) {
  if (nrow(bin_summary) == 0)
    stop("bin_summary is empty for: ", feature_name)

  groups   <- sort(unique(bin_summary$group))
  n_groups <- length(groups)
  pal      <- make_group_palette(n_groups)
  shapes   <- rep(15:25, length.out = n_groups)

  p <- ggplot2::ggplot(
    bin_summary,
    ggplot2::aes(x = age_bin_center, y = mean_value,
                 colour = group, shape = group)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mean_value - se_value,
                   ymax = mean_value + se_value, fill = group),
      alpha = 0.10, colour = NA
    ) +
    ggplot2::geom_point(ggplot2::aes(size = n_obs), alpha = 0.8) +
    ggplot2::geom_smooth(ggplot2::aes(group = group), method = smooth_method,
                         se = FALSE, alpha = 0.7, linewidth = 0.6, formula = y ~ x) +
    ggplot2::facet_grid(sex ~ data_type, scales = "free_y",
                        labeller = ggplot2::labeller(
                          sex = function(x) paste("Sex:", x))) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal,    groups)) +
    ggplot2::scale_fill_manual(values   = stats::setNames(pal,    groups)) +
    ggplot2::scale_shape_manual(values  = stats::setNames(shapes, groups)) +
    ggplot2::scale_size_continuous(range = c(1, 4), name = "n") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::labs(
      title    = paste("Age Trajectories:", feature_name),
      subtitle = paste("Grouped by", group_col, "|",
                       toupper(smooth_method), "smooth | Ribbons = +/-1 SE"),
      x = "Age", y = paste("Mean", feature_name),
      colour = group_col, fill = group_col, shape = group_col
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = if (n_groups > 20) "none" else "bottom",
      plot.title       = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 11),
      strip.text       = ggplot2::element_text(size = 10, face = "bold"),
      axis.text        = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = min(6, n_groups), title.position = "top",
                                     override.aes = list(size = 3)),
      fill   = ggplot2::guide_legend(ncol = min(6, n_groups), title.position = "top"),
      shape  = ggplot2::guide_legend(ncol = min(6, n_groups), title.position = "top"),
      size   = ggplot2::guide_legend(title.position = "top")
    )

  if (!is.null(value_limits) &&
      is.finite(value_limits$lo) && is.finite(value_limits$hi))
    p <- p + ggplot2::coord_cartesian(ylim = c(value_limits$lo, value_limits$hi))

  out_file <- file.path(output_dir, paste0(feature_name, "_trajectory.png"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p, width = 16, height = 10, dpi = 300)
  message("  Saved: ", out_file)
  invisible(p)
}

plot_trajectories <- function(features, pre_data = NULL, post_data = NULL,
                              output_dir, batch_var = "cohort", id_var = "id",
                              group_col = "cohort", smooth_method = "loess",
                              age_bin_width = 5, fix_y_limits = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (feat in features) {
    message("Plotting: ", feat)
    tryCatch({
      traj_data    <- prepare_trajectory_data(pre_data = pre_data, post_data = post_data,
                                              feature_name = feat,
                                              batch_var = batch_var, id_var = id_var)
      bin_summary  <- create_age_bin_summary(traj_data, age_bin_width)
      value_limits <- if (fix_y_limits) compute_value_limits(traj_data$value) else NULL
      create_trajectory_plot(bin_summary = bin_summary, feature_name = feat,
                             output_dir = output_dir, group_col = group_col,
                             smooth_method = smooth_method, value_limits = value_limits)
    }, error = function(e) message("  ERROR: ", e$message))
  }
  invisible(NULL)
}
