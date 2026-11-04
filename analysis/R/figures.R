figure_blue <- "#0072B2"
figure_orange <- "#D55E00"
display_months <- c(0, 36)
figure_width <- 7.4

size_palette <- setNames(c(figure_blue, figure_orange, "#009E73", "#CC79A7"), size_levels)
stage_palette <- setNames(c(figure_blue, figure_orange, "#666666"), stage_levels)

main_theme <- function() {
  theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 8, 0, 5.5)
    )
}

probability_axes <- function(y_label) {
  list(
    scale_x_continuous(breaks = risk_table_months, expand = expansion(mult = c(0, 0.015))),
    scale_y_continuous(labels = percent_format(accuracy = 1)),
    coord_cartesian(xlim = display_months, ylim = c(0, 1)),
    labs(x = NULL, y = y_label),
    main_theme()
  )
}

cif_curves <- function(fit, cause, stratum_label = NULL) {
  suffix <- paste0(" ", cause, "$")
  curves <- setdiff(names(fit), "Tests")
  map_dfr(curves[grepl(suffix, curves)], function(curve) {
    points <- fit[[curve]]
    se <- sqrt(pmax(0, points$var))
    tibble(
      time = points$time,
      estimate = points$est,
      lower = pmax(0, points$est - normal_95 * se),
      upper = pmin(1, points$est + normal_95 * se),
      stratum = if (is.null(stratum_label)) sub(suffix, "", curve) else stratum_label
    )
  })
}

km_curves <- function(fit, strata = NULL) {
  estimates <- broom::tidy(fit)
  start <- tibble(time = 0, estimate = 1, conf.low = 1, conf.high = 1, n.censor = 0)
  if (is.null(strata)) {
    return(bind_rows(start, estimates))
  }
  bind_rows(
    tibble(start, stratum = strata),
    mutate(estimates, stratum = sub("^group=", "", strata))
  )
}

risk_counts <- function(times, groups, levels) {
  map_dfr(levels, function(level) {
    tibble(
      stratum = level,
      month = risk_table_months,
      number_at_risk = number_at_risk(times[groups == level], risk_table_months)
    )
  })
}

overall_risk_counts <- function(times) {
  risk_counts(times, rep("Overall", length(times)), "Overall")
}

risk_table_plot <- function(risk, levels) {
  risk$stratum <- factor(risk$stratum, levels = rev(levels))
  ggplot(risk, aes(month, stratum, label = number_at_risk)) +
    geom_text(size = 3.6, color = "black") +
    scale_x_continuous(breaks = risk_table_months, limits = display_months, expand = expansion(mult = c(0, 0.015))) +
    coord_cartesian(xlim = display_months, clip = "off") +
    labs(x = "Months after cryoablation", y = NULL, title = "Number at risk") +
    theme_classic(base_size = 10) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(margin = margin(r = 18)),
      axis.title.x = element_text(face = "bold"),
      plot.title = element_text(size = 10),
      plot.margin = margin(0, 8, 5.5, 5.5)
    )
}

with_risk_table <- function(curve_plot, risk, levels, table_height) {
  curve_plot / risk_table_plot(risk, levels) + plot_layout(heights = c(4, table_height))
}

stratified_plot <- function(curves, palette, y_label, interval, interval_alpha, censor_size = NULL, annotation = NULL) {
  curves$stratum <- factor(curves$stratum, levels = names(palette))
  plot <- ggplot(curves, aes(time, estimate, color = stratum, fill = stratum)) +
    geom_ribbon(aes(ymin = .data[[interval[1]]], ymax = .data[[interval[2]]]), alpha = interval_alpha, color = NA) +
    geom_step(linewidth = 0.95, linetype = "solid")
  if (!is.null(censor_size)) {
    plot <- plot + geom_point(data = filter(curves, n.censor > 0), shape = 3, size = censor_size)
  }
  plot <- plot +
    scale_color_manual(values = palette, breaks = names(palette)) +
    scale_fill_manual(values = palette, breaks = names(palette)) +
    probability_axes(y_label)
  if (!is.null(annotation)) {
    plot <- plot + annotate("text", x = 25, y = 0.92, label = annotation, size = 3.5, color = "black")
  }
  plot
}

export_figure <- function(plot, figure_dir, stem, height) {
  path <- function(suffix) file.path(figure_dir, paste0(stem, suffix))
  ggsave(path(".svg"), plot, width = figure_width, height = height, device = svglite)
  rsvg::rsvg_pdf(path(".svg"), path(".pdf"))
  ggsave(
    path("_preview.png"), plot,
    width = figure_width, height = height, dpi = 180, device = ragg::agg_png, bg = "white"
  )
  ggsave(
    path("_600dpi.tiff"), plot,
    width = figure_width, height = height, dpi = 600, device = ragg::agg_tiff, bg = "white", compression = "lzw"
  )
}

write_figure_data <- function(table, figure_dir, name) {
  write_table(table, file.path(figure_dir, paste0(name, ".csv")), na = "NA")
}

overall_incidence_figure <- function(data, figure_dir) {
  fit <- cuminc(data$cr_time_months, data$cr_status, cencode = 0)
  curves <- bind_rows(
    cif_curves(fit, 1, "Local progression"),
    cif_curves(fit, 2, "Death before progression")
  )
  outcome_colors <- c(`Local progression` = figure_blue, `Death before progression` = figure_orange)
  plot <- ggplot(curves, aes(time, estimate, color = stratum, fill = stratum, linetype = stratum)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = NA) +
    geom_step(linewidth = 0.95) +
    scale_color_manual(values = outcome_colors) +
    scale_fill_manual(values = outcome_colors) +
    scale_linetype_manual(values = c(`Local progression` = "solid", `Death before progression` = "dashed")) +
    probability_axes("Cumulative incidence")
  risk <- overall_risk_counts(data$cr_time_months)
  export_figure(with_risk_table(plot, risk, "Overall", 1), figure_dir, "figure1_overall_cif", 5.4)
  write_figure_data(risk, figure_dir, "figure1_risk_counts")
  write_figure_data(curves, figure_dir, "figure1_curve_data")
}

overall_survival_figure <- function(data, figure_dir) {
  fit <- survfit(Surv(os_time_months, os_event) ~ 1, data = data, conf.type = "log-log")
  curves <- km_curves(fit)
  plot <- ggplot(curves, aes(time, estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = figure_blue, alpha = 0.12) +
    geom_step(color = figure_blue, linewidth = 0.95) +
    geom_point(data = filter(curves, n.censor > 0), shape = 3, size = 1.7, color = figure_blue) +
    probability_axes("Overall survival probability")
  risk <- overall_risk_counts(data$os_time_months)
  export_figure(with_risk_table(plot, risk, "Overall", 1), figure_dir, "figure3_overall_os", 5.4)
  write_figure_data(risk, figure_dir, "figure3_risk_counts")
  write_figure_data(curves, figure_dir, "figure3_curve_data")
}

incidence_by_group_figure <- function(data, grouping, palette, stem, data_prefix, interval_alpha, table_height, height,
                                      figure_dir, annotation = NULL) {
  groups <- data[[grouping]]
  fit <- cuminc(data$cr_time_months, data$cr_status, group = groups, cencode = 0)
  curves <- cif_curves(fit, 1)
  plot <- stratified_plot(
    curves, palette, "Cumulative incidence of local progression", c("lower", "upper"), interval_alpha,
    annotation = annotation
  )
  risk <- risk_counts(data$cr_time_months, groups, names(palette))
  export_figure(with_risk_table(plot, risk, names(palette), table_height), figure_dir, stem, height)
  write_figure_data(risk, figure_dir, paste0(data_prefix, "_risk_counts"))
  write_figure_data(curves, figure_dir, paste0(data_prefix, "_curve_data"))
}

survival_by_size_figure <- function(data, endpoint, stem, y_label, annotation, figure_dir) {
  columns <- survival_endpoints[[endpoint]]
  plotted <- tibble(time = data[[columns$time]], event = data[[columns$event]], group = data$tumor_size_group)
  fit <- survfit(Surv(time, event) ~ group, data = plotted, conf.type = "log-log")
  curves <- km_curves(fit, size_levels)
  plot <- stratified_plot(
    curves, size_palette, y_label, c("conf.low", "conf.high"), 0.05,
    censor_size = 1.5, annotation = annotation
  )
  risk <- risk_counts(plotted$time, plotted$group, size_levels)
  export_figure(with_risk_table(plot, risk, size_levels, 1.8), figure_dir, stem, 6.1)
  write_figure_data(risk, figure_dir, paste0(stem, "_risk_counts"))
  write_figure_data(curves, figure_dir, paste0(stem, "_curve_data"))
}

size_group_p_label <- function(group_tests, endpoint, test) {
  p <- group_tests$p_value[
    group_tests$grouping == "tumor_size_group" & group_tests$endpoint == endpoint & group_tests$test == test
  ]
  formatted <- if (p < 0.001) "<.001" else sub("^0", "", sprintf("%.3f", p))
  paste0(test, if (test == "Gray") " test", " P = ", formatted)
}

# Reads the exported CSV because the reported figures were drawn from it.
run_figures <- function(out_dir) {
  figure_dir <- file.path(out_dir, "publication_figures")
  dir.create(figure_dir, showWarnings = FALSE)
  data <- read_csv(file.path(out_dir, "analysis_data", "reanalysis_dataset.csv"), show_col_types = FALSE)
  data$tumor_size_group <- factor(data$tumor_size_group, levels = size_levels)
  data$stage_group <- factor(data$stage_group, levels = stage_levels)
  group_tests <- read_csv(file.path(out_dir, "model_results", "group_tests.csv"), show_col_types = FALSE)

  overall_incidence_figure(data, figure_dir)
  overall_survival_figure(data, figure_dir)

  incidence_by_group_figure(
    data, "tumor_size_group", size_palette, "supplement_figure1_ltp_size", "supplement_figure1_ltp_size",
    0.05, 1.8, 6.1, figure_dir,
    annotation = size_group_p_label(group_tests, "LTP", "Gray")
  )
  incidence_by_group_figure(
    data, "stage_group", stage_palette, "figure2_ltp_cif_stage", "figure2_stage", 0.07, 1.5, 5.9, figure_dir
  )
  survival_by_size_figure(
    data, "OS", "supplement_figure2_os_size", "Overall survival probability",
    size_group_p_label(group_tests, "OS", "Log-rank"), figure_dir
  )
  survival_by_size_figure(
    data, "Composite", "supplement_figure3_composite_size", "Progression-or-death-free probability",
    size_group_p_label(group_tests, "Composite", "Log-rank"), figure_dir
  )
}
