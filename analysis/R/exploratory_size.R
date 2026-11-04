# Exploratory, unadjusted. Separate from the four primary diameter groups.

size3_levels <- c("<3 cm", "3 to <5 cm", ">=5 cm")

read_exploratory_data <- function(input, cohort) {
  data <- read_csv(input, show_col_types = FALSE) %>%
    mutate(size3 = factor(
      case_when(tumor_size_cm < 3 ~ "<3 cm", tumor_size_cm < 5 ~ "3 to <5 cm", TRUE ~ ">=5 cm"),
      levels = size3_levels
    ))
  stopifnot(
    nrow(data) == nrow(cohort), !anyNA(data$size3),
    sum(data$cr_status == 1) == sum(cohort$cr_status == 1),
    sum(data$cr_status == 2) == sum(cohort$cr_status == 2),
    sum(data$lpfs_event) == sum(cohort$lpfs_event)
  )
  data
}

local_control_milestones <- function(subgroup, cohort, size) {
  fit <- survfit(Surv(cr_time_months, cr_status == 1) ~ 1, data = subgroup, conf.type = "log-log")
  at_times <- summary(fit, times = milestone_months, extend = TRUE)
  tibble(
    cohort = cohort,
    stratum = size,
    time_months = milestone_months,
    n = nrow(subgroup),
    estimate = at_times$surv,
    lower = at_times$lower,
    upper = at_times$upper,
    number_at_risk = number_at_risk(subgroup$cr_time_months, milestone_months)
  ) %>%
    mutate(across(c(estimate, lower, upper), ~ if_else(number_at_risk == 0, NA_real_, .x)))
}

summarize_exploratory_size <- function(data, cohort) {
  by_size <- map(set_names(levels(data$size3)), ~ filter(data, size3 == .x))
  gray <- cuminc(data$cr_time_months, data$cr_status, group = data$size3, cencode = 0)
  composite_test <- survdiff(Surv(cr_time_months, lpfs_event) ~ size3, data = data)
  control_test <- survdiff(Surv(cr_time_months, cr_status == 1) ~ size3, data = data)
  list(
    cif = imap_dfr(by_size, ~ cif_estimates(.x, cohort, .y)),
    km = imap_dfr(by_size, ~ km_estimates(.x, cohort, "Composite", .y)$milestones),
    lc = imap_dfr(by_size, ~ local_control_milestones(.x, cohort, .y)),
    tests = tibble(
      cohort = cohort,
      endpoint = c("LTP with competing death", "Local progression or death", "Death-censored local control"),
      test = c("Gray", "Log-rank", "Log-rank"),
      p_value = c(
        gray$Tests["1", "pv"],
        pchisq(composite_test$chisq, 2, lower.tail = FALSE),
        pchisq(control_test$chisq, 2, lower.tail = FALSE)
      )
    )
  )
}

exploratory_curves <- function(data) {
  map_dfr(levels(data$size3), function(size) {
    subgroup <- filter(data, size3 == size)
    cif <- cuminc(subgroup$cr_time_months, subgroup$cr_status, cencode = 0)[["1 1"]]
    composite <- survfit(Surv(cr_time_months, lpfs_event) ~ 1, data = subgroup)
    bind_rows(
      tibble(time = cif$time, estimate = cif$est, endpoint = "Local tumor progression", group = size),
      tibble(
        time = c(0, composite$time), estimate = c(1, composite$surv),
        endpoint = "Local progression-or-death-free survival", group = size
      )
    )
  })
}

plot_exploratory_curves <- function(curves, n_records, path) {
  # Colour order from the reported figure.
  curves$group <- factor(curves$group, levels = c("<3 cm", ">=5 cm", "3 to <5 cm"))
  plot <- ggplot(curves, aes(time, estimate, color = group)) +
    geom_step(linewidth = 0.85) +
    facet_wrap(~endpoint) +
    coord_cartesian(xlim = c(0, 36), ylim = c(0, 1)) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = seq(0, 36, 12)) +
    scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
    labs(
      x = "Months after cryoablation", y = "Estimated probability", color = "Tumor diameter",
      title = "Outcomes by tumor diameter before salvage cryoablation",
      subtitle = sprintf("Exploratory, unadjusted analysis; %d records; updated supplied follow-up", n_records),
      caption = paste0(
        "LTP accounts for death as a competing event. Survival counts local progression or death.\n",
        "Subgroup follow-up becomes sparse; see milestone tables for 95% CIs and numbers at risk."
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(path, plot, width = 11, height = 5.5, dpi = 180)
}

run_exploratory_size <- function(cohort, out_dir, make_figures) {
  exploratory_dir <- file.path(out_dir, "exploratory_size")
  dir.create(exploratory_dir, showWarnings = FALSE)
  input <- file.path(out_dir, "analysis_data", "reanalysis_dataset.csv")
  data <- read_exploratory_data(input, cohort)

  counts <- data %>%
    group_by(size3) %>%
    summarise(
      n = n(), min_cm = min(tumor_size_cm), max_cm = max(tumor_size_cm),
      local_progressions = sum(cr_status == 1), deaths_before_progression = sum(cr_status == 2),
      .groups = "drop"
    )
  write_table(counts, file.path(exploratory_dir, "group_counts.csv"), na = "NA")

  full_cohort <- summarize_exploratory_size(data, "Full107")
  excluding_subcentimeter <- summarize_exploratory_size(filter(data, tumor_size_cm >= 1), "Exclude_subcentimeter104")
  for (name in names(full_cohort)) {
    write_table(
      bind_rows(full_cohort[[name]], excluding_subcentimeter[[name]]),
      file.path(exploratory_dir, paste0(name, ".csv")),
      na = "NA"
    )
  }
  curves <- exploratory_curves(data)
  write_table(curves, file.path(exploratory_dir, "curve_data.csv"), na = "NA")
  if (make_figures) {
    plot_exploratory_curves(curves, nrow(data), file.path(exploratory_dir, "size_group_curves.png"))
  }
  writeLines(
    c(paste("Input MD5:", tools::md5sum(input)), capture.output(sessionInfo())),
    file.path(exploratory_dir, "provenance.txt")
  )
}
