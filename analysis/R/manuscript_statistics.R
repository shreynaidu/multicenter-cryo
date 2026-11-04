write_cohort_flow <- function(dat, out_dir) {
  flow <- tibble(
    analysis = c(
      "Primary", "Technical-success sensitivity", "Known-stage sensitivity",
      "Unknown-stage category sensitivity", "Three-group size sensitivity"
    ),
    criterion = c(
      "Retain every source row", "technical_success == Yes", "stage_raw is nonmissing",
      "Retain unknown as a stage category", "tumor_size_cm >= 1"
    ),
    records_before = nrow(dat),
    records_after = c(
      nrow(dat), sum(dat$technical_success == "Yes"), sum(!is.na(dat$stage_raw)), nrow(dat), sum(dat$tumor_size_cm >= 1)
    )
  ) %>%
    mutate(records_excluded = records_before - records_after)
  write_table(flow, file.path(out_dir, "tables", "cohort_flow.csv"))
}

descriptive_counts <- function(dat) {
  counts <- c(
    patients = nrow(dat),
    copd = sum(dat$copd == "Yes"),
    home_oxygen = sum(dat$oxygen_any == "Any home oxygen"),
    diameter_at_least_3_cm = sum(dat$tumor_size_cm >= 3),
    central_location = sum(dat$central_location == "Central"),
    known_stage = sum(!is.na(dat$stage_raw)),
    stage_i_ii = sum(dat$stage_group == "Stage I-II"),
    stage_iii_iv = sum(dat$stage_group == "Stage III-IV"),
    systemic_therapy = sum(dat$systemic_therapy == "Yes"),
    technical_success = sum(dat$technical_success == "Yes"),
    pneumothorax = sum(dat$pneumothorax),
    chest_tube = sum(dat$chest_tube),
    grade_3 = sum(dat$sir_comp_grade == 3),
    grade_3_with_tube = sum(dat$sir_comp_grade == 3 & dat$chest_tube == 1),
    grade_3_without_tube = sum(dat$sir_comp_grade == 3 & dat$chest_tube == 0),
    grade_2_with_tube = sum(dat$sir_comp_grade == 2 & dat$chest_tube == 1),
    death_30d = sum(dat$death_30d),
    death_90d = sum(dat$death_90d),
    local_progressions = sum(dat$cr_status == 1),
    competing_deaths = sum(dat$cr_status == 2),
    composite_events = sum(dat$lpfs_event),
    all_cause_deaths = sum(dat$os_event),
    progression_before_6_months = sum(dat$cr_status == 1 & dat$cr_time_months < 6),
    progression_at_or_after_6_months = sum(dat$cr_status == 1 & dat$cr_time_months >= 6),
    repeat_ablation = sum(dat$repeat_ablation),
    repeat_ablation_with_progression = sum(dat$repeat_ablation == 1 & dat$local_progression_corrected == 1),
    unknown_stage = sum(is.na(dat$stage_raw)),
    progression_status_recodes = sum(dat$local_progression_correction_flag)
  )
  over_cohort <- tibble(
    statistic = names(counts),
    numerator = as.integer(counts),
    denominator = nrow(dat),
    percent = 100 * numerator / denominator
  )
  known_stage <- over_cohort %>%
    filter(statistic %in% c("stage_i_ii", "stage_iii_iv")) %>%
    mutate(
      statistic = paste0(statistic, "_among_known"),
      denominator = sum(!is.na(dat$stage_raw)),
      percent = 100 * numerator / denominator
    )
  among_progressions <- over_cohort %>%
    filter(statistic == "repeat_ablation_with_progression") %>%
    mutate(
      statistic = "repeat_ablation_among_progressions",
      denominator = sum(dat$local_progression_corrected),
      percent = 100 * numerator / denominator
    )
  bind_rows(over_cohort, known_stage, among_progressions)
}

write_derived_probability <- function(statistics_dir, out_dir) {
  ltp <- read_csv(file.path(out_dir, "tables", "cif_all_cohorts.csv"), show_col_types = FALSE) %>%
    filter(cohort == "Full107", event == "Local tumor progression", time_months == 24)
  complement <- tibble(
    statistic = "one_minus_ltp_cif_24_months",
    estimate = 1 - ltp$estimate,
    percent = 100 * estimate,
    source_file = "tables/cif_all_cohorts.csv",
    selection = "Full107; Local tumor progression; 24 months",
    derivation = "1 - estimate"
  )
  write_table(complement, file.path(statistics_dir, "derived_probabilities.csv"))
}

source_cell_registry <- function(out_dir) {
  folders <- c("tables", "model_results", "diagnostics", "supplemental", "reviewer_sensitivities")
  label_columns <- c(
    "characteristic", "level", "outcome", "category", "cohort", "stratum", "endpoint",
    "model", "term", "time_months", "statistic", "population", "stage_group"
  )
  files <- unlist(lapply(folders, function(folder) {
    paths <- list.files(file.path(out_dir, folder), pattern = "\\.csv$", full.names = TRUE)
    # Keeps x_summary.csv ahead of x.csv in every locale.
    paths[order(chartr(".", "~", paths))]
  }))
  map_dfr(files, function(path) {
    table <- read.csv(path, colClasses = "character", check.names = FALSE, na.strings = character())
    if (!nrow(table)) {
      return(tibble())
    }
    labels <- intersect(label_columns, names(table))
    item <- if (length(labels)) {
      apply(table[, labels, drop = FALSE], 1, paste, collapse = " | ")
    } else {
      as.character(seq_len(nrow(table)))
    }
    map_dfr(setdiff(names(table), labels), function(column) {
      tibble(
        source_file = substring(path, nchar(out_dir) + 2L),
        source_row = seq_len(nrow(table)),
        item = item,
        field = column,
        value = table[[column]]
      )
    })
  })
}

run_manuscript_statistics <- function(dat, out_dir) {
  statistics_dir <- file.path(out_dir, "manuscript_statistics")
  dir.create(statistics_dir, showWarnings = FALSE)
  write_cohort_flow(dat, out_dir)
  write_table(descriptive_counts(dat), file.path(statistics_dir, "descriptive_counts.csv"))
  write_derived_probability(statistics_dir, out_dir)
  write_table(source_cell_registry(out_dir), file.path(statistics_dir, "source_cells.csv"))
}
