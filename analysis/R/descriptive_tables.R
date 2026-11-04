summarize_continuous <- function(values) {
  sprintf(
    "%.1f (IQR %.1f–%.1f; range %.1f–%.1f)",
    median(values, na.rm = TRUE), sample_quantile(values, 0.25), sample_quantile(values, 0.75),
    min(values, na.rm = TRUE), max(values, na.rm = TRUE)
  )
}

format_fraction <- function(count, denominator) {
  sprintf("%d/%d (%.1f%%)", count, denominator, 100 * count / denominator)
}

summarize_category <- function(values, level) {
  format_fraction(sum(values == level, na.rm = TRUE), sum(!is.na(values)))
}

summarize_events <- function(events) {
  format_fraction(sum(events), length(events))
}

row_by_size <- function(dat, characteristic, level, values, summarize) {
  row <- tibble(characteristic = characteristic, level = level, overall = summarize(values))
  for (size in size_levels) row[[size]] <- summarize(values[dat$tumor_size_group == size])
  row
}

summarize_baseline <- function(dat, variable, characteristic) {
  values <- dat[[variable]]
  if (!is.factor(values)) {
    return(row_by_size(dat, characteristic, "Median (IQR; range)", values, summarize_continuous))
  }
  map_dfr(levels(values), function(level) {
    row_by_size(dat, characteristic, level, values, function(x) summarize_category(x, level))
  })
}

baseline_table <- function(dat) {
  stage_characteristic <- "Disease stage (descriptive only)"
  bind_rows(
    summarize_baseline(dat, "age_years", "Age, years"),
    summarize_baseline(dat, "sex", "Sex"),
    summarize_baseline(dat, "smoking_history", "Smoking history"),
    summarize_baseline(dat, "pack_years", "Pack-years"),
    summarize_baseline(dat, "copd", "COPD"),
    summarize_baseline(dat, "oxygen_any", "Baseline oxygen"),
    summarize_baseline(dat, "pathology_detail", "Pathology"),
    summarize_baseline(dat, "pd_l1", "PD-L1"),
    summarize_baseline(dat, "systemic_therapy", "Concurrent or post-ablation systemic therapy"),
    summarize_baseline(dat, "tumor_size_cm", "Tumor size, cm"),
    summarize_baseline(dat, "tumor_size_group", "Tumor-size group"),
    summarize_baseline(dat, "central_location", "Tumor location"),
    map_dfr(1:4, function(stage) {
      row_by_size(
        dat, stage_characteristic, paste("Stage", as.character(as.roman(stage))),
        dat$stage_raw, function(x) summarize_category(x, stage)
      )
    }),
    row_by_size(
      dat, stage_characteristic, "Missing", dat$stage_raw,
      function(x) format_fraction(sum(is.na(x)), length(x))
    )
  )
}

procedure_safety_table <- function(dat) {
  repeat_after_progression <- sum(dat$repeat_ablation == 1 & dat$local_progression_corrected == 1)
  tribble(
    ~outcome, ~category, ~result,
    "Probe count", "Median (IQR; range)", summarize_continuous(dat$probe_count),
    "Technical success", "Complete ice-ball coverage with 1-cm margin and no premature termination",
    summarize_category(dat$technical_success, "Yes"),
    "Hospital stay, days", "Median (IQR; range)", summarize_continuous(dat$hospital_stay_days),
    "Pneumothorax", "Yes", summarize_events(dat$pneumothorax),
    "Hemothorax/hemorrhage", "Yes", summarize_events(dat$hemothorax_hemorrhage),
    "Pleural effusion", "Yes", summarize_events(dat$pleural_effusion),
    "Chest-tube intervention", "Yes", summarize_events(dat$chest_tube),
    "Any recorded complication", "Yes", summarize_events(dat$any_recorded_complication),
    "SIR grade", "0", summarize_category(dat$sir_comp_grade, 0),
    "SIR grade", "1", summarize_category(dat$sir_comp_grade, 1),
    "SIR grade", "2", summarize_category(dat$sir_comp_grade, 2),
    "SIR grade", "3", summarize_category(dat$sir_comp_grade, 3),
    "Major complication", "SIR grade ≥3", summarize_events(dat$major_complication),
    "30-day all-cause mortality", "Yes", summarize_events(dat$death_30d),
    "90-day all-cause mortality", "Yes", summarize_events(dat$death_90d),
    "Procedure-related mortality", "Causality", "Not available in source data",
    "Repeat ablation", "Overall", summarize_events(dat$repeat_ablation),
    "Repeat ablation", "Among corrected local progressions",
    format_fraction(repeat_after_progression, sum(dat$local_progression_corrected))
  )
}

run_descriptive_tables <- function(dat, raw, out_dir) {
  # Table 1 was reported with sex coded the opposite way. Still unresolved.
  dat$sex <- factor(raw$sex, c(0, 1), c("Female", "Male"))
  write_table(baseline_table(dat), file.path(out_dir, "tables", "table1_baseline_characteristics.csv"))
  write_table(procedure_safety_table(dat), file.path(out_dir, "tables", "table2_procedure_safety.csv"))
}
