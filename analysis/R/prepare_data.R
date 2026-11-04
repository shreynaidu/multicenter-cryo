required_source_columns <- c(
  "age_at_procedure", "sex", "smoking_history", "pack_years", "copd_status",
  "o2_requirement_baseline", "pathology_result", "pd_l1_status", "disease_stage",
  "systemic_therapy", "tumor_size_cm", "tumor_location_cp", "probe_count",
  "technical_success", "hospital_stay_days", "comp_pneumothorax", "comp_hemothorax",
  "comp_effusion", "intervention_chest_tube", "cirse_comp_grade", "local_recurrence",
  "repeat_ablation", "death_status", "time_to_recurrence_days", "time_to_death_days",
  "time_to_lastfu_days"
)

binary_source_columns <- c(
  "sex", "smoking_history", "copd_status", "systemic_therapy", "tumor_location_cp",
  "technical_success", "comp_pneumothorax", "comp_hemothorax", "comp_effusion",
  "intervention_chest_tube", "local_recurrence", "repeat_ablation", "death_status"
)

numeric_source_columns <- c(
  "age_at_procedure", "pack_years", "tumor_size_cm", "probe_count", "hospital_stay_days",
  "time_to_recurrence_days", "time_to_death_days", "time_to_lastfu_days"
)

coded_source_columns <- list(
  o2_requirement_baseline = 0:3,
  pathology_result = c(1:4, 99),
  pd_l1_status = c(0, 1, 99)
)

# No source IDs or dates here.
analysis_columns <- c(
  "record_id", "age_years", "age_per_10y", "sex", "smoking_history", "pack_years",
  "pack_years_per_10", "copd", "oxygen_category", "oxygen_any", "pathology_original",
  "pathology_analysis", "pd_l1", "systemic_therapy", "tumor_size_cm", "tumor_size_group",
  "central_location", "probe_count", "technical_success", "hospital_stay_days", "pneumothorax",
  "hemothorax_hemorrhage", "pleural_effusion", "chest_tube", "sir_comp_grade",
  "major_complication", "any_recorded_complication", "local_progression_original",
  "local_progression_corrected", "local_progression_correction_flag", "repeat_ablation",
  "death_status", "supplied_recurrence_days", "supplied_death_days", "supplied_followup_days",
  "corrected_followup_days", "administrative_followup_flag", "os_time_days", "os_time_months",
  "os_event", "cr_time_days", "cr_time_months", "cr_status", "lpfs_time_days",
  "lpfs_time_months", "lpfs_event", "death_30d", "death_90d"
)

read_source_data <- function(path) {
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    return(readxl::read_excel(path, sheet = input_sheet, guess_max = 10000))
  }
  source_table <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_double()))
  tibble::as_tibble(as.data.frame(source_table))
}

read_expected_source_checks <- function(path) {
  expected <- read.csv(path, colClasses = "character", na.strings = character(), encoding = "UTF-8")
  if (!identical(names(expected), c("check", "expected", "note")) || anyDuplicated(expected$check)) {
    stop("Unexpected structure in ", path)
  }
  expected
}

invalid_source_values <- function(raw) {
  invalid_binary <- map_dfr(binary_source_columns, function(column) {
    values <- raw[[column]]
    invalid <- which(is.na(values) | !values %in% c(0, 1))
    tibble(
      record_id = invalid, variable = column,
      issue = "Missing or invalid binary code", value = as.character(values[invalid])
    )
  })
  invalid_numeric <- map_dfr(numeric_source_columns, function(column) {
    values <- raw[[column]]
    numbers <- suppressWarnings(as.numeric(values))
    invalid <- which(!is.na(values) & (!is.finite(numbers) | numbers < 0))
    tibble(
      record_id = invalid, variable = column,
      issue = "Negative, infinite, or nonnumeric value", value = as.character(values[invalid])
    )
  })
  bind_rows(invalid_binary, invalid_numeric)
}

validate_source <- function(raw, expected_rows, out_dir) {
  missing_columns <- setdiff(required_source_columns, names(raw))
  if (length(missing_columns)) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  if (nrow(raw) != expected_rows) {
    stop("Expected ", expected_rows, " rows; found ", nrow(raw))
  }
  invalid <- invalid_source_values(raw)
  if (nrow(invalid)) {
    write_table(invalid, file.path(out_dir, "analysis_data", "invalid_source_values.csv"))
    stop("Invalid source values; see analysis_data/invalid_source_values.csv.")
  }
  for (column in names(coded_source_columns)) {
    if (any(is.na(raw[[column]]) | !raw[[column]] %in% coded_source_columns[[column]])) {
      stop("Missing or invalid source code in ", column)
    }
  }
  if (any(!is.na(raw$disease_stage) & !raw$disease_stage %in% 1:4)) {
    stop("Invalid disease stage code.")
  }
  followup <- as.numeric(raw$time_to_lastfu_days)
  if (any(!is.finite(followup) | followup <= 0)) {
    stop("Every record requires a finite, positive supplied follow-up duration.")
  }
}

derive_cohort <- function(raw) {
  raw %>%
    mutate(
      # pt_id repeats across rows; each row was adjudicated as a separate patient.
      record_id = row_number(),
      age_years = as.numeric(age_at_procedure),
      age_per_10y = age_years / 10,
      sex = factor(sex, c(0, 1), c("Male", "Female")),
      smoking_history = factor(smoking_history, c(0, 1), c("Never", "Ever")),
      pack_years = as.numeric(pack_years),
      pack_years_per_10 = pack_years / 10,
      copd = factor(copd_status, c(0, 1), c("No", "Yes")),
      oxygen_category = factor(
        o2_requirement_baseline, c(0, 1, 2, 3), c("Room air", "2 L/min", "3 L/min", ">3 L/min")
      ),
      oxygen_any = factor(
        if_else(o2_requirement_baseline == 0, 0, 1), c(0, 1), c("Room air", "Any home oxygen")
      ),
      pathology_original = as.integer(pathology_result),
      pathology_detail = factor(
        pathology_result, c(1, 2, 3, 4, 99),
        c("Adenocarcinoma", "Squamous", "Adenosquamous", "Other", "Unknown/not available")
      ),
      pathology_analysis = factor(
        case_when(
          pathology_result == 1 ~ "Adenocarcinoma",
          pathology_result == 2 ~ "Squamous",
          TRUE ~ "Other/unknown"
        ),
        c("Adenocarcinoma", "Squamous", "Other/unknown")
      ),
      pd_l1 = factor(pd_l1_status, c(0, 1, 99), c("Negative", "Positive", "Unknown")),
      systemic_therapy = factor(systemic_therapy, c(0, 1), c("No", "Yes")),
      tumor_size_cm = as.numeric(tumor_size_cm),
      tumor_size_group = factor(
        cut(tumor_size_cm, c(-Inf, 2, 3, 4, Inf), right = FALSE, labels = size_levels),
        levels = size_levels
      ),
      central_location = factor(tumor_location_cp, c(0, 1), c("Peripheral", "Central")),
      # Missing stage is its own category, not imputed.
      stage_raw = as.integer(disease_stage),
      stage_group = factor(
        case_when(
          stage_raw %in% c(1, 2) ~ "Stage I-II",
          stage_raw %in% c(3, 4) ~ "Stage III-IV",
          TRUE ~ "Unknown stage"
        ),
        levels = stage_levels
      ),
      probe_count = as.numeric(probe_count),
      technical_success = factor(technical_success, c(0, 1), c("No", "Yes")),
      hospital_stay_days = as.numeric(hospital_stay_days),
      pneumothorax = as.integer(comp_pneumothorax),
      hemothorax_hemorrhage = as.integer(comp_hemothorax),
      pleural_effusion = as.integer(comp_effusion),
      chest_tube = as.integer(intervention_chest_tube),
      # Column is labelled CIRSE; investigators adjudicated these as SIR grades.
      sir_comp_grade = as.integer(cirse_comp_grade),
      major_complication = as.integer(sir_comp_grade >= 3),
      any_recorded_complication = as.integer(
        pneumothorax == 1 | hemothorax_hemorrhage == 1 | pleural_effusion == 1 | sir_comp_grade > 0
      ),
      # Investigators count any recorded recurrence time as progression, even with status 0.
      local_progression_original = as.integer(local_recurrence),
      local_progression_correction_flag = as.integer(
        local_recurrence == 0 & !is.na(time_to_recurrence_days)
      ),
      local_progression_corrected = as.integer(
        local_recurrence == 1 | !is.na(time_to_recurrence_days)
      ),
      repeat_ablation = as.integer(repeat_ablation),
      death_status = as.integer(death_status),
      supplied_recurrence_days = as.numeric(time_to_recurrence_days),
      supplied_death_days = as.numeric(time_to_death_days),
      supplied_followup_days = as.numeric(time_to_lastfu_days),
      # Always 0 now that the workbook has every duration. Kept as an export column.
      administrative_followup_flag = 0L,
      # Follow-up can't end before a documented event.
      corrected_followup_days = pmax(
        supplied_followup_days,
        if_else(local_progression_corrected == 1L, supplied_recurrence_days, NA_real_),
        if_else(death_status == 1L, supplied_death_days, NA_real_),
        na.rm = TRUE
      ),
      recurrence_after_death_flag = as.integer(
        local_progression_corrected == 1L & death_status == 1L &
          supplied_recurrence_days > supplied_death_days
      ),
      # Same-day progression and death counts as progression.
      cr_status = case_when(
        local_progression_corrected == 1L &
          (death_status == 0L | supplied_recurrence_days <= supplied_death_days) ~ 1L,
        death_status == 1L ~ 2L,
        TRUE ~ 0L
      ),
      cr_time_days = case_when(
        cr_status == 1L ~ supplied_recurrence_days,
        cr_status == 2L ~ supplied_death_days,
        TRUE ~ corrected_followup_days
      ),
      cr_time_months = cr_time_days / days_per_month,
      lpfs_event = as.integer(cr_status != 0L),
      lpfs_time_days = cr_time_days,
      lpfs_time_months = cr_time_months,
      os_event = death_status,
      os_time_days = if_else(os_event == 1L, supplied_death_days, corrected_followup_days),
      os_time_months = os_time_days / days_per_month,
      death_30d = as.integer(death_status == 1L & supplied_death_days <= 30),
      death_90d = as.integer(death_status == 1L & supplied_death_days <= 90)
    )
}

assert_valid_endpoints <- function(dat) {
  if (any(is.na(dat$cr_time_days) | dat$cr_time_days <= 0 | is.na(dat$os_time_days) | dat$os_time_days <= 0)) {
    stop("Invalid endpoint time remains after adjudication.")
  }
  if (any(dat$local_progression_corrected == 1L & is.na(dat$supplied_recurrence_days))) {
    stop("Progression without recurrence time.")
  }
  if (any(dat$death_status == 1L & is.na(dat$supplied_death_days))) {
    stop("Death without death time.")
  }
}

correction_log <- function(dat) {
  bind_rows(
    dat %>%
      filter(local_progression_correction_flag == 1L) %>%
      transmute(
        record_id,
        variable = "local_progression",
        original_value = as.character(local_progression_original),
        corrected_value = as.character(local_progression_corrected),
        reason = "Nonmissing recurrence duration with original status 0; investigator-adjudicated recode"
      ),
    dat %>%
      filter(corrected_followup_days > supplied_followup_days) %>%
      transmute(
        record_id,
        variable = "analysis_followup_duration_days",
        original_value = as.character(supplied_followup_days),
        corrected_value = as.character(corrected_followup_days),
        reason = case_when(
          local_progression_corrected == 1L & supplied_recurrence_days > supplied_followup_days ~
            "Extended analytic observation through documented progression",
          death_status == 1L & supplied_death_days > supplied_followup_days ~
            "Extended analytic observation through documented death",
          TRUE ~ "Extended analytic observation through documented event"
        )
      ),
    tibble(
      record_id = NA_integer_,
      variable = "cirse_comp_grade",
      original_value = "Workbook column label",
      corrected_value = "sir_comp_grade",
      reason = "Investigator adjudicated supplied grade as SIR grade"
    )
  )
}

flag_records <- function(dat, condition, variable, issue, value) {
  selected <- which(condition)
  tibble(
    record_id = dat$record_id[selected],
    variable = variable,
    issue = issue,
    value = as.character(value[selected])
  )
}

# Flagged for review, not excluded.
quality_flags <- function(dat) {
  implausible <- function(condition, variable) {
    flag_records(dat, condition, variable, "Implausible value", dat[[variable]])
  }
  recurrence_after_followup <- dat$local_progression_corrected == 1L &
    dat$supplied_recurrence_days > dat$supplied_followup_days
  death_after_followup <- dat$death_status == 1L & dat$supplied_death_days > dat$supplied_followup_days
  bind_rows(
    implausible(dat$age_years < 18 | dat$age_years > 110, "age_years"),
    implausible(dat$tumor_size_cm <= 0 | dat$tumor_size_cm > 20, "tumor_size_cm"),
    implausible(dat$probe_count < 1 | dat$probe_count > 20, "probe_count"),
    implausible(dat$hospital_stay_days < 0 | dat$hospital_stay_days > 90, "hospital_stay_days"),
    flag_records(
      dat, dat$recurrence_after_death_flag == 1L, "endpoint_order", "Recurrence after death",
      paste(dat$supplied_recurrence_days, dat$supplied_death_days, sep = "/")
    ),
    flag_records(
      dat, recurrence_after_followup, "followup_order",
      "Supplied follow-up precedes documented recurrence; analytic follow-up extended",
      paste(dat$supplied_followup_days, dat$supplied_recurrence_days, sep = "/")
    ),
    flag_records(
      dat, death_after_followup, "followup_order",
      "Supplied imaging follow-up precedes documented death; endpoint uses death time",
      paste(dat$supplied_followup_days, dat$supplied_death_days, sep = "/")
    )
  )
}

observed_source_checks <- function(raw, dat) {
  duplicate_rows <- duplicated(as.data.frame(raw[, setdiff(names(raw), "pt_id")]))
  c(
    "Source rows" = nrow(raw),
    "Analytic patients" = nrow(dat),
    "Deaths" = sum(dat$death_status),
    "Originally coded local progressions" = sum(dat$local_progression_original),
    "Status-zero records with recurrence duration" = sum(dat$local_progression_correction_flag),
    "Corrected local progressions" = sum(dat$local_progression_corrected),
    "Administrative follow-up derivations" = sum(dat$administrative_followup_flag),
    "Tumor-size group counts" = paste(as.integer(table(dat$tumor_size_group)), collapse = "|"),
    "Technical successes" = sum(dat$technical_success == "Yes"),
    "Pneumothoraces" = sum(dat$pneumothorax),
    "Hemothorax/hemorrhage" = sum(dat$hemothorax_hemorrhage),
    "Pleural effusions" = sum(dat$pleural_effusion),
    "Chest tubes" = sum(dat$chest_tube),
    "Repeat ablations" = sum(dat$repeat_ablation),
    "Unknown pathology" = sum(dat$pathology_original == 99),
    "30-day deaths" = sum(dat$death_30d),
    "90-day deaths" = sum(dat$death_90d),
    "Recurrence after death" = sum(dat$recurrence_after_death_flag),
    "Exact duplicate rows excluding pt_id" = sum(duplicate_rows)
  )
}

readiness_report <- function(raw, dat, expected_checks) {
  calculated <- observed_source_checks(raw, dat)
  if (!setequal(names(calculated), expected_checks$check)) {
    stop("Expected source checks do not match those calculated.")
  }
  expected_checks %>%
    mutate(
      observed = unname(calculated[check]),
      status = ifelse(observed == expected, "PASS", "FAIL")
    ) %>%
    select(check, expected, observed, status, note)
}

prepare_cohort <- function(raw, expected_checks, out_dir) {
  expected_rows <- as.integer(expected_checks$expected[expected_checks$check == "Source rows"])
  validate_source(raw, expected_rows, out_dir)
  dat <- derive_cohort(raw)
  assert_valid_endpoints(dat)

  analysis_data_dir <- file.path(out_dir, "analysis_data")
  write_table(correction_log(dat), file.path(analysis_data_dir, "correction_log.csv"))
  write_table(
    dat[, analysis_columns],
    file.path(analysis_data_dir, "deidentified_corrected_analysis_dataset.csv")
  )
  write_table(quality_flags(dat), file.path(analysis_data_dir, "data_quality_flags.csv"))

  readiness <- readiness_report(raw, dat, expected_checks)
  write_table(readiness, file.path(analysis_data_dir, "data_readiness_report.csv"))
  if (any(readiness$status != "PASS")) {
    stop("Source checks failed; see analysis_data/data_readiness_report.csv.")
  }
  dat
}

analysis_cohort <- function(dat) {
  dat %>% select(all_of(analysis_columns), stage_raw, stage_group)
}
