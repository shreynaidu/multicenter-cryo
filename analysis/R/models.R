competing_risk_counts <- function(data) {
  tibble(
    n = nrow(data),
    local_failures = sum(data$cr_status == 1),
    competing_deaths = sum(data$cr_status == 2)
  )
}

write_cohort_summaries <- function(sets, out_dir) {
  counts <- imap_dfr(sets, function(data, cohort) {
    tibble(cohort = cohort, competing_risk_counts(data)) %>%
      mutate(
        composite_events = sum(data$lpfs_event),
        os_deaths = sum(data$os_event),
        technical_failures = sum(data$technical_success == "No"),
        early_local_failures = sum(data$cr_status == 1 & data$cr_time_months < 6),
        stage_I_II = sum(data$stage_group == "Stage I-II"),
        stage_III_IV = sum(data$stage_group == "Stage III-IV"),
        stage_unknown = sum(data$stage_group == "Unknown stage")
      )
  })
  write_table(counts, file.path(out_dir, "tables", "cohort_counts.csv"))
  write_table(imap_dfr(sets, ~ cif_estimates(.x, .y)), file.path(out_dir, "tables", "cif_all_cohorts.csv"))

  survival <- unlist(imap(sets, function(data, cohort) {
    map(c("Composite", "OS"), ~ km_estimates(data, cohort, .x))
  }), recursive = FALSE)
  write_table(map_dfr(survival, "milestones"), file.path(out_dir, "tables", "km_all_cohorts.csv"))
  write_table(map_dfr(survival, "median"), file.path(out_dir, "tables", "survival_medians.csv"))
}

# Recorded stage only. It doesn't tell us about isolated recurrence or treatment intent.
write_stage_summaries <- function(dat, out_dir) {
  cohort <- "Full107"
  by_stage <- map(levels(dat$stage_group), function(stage) {
    list(stage = stage, data = filter(dat, stage_group == stage))
  })
  stage_counts <- map_dfr(by_stage, function(group) {
    tibble(cohort = cohort, stage_group = group$stage, competing_risk_counts(group$data)) %>%
      mutate(os_deaths = sum(group$data$os_event))
  })
  write_table(stage_counts, file.path(out_dir, "tables", "stage_group_counts.csv"))
  write_table(
    map_dfr(by_stage, ~ cif_estimates(.x$data, cohort, .x$stage, times = 24)),
    file.path(out_dir, "tables", "stage_cif_24months.csv")
  )
  write_table(
    map_dfr(by_stage, ~ km_estimates(.x$data, cohort, "OS", .x$stage, times = 24)$milestones),
    file.path(out_dir, "tables", "stage_os_24months.csv")
  )
}

fit_primary_models <- function(dat) {
  list(
    FG107 = fit_fine_gray(dat, "Full107 Fine-Gray"),
    Composite107 = fit_cox(dat, "Full107 Composite Cox", "Composite"),
    OS107 = fit_cox(dat, "Full107 OS Cox", "OS")
  )
}

fit_continuous_models <- function(sets, primary) {
  fits <- imap(sets, function(data, cohort) {
    full_cohort <- cohort == "Full107"
    fine_gray <- if (full_cohort) primary$FG107 else fit_fine_gray(data, paste(cohort, "Fine-Gray"))
    composite <- if (full_cohort) primary$Composite107 else fit_cox(data, paste(cohort, "Composite Cox"), "Composite")
    overall <- if (full_cohort) primary$OS107 else fit_cox(data, paste(cohort, "OS Cox"), "OS")
    local_failure <- fit_cox(data, paste(cohort, "LocalFailure Cox"), "LocalFailure")
    list(fine_gray, composite, overall, local_failure)
  })
  unlist(fits, recursive = FALSE)
}

write_categorical_models <- function(dat, out_dir) {
  covariates <- c("tumor_size_group", setdiff(model_covariates, "tumor_size_cm"))
  fine_gray <- fit_fine_gray(dat, "Full107 categorical Fine-Gray", covariates)
  size_terms <- grep("tumor_size_group", names(fine_gray$fit$coef))
  cox <- map(c("Composite", "OS"), function(endpoint) {
    full <- fit_cox(dat, paste("Full107 categorical", endpoint), endpoint, covariates)
    reduced <- fit_cox(dat, paste("Full107 reduced", endpoint), endpoint, setdiff(covariates, "tumor_size_group"))
    likelihood_ratio <- 2 * (as.numeric(logLik(full$fit)) - as.numeric(logLik(reduced$fit)))
    list(
      table = full$table,
      global = tibble(
        cohort = "Full107",
        endpoint = endpoint,
        global_size_p = pchisq(likelihood_ratio, 3, lower.tail = FALSE)
      )
    )
  })
  write_table(
    bind_rows(fine_gray$table, map_dfr(cox, "table")),
    file.path(out_dir, "model_results", "categorical_models_with_oxygen.csv")
  )
  fine_gray_global <- tibble(
    cohort = "Full107",
    endpoint = "Fine-Gray",
    global_size_p = wald_p_value(
      fine_gray$fit$coef[size_terms],
      fine_gray$fit$var[size_terms, size_terms, drop = FALSE]
    )
  )
  write_table(
    bind_rows(fine_gray_global, map_dfr(cox, "global")),
    file.path(out_dir, "model_results", "categorical_global_tests.csv")
  )
}

# Unknown stage as its own category, then known-stage records only.
fit_stage_adjusted_models <- function(dat) {
  known_stage <- droplevels(filter(dat, !is.na(stage_raw)))
  covariates <- c(model_covariates, "stage_group")
  fits <- map(c("OS", "Composite"), function(endpoint) {
    list(
      fit_cox(dat, paste(endpoint, "stage sensitivity: Unknown category"), endpoint, covariates),
      fit_cox(known_stage, paste(endpoint, "stage sensitivity: Known-stage only"), endpoint, covariates)
    )
  })
  unlist(fits, recursive = FALSE)
}

# Too few events for more than three covariates.
write_competing_death_model <- function(dat, out_dir) {
  fit <- fit_cox(
    dat, "Full107 exploratory competing-death Cox", "CompetingDeath",
    setdiff(model_covariates, "tumor_size_cm")
  )
  write_table(fit$table, file.path(out_dir, "model_results", "competing_death_exploratory.csv"))
  write_table(fit$ph, file.path(out_dir, "diagnostics", "competing_death_ph_tests.csv"))
}

# PH check: one covariate-by-log(time) term at a time.
fine_gray_time_interactions <- function(dat, primary_fine_gray) {
  design <- primary_fine_gray$mm
  map_dfr(seq_len(ncol(design)), function(column) {
    fit <- tryCatch(
      cmprsk::crr(
        dat$cr_time_months, dat$cr_status,
        cov1 = design,
        cov2 = design[, column, drop = FALSE],
        tf = function(t) matrix(log(pmax(t, 0.001)), ncol = 1),
        failcode = 1,
        cencode = 0
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(tibble(term = colnames(design)[column], p_value = NA_real_, converged = FALSE))
    }
    last <- length(fit$coef)
    tibble(
      term = colnames(design)[column],
      p_value = 2 * pnorm(abs(fit$coef[last] / sqrt(fit$var[last, last])), lower.tail = FALSE),
      converged = isTRUE(fit$converged)
    )
  })
}

location_leave_one_out <- function(dat) {
  map_dfr(c("FG", "Composite", "OS"), function(endpoint) {
    map_dfr(seq_len(nrow(dat)), function(record) {
      remaining <- dat[-record, ]
      fit <- tryCatch(
        if (endpoint == "FG") fit_fine_gray(remaining, "LOO") else fit_cox(remaining, "LOO", endpoint),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(tibble(
          endpoint = endpoint, omitted_record = record, estimate = NA_real_, p_value = NA_real_, converged = FALSE
        ))
      }
      fit$table %>%
        filter(term == "central_locationCentral") %>%
        transmute(endpoint = .env$endpoint, omitted_record = record, estimate, p_value, converged)
    })
  })
}

write_location_influence <- function(dat, out_dir) {
  influence <- location_leave_one_out(dat)
  write_table(influence, file.path(out_dir, "diagnostics", "location_leave_one_out.csv"))
  summary <- influence %>%
    filter(converged) %>%
    group_by(endpoint) %>%
    summarise(
      successful = n(),
      min_hr = min(estimate),
      max_hr = max(estimate),
      min_p = min(p_value),
      max_p = max(p_value),
      omissions_p_ge_05 = sum(p_value >= 0.05),
      .groups = "drop"
    )
  write_table(summary, file.path(out_dir, "diagnostics", "location_leave_one_out_summary.csv"))
}

write_pack_year_summary <- function(dat, out_dir) {
  populations <- list(All107 = dat$pack_years, EverSmokers94 = dat$pack_years[dat$smoking_history == "Ever"])
  summary <- imap_dfr(populations, function(values, population) {
    tibble(
      population = population,
      n_nonmissing = sum(!is.na(values)),
      median = median(values),
      q1 = sample_quantile(values, 0.25),
      q3 = sample_quantile(values, 0.75),
      minimum = min(values),
      maximum = max(values),
      quantile_method = "R quantile type 2; median uses standard median"
    )
  })
  write_table(summary, file.path(out_dir, "tables", "pack_years_denominators.csv"))
}

# Reverse KM: death counts as censoring.
write_followup_summary <- function(sets, out_dir) {
  summary <- imap_dfr(sets, function(data, cohort) {
    reverse_km <- survfit(Surv(os_time_months, 1 - os_event) ~ 1, data = data, conf.type = "log-log")
    overall <- summary(reverse_km)$table
    quartiles <- quantile(reverse_km, probs = c(0.25, 0.5, 0.75))$quantile
    observed <- data$corrected_followup_days / days_per_month
    tibble(
      cohort = cohort,
      reverse_km_q1 = unname(quartiles[1]),
      reverse_km_median = unname(overall["median"]),
      reverse_km_q3 = unname(quartiles[3]),
      reverse_km_lower = unname(overall["0.95LCL"]),
      reverse_km_upper = unname(overall["0.95UCL"]),
      observed_followup_median = median(observed),
      observed_followup_q1 = sample_quantile(observed, 0.25),
      observed_followup_q3 = sample_quantile(observed, 0.75),
      observed_followup_min = min(observed),
      observed_followup_max = max(observed)
    )
  })
  write_table(summary, file.path(out_dir, "tables", "followup_summary.csv"))
}

write_overlap_tables <- function(dat, out_dir) {
  write_table(
    count(dat, sir_comp_grade, chest_tube, pneumothorax, hemothorax_hemorrhage, pleural_effusion, name = "n"),
    file.path(out_dir, "tables", "recorded_grade_tube_event_overlap.csv")
  )
  technical <- dat %>%
    group_by(technical_success, administrative_followup_flag) %>%
    summarise(
      n = n(),
      documented_local_failures = sum(cr_status == 1),
      competing_deaths = sum(cr_status == 2),
      os_deaths = sum(os_event),
      deaths_30d = sum(death_30d),
      deaths_90d = sum(death_90d),
      .groups = "drop"
    )
  write_table(technical, file.path(out_dir, "tables", "technical_followup_overlap.csv"))
  write_table(
    count(dat, stage_group, systemic_therapy, name = "n"),
    file.path(out_dir, "tables", "stage_systemic_therapy_overlap.csv"),
    na = "NA"
  )
}

# Clear formula environments or the raw data gets saved along with the models.
save_analysis_objects <- function(dat, sets, primary, out_dir) {
  fitted_models <- map(primary, function(model) {
    fit <- model$fit
    if (inherits(fit, "coxph")) {
      attr(fit$terms, ".Environment") <- baseenv()
      if (!is.null(fit$formula)) environment(fit$formula) <- baseenv()
    }
    fit
  })
  saveRDS(
    list(data = dat, cohorts = sets, main_models = fitted_models),
    file.path(out_dir, "analysis_objects.rds")
  )
}

run_models <- function(dat, out_dir) {
  sets <- list(Full107 = dat, Success103 = filter(dat, technical_success == "Yes"))
  write_cohort_summaries(sets, out_dir)
  write_stage_summaries(dat, out_dir)

  primary <- fit_primary_models(dat)
  write_table(
    map_dfr(primary, "table"),
    file.path(out_dir, "model_results", "main_continuous_models.csv")
  )
  continuous <- fit_continuous_models(sets, primary)
  write_table(
    map_dfr(continuous, "table"),
    file.path(out_dir, "model_results", "continuous_models_all_cohorts.csv")
  )
  write_categorical_models(dat, out_dir)

  stage_adjusted <- fit_stage_adjusted_models(dat)
  write_table(
    map_dfr(stage_adjusted, "table"),
    file.path(out_dir, "model_results", "stage_adjusted_sensitivity.csv")
  )
  schoenfeld_models <- c(
    Filter(function(model) !is.null(model$ph), continuous), stage_adjusted, primary[c("Composite107", "OS107")]
  )
  write_table(
    map_dfr(schoenfeld_models, "ph"),
    file.path(out_dir, "diagnostics", "cox_ph_tests.csv")
  )
  write_competing_death_model(dat, out_dir)
  write_table(
    fine_gray_time_interactions(dat, primary$FG107),
    file.path(out_dir, "diagnostics", "primary_fine_gray_time_interactions.csv")
  )
  write_location_influence(dat, out_dir)
  write_pack_year_summary(dat, out_dir)
  write_followup_summary(sets, out_dir)
  write_overlap_tables(dat, out_dir)

  reanalysis <- dat %>%
    transmute(record_id, stage_raw, stage_group, across(all_of(analysis_columns[-1])))
  write_table(reanalysis, file.path(out_dir, "analysis_data", "reanalysis_dataset.csv"))
  save_analysis_objects(dat, sets, primary, out_dir)
  primary
}
