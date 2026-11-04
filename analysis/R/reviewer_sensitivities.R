# The source has no center column. These row blocks were inferred from workbook fill
# colours and haven't been confirmed, so treat center results as exploratory.
read_center_assignments <- function(path, n_records) {
  blocks <- read.csv(path, stringsAsFactors = FALSE, encoding = "UTF-8")
  if (!identical(names(blocks), c("center", "first_record", "last_record"))) {
    stop("Unexpected columns in ", path)
  }
  blocks <- blocks[order(blocks$first_record), ]
  covered <- unlist(Map(seq, blocks$first_record, blocks$last_record))
  if (!identical(as.integer(covered), seq_len(n_records))) {
    stop("Center blocks must cover every record exactly once: ", path)
  }
  factor(rep(blocks$center, blocks$last_record - blocks$first_record + 1L))
}

leave_one_center_out <- function(dat, centers) {
  map_dfr(levels(centers), function(center) {
    remaining <- dat[centers != center, ]
    bind_rows(
      fit_fine_gray(remaining, "Fine-Gray LTP")$table,
      fit_cox(remaining, "Cox OS", "OS")$table
    ) %>%
      filter(term %in% c("tumor_size_cm", "central_locationCentral")) %>%
      transmute(
        term, estimate, lower, upper,
        p = p_value, model,
        omitted_center = center, n = nrow(remaining)
      )
  })
}

center_counts <- function(dat, centers) {
  map_dfr(levels(centers), function(center) {
    records <- dat[centers == center, ]
    data.frame(
      center = center,
      n = nrow(records),
      ltp = sum(records$cr_status == 1),
      death_before_ltp = sum(records$cr_status == 2),
      deaths = sum(records$os_event),
      pneumothorax = sum(records$pneumothorax),
      chest_tube = sum(records$chest_tube),
      recorded_grade_ge3 = sum(records$major_complication)
    )
  })
}

firth_probe_sensitivity <- function(safety) {
  map_dfr(names(safety_outcomes), function(outcome) {
    map_dfr(c(FALSE, TRUE), function(probe_adjusted) {
      estimates <- if (probe_adjusted) safety$adjusted[[outcome]] else safety$probe_omitted[[outcome]]
      estimates %>%
        filter(grepl("^tumor_size_group", term)) %>%
        transmute(
          outcome = .env$outcome, probe_adjusted = .env$probe_adjusted,
          term, OR = estimate, lower, upper, p = p_value, global_p = global_size_p
        )
    })
  })
}

# Death is censored here, so this isn't 1 minus the LTP CIF.
death_censored_local_control <- function(dat) {
  fit <- survfit(Surv(cr_time_months, cr_status == 1) ~ 1, data = dat, conf.type = "log-log")
  at_times <- summary(fit, times = milestone_months)
  data.frame(
    month = at_times$time, estimate = at_times$surv,
    lower = at_times$lower, upper = at_times$upper, n_risk = at_times$n.risk
  )
}

run_reviewer_sensitivities <- function(dat, safety, center_file, out_dir) {
  reviewer_dir <- file.path(out_dir, "reviewer_sensitivities")
  dir.create(reviewer_dir, showWarnings = FALSE)
  write_base_csv(firth_probe_sensitivity(safety), file.path(reviewer_dir, "firth_probe_sensitivity.csv"))
  write_base_csv(death_censored_local_control(dat), file.path(reviewer_dir, "death_censored_local_control.csv"))
  if (!file.exists(center_file)) {
    message("Center assignments are not available; skipping the center analyses.")
    return(invisible(NULL))
  }
  centers <- read_center_assignments(center_file, nrow(dat))
  write_base_csv(leave_one_center_out(dat, centers), file.path(reviewer_dir, "leave_one_center_out.csv"))
  write_base_csv(center_counts(dat, centers), file.path(reviewer_dir, "center_counts.csv"))
}
