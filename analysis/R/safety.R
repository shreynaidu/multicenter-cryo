# Probe count may sit between size and complications, so fit with and without it.

safety_outcomes <- c(
  pneumothorax = "Pneumothorax",
  chest_tube = "Chest-tube intervention",
  major_complication = "SIR grade ≥3 adverse event"
)

run_safety_models <- function(dat, out_dir) {
  outcomes <- set_names(names(safety_outcomes))
  adjusted <- map(outcomes, ~ fit_firth_safety(dat, .x))
  probe_omitted <- map(outcomes, ~ fit_firth_safety(dat, .x, include_probe = FALSE))

  summary <- imap_dfr(adjusted, function(estimates, outcome) {
    estimates %>%
      transmute(
        outcome = unname(safety_outcomes[[.env$outcome]]), term, estimate, lower, upper, p_value,
        model_n = nrow(dat), events = sum(dat[[.env$outcome]]), global_size_p
      )
  })
  write_table(summary, file.path(out_dir, "model_results", "safety_models.csv"))

  dir.create(file.path(out_dir, "supplemental"), showWarnings = FALSE)
  omitted_summary <- imap_dfr(probe_omitted, function(estimates, outcome) {
    estimates %>%
      transmute(outcome = .env$outcome, term, OR = estimate, lower, upper, p = p_value, global_p = global_size_p)
  })
  write_base_csv(omitted_summary, file.path(out_dir, "supplemental", "refitted_probe_omitted_models.csv"))
  list(adjusted = adjusted, probe_omitted = probe_omitted)
}
