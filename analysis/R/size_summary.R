gray_tests <- function(dat, grouping) {
  tests <- cuminc(dat$cr_time_months, dat$cr_status, group = dat[[grouping]], cencode = 0)$Tests
  as.data.frame(tests) %>%
    rownames_to_column("cause") %>%
    transmute(
      grouping = grouping,
      endpoint = if_else(cause == "1", "LTP", "Competing death"),
      test = "Gray",
      statistic = stat,
      df,
      p_value = pv
    )
}

logrank_test <- function(dat, grouping, endpoint) {
  fit <- survdiff(reformulate(grouping, response = survival_response(endpoint)), data = dat)
  df <- length(fit$n) - 1
  tibble(
    grouping = grouping,
    endpoint = endpoint,
    test = "Log-rank",
    statistic = fit$chisq,
    df = df,
    p_value = pchisq(fit$chisq, df, lower.tail = FALSE)
  )
}

run_size_summaries <- function(dat, out_dir) {
  by_size <- map(set_names(levels(dat$tumor_size_group)), ~ filter(dat, tumor_size_group == .x))
  write_table(
    imap_dfr(by_size, ~ cif_estimates(.x, "Full107", .y)),
    file.path(out_dir, "tables", "size_cif_milestones.csv"),
    na = "NA"
  )
  write_table(
    imap_dfr(by_size, function(data, size) {
      map_dfr(c("OS", "Composite"), ~ km_estimates(data, "Full107", .x, size)$milestones)
    }),
    file.path(out_dir, "tables", "size_km_milestones.csv"),
    na = "NA"
  )
  group_tests <- map_dfr(c("tumor_size_group", "stage_group"), function(grouping) {
    bind_rows(
      gray_tests(dat, grouping),
      logrank_test(dat, grouping, "OS"),
      logrank_test(dat, grouping, "Composite")
    )
  })
  write_table(group_tests, file.path(out_dir, "model_results", "group_tests.csv"), na = "NA")
}
