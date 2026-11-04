#!/usr/bin/env Rscript
# Check the Python figures' numbers against R. Calculated labels have to agree. Values
# kept from the supplied figures are only reported, unless --strict is set.
script_arg <- gsub("~+~", " ", grep("^--file=", commandArgs(), value = TRUE), fixed = TRUE)
project_dir <- dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg)))))
source(file.path(project_dir, "analysis", "cli.R"))
opts <- parse_cli(
  paste(
    "Usage: Rscript --vanilla analysis/scripts/audit_python_figure_statistics.R",
    "[--figures DIRECTORY] [--analysis DIRECTORY] [--library DIRECTORY] [--strict]"
  ),
  options = list(
    figures = file.path(project_dir, "outputs/python_publication_figures"),
    analysis = file.path(project_dir, "outputs/current"),
    library = file.path(project_dir, "renv", "library")
  ),
  switches = "--strict"
)
configure_locale()
use_project_library(project_dir, opts$library)
source(file.path(project_dir, "analysis", "config.R"))
suppressPackageStartupMessages({
  library(survival)
  library(cmprsk)
})

agreement_tolerance <- 1e-12
figure_dir <- opts$figures
numerical_file <- function(name) file.path(figure_dir, "numerical", name)

records <- read.csv(numerical_file("figure_analysis_dataset.csv"), check.names = FALSE)
records$tumor_size_group <- factor(records$tumor_size_group, levels = size_levels)
annotations <- jsonlite::read_json(numerical_file("computed_annotations.json"), simplifyVector = TRUE)

agreement_check <- function(check, figure_value, independent_value, disagreement_status, note) {
  difference <- abs(figure_value - independent_value)
  data.frame(
    check = check, figure_value = figure_value, independent_value = independent_value,
    absolute_difference = difference,
    status = ifelse(difference <= agreement_tolerance, "PASS", disagreement_status),
    note = note
  )
}

calculated_label_note <- "Calculated from the actual retained Python plotting data"

logrank_checks <- lapply(c("os", "lpfs"), function(endpoint) {
  fit <- survdiff(
    reformulate("tumor_size_group", paste0("Surv(", endpoint, "_time_months,", endpoint, "_event)")),
    data = records
  )
  agreement_check(
    paste(endpoint, "log-rank P"), annotations$logrank_p[[endpoint]],
    pchisq(fit$chisq, 3, lower.tail = FALSE), "DIFFERENCE", calculated_label_note
  )
})

os_fit <- survfit(Surv(os_time_months, os_event) ~ 1, data = records, conf.type = "log-log")
median_check <- agreement_check(
  c("OS median", "OS median lower CI", "OS median upper CI"),
  unlist(annotations$median_os_months),
  unname(summary(os_fit)$table[c("median", "0.95LCL", "0.95UCL")]),
  "DIFFERENCE", calculated_label_note
)

gray_p <- cuminc(records$cr_time_months, records$cr_status, group = records$tumor_size_group)$Tests[1, "pv"]
reported_gray <- annotations$retained_reported_gray_p
gray_check <- data.frame(
  check = "Reported Gray P", figure_value = reported_gray, independent_value = gray_p,
  absolute_difference = abs(reported_gray - gray_p),
  status = ifelse(sprintf("%.3f", reported_gray) == sprintf("%.3f", gray_p), "PASS", "OWNER_DECISION"),
  note = "Original reported label retained; it does not agree at displayed precision with the plotted data"
)

incidence_checks <- lapply(0:length(size_levels), function(index) {
  subgroup <- if (index == 0) records else records[records$tumor_size_group == size_levels[index], ]
  prefix <- if (index == 0) "overall" else paste0("size_", index)
  incidence <- cuminc(subgroup$cr_time_months, subgroup$cr_status, cencode = 0)
  at_times <- timepoints(incidence, milestone_months)
  do.call(rbind, lapply(1:2, function(cause) {
    curve <- read.csv(numerical_file(paste0(prefix, "_cause_", cause, "_cif.csv")))
    row <- grep(paste0(" ", cause, "$"), rownames(at_times$est))
    do.call(rbind, lapply(seq_along(milestone_months), function(column) {
      month <- milestone_months[column]
      plotted <- tail(curve[curve$time <= month, ], 1)
      agreement_check(
        paste(prefix, "cause", cause, month, c("estimate", "standard error")),
        c(plotted$estimate, plotted$se),
        c(at_times$est[row, column], sqrt(at_times$var[row, column])),
        "OWNER_DECISION",
        "Point estimates/SE checked independently; original log-log CIF limits remain unchanged"
      )
    }))
  }))
})

# The stage curves came from an image, so compare them with estimates from the records.
canonical <- read.csv(file.path(opts$analysis, "analysis_data/reanalysis_dataset.csv"))
stopifnot(!anyDuplicated(canonical$record_id), identical(records$record_id, canonical$record_id))
records$stage_group <- canonical$stage_group
reference_risks <- read.csv(file.path(project_dir, "analysis/figure_reference/retained_stage_risk_counts.csv"))

stage_checks <- lapply(stage_levels, function(stage) {
  subgroup <- records[records$stage_group == stage, ]
  retained <- read.csv(numerical_file(paste0(gsub(" ", "_", stage), "_stage_steps.csv")))
  fit <- cuminc(subgroup$cr_time_months, subgroup$cr_status, cencode = 0)
  estimated <- timepoints(fit, retained$time)$est["1 1", ] * 100
  difference <- max(abs(retained$estimate_pct - estimated), na.rm = TRUE)
  curve_check <- data.frame(
    check = paste(stage, "maximum curve-coordinate difference (percentage points)"),
    figure_value = NA_real_, independent_value = NA_real_, absolute_difference = difference,
    status = ifelse(difference <= 1e-10, "PASS", "OWNER_DECISION"),
    note = "Retained image-derived coordinates differ from patient-derived estimates"
  )
  independent_risks <- sapply(risk_table_months, function(month) sum(subgroup$cr_time_months >= month))
  retained_risks <- reference_risks$number_at_risk[reference_risks$stage == stage]
  risk_check <- data.frame(
    check = paste(stage, "risk count", risk_table_months), figure_value = retained_risks,
    independent_value = independent_risks, absolute_difference = abs(retained_risks - independent_risks),
    status = ifelse(retained_risks == independent_risks, "PASS", "OWNER_DECISION"),
    note = "Counts checked from retained plotting times and canonical stage membership"
  )
  rbind(curve_check, risk_check)
})

report <- do.call(rbind, c(logrank_checks, list(median_check, gray_check), incidence_checks, stage_checks))
write.csv(report, file.path(figure_dir, "statistical_consistency.csv"), row.names = FALSE)
summary <- list(
  reproducibility_is_separate_from_statistical_consistency = TRUE,
  checks = nrow(report),
  passed = sum(report$status == "PASS"),
  owner_decisions = sum(report$status == "OWNER_DECISION"),
  implementation_differences = sum(report$status == "DIFFERENCE"),
  reported_gray_p = reported_gray,
  plotted_data_gray_p = unname(gray_p)
)
jsonlite::write_json(
  summary, file.path(figure_dir, "statistical_consistency.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = 17
)
if (any(report$status == "DIFFERENCE")) stop("Calculated Python annotations do not match direct R calculations")
if (opts$strict && any(report$status != "PASS")) stop("Retained historical figures require statistical owner decisions")
cat(
  "Python figure audit:", summary$passed, "checks pass,", summary$owner_decisions,
  "retained differences.\n"
)
