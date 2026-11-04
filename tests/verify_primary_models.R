#!/usr/bin/env Rscript
# Refit the primary models from the source file without using any analysis code.
script_arg <- gsub("~+~", " ", grep("^--file=", commandArgs(), value = TRUE), fixed = TRUE)
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(project_dir, "analysis", "cli.R"))
opts <- parse_cli(
  "Usage: Rscript --vanilla tests/verify_primary_models.R [--input SOURCE] [--results DIRECTORY] [--library DIRECTORY]",
  options = list(
    input = default_source(project_dir),
    results = file.path(project_dir, "outputs/current"),
    library = file.path(project_dir, "renv", "library")
  )
)
configure_locale()
use_project_library(project_dir, opts$library)
suppressPackageStartupMessages(library(survival))

input_file <- opts$input
result_dir <- opts$results
source_hash <- digest::digest(file = input_file, algo = "sha256")
source_manifest <- jsonlite::read_json(file.path(result_dir, "source_manifest.json"))
stopifnot(identical(source_hash, source_manifest$source_sha256))

raw <- if (grepl("\\.xlsx$", input_file, ignore.case = TRUE)) {
  as.data.frame(readxl::read_excel(input_file, sheet = "Sheet2"))
} else {
  read.csv(input_file)
}
progression <- raw$local_recurrence == 1 | !is.na(raw$time_to_recurrence_days)
death <- raw$death_status == 1
first_status <- ifelse(
  progression & (!death | raw$time_to_recurrence_days <= raw$time_to_death_days), 1,
  ifelse(death, 2, 0)
)
followup <- pmax(
  raw$time_to_lastfu_days, ifelse(progression, raw$time_to_recurrence_days, NA),
  ifelse(death, raw$time_to_death_days, NA),
  na.rm = TRUE
)
days_per_month <- 365.25 / 12
model_data <- data.frame(
  diameter = raw$tumor_size_cm,
  age10 = raw$age_at_procedure / 10,
  central = raw$tumor_location_cp,
  oxygen = as.integer(raw$o2_requirement_baseline > 0),
  first_time = ifelse(
    first_status == 1, raw$time_to_recurrence_days,
    ifelse(first_status == 2, raw$time_to_death_days, followup)
  ) / days_per_month,
  status = first_status,
  death_time = ifelse(death, raw$time_to_death_days, followup) / days_per_month,
  death = as.integer(death)
)

expected <- read.csv(file.path(project_dir, "tests", "expected", "manuscript_counts.csv"))
expected_count <- function(statistic) expected$numerator[expected$statistic == statistic]
stopifnot(
  nrow(model_data) == expected_count("patients"),
  sum(first_status == 1) == expected_count("local_progressions"),
  sum(first_status == 2) == expected_count("competing_deaths"),
  sum(death) == expected_count("all_cause_deaths")
)

models <- list(
  `Full107 Fine-Gray` = cmprsk::crr(
    model_data$first_time, model_data$status,
    cov1 = as.matrix(model_data[, c("diameter", "age10", "central", "oxygen")]),
    failcode = 1, cencode = 0
  ),
  `Full107 Composite Cox` = coxph(
    Surv(first_time, status != 0) ~ diameter + age10 + central + oxygen,
    data = model_data, ties = "efron"
  ),
  `Full107 OS Cox` = coxph(
    Surv(death_time, death) ~ diameter + age10 + central + oxygen,
    data = model_data, ties = "efron"
  )
)

production <- read.csv(file.path(result_dir, "model_results/main_continuous_models.csv"))
terms <- c("tumor_size_cm", "age_per_10y", "central_locationCentral", "oxygen_anyAny home oxygen")
checks <- lapply(names(models), function(name) {
  fit <- models[[name]]
  coefficients <- if (inherits(fit, "crr")) fit$coef else coef(fit)
  standard_errors <- sqrt(diag(if (inherits(fit, "crr")) fit$var else vcov(fit)))
  direct <- cbind(
    estimate = exp(coefficients),
    lower = exp(coefficients - qnorm(0.975) * standard_errors),
    upper = exp(coefficients + qnorm(0.975) * standard_errors),
    p_value = 2 * pnorm(-abs(coefficients / standard_errors))
  )
  reference <- production[production$model == name, ]
  reference <- reference[match(terms, reference$term), ]
  differences <- abs(direct - as.matrix(reference[, colnames(direct)]))
  if (anyNA(differences) || any(differences > 1e-12)) stop("Direct-source model differs: ", name)
  data.frame(model = name, term = terms, maximum_absolute_difference = apply(differences, 1, max))
})
stopifnot(identical(digest::digest(file = input_file, algo = "sha256"), source_hash))
dir.create(file.path(result_dir, "validation"), showWarnings = FALSE)
write.csv(do.call(rbind, checks), file.path(result_dir, "validation/primary_model_check.csv"), row.names = FALSE)
cat("Primary model refits match.\n")
