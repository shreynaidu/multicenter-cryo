accepted_source_sha256 <- c(
  workbook = "9a30edb9325e20dd08a60424b56a95281ac5b9eacdcfe59cf59a3c4fb6c0930a",
  release = "47e5c158fe2270eac5e5da280f34fff734bfac376cbc930fce69116a14e0cff5"
)
input_sheet <- "Sheet2"

days_per_month <- 365.25 / 12
milestone_months <- c(12, 24, 36)
risk_table_months <- c(0, 12, 24, 36)

size_levels <- c("<2.0 cm", "2.0 to <3.0 cm", "3.0 to <4.0 cm", "≥4.0 cm")
stage_levels <- c("Stage I-II", "Stage III-IV", "Unknown stage")

model_covariates <- c("tumor_size_cm", "age_per_10y", "central_location", "oxygen_any")

# Nothing in the analysis is random; this is only recorded in the run manifest.
seed <- 20260829L

packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "readr", "tibble", "survival", "cmprsk", "logistf",
  "ggplot2", "broom", "patchwork", "ragg", "svglite", "rsvg", "scales", "digest", "jsonlite"
)
