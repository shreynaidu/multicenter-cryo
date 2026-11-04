#!/usr/bin/env Rscript
# Rerun the analysis and compare it with the frozen baseline. Numbers have to match
# everywhere; figure files only on the platform the baseline came from.
# Don't regenerate the baseline to get a pass.
script_arg <- gsub("~+~", " ", grep("^--file=", commandArgs(), value = TRUE), fixed = TRUE)
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(project_dir, "analysis", "cli.R"))

numeric_tolerance <- 1e-12

default_baseline <- function() {
  private <- file.path(project_dir, "tests", "baseline")
  if (dir.exists(private)) private else file.path(project_dir, "tests", "public_baseline")
}

on_reference_platform <- function() {
  info <- Sys.info()
  info[["sysname"]] == "Darwin" && info[["machine"]] == "arm64"
}

run_pipeline <- function(opts, result_dir) {
  dir.create(dirname(result_dir), recursive = TRUE, showWarnings = FALSE)
  log_file <- tempfile("cryo-test-", tmpdir = dirname(result_dir), fileext = ".log")
  arguments <- c(
    "--vanilla", file.path(project_dir, "run_all.R"),
    "--input", opts$input, "--output", result_dir, "--library", opts$library
  )
  if (opts$python_figures) {
    arguments <- c(arguments, "--python-figures", "--figure-python", opts$figure_python)
  }
  started <- proc.time()[["elapsed"]]
  run_command(rscript_path(), arguments, stdout = log_file, stderr = log_file)
  list(log = log_file, runtime = proc.time()[["elapsed"]] - started)
}

verify_baseline_checksums <- function(baseline_dir) {
  checksums <- jsonlite::read_json(file.path(baseline_dir, "checksums.json"))
  for (entry in checksums) {
    path <- file.path(baseline_dir, entry$path)
    if (!file.exists(path) || !identical(digest::digest(file = path, algo = "sha256"), entry$sha256)) {
      stop("Baseline integrity failure: ", entry$path)
    }
  }
  length(checksums)
}

read_cells <- function(path) {
  read.csv(path, colClasses = "character", check.names = FALSE, na.strings = character())
}

csv_comparison <- function(file, rows = NA_integer_, numeric_cells = 0L, maximum = NA_real_, status) {
  data.frame(
    file = file, rows = rows, numeric_cells = numeric_cells,
    max_absolute_difference = maximum, status = status
  )
}

compare_csv <- function(before, after, label, ignore_row_order = FALSE) {
  if (!file.exists(after)) {
    return(csv_comparison(label, status = "missing"))
  }
  expected <- read_cells(before)
  actual <- read_cells(after)
  if (!identical(names(expected), names(actual)) || !identical(dim(expected), dim(actual))) {
    return(csv_comparison(label, status = "different shape"))
  }
  if (ignore_row_order) {
    expected <- expected[do.call(order, expected), , drop = FALSE]
    actual <- actual[do.call(order, actual), , drop = FALSE]
  }
  numeric_cells <- 0L
  maximum <- 0
  differing <- character()
  for (column in names(expected)) {
    x <- suppressWarnings(as.numeric(expected[[column]]))
    y <- suppressWarnings(as.numeric(actual[[column]]))
    finite <- is.finite(x) & is.finite(y)
    numeric_cells <- numeric_cells + sum(finite)
    difference <- abs(x[finite] - y[finite])
    maximum <- max(c(maximum, difference))
    same <- expected[[column]] == actual[[column]]
    same[finite] <- difference <= numeric_tolerance + numeric_tolerance * abs(x[finite])
    if (any(!same)) differing <- c(differing, column)
  }
  status <- if (length(differing)) paste("differs in", paste(differing, collapse = ", ")) else "match"
  csv_comparison(label, nrow(actual), numeric_cells, maximum, status)
}

compare_all_csv <- function(baseline_dir, result_dir) {
  results_dir <- file.path(baseline_dir, "results")
  primary <- lapply(list.files(results_dir, pattern = "\\.csv$", recursive = TRUE), function(relative) {
    # Row order follows library folders, which vary by install.
    compare_csv(
      file.path(results_dir, relative), file.path(result_dir, relative), relative,
      ignore_row_order = relative == "package_versions.csv"
    )
  })
  # Analyses added after the original run have their own baselines.
  additional_output <- c(
    descriptive = "tables", reviewer = "reviewer_sensitivities",
    probe_omitted = "supplemental", size_three_groups = "exploratory_size"
  )
  additional <- lapply(names(additional_output), function(folder) {
    frozen <- file.path(baseline_dir, "additional", folder)
    lapply(list.files(frozen, pattern = "\\.csv$"), function(relative) {
      label <- file.path(additional_output[[folder]], relative)
      compare_csv(file.path(frozen, relative), file.path(result_dir, label), label)
    })
  })
  do.call(rbind, c(primary, unlist(additional, recursive = FALSE)))
}

# Only the private baseline has this file.
compare_model_objects <- function(baseline_dir, result_dir) {
  frozen <- file.path(baseline_dir, "results", "analysis_objects.rds")
  if (!file.exists(frozen)) {
    return(0L)
  }
  expected <- readRDS(frozen)
  actual <- readRDS(file.path(result_dir, "analysis_objects.rds"))
  for (component in c("data", "cohorts")) {
    if (!isTRUE(all.equal(expected[[component]], actual[[component]], tolerance = numeric_tolerance))) {
      stop("Serialized analysis data differ: ", component)
    }
  }
  fields <- c(
    "coef", "coefficients", "var", "loglik", "score", "linear.predictors", "residuals", "n", "nevent", "converged"
  )
  compared <- 0L
  for (name in names(expected$main_models)) {
    left <- expected$main_models[[name]]
    right <- actual$main_models[[name]]
    if (!identical(class(left), class(right))) stop("Model class differs: ", name)
    for (field in fields) {
      if (!isTRUE(all.equal(left[[field]], right[[field]], tolerance = numeric_tolerance))) {
        stop("Model component differs: ", name, "/", field)
      }
      if (!is.null(left[[field]])) compared <- compared + 1L
    }
  }
  compared
}

compare_figures <- function(baseline_dir, result_dir, python, report_only) {
  log <- tempfile("cryo-figures-", fileext = ".json")
  run_command(
    python,
    c(file.path(project_dir, "tests", "compare_artifacts.py"), file.path(baseline_dir, "results"), result_dir),
    stdout = log
  )
  figures <- jsonlite::read_json(log)
  differences <- unlist(figures$differences)
  frozen_png <- file.path(baseline_dir, "additional", "size_three_groups", "size_group_curves.png")
  exploratory_png <- file.path(result_dir, "exploratory_size", "size_group_curves.png")
  if (!identical(unname(tools::md5sum(frozen_png)), unname(tools::md5sum(exploratory_png)))) {
    differences <- c(differences, "exploratory_size/size_group_curves.png")
  }
  if (length(differences)) {
    message <- paste("Figures differ from the baseline:", paste(differences, collapse = ", "))
    if (report_only) warning(message, call. = FALSE) else stop(message)
  }
  list(compared = figures$figures_compared, differences = differences)
}

check_manuscript_counts <- function(result_dir) {
  expected_dir <- file.path(project_dir, "tests", "expected")
  expected_flow <- read.csv(file.path(expected_dir, "cohort_flow.csv"))
  flow <- read.csv(file.path(result_dir, "tables", "cohort_flow.csv"))
  stopifnot(
    identical(flow$analysis, expected_flow$analysis),
    identical(flow$records_after, expected_flow$records_after),
    identical(flow$records_excluded, expected_flow$records_excluded)
  )
  expected <- read.csv(file.path(expected_dir, "manuscript_counts.csv"))
  counts <- read.csv(file.path(result_dir, "manuscript_statistics", "descriptive_counts.csv"))
  if (!identical(counts$statistic, expected$statistic) || !identical(counts$numerator, expected$numerator)) {
    stop("Manuscript descriptive counts differ")
  }
  if (!identical(counts$denominator, expected$denominator)) stop("Manuscript denominators differ")
  stopifnot(max(abs(counts$percent - 100 * expected$numerator / expected$denominator)) < numeric_tolerance)
  nrow(expected)
}

check_publication_tables <- function(baseline_dir, result_dir) {
  descriptive <- file.path(baseline_dir, "additional", "descriptive")
  for (pair in list(
    c("table1_baseline_characteristics.csv", "MainTable1.csv"),
    c("table2_procedure_safety.csv", "MainTable2.csv")
  )) {
    stopifnot(identical(
      unname(tools::md5sum(file.path(descriptive, pair[1]))),
      unname(tools::md5sum(file.path(result_dir, "publication_tables", pair[2])))
    ))
  }
  omitted <- read.csv(file.path(baseline_dir, "additional", "probe_omitted", "refitted_probe_omitted_models.csv"))
  panel_b <- read.csv(
    file.path(result_dir, "publication_tables", "SupplementaryTable5_PanelB.csv"),
    colClasses = "character", check.names = FALSE
  )
  intervals <- sprintf("%.2f (%.2f–%.2f)", omitted$OR, omitted$lower, omitted$upper)
  p_values <- ifelse(omitted$p < 0.001, "<.001", sub("^0", "", sprintf("%.3f", omitted$p)))
  stopifnot(identical(panel_b[["OR (95% CI)"]], intervals), identical(panel_b$P, p_values))

  frozen <- read.csv(file.path(baseline_dir, "results", "tables", "cif_all_cohorts.csv"))
  frozen_ltp <- frozen$estimate[
    frozen$cohort == "Full107" & frozen$event == "Local tumor progression" & frozen$time_months == 24
  ]
  derived <- read.csv(file.path(result_dir, "manuscript_statistics", "derived_probabilities.csv"))
  stopifnot(abs(derived$estimate - (1 - frozen_ltp)) < numeric_tolerance)
}

check_source_cell_registry <- function(result_dir) {
  registry <- read_cells(file.path(result_dir, "manuscript_statistics", "source_cells.csv"))
  for (source_file in unique(registry$source_file)) {
    cells <- registry[registry$source_file == source_file, ]
    source_table <- read_cells(file.path(result_dir, source_file))
    for (field in unique(cells$field)) {
      records <- cells[cells$field == field, ]
      stopifnot(identical(records$value, source_table[[field]][as.integer(records$source_row)]))
    }
  }
}

main <- function() {
  opts <- parse_cli(
    paste(
      "Usage: Rscript --vanilla tests/test_reproducibility.R [--input SOURCE]",
      "[--output DIRECTORY | --compare-only DIRECTORY] [--library DIRECTORY] [--baseline DIRECTORY]",
      "[--python-figures] [--figure-python EXECUTABLE] [--figures-report-only]"
    ),
    options = list(
      input = default_source(project_dir),
      output = file.path(project_dir, "outputs", "reproducibility"),
      compare_only = NULL,
      library = file.path(project_dir, "renv", "library"),
      baseline = default_baseline(),
      figure_python = venv_python(project_dir)
    ),
    switches = c("--python-figures", "--figures-report-only")
  )
  configure_locale()
  use_project_library(project_dir, opts$library)

  generated <- is.null(opts$compare_only)
  run <- if (generated) run_pipeline(opts, opts$output) else list(log = NULL, runtime = NA_real_)
  result_dir <- normalizePath(if (generated) opts$output else opts$compare_only)
  if (!file.exists(file.path(result_dir, "run_manifest.json"))) stop("Incomplete analysis run")
  report_dir <- file.path(result_dir, "validation")
  dir.create(report_dir, showWarnings = FALSE)

  baseline_dir <- normalizePath(opts$baseline)
  baseline_files <- verify_baseline_checksums(baseline_dir)
  comparisons <- compare_all_csv(baseline_dir, result_dir)
  write.csv(comparisons, file.path(report_dir, "baseline_comparison.csv"), row.names = FALSE)
  failed <- comparisons[comparisons$status != "match", ]
  if (nrow(failed)) {
    stop("Tables differ from the baseline:\n", paste(failed$file, failed$status, sep = ": ", collapse = "\n"))
  }
  model_components <- compare_model_objects(baseline_dir, result_dir)
  report_only <- opts$figures_report_only || !on_reference_platform()
  figures <- compare_figures(baseline_dir, result_dir, standard_python(opts$figure_python), report_only)
  manuscript_counts <- check_manuscript_counts(result_dir)
  check_publication_tables(baseline_dir, result_dir)
  check_source_cell_registry(result_dir)

  run_command(
    rscript_path(),
    c(
      "--vanilla", file.path(project_dir, "tests", "verify_primary_models.R"),
      "--input", opts$input, "--results", result_dir, "--library", opts$library
    ),
    stdout = file.path(report_dir, "primary_model_check.log"),
    stderr = file.path(report_dir, "primary_model_check.log")
  )
  python_report <- NULL
  if (opts$python_figures) {
    run_command(opts$figure_python, c(
      file.path(project_dir, "analysis", "validation", "independent_qc.py"),
      "--input", opts$input, "--output", result_dir
    ))
    run_command(opts$figure_python, c(
      file.path(project_dir, "tests", "test_python_figures.py"),
      "--baseline", file.path(baseline_dir, "python_figures"),
      "--compare-only", file.path(result_dir, "python_publication_figures"),
      "--analysis", result_dir, "--library", opts$library, "--rscript", rscript_path(),
      if (report_only) "--figures-report-only"
    ))
    python_report <- jsonlite::read_json(file.path(result_dir, "python_publication_figures", "reproducibility.json"))
  }
  report <- list(
    python_figure_verification = python_report,
    status = "passed",
    direct_source_refit = TRUE,
    baseline_files_verified = baseline_files,
    csv_files_compared = nrow(comparisons),
    numeric_cells_compared = sum(comparisons$numeric_cells),
    maximum_absolute_difference = max(comparisons$max_absolute_difference),
    model_components_compared = model_components,
    figures_compared = figures$compared,
    figure_differences_reported = as.list(figures$differences),
    figures_required_to_match = !report_only,
    additional_figures_compared = 1L,
    manuscript_counts_compared = manuscript_counts,
    runtime_seconds = run$runtime,
    tolerance = "absolute 1e-12 plus relative 1e-12; labels and dimensions exact",
    figure_rule = "PDF pixels identical at 300 dpi; other figure bytes identical"
  )
  jsonlite::write_json(report, file.path(report_dir, "reproducibility.json"), pretty = TRUE, auto_unbox = TRUE)
  if (!is.null(run$log)) invisible(file.rename(run$log, file.path(report_dir, "run.log")))
  cat(
    "PASS:", nrow(comparisons), "CSV files;", report$numeric_cells_compared, "numeric cells;",
    model_components, "model components;", figures$compared, "figures. Maximum numeric difference:",
    report$maximum_absolute_difference, "\n"
  )
}

main()
