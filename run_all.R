#!/usr/bin/env Rscript
script_arg <- gsub("~+~", " ", grep("^--file=", commandArgs(), value = TRUE), fixed = TRUE)
project_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(project_dir, "analysis", "cli.R"))
source(file.path(project_dir, "analysis", "config.R"))
for (module in c(
  "survival_estimates", "prepare_data", "descriptive_tables", "models", "diagnostics", "size_summary",
  "safety", "reviewer_sensitivities", "exploratory_size", "manuscript_statistics", "figures"
)) {
  source(file.path(project_dir, "analysis", "R", paste0(module, ".R")))
}

usage <- paste(
  "Usage: Rscript --vanilla run_all.R [--input SOURCE] [--output DIRECTORY]",
  "[--library DIRECTORY] [--no-figures] [--no-package] [--python-figures]",
  "[--figure-python EXECUTABLE]"
)

check_environment <- function(opts) {
  mismatches <- locked_package_mismatches(project_dir)
  if (length(mismatches)) {
    stop(
      "Package versions differ from renv.lock: ", paste(mismatches, collapse = ", "),
      ". Run analysis/restore_environment.R to restore the recorded versions."
    )
  }
  if (getRversion() != package_version("4.5.0")) {
    warning("The recorded environment used R 4.5.0; this run uses ", getRversion())
  }
  if (opts$no_figures && !opts$no_package) stop("Use --no-package with --no-figures.")
  if (opts$python_figures) {
    if (opts$no_figures) stop("--python-figures conflicts with --no-figures.")
    if (!is_runnable(opts$figure_python)) {
      stop("Set up .venv-figures or supply --figure-python.")
    }
    run_command(opts$figure_python, c(python_renderer(), "--check-environment"))
  }
}

python_renderer <- function() {
  file.path(project_dir, "analysis", "scripts", "generate_publication_figures.py")
}

# Only replaces directories this pipeline created. The last run is kept as .previous.
prepare_output_directory <- function(out_dir, input_file) {
  marker <- file.path(out_dir, ".cryo-analysis")
  populated <- dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))
  if (populated) {
    if (!file.exists(marker)) {
      stop("Output directory is not empty and has no analysis marker; choose a new --output directory.")
    }
    inside_output <- startsWith(c(input_file, project_dir), paste0(out_dir, .Platform$file.sep))
    if (any(inside_output) || identical(out_dir, project_dir)) {
      stop("Output directory overlaps source data or project code.")
    }
    if (nzchar(Sys.readlink(out_dir))) stop("Output directory must not be a symbolic link.")
    previous <- paste0(out_dir, ".previous")
    if (dir.exists(previous)) {
      if (!file.exists(file.path(previous, ".cryo-analysis"))) {
        stop("Refusing to replace an unmarked directory: ", previous)
      }
      unlink(previous, recursive = TRUE)
    }
    if (!file.rename(out_dir, previous)) stop("Could not set aside the previous run.")
  }
  for (folder in c("analysis_data", "model_results", "tables", "diagnostics")) {
    dir.create(file.path(out_dir, folder), recursive = TRUE, showWarnings = FALSE)
  }
  writeLines("Created by run_all.R.", marker)
}

run_python_figures <- function(opts, out_dir) {
  message("Generating the Python figures")
  figure_dir <- file.path(out_dir, "python_publication_figures")
  run_command(opts$figure_python, c(
    python_renderer(),
    "--input", file.path(out_dir, "analysis_data", "deidentified_corrected_analysis_dataset.csv"),
    "--output", figure_dir
  ))
  run_command(rscript_path(), c(
    "--vanilla", file.path(project_dir, "analysis", "scripts", "audit_python_figure_statistics.R"),
    "--figures", figure_dir, "--analysis", out_dir, "--library", opts$library
  ))
}

file_digest <- function(path, root) {
  list(
    path = substring(path, nchar(root) + 2L),
    sha256 = digest::digest(file = path, algo = "sha256")
  )
}

write_run_records <- function(opts, input_file, out_dir, source_sha256, records) {
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
  installed <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")])
  loaded <- installed[installed$Package %in% loadedNamespaces(), ]
  loaded <- loaded[order(match(loaded$LibPath, .libPaths()), tolower(loaded$Package)), ]
  write_table(loaded[, c("Package", "Version")], file.path(out_dir, "package_versions.csv"))

  code_dir <- file.path(project_dir, "analysis")
  code_files <- c(
    file.path(code_dir, c("cli.R", "config.R")),
    file.path(project_dir, c("run_all.R", "renv.lock")),
    list.files(file.path(code_dir, "R"), full.names = TRUE),
    file.path(code_dir, "scripts", "package_assets.py")
  )
  if (opts$python_figures) {
    code_files <- c(
      code_files,
      list.files(file.path(code_dir, "scripts"), pattern = "figure|generate_publication", full.names = TRUE),
      list.files(file.path(code_dir, "figure_reference"), full.names = TRUE),
      file.path(project_dir, "requirements-figures.lock")
    )
  }
  outputs <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  outputs <- outputs[!dir.exists(outputs)]
  jsonlite::write_json(
    list(
      code = lapply(code_files, file_digest, root = project_dir),
      source_name = basename(input_file),
      source_sheet = input_sheet,
      source_sha256 = source_sha256,
      seed = seed,
      records = records,
      figures_generated = !opts$no_figures,
      python_figures_generated = opts$python_figures,
      package_generated = !opts$no_package,
      R_version = as.character(getRversion()),
      files = lapply(outputs, file_digest, root = out_dir)
    ),
    file.path(out_dir, "run_manifest.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

main <- function() {
  opts <- parse_cli(
    usage,
    options = list(
      input = default_source(project_dir),
      output = file.path(project_dir, "outputs", "current"),
      library = file.path(project_dir, "renv", "library"),
      figure_python = venv_python(project_dir)
    ),
    switches = c("--no-figures", "--no-package", "--python-figures")
  )
  configure_locale()
  if (!file.exists(opts$input)) {
    stop("Source data not found: ", opts$input, "\nThe repository contains no patient-level data; see README.md.")
  }
  input_file <- normalizePath(opts$input)
  out_dir <- normalizePath(opts$output, mustWork = FALSE)
  use_project_library(project_dir, opts$library)
  code_dir <- file.path(project_dir, "analysis")
  check_environment(opts)

  source_sha256 <- digest::digest(file = input_file, algo = "sha256")
  if (!source_sha256 %in% accepted_source_sha256) {
    stop("Source SHA-256 differs from the recorded versions. Review the input version before running.")
  }
  packaging_python <- if (opts$no_package) NULL else standard_python(opts$figure_python)
  prepare_output_directory(out_dir, input_file)

  # Stop on missing covariates rather than dropping patients.
  old_options <- options(
    stringsAsFactors = FALSE,
    warn = 1,
    na.action = "na.fail",
    contrasts = c("contr.treatment", "contr.poly")
  )
  on.exit(options(old_options), add = TRUE)
  Sys.setenv(TZ = "UTC")
  set.seed(seed)
  suppressPackageStartupMessages(
    for (package in packages) library(package, character.only = TRUE)
  )
  make_figures <- !opts$no_figures

  message("Preparing the analysis cohort")
  raw <- read_source_data(input_file)
  expected_checks <- read_expected_source_checks(file.path(code_dir, "reference", "expected_source_checks.csv"))
  full_cohort <- prepare_cohort(raw, expected_checks, out_dir)
  run_descriptive_tables(full_cohort, raw, out_dir)

  message("Fitting models and summarizing outcomes")
  dat <- analysis_cohort(full_cohort)
  primary_models <- run_models(dat, out_dir)
  run_diagnostics(dat, primary_models, out_dir)
  run_size_summaries(dat, out_dir)
  safety <- run_safety_models(dat, out_dir)
  center_file <- file.path(project_dir, "data", "reference", "center_assignments.csv")
  run_reviewer_sensitivities(dat, safety, center_file, out_dir)
  run_exploratory_size(dat, out_dir, make_figures)
  run_manuscript_statistics(dat, out_dir)

  if (make_figures) {
    message("Exporting figures")
    run_figures(out_dir)
  }
  if (!identical(digest::digest(file = input_file, algo = "sha256"), source_sha256)) {
    stop("The source workbook changed during the run.")
  }
  jsonlite::write_json(
    list(source_name = basename(input_file), source_sheet = input_sheet, source_sha256 = source_sha256),
    file.path(out_dir, "source_manifest.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
  if (!opts$no_package) {
    message("Preparing publication tables and archive")
    run_command(packaging_python, c(file.path(code_dir, "scripts", "package_assets.py"), "--output", out_dir))
  }
  if (opts$python_figures) run_python_figures(opts, out_dir)
  write_run_records(opts, input_file, out_dir, source_sha256, nrow(dat))
  message("Analysis complete: ", out_dir)
}

main()
