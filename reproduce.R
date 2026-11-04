#!/usr/bin/env Rscript
# Set up dependencies, run everything, and optionally verify.
script_arg <- gsub("~+~", " ", grep("^--file=", commandArgs(), value = TRUE), fixed = TRUE)
project_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(project_dir, "analysis", "cli.R"))

usage <- paste(
  "Usage: Rscript --vanilla reproduce.R [--verify] [--skip-setup]",
  "[--input SOURCE] [--output DIRECTORY] [--library DIRECTORY]",
  "[--figure-python EXECUTABLE] [--python PYTHON_3_12]"
)

restore_r_packages <- function(opts, default_library) {
  use_project_library(project_dir, opts$library)
  if (!length(locked_package_mismatches(project_dir))) {
    return(invisible(NULL))
  }
  if (opts$skip_setup) {
    stop("Recorded R packages are unavailable; run without --skip-setup to restore locally.")
  }
  if (opts$library != default_library) {
    stop("The supplied --library does not match renv.lock. Restore the project library or supply a matching one.")
  }
  run_command(rscript_path(), c("--vanilla", file.path(project_dir, "analysis", "restore_environment.R")))
}

python_environment_ready <- function(python, renderer) {
  is_runnable(python) && system2(
    python, shQuote(c(renderer, "--check-environment")),
    stdout = FALSE, stderr = FALSE
  ) == 0L
}

is_python_3_12 <- function(python) {
  check <- "import sys; sys.exit(sys.version_info[:2] != (3, 12))"
  is_runnable(python) && system2(python, shQuote(c("-c", check)), stdout = FALSE, stderr = FALSE) == 0L
}

create_figure_environment <- function(opts, renderer) {
  candidates <- if (is.null(opts$python)) c("python3.12", "python3", "python") else opts$python
  suitable <- candidates[vapply(candidates, is_python_3_12, logical(1))]
  if (!length(suitable)) {
    stop("Install Python 3.12, or supply --python with its executable. This project installs no system software.")
  }
  run_command(suitable[1], c("-m", "venv", file.path(project_dir, ".venv-figures")))
  Sys.setenv(PIP_CACHE_DIR = file.path(project_dir, ".cache", "pip"))
  run_command(opts$figure_python, c("-m", "pip", "install", "-r", file.path(project_dir, "requirements-figures.lock")))
  run_command(opts$figure_python, c(renderer, "--check-environment"))
}

main <- function() {
  default_library <- file.path(project_dir, "renv", "library")
  default_python <- venv_python(project_dir)
  opts <- parse_cli(
    usage,
    options = list(
      input = default_source(project_dir),
      output = file.path(project_dir, "outputs", "current"),
      library = default_library,
      figure_python = default_python,
      python = NULL
    ),
    switches = c("--verify", "--skip-setup")
  )
  configure_locale()
  if (!file.exists(opts$input)) {
    stop("Source data not found: ", opts$input, "\nThe repository contains no patient-level data; see README.md.")
  }
  if (opts$verify && !nzchar(Sys.which("pdftoppm"))) {
    stop("Verification requires Poppler's pdftoppm on PATH; generation does not.")
  }
  restore_r_packages(opts, default_library)

  renderer <- file.path(project_dir, "analysis", "scripts", "generate_publication_figures.py")
  if (!python_environment_ready(opts$figure_python, renderer)) {
    if (opts$skip_setup || opts$figure_python != default_python) {
      stop(
        "Python plotting dependencies are unavailable; ",
        "restore the local environment or supply a matching --figure-python."
      )
    }
    create_figure_environment(opts, renderer)
  }

  entry <- if (opts$verify) "tests/test_reproducibility.R" else "run_all.R"
  run_command(rscript_path(), c(
    "--vanilla", file.path(project_dir, entry),
    "--input", opts$input, "--output", opts$output, "--library", opts$library,
    "--python-figures", "--figure-python", opts$figure_python
  ))
  cat("Results:", normalizePath(opts$output), "\n")
}

main()
