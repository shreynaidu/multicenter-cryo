parse_cli <- function(usage, options = list(), switches = character()) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat(usage, "\n")
    quit(status = 0L)
  }
  key <- function(flag) gsub("-", "_", sub("^--", "", flag))
  values <- options
  values[key(switches)] <- FALSE
  index <- 1L
  while (index <= length(args)) {
    flag <- args[index]
    if (flag %in% switches) {
      values[[key(flag)]] <- TRUE
      index <- index + 1L
    } else if (key(flag) %in% names(options) && startsWith(flag, "--")) {
      if (index == length(args) || startsWith(args[index + 1L], "--")) {
        stop("Missing value for ", flag)
      }
      values[[key(flag)]] <- path.expand(args[index + 1L])
      index <- index + 2L
    } else {
      stop("Unknown argument: ", flag, "\n", usage)
    }
  }
  values
}

# Size labels contain ≥, which R misreads outside UTF-8.
configure_locale <- function() {
  Sys.setlocale("LC_COLLATE", "C")
  if (isTRUE(l10n_info()[["UTF-8"]])) {
    return(invisible(NULL))
  }
  for (locale in c("C.UTF-8", "en_US.UTF-8", "UTF-8")) {
    suppressWarnings(Sys.setlocale("LC_CTYPE", locale))
    if (isTRUE(l10n_info()[["UTF-8"]])) {
      return(invisible(NULL))
    }
  }
  stop("A UTF-8 locale is required; set LC_ALL to one such as C.UTF-8.")
}

default_source <- function(project_dir) {
  candidates <- file.path(
    project_dir, "data",
    c("raw/Combined Multicenter.xlsx", "release/multicenter_cryo_release.csv")
  )
  available <- candidates[file.exists(candidates)]
  if (length(available)) available[1] else candidates[1]
}

use_project_library <- function(project_dir, library_dir = NULL) {
  candidates <- c(library_dir, file.path(project_dir, "renv", "library"))
  .libPaths(c(candidates[dir.exists(candidates)], .libPaths()))
}

locked_package_mismatches <- function(project_dir) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return("jsonlite")
  }
  locked <- jsonlite::read_json(file.path(project_dir, "renv.lock"))$Packages
  differs <- vapply(locked, function(record) {
    !requireNamespace(record$Package, quietly = TRUE) ||
      packageVersion(record$Package) != package_version(record$Version)
  }, logical(1))
  names(locked)[differs]
}

venv_python <- function(project_dir) {
  relative <- if (.Platform$OS.type == "windows") "Scripts/python.exe" else "bin/python"
  file.path(project_dir, ".venv-figures", relative)
}

is_runnable <- function(command) {
  file.exists(command) || nzchar(Sys.which(command))
}

# Packaging and comparison only need the standard library, so any Python 3 works.
standard_python <- function(figure_python) {
  candidates <- c(figure_python, "python3", "python")
  available <- candidates[vapply(candidates, is_runnable, logical(1))]
  if (!length(available)) {
    stop("Python 3 is required for publication packaging and comparison.")
  }
  available[1]
}

rscript_path <- function() {
  file.path(R.home("bin"), "Rscript")
}

run_command <- function(command, arguments = character(), ...) {
  status <- system2(command, vapply(arguments, shQuote, character(1)), ...)
  if (status != 0L) {
    stop("Command failed: ", basename(command), " (exit ", status, ")")
  }
}
