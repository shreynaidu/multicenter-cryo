#!/usr/bin/env Rscript
# Install the renv.lock versions into renv/library.
script_arg <- gsub("~+~", " ", grep("^--file=", commandArgs(), value = TRUE), fixed = TRUE)
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
library_dir <- file.path(project_dir, "renv", "library")
Sys.setenv(
  RENV_PATHS_ROOT = file.path(project_dir, "renv", "state"),
  RENV_PATHS_CACHE = file.path(project_dir, "renv", "cache"),
  RENV_CONFIG_CACHE_SYMLINKS = "FALSE"
)
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))

# renv can't lock itself. Older versions move to the CRAN archive.
renv_version <- "1.3.1"
has_locked_renv <- function() {
  requireNamespace("renv", quietly = TRUE) && packageVersion("renv") == renv_version
}
cran <- "https://cloud.r-project.org/src/contrib"
for (location in c(cran, file.path(cran, "Archive", "renv"))) {
  if (!has_locked_renv()) {
    url <- sprintf("%s/renv_%s.tar.gz", location, renv_version)
    try(install.packages(url, repos = NULL, type = "source", lib = library_dir))
  }
}
if (!has_locked_renv()) {
  stop("renv ", renv_version, " could not be installed. Restore requires package-repository access.")
}
renv::restore(
  project = project_dir, lockfile = file.path(project_dir, "renv.lock"),
  library = library_dir, prompt = FALSE
)
