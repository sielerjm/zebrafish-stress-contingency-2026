# source_project_init.R
# Created by: Michael Sieler
# Date last updated: 2026-03-29
#
# Description: Bootstrap for R Markdown and interactive sessions — finds the Sieler2026 repo
#   root, sets working directory, optionally knitr root.dir, then sources
#   00__InitializeEnvironment.R. Keep this file in Code/00__Setup/ (not in 03__HelperFunctions.R:
#   helpers load after init, so they cannot define the step that sources init).
#
# Expected input:  none (uses knitr::current_input() when knitting, else getwd()).
# Expected output:  proj.path and other objects from 00__InitializeEnvironment.R in .GlobalEnv;
#   invisibly returns proj.path.

# Find repo root (directory containing Code/00__Setup/00__InitializeEnvironment.R).
# start: directory to walk up from (default: knitr input dir when knitting, else getwd()).
sieler2026_find_project_root <- function(start = NULL) {
  if (is.null(start)) {
    start <- if (requireNamespace("knitr", quietly = TRUE) && !is.null(knitr::current_input())) {
      dirname(knitr::current_input(dir = TRUE))
    } else {
      getwd()
    }
  }
  init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
  root <- normalizePath(start, mustWork = TRUE)
  while (!file.exists(file.path(root, init_rel))) {
    parent <- dirname(root)
    if (identical(parent, root)) {
      stop("Could not find Sieler2026 project root (", init_rel, ") above ", start)
    }
    root <- parent
  }
  root
}

# Set WD to repo root, knitr root.dir, source full environment init.
# root: optional path to repo root (e.g. after a single walk in an Rmd setup chunk); if NULL, walks from knitr/getwd().
sieler2026_source_initialize_project <- function(root = NULL) {
  if (is.null(root)) {
    root <- sieler2026_find_project_root()
  } else {
    root <- normalizePath(root, mustWork = TRUE)
  }
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::opts_knit$set(root.dir = root)
  }
  init_path <- file.path(root, "Code", "00__Setup", "00__InitializeEnvironment.R")
  source(init_path, local = FALSE)
  # Init assigns into .GlobalEnv; check there (not inherits = FALSE in this function frame).
  stopifnot(
    exists("proj.path", envir = .GlobalEnv, inherits = FALSE),
    identical(normalizePath(root), normalizePath(get("proj.path", envir = .GlobalEnv)))
  )
  invisible(get("proj.path", envir = .GlobalEnv))
}
