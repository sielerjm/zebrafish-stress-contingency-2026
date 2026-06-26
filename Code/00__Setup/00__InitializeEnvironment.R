# 00__InitializeEnvironment.R
# Created by: Michael Sieler
# Date last updated: 2026-03-29
#
# Expected input: none (uses here::here() for project root)
# Expected output: loads libraries, plot settings, and HelperFunctions; defines path.*; loads
#   ps.list / data.list when dated RDS exist under Data/r_objects/.
#
# Loads `Code/00__Setup/` in dependency order: packages → ggplot theme/palettes → helpers.
# Does NOT source 04__DataPreProcess.R — that script is a batch driver (rebuild phyloseq + r_objects);
#   run it manually when DADA2 inputs or filters change, then reload the session or re-knit.
# Archive pipeline chunks under `Code/99__Archive/Functions/AnalysisFunctions/` (they mutate ps.list);
# source those from analysis drivers or Rmd chunks, not from this init script.
#
# Phyloseq objects: Data/DADA2/ holds pseq_uncleaned_*.rds and dated pseq_cleaned_filtered_*.rds;
# path.pseq.* point at the newest file of each pattern. Optional list objects (ps-list__*.rds, …)
# belong in Data/r_objects/ once built by preprocessing drivers.

suppressPackageStartupMessages(library(here))

# Project root must be this repository (Sieler2026/), not a parent multi-root workspace.
# When run via Rscript, anchor from --file=; otherwise use here::i_am() with the path below.
script_path_from_cli <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[startsWith(ca, "--file=")])
  if (length(f) == 1L && nzchar(f)) {
    return(normalizePath(f, winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

.sp <- script_path_from_cli()
if (!is.na(.sp)) {
  proj.path <- normalizePath(file.path(dirname(.sp), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  here::i_am("Code/00__Setup/00__InitializeEnvironment.R")
  proj.path <- here::here()
}

# --- Shared R utilities: attach packages first, then session plot defaults, then helper functions
path.setup <- file.path(proj.path, "Code", "00__Setup")
source(file.path(path.setup, "01__Libraries.R"))
source(file.path(path.setup, "02__PlotSettings.R"))
source(file.path(path.setup, "03__HelperFunctions.R"))

# --- Canonical data / results paths (see repository README)
path.code <- file.path(proj.path, "Code")
path.data <- file.path(proj.path, "Data")
path.objects <- file.path(path.data, "r_objects")
path.input <- file.path(path.data, "input")
path.results <- file.path(proj.path, "Results")
path.stats.input <- file.path(path.data, "stats_intermediate")
path.dada2 <- file.path(path.data, "DADA2")
path.deg <- file.path(path.data, "DEG")
path.metadata <- file.path(path.data, "Metadata", "metadata.tsv")
path.context <- file.path(path.data, "Context", "ExperimentalDesignContext.md")

# Primary phyloseq RDS paths (resolve_latest_rds is defined in 03, sourced above).
path.pseq.cleaned <- resolve_latest_rds(path.dada2, "^pseq_cleaned_filtered_.*\\.rds$", error_if_empty = FALSE)
path.pseq.uncleaned <- resolve_latest_rds(path.dada2, "^pseq_uncleaned_.*\\.rds$", error_if_empty = FALSE)

path.data.legacy <- path.data
path.r.objects.legacy <- path.objects

analysis.ID <- paste0("Sieler_2026__", Sys.Date())

# --- Flag missing paths once (directories or files expected by the project layout)
flag_followup_missing <- function(path, description) {
  exists <- if (grepl("\\.(rds|tsv|md)$", path, ignore.case = TRUE)) {
    file.exists(path)
  } else {
    dir.exists(path)
  }
  if (!exists) {
    message("FOLLOW-UP: ", description, " not found: ", path)
  }
}

flag_followup_missing(path.objects, "Data/r_objects (for ps-list / data-list RDS)")
flag_followup_missing(path.input, "Data/input")
flag_followup_missing(path.stats.input, "Data/stats_intermediate")
flag_followup_missing(path.metadata, "Sample metadata TSV")
flag_followup_missing(path.context, "Experimental design context")
flag_followup_missing(path.dada2, "Data/DADA2")
flag_followup_missing(path.deg, "Data/DEG")

if (dir.exists(path.dada2)) {
  if (is.na(path.pseq.cleaned) || !file.exists(path.pseq.cleaned)) {
    message("FOLLOW-UP: No pseq_cleaned_filtered_*.rds in ", path.dada2, " (run 04__DataPreProcess.R).")
  } else {
    message("path.pseq.cleaned -> ", basename(path.pseq.cleaned))
  }
  if (is.na(path.pseq.uncleaned) || !file.exists(path.pseq.uncleaned)) {
    message("FOLLOW-UP: No pseq_uncleaned_*.rds in ", path.dada2, ".")
  } else {
    message("path.pseq.uncleaned -> ", basename(path.pseq.uncleaned))
  }
}

# Load ps.list and data.list from dated RDS files in Data/r_objects/
if (dir.exists(path.objects)) {
  ps_list_files <- list.files(path.objects, pattern = "^ps-list__.*\\.rds$", full.names = TRUE)
  data_list_files <- list.files(path.objects, pattern = "^data-list__.*\\.rds$", full.names = TRUE)

  if (length(ps_list_files) > 0) {
    ps_list_file <- ps_list_files[which.max(file.mtime(ps_list_files))]
    ps.list <- readRDS(ps_list_file)
    message("Loaded ps.list from: ", basename(ps_list_file))
  } else {
    message("No ps-list__*.rds files in ", path.objects, " (optional until preprocessing exports them).")
  }

  if (length(data_list_files) > 0) {
    data_list_file <- data_list_files[which.max(file.mtime(data_list_files))]
    data.list <- readRDS(data_list_file)
    message("Loaded data.list from: ", basename(data_list_file))
  } else {
    message("No data-list__*.rds files in ", path.objects, " (optional until preprocessing exports them).")
  }
} else {
  message(
    "FOLLOW-UP: Create ", path.objects,
    " and add ps-list__*.rds / data-list__*.rds when list-based workflows are ready."
  )
}

message("Environment setup complete. proj.path = ", proj.path)
message("Analysis ID: ", analysis.ID)
