# 05__WriteAuditManifest.R
#
# Created by: Michael Sieler
# Date last updated: 2026-04-22
#
# Description:
#   Writes a dated audit manifest (RDS + JSON) under Data/r_objects/ for before/after
#   comparisons when preprocessing or analysis outputs are regenerated.
#
# Expected input:
#   - Optional first CLI argument: label (e.g. "pre", "post"); default "snapshot".
#
# Expected output:
#   - Data/r_objects/audit__<label>__YYYY-MM-DD__HHMMSS.rds
#   - Data/r_objects/audit__<label>__YYYY-MM-DD__HHMMSS.json
#
# Usage:
#   Rscript Code/00__Setup/05__WriteAuditManifest.R pre
#   Rscript Code/00__Setup/05__WriteAuditManifest.R post

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(here))

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
  here::i_am("Code/00__Setup/05__WriteAuditManifest.R")
  proj.path <- here::here()
}

old_wd <- setwd(proj.path)
on.exit(setwd(old_wd), add = TRUE)

# Canonical init provides path.* plus the resolve_latest_rds helper.
source(file.path(proj.path, "Code", "00__Setup", "00__InitializeEnvironment.R"), local = FALSE)

label_arg <- commandArgs(trailingOnly = TRUE)
label <- if (length(label_arg) >= 1L && nzchar(label_arg[[1]])) {
  label_arg[[1]]
} else {
  "snapshot"
}
# Restrict label to safe filename fragment
label <- gsub("[^A-Za-z0-9._-]", "_", label)

stamp <- format(Sys.time(), "%Y-%m-%d__%H%M%S")
base_name <- paste0("audit__", label, "__", stamp)
out_dir <- path.objects
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

file_meta <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(list(path = path, exists = FALSE))
  }
  info <- file.info(path)
  list(
    path = path,
    exists = TRUE,
    basename = basename(path),
    size_bytes = unname(info$size),
    mtime_utc = format(info$mtime, tz = "UTC", usetz = TRUE)
  )
}

phyloseq_dims <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(list(loaded = FALSE, n_samples = NA_integer_, n_taxa = NA_integer_))
  }
  ps <- readRDS(path)
  if (!inherits(ps, "phyloseq")) {
    return(list(loaded = FALSE, n_samples = NA_integer_, n_taxa = NA_integer_, note = "not_phyloseq"))
  }
  list(
    loaded = TRUE,
    n_samples = phyloseq::nsamples(ps),
    n_taxa = phyloseq::ntaxa(ps)
  )
}

bundle_inventory <- function(results_root) {
  if (!dir.exists(results_root)) {
    return(tibble::tibble())
  }
  hits <- list.files(
    results_root,
    pattern = "__bundle\\.rds$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = FALSE
  )
  if (length(hits) == 0L) {
    return(tibble::tibble())
  }
  info <- file.info(hits)
  tibble::tibble(
    rel_path = gsub(
      paste0("^", normalizePath(results_root, winslash = "/", mustWork = TRUE), "/"),
      "",
      normalizePath(hits, winslash = "/", mustWork = TRUE)
    ),
    mtime_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    size_bytes = unname(info$size)
  ) |>
    dplyr::arrange(.data$rel_path)
}

rds_inventory <- function(dir, pattern) {
  if (!dir.exists(dir)) {
    return(tibble::tibble())
  }
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0L) {
    return(tibble::tibble())
  }
  info <- file.info(f)
  tibble::tibble(
    basename = basename(f),
    mtime_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
    size_bytes = unname(info$size)
  ) |>
    dplyr::arrange(.data$basename)
}

session_txt <- paste(capture.output(print(sessionInfo())), collapse = "\n")

manifest <- list(
  label = label,
  stamp = stamp,
  proj_path = normalizePath(proj.path, winslash = "/", mustWork = TRUE),
  R_version = R.version.string,
  session_info_text = session_txt,
  paths = list(
    path_dada2 = path.dada2,
    path_objects = path.objects,
    path_results = path.results
  ),
  phyloseq_rds = list(
    uncleaned = file_meta(path.pseq.uncleaned),
    cleaned = file_meta(path.pseq.cleaned),
    uncleaned_dims = phyloseq_dims(path.pseq.uncleaned),
    cleaned_dims = phyloseq_dims(path.pseq.cleaned)
  ),
  r_objects_inventory = list(
    ps_list = rds_inventory(path.objects, "^ps-list__.*\\.rds$"),
    data_list = rds_inventory(path.objects, "^data-list__.*\\.rds$"),
    decontam = rds_inventory(path.objects, "^decontam__.*\\.rds$"),
    beta_dist = rds_inventory(file.path(path.objects, "Rds"), "^beta\\.dist\\.mat__.*\\.rds$")
  ),
  results_bundle_inventory = bundle_inventory(path.results)
)

out_rds <- file.path(out_dir, paste0(base_name, ".rds"))
out_json <- file.path(out_dir, paste0(base_name, ".json"))
saveRDS(manifest, out_rds)

# JSON-friendly flatten: tibbles -> data.frames; keep session as one string.
manifest_json <- manifest
manifest_json$results_bundle_inventory <- as.data.frame(manifest$results_bundle_inventory)
for (nm in names(manifest_json$r_objects_inventory)) {
  manifest_json$r_objects_inventory[[nm]] <- as.data.frame(manifest_json$r_objects_inventory[[nm]])
}
jsonlite::write_json(
  manifest_json,
  out_json,
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null"
)

message("Wrote audit manifest:")
message("  RDS:  ", out_rds)
message("  JSON: ", out_json)

