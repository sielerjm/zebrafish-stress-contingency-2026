# 99__SupplementBuild.R
#
# Created by: Michael Sieler
# Date last updated: 2026-06-25
#
# Description:
#   Single documented entry point for supplementary table/figure assembly. Dispatches to the
#   existing `99__BuildSupplement*.R` scripts (same behavior as calling them directly).
#
# Expected input: CLI subcommand and arguments (see Usage below). Run from repo root or via Rscript
#   with `--file=` so the project root resolves.
#
# Expected output: Submission zips + INDEX xlsx at `Manuscript/Supplementary/` root;
#   loose artifacts under `Manuscript/Supplementary/_build/<mode>__<date>/`.

options(stringsAsFactors = FALSE)

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
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Install here: install.packages('here')")
  }
  here::i_am("Code/01__Analysis/99__SupplementBuild.R")
  proj.path <- as.character(here::here())
}

script_dir <- file.path(proj.path, "Code", "01__Analysis")

rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) {
  rscript <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows") {
    rscript <- paste0(rscript, ".exe")
  }
}
if (!nzchar(rscript) || !file.exists(rscript)) {
  stop("Could not find Rscript executable; install R or add Rscript to PATH.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(
    "Usage:\n",
    "  Rscript Code/01__Analysis/99__SupplementBuild.R combined <mode> [--tables-only] [--figures-only]\n",
    "  Rscript Code/01__Analysis/99__SupplementBuild.R all <mode> [--generate-missing-maps] [--sync-figure-maps] [--allow-map-writes] [--force-map-overwrite] [--combined] [--combined-tables] [--combined-figures]\n",
    "  Rscript Code/01__Analysis/99__SupplementBuild.R excel <module_dir> <mode>\n",
    "  Rscript Code/01__Analysis/99__SupplementBuild.R figures-pdf <module_dir> <mode>\n",
    "Default `all`: per-module Excel + figure PDFs in `_build/`, INDEX xlsx + zip bundles at supplementary root. Legacy __ALL__ via --combined* flags.\n",
    "Legacy equivalents: 99__BuildSupplementCombined.R, 99__BuildSupplementAll.R, ...",
    call. = FALSE
  )
}

cmd <- args[[1L]]
tail_args <- args[-1L]

legacy <- switch(
  cmd,
  combined = "99__BuildSupplementCombined.R",
  all = "99__BuildSupplementAll.R",
  excel = "99__BuildSupplementExcel.R",
  `figures-pdf` = "99__BuildSupplementFiguresPdf.R",
  stop("Unknown subcommand: ", encodeString(cmd), call. = FALSE)
)

target <- normalizePath(file.path(script_dir, legacy), winslash = "/", mustWork = TRUE)

old_wd <- setwd(proj.path)
on.exit(setwd(old_wd), add = TRUE)

st <- system2(rscript, c(target, tail_args), stdout = "", stderr = "", wait = TRUE)
if (!is.numeric(st) || length(st) != 1L) {
  stop("system2 returned unexpected status for ", legacy)
}
if (st != 0L) {
  stop(legacy, " exited with code ", st)
}
