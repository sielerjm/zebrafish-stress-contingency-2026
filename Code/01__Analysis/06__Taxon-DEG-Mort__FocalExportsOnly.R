# 06__Taxon-DEG-Mort__FocalExportsOnly.R
# Created by: Michael Sieler
# Date last updated: 2026-04-26
#
# Description:
#   Thin wrapper around `06__Taxon-DEG-Mort.R --focal-exports-only` for backwards-compatible paths
#   and documentation. Refreshes focal-genus edge CSVs and focal-four gene-set Venns from existing
#   `combined_sig_partial_correlations.csv` without rerunning partial-correlation inference.
#
# Expected input: Run from Sieler2026 repository root (same as main module 06 driver).
# Expected output: Same as `Rscript Code/01__Analysis/06__Taxon-DEG-Mort.R --focal-exports-only`.

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
if (is.na(.sp)) {
  stop(
    "Run via:\n",
    "  Rscript Code/01__Analysis/06__Taxon-DEG-Mort__FocalExportsOnly.R\n",
    "from the project root (so --file= resolves), or use:\n",
    "  Rscript Code/01__Analysis/06__Taxon-DEG-Mort.R --focal-exports-only"
  )
}

proj.path <- normalizePath(file.path(dirname(.sp), "..", ".."), winslash = "/", mustWork = TRUE)
main <- normalizePath(
  file.path(proj.path, "Code", "01__Analysis", "06__Taxon-DEG-Mort.R"),
  winslash = "/",
  mustWork = TRUE
)

rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) {
  rscript <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows") {
    rscript <- paste0(rscript, ".exe")
  }
}
if (!nzchar(rscript) || !file.exists(rscript)) {
  stop("Could not find Rscript executable.")
}

old_wd <- setwd(proj.path)
on.exit(setwd(old_wd), add = TRUE)

st <- system2(rscript, c(main, "--focal-exports-only"), stdout = "", stderr = "", wait = TRUE)
if (!is.numeric(st) || length(st) != 1L) {
  stop("system2 returned unexpected status forwarding to 06__Taxon-DEG-Mort.R")
}
if (st != 0L) {
  stop("06__Taxon-DEG-Mort.R --focal-exports-only exited with code ", st)
}
