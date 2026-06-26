# 99__BuildSupplementExcel.R
#
# Created by: Michael Sieler
# Date last updated: 2026-06-25
#
# Description:
#   Assemble a supplementary Excel workbook from `Results/<module>/Tables/*.csv` with:
#   - META sheet (run provenance)
#   - TOC sheet (stable `supp_id` + `Table S<supp_id>` citation label)
#   - one sheet per CSV (worksheet name = sanitized `supp_id`, e.g. `1.2.3`; ID-like columns as text)
#   - manifest CSV alongside workbook for grep/search
#
# Expected input:
#   - CLI: `<module_dir>` `<mode>` e.g. `01__Diversity submission`
#   - Optional overrides: `Manuscript/Supplementary/supplement_overrides__<module>.yml`
#
# Expected output:
#   - `Manuscript/Supplementary/_build/<mode>__<date>/Supplementary_Tables__<module>__<mode>__YYYY-MM-DD.xlsx`
#   - `Manuscript/Supplementary/_build/<mode>__<date>/Supplementary_Tables__<module>__<mode>__YYYY-MM-DD__manifest.csv`

options(stringsAsFactors = FALSE)

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install here: install.packages('here')")
}
if (!requireNamespace("tibble", quietly = TRUE)) {
  stop("Install tibble: install.packages('tibble')")
}

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
  here::i_am("Code/01__Analysis/99__BuildSupplementExcel.R")
  proj.path <- as.character(here::here())
}

old_wd <- setwd(proj.path)
on.exit(setwd(old_wd), add = TRUE)

shared <- file.path(proj.path, "Code", "01__Analysis", "98__SupplementShared.R")
source(shared, local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript Code/01__Analysis/99__BuildSupplementExcel.R <module_dir> <mode>\n",
    "Example: Rscript Code/01__Analysis/99__BuildSupplementExcel.R 01__Diversity submission"
  )
}

module_dir <- args[[1]]
mode <- args[[2]]

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Install openxlsx: install.packages('openxlsx')")
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Install yaml: install.packages('yaml')")
}

path.manuscript_supp <- file.path(proj.path, "Manuscript", "Supplementary")
dir.create(path.manuscript_supp, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.Date(), "%Y-%m-%d")
paths <- sieler2026_supp_paths(proj.path, mode, stamp)

tr <- sieler2026_supp_build_table_registry(
  proj.path,
  module_dir,
  path.manuscript_supp,
  require_csv = TRUE
)
base_tbl <- tr$registry

meta_tbl <- tibble::tibble(
  key = c(
    "build_date_utc",
    "git_commit",
    "R_version",
    "openxlsx",
    "readr",
    "yaml",
    "module",
    "mode",
    "tables_dir"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    sieler2026_supp_git_sha(proj.path),
    R.version.string,
    sieler2026_supp_pkg_version("openxlsx"),
    sieler2026_supp_pkg_version("readr"),
    sieler2026_supp_pkg_version("yaml"),
    module_dir,
    mode,
    sieler2026_supp_path_for_workbook(tr$path_tables, proj.path)
  )
)

out_xlsx <- file.path(
  paths$build,
  paste0("Supplementary_Tables__", module_dir, "__", mode, "__", stamp, ".xlsx")
)
out_manifest <- file.path(
  paths$build,
  paste0("Supplementary_Tables__", module_dir, "__", mode, "__", stamp, "__manifest.csv")
)

res <- sieler2026_supp_write_excel_workbook(
  base_tbl,
  meta_tbl,
  out_xlsx,
  out_manifest,
  proj.path
)

message("Wrote supplement Excel: ", res$out_xlsx)
message("Wrote manifest CSV:    ", res$out_manifest)
