# 99__SessionInfoExport.R
# Created by: Michael Sieler
# Date last updated: 2026-06-26
#
# Description: Export unified software/session provenance CSV for supplementary materials.
#   Captures environment metadata, all attached/loaded packages, and analysis-driver rows
#   for modules 01–08. Also invoked automatically by sieler2026_supp_finalize_packages().
#
# Expected input:
#   - CLI: `<mode>` (e.g. `submission`); optional `--stamp=YYYY-MM-DD`
#   - Project root via Rscript --file= resolution
#
# Expected output:
#   - Manuscript/Supplementary/_build/<mode>__<stamp>/Software_SessionInfo__<mode>__<stamp>.csv
#   - Copy at Manuscript/Supplementary/Software_SessionInfo__<mode>__<stamp>.csv

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
  here::i_am("Code/01__Analysis/99__SessionInfoExport.R")
  proj.path <- here::here()
}

old_wd <- setwd(proj.path)
on.exit(setwd(old_wd), add = TRUE)

shared <- file.path(proj.path, "Code", "01__Analysis", "98__SupplementShared.R")
source(shared, local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
stamp_arg <- grep("^--stamp=", args, value = TRUE)
stamp <- if (length(stamp_arg) >= 1L) {
  sub("^--stamp=", "", stamp_arg[[1L]])
} else {
  format(Sys.Date(), "%Y-%m-%d")
}
mode_args <- args[!grepl("^--", args)]
if (length(mode_args) < 1L) {
  stop(
    "Usage: Rscript Code/01__Analysis/99__SessionInfoExport.R <mode> [--stamp=YYYY-MM-DD]\n",
    "Example: Rscript Code/01__Analysis/99__SessionInfoExport.R submission"
  )
}
mode <- mode_args[[1L]]

out_csv <- sieler2026_export_session_info_csv(proj.path, mode, stamp)

if (requireNamespace("gt", quietly = TRUE)) {
  tbl <- readr::read_csv(out_csv, show_col_types = FALSE)
  key_pkgs <- c(
    "phyloseq", "vegan", "glmmTMB", "DESeq2", "maaslin3",
    "clusterProfiler", "nptest", "microViz", "emmeans"
  )
  pkg_preview <- tbl |>
    dplyr::filter(.data$record_type == "package", .data$key %in% key_pkgs) |>
    dplyr::arrange(.data$key)
  if (nrow(pkg_preview) >= 1L) {
    print(
      gt::gt(pkg_preview) |>
        gt::tab_header(
          title = "Key analysis packages (preview)",
          subtitle = basename(out_csv)
        )
    )
  }
}

message("99__SessionInfoExport.R complete.")
