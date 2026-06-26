# 99__BuildSupplementAll.R
#
# Created by: Michael Sieler
# Date last updated: 2026-06-25
#
# Description:
#   Optional: create missing per-module `supplement_map__NN__*.yml` from `Figures/` stems
#   (`--generate-missing-maps`) OR sync/overwrite all maps (`--sync-figure-maps`).
#   Default: build **per-module** supplementary Excel workbooks + per-module figure PDFs,
#   then a master INDEX manifest (and slim INDEX xlsx). Legacy combined `__ALL__` outputs
#   are opt-in via `--combined`, `--combined-tables`, or `--combined-figures`.
#
# Expected input:
#   - CLI: `<mode>` and optional flags (see Usage)
#
# Expected output:
#   - `Manuscript/Supplementary/_build/<mode>__<date>/` — per-module xlsx/pdf + manifests
#   - `Manuscript/Supplementary/Supplementary_Tables__INDEX__<mode>__*.xlsx` (root)
#   - `Manuscript/Supplementary/Sieler2026_SupplementaryTables__<mode>__<date>.zip` (+ figures zip)
#   - Optional: `Supplementary_Tables__ALL__*` / `Supplementary_Figures__ALL__*` under `_build/`

options(stringsAsFactors = FALSE)

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install here: install.packages('here')")
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Install yaml: install.packages('yaml')")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Install readr: install.packages('readr')")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Install dplyr: install.packages('dplyr')")
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
  here::i_am("Code/01__Analysis/99__BuildSupplementAll.R")
  proj.path <- as.character(here::here())
}

old_wd <- setwd(proj.path)
on.exit(setwd(old_wd), add = TRUE)

shared <- file.path(proj.path, "Code", "01__Analysis", "98__SupplementShared.R")
source(shared, local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
flag_gen <- "--generate-missing-maps" %in% args || "--generate-missing-map" %in% args
flag_sync <- "--sync-figure-maps" %in% args || "--sync-figure-map" %in% args
flag_allow_map_writes <- "--allow-map-writes" %in% args
flag_force_map_overwrite <- "--force-map-overwrite" %in% args
flag_combined <- "--combined" %in% args
flag_combined_tables <- "--combined-tables" %in% args || isTRUE(flag_combined)
flag_combined_figures <- "--combined-figures" %in% args || isTRUE(flag_combined)
args_mode <- args[!args %in% c(
  "--generate-missing-maps",
  "--generate-missing-map",
  "--sync-figure-maps",
  "--sync-figure-map",
  "--allow-map-writes",
  "--force-map-overwrite",
  "--combined",
  "--combined-tables",
  "--combined-figures"
)]
if (length(args_mode) < 1L) {
  stop(
    "Usage: Rscript Code/01__Analysis/99__BuildSupplementAll.R <mode> [flags]\n",
    "Flags:\n",
    "  --generate-missing-maps --allow-map-writes\n",
    "  --sync-figure-maps --allow-map-writes --force-map-overwrite\n",
    "  --combined-tables   (also build legacy __ALL__ Excel)\n",
    "  --combined-figures  (also build legacy __ALL__ figure PDF)\n",
    "  --combined          (both legacy __ALL__ outputs)\n",
    "Examples:\n",
    "  Rscript Code/01__Analysis/99__BuildSupplementAll.R submission\n",
    "  Rscript Code/01__Analysis/99__BuildSupplementAll.R submission --combined"
  )
}
mode <- args_mode[[1]]

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

path.results <- file.path(proj.path, "Results")
path.supp <- file.path(proj.path, "Manuscript", "Supplementary")
dir.create(path.supp, recursive = TRUE, showWarnings = FALSE)

mods <- sieler2026_supp_list_result_modules(proj.path)
stamp <- format(Sys.Date(), "%Y-%m-%d")
paths <- sieler2026_supp_paths(proj.path, mode, stamp)
sieler2026_supp_relocate_stray_build_artifacts(paths)

n_maps_written <- 0L

if (isTRUE(flag_gen)) {
  if (!isTRUE(flag_allow_map_writes)) {
    warning(
      "Requested `--generate-missing-maps` but map writes are disabled. ",
      "Skipping YAML generation. To enable: add `--allow-map-writes`.",
      call. = FALSE
    )
  } else {
    for (m in mods) {
      path.fig <- file.path(path.results, m, "Figures")
      map_path <- file.path(path.supp, paste0("supplement_map__", m, ".yml"))
      if (file.exists(map_path) || !dir.exists(path.fig)) {
        next
      }
      stems <- sieler2026_supp_list_figure_basenames(path.fig)
      if (length(stems) == 0L) {
        next
      }
      header <- c(
        paste0("# supplement_map__", m, ".yml"),
        "# Auto-generated by 99__BuildSupplementAll.R (--generate-missing-maps).",
        "# Figures are PDF-first stems under Results/<module>/Figures/; edit order or add per-figure `group:` as needed.",
        paste0("# Date: ", Sys.Date()),
        ""
      )
      body <- yaml::as.yaml(list(
        include_cover = TRUE,
        figures = as.list(stems)
      ))
      writeLines(c(header, body), map_path)
      n_maps_written <- n_maps_written + 1L
      message("Wrote default figure map: ", map_path)
    }
  }
}

if (isTRUE(flag_sync)) {
  if (!isTRUE(flag_allow_map_writes) || !isTRUE(flag_force_map_overwrite)) {
    warning(
      "Requested `--sync-figure-maps` but overwrite is not explicitly enabled. ",
      "Skipping YAML sync. To enable: add `--allow-map-writes --force-map-overwrite`.",
      call. = FALSE
    )
  } else {
    for (m in mods) {
      path.fig <- file.path(path.results, m, "Figures")
      map_path <- file.path(path.supp, paste0("supplement_map__", m, ".yml"))
      if (!dir.exists(path.fig)) {
        next
      }
      stems <- sieler2026_supp_list_figure_basenames(path.fig)
      if (length(stems) == 0L) {
        next
      }

      include_cover <- TRUE
      if (file.exists(map_path)) {
        old <- tryCatch(yaml::read_yaml(map_path), error = function(e) NULL)
        if (!is.null(old)) {
          include_cover <- sieler2026_supp_yaml_truthy(old$include_cover, default = TRUE)
        }
      }

      header <- c(
        paste0("# supplement_map__", m, ".yml"),
        "# Synced by 99__BuildSupplementAll.R (--sync-figure-maps).",
        "# Figures are PDF-first stems under Results/<module>/Figures/; edit order or add per-figure `group:` as needed.",
        paste0("# Date: ", Sys.Date()),
        ""
      )
      body <- yaml::as.yaml(list(
        include_cover = include_cover,
        figures = as.list(stems)
      ))
      writeLines(c(header, body), map_path)
      n_maps_written <- n_maps_written + 1L
      message("Synced figure map: ", map_path)
    }
  }
}

build_script <- normalizePath(
  file.path(proj.path, "Code", "01__Analysis", "99__SupplementBuild.R"),
  winslash = "/",
  mustWork = TRUE
)

run_delegate <- function(subcmd, ...) {
  st <- system2(rscript, c(build_script, subcmd, ...))
  if (!is.numeric(st) || length(st) != 1L || st != 0L) {
    stop(
      "99__SupplementBuild.R ",
      subcmd,
      " failed (exit ",
      st,
      "): ",
      paste(c(subcmd, ...), collapse = " ")
    )
  }
  invisible(TRUE)
}

# ----- Per-module tables (default) -------------------------------------------

modules_with_tables <- sieler2026_supp_modules_with_tables(proj.path)
if (length(modules_with_tables) < 1L) {
  stop("No module with Tables/*.csv found; nothing to build.")
}

message(
  "Building per-module supplement tables (mode=",
  mode,
  ", ",
  length(modules_with_tables),
  " module(s)) ..."
)
for (m in modules_with_tables) {
  message("  excel: ", m)
  run_delegate("excel", m, mode)
}

# ----- Per-module figures (default) ------------------------------------------

message("Building per-module supplement figures where maps exist ...")
for (m in mods) {
  map_path <- file.path(path.supp, paste0("supplement_map__", m, ".yml"))
  if (!file.exists(map_path)) {
    next
  }
  message("  figures-pdf: ", m)
  run_delegate("figures-pdf", m, mode)
}

# ----- Master INDEX (manifest + slim xlsx) -----------------------------------

index_manifest_paths <- file.path(
  paths$build,
  paste0(
    "Supplementary_Tables__",
    modules_with_tables,
    "__",
    mode,
    "__",
    stamp,
    "__manifest.csv"
  )
)
names(index_manifest_paths) <- modules_with_tables

missing_manifests <- index_manifest_paths[!file.exists(index_manifest_paths)]
if (length(missing_manifests) > 0L) {
  stop(
    "Expected per-module manifest(s) missing after excel build:\n",
    paste(missing_manifests, collapse = "\n")
  )
}

index_parts <- lapply(modules_with_tables, function(m) {
  mp <- index_manifest_paths[[m]]
  wb_file <- paste0(
    "Supplementary_Tables__",
    m,
    "__",
    mode,
    "__",
    stamp,
    ".xlsx"
  )
  df <- readr::read_csv(mp, show_col_types = FALSE)
  df$workbook_file <- wb_file
  df
})
index_tbl <- dplyr::bind_rows(index_parts) |>
  dplyr::arrange(.data$module_index, .data$group_index, .data$within_group_idx)

out_index_manifest <- file.path(
  paths$build,
  paste0("Supplementary_Tables__INDEX__", mode, "__", stamp, "__manifest.csv")
)
readr::write_csv(index_tbl, out_index_manifest)

meta_index <- tibble::tibble(
  key = c(
    "build_date_utc",
    "git_commit",
    "R_version",
    "openxlsx",
    "readr",
    "yaml",
    "bundle",
    "mode",
    "modules_included",
    "tables_scope",
    "software_session_info_csv"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    sieler2026_supp_git_sha(proj.path),
    R.version.string,
    sieler2026_supp_pkg_version("openxlsx"),
    sieler2026_supp_pkg_version("readr"),
    sieler2026_supp_pkg_version("yaml"),
    "INDEX",
    mode,
    paste(modules_with_tables, collapse = ", "),
    "Master TOC across per-module supplementary table workbooks",
    sieler2026_session_info_csv_basename(mode, stamp)
  )
)

out_index_xlsx <- file.path(
  paths$build,
  paste0("Supplementary_Tables__INDEX__", mode, "__", stamp, ".xlsx")
)
sieler2026_supp_write_index_xlsx(index_tbl, meta_index, out_index_xlsx)

message("Wrote master INDEX manifest: ", out_index_manifest)
message("Wrote master INDEX xlsx:     ", out_index_xlsx)

# ----- Optional legacy combined outputs --------------------------------------

if (isTRUE(flag_combined_tables) || isTRUE(flag_combined_figures)) {
  combined_args <- c("combined", mode)
  if (isTRUE(flag_combined_tables) && !isTRUE(flag_combined_figures)) {
    combined_args <- c(combined_args, "--tables-only")
  }
  if (isTRUE(flag_combined_figures) && !isTRUE(flag_combined_tables)) {
    combined_args <- c(combined_args, "--figures-only")
  }
  message("Running legacy combined supplement build ...")
  st <- system2(rscript, c(build_script, combined_args))
  if (st != 0L) {
    stop("99__SupplementBuild.R combined exited with code ", st)
  }
}

pkg_res <- sieler2026_supp_finalize_packages(
  proj.path,
  mode,
  stamp,
  combined = isTRUE(flag_combined_tables) || isTRUE(flag_combined_figures)
)

message(
  "Done. Submission packages at ",
  path.supp,
  " (tables zip: ",
  basename(pkg_res$tables_zip),
  if (!is.null(pkg_res$figures_zip)) {
    paste0("; figures zip: ", basename(pkg_res$figures_zip))
  } else {
    ""
  },
  "). Loose artifacts under ",
  paths$build,
  ". Default maps written this run: ",
  n_maps_written
)
