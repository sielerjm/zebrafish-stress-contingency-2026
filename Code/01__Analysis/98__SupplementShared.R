# 98__SupplementShared.R
#
# Created by: Michael Sieler
# Date last updated: 2026-06-26
#
# Description:
#   Small shared helpers for supplementary Excel / PDF assembly (sourced by 99__* scripts).
#
# Expected input:  none (sourced).
# Expected output: function definitions in the global environment.

#' Paths for supplementary workbooks: `/<project_basename>/Results/...` (no home directory).
#'
#' @param abs_path Absolute file or directory path.
#' @param proj_path Absolute project root (e.g. `here::here()`).
sieler2026_supp_path_for_workbook <- function(abs_path, proj_path) {
  proj_path <- normalizePath(proj_path, winslash = "/", mustWork = FALSE)
  ap <- normalizePath(abs_path, winslash = "/", mustWork = FALSE)
  proj_nm <- basename(proj_path)
  ncp <- nchar(proj_path, type = "chars", allowNA = TRUE)
  if (is.na(ncp) || ncp < 1L) {
    return(as.character(abs_path))
  }
  if (!startsWith(ap, proj_path)) {
    return(ap)
  }
  rel <- if (nchar(ap, type = "chars") <= ncp) {
    ""
  } else {
    substring(ap, ncp + 2L)
  }
  if (!nzchar(rel)) {
    return(paste0("/", proj_nm))
  }
  paste0("/", proj_nm, "/", rel)
}

#' Leading two-digit module code from folder name (e.g. `01__Diversity` -> 1L).
sieler2026_supp_module_index <- function(module_dir) {
  m <- regmatches(module_dir, regexpr("^\\d{2}", module_dir, perl = TRUE))
  if (length(m) != 1L || !nzchar(m)) {
    return(NA_integer_)
  }
  as.integer(m)
}

#' Infer a stable "group" string from a CSV basename without extension.
sieler2026_supp_infer_group <- function(basename_no_ext) {
  parts <- strsplit(basename_no_ext, "__", fixed = TRUE)[[1]]
  if (length(parts) >= 2L) {
    paste(parts[-length(parts)], collapse = "__")
  } else {
    basename_no_ext
  }
}

#' Sanitize Excel worksheet names (31 chars max; illegal chars removed).
sieler2026_supp_sanitize_sheet_name <- function(name, used = character()) {
  x <- name
  x <- gsub("[:\\\\/?*\\[\\]]", "_", x)
  x <- gsub("\\s+", "_", x, perl = TRUE)
  x <- substring(x, 1L, 31L)
  if (!nzchar(x)) {
    x <- "SHEET"
  }
  base <- x
  k <- 1L
  while (x %in% used) {
    suffix <- paste0("_", k)
    max_base <- 31L - nchar(suffix, type = "chars", allowNA = TRUE)
    if (is.na(max_base) || max_base < 1L) {
      max_base <- 1L
    }
    x <- paste0(substring(base, 1L, max_base), suffix)
    k <- k + 1L
  }
  x
}

#' Column names likely to be corrupted by Excel if treated as numeric.
sieler2026_supp_idish_col <- function(col_name) {
  cn <- tolower(col_name)
  grepl(
    "gene|ensembl|entrez|symbol|taxon|asv|otu|feature|^id$|_id$|sample|worm|tank",
    cn,
    ignore.case = TRUE,
    perl = TRUE
  )
}

#' Coerce YAML-parsed values to logical (handles `yes` / `true` from `yaml::as.yaml()`).
sieler2026_supp_yaml_truthy <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  if (is.logical(x) && length(x) == 1L) {
    return(isTRUE(x))
  }
  if (is.character(x) && length(x) == 1L) {
    v <- tolower(trimws(x))
    return(v %in% c("true", "yes", "y", "1"))
  }
  if (is.numeric(x) && length(x) == 1L) {
    return(!is.na(x) && x != 0)
  }
  default
}

#' Basenames for supplement figure maps: prefer `.pdf` stems; include `.png` only when no PDF of same name.
#'
#' @param fig_dir `Results/<module>/Figures` (must exist).
sieler2026_supp_list_figure_basenames <- function(fig_dir) {
  if (!dir.exists(fig_dir)) {
    return(character(0))
  }
  pdfs <- list.files(fig_dir, pattern = "\\.pdf$", ignore.case = TRUE)
  stems_pdf <- unique(tools::file_path_sans_ext(pdfs))
  pngs <- list.files(fig_dir, pattern = "\\.png$", ignore.case = TRUE)
  stems_png <- unique(tools::file_path_sans_ext(pngs))
  only_png <- setdiff(stems_png, stems_pdf)
  sort(unique(c(stems_pdf, only_png)))
}

#' Resolve a figure basename (no extension) to a PDF path, else PNG path.
sieler2026_supp_resolve_figure_path_pdf_first <- function(fig_dir, basename_no_ext) {
  pdf_path <- file.path(fig_dir, paste0(basename_no_ext, ".pdf"))
  png_path <- file.path(fig_dir, paste0(basename_no_ext, ".png"))
  if (file.exists(pdf_path)) {
    return(list(path = pdf_path, ext = "pdf"))
  }
  if (file.exists(png_path)) {
    return(list(path = png_path, ext = "png"))
  }
  list(path = NA_character_, ext = NA_character_)
}

#' Normalize `figures:` entries from `supplement_map__*.yml` to a tibble.
sieler2026_supp_normalize_map_figures <- function(figures) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Install tibble: install.packages('tibble')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install dplyr: install.packages('dplyr')")
  }
  if (is.null(figures) || length(figures) == 0L) {
    return(tibble::tibble())
  }
  out <- list()
  for (item in figures) {
    if (is.character(item) && length(item) == 1L && nzchar(item)) {
      out[[length(out) + 1L]] <- tibble::tibble(
        basename = item,
        caption = NA_character_,
        yaml_group = NA_character_,
        order_within_group = NA_integer_,
        supp_id_override = NA_character_
      )
    } else if (is.list(item)) {
      bn <- item$basename
      if (is.null(bn) || !nzchar(as.character(bn))) {
        stop("Each `figures:` list entry must be a string or a list with `basename:`.")
      }
      cap <- item$caption
      yg <- item$group
      owg <- item$order_within_group
      sid <- item$supp_id_override
      out[[length(out) + 1L]] <- tibble::tibble(
        basename = as.character(bn),
        caption = if (is.null(cap)) NA_character_ else as.character(cap),
        yaml_group = if (is.null(yg)) NA_character_ else as.character(yg),
        order_within_group = if (is.null(owg)) {
          NA_integer_
        } else {
          suppressWarnings(as.integer(owg))
        },
        supp_id_override = if (is.null(sid)) {
          NA_character_
        } else {
          as.character(sid)
        }
      )
    } else {
      stop("Unsupported `figures:` entry type.")
    }
  }
  dplyr::bind_rows(out)
}

#' Read optional `supplement_overrides__<module>.yml` (tables section only).
sieler2026_supp_read_table_overrides <- function(path_manuscript_supp, module_dir) {
  overrides_path <- file.path(
    path_manuscript_supp,
    paste0("supplement_overrides__", module_dir, ".yml")
  )
  if (file.exists(overrides_path)) {
    yaml::read_yaml(overrides_path)
  } else {
    list()
  }
}

#' Read a CSV/TSV for Excel export (ID-like and very large numeric columns as text).
sieler2026_supp_read_csv_for_excel <- function(path) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Install readr: install.packages('readr')")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Install tibble: install.packages('tibble')")
  }
  ext <- tolower(tools::file_ext(path))
  df <- if (identical(ext, "tsv")) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }
  df <- tibble::as_tibble(df)
  nm <- colnames(df)
  for (j in seq_along(nm)) {
    col <- df[[j]]
    if (sieler2026_supp_idish_col(nm[[j]])) {
      df[[j]] <- as.character(col)
      next
    }
    if (is.numeric(col)) {
      mx <- suppressWarnings(max(abs(col), na.rm = TRUE))
      if (is.finite(mx) && mx >= 1e9) {
        df[[j]] <- as.character(col)
      }
    }
  }
  df
}

#' Add workbook-relative paths and unique Excel sheet names to a table registry tibble.
sieler2026_supp_prepare_registry_workbook_cols <- function(base_tbl, proj_path) {
  if (nrow(base_tbl) < 1L) {
    return(base_tbl)
  }
  base_tbl$csv_path_workbook <- vapply(
    base_tbl$csv_path,
    function(p) sieler2026_supp_path_for_workbook(p, proj_path),
    character(1L)
  )
  used <- c("META", "TOC")
  sheet_names <- character(nrow(base_tbl))
  for (i in seq_len(nrow(base_tbl))) {
    sn <- sieler2026_supp_sanitize_sheet_name(as.character(base_tbl$supp_id[[i]]), used = used)
    used <- c(used, sn)
    sheet_names[[i]] <- sn
  }
  base_tbl$sheet_name <- sheet_names
  base_tbl
}

#' Write META + TOC + one sheet per registry row; export manifest CSV.
#'
#' @return Named list: `out_xlsx`, `out_manifest`, `manifest_tbl` (paths on disk).
sieler2026_supp_write_excel_workbook <- function(
    base_tbl,
    meta_tbl,
    out_xlsx,
    out_manifest,
    proj_path
) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Install openxlsx: install.packages('openxlsx')")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Install readr: install.packages('readr')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install dplyr: install.packages('dplyr')")
  }
  if (nrow(base_tbl) < 1L) {
    stop("Cannot write supplement Excel: registry has zero rows.")
  }

  base_tbl <- sieler2026_supp_prepare_registry_workbook_cols(base_tbl, proj_path)

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "META")
  openxlsx::writeData(wb, "META", meta_tbl, withFilter = TRUE)
  openxlsx::freezePane(wb, "META", firstActiveRow = 2, firstActiveCol = 1)

  toc <- base_tbl |>
    dplyr::transmute(
      supp_id = .data$supp_id,
      citation_table = .data$citation_table,
      module = .data$module_dir,
      group = .data$group,
      csv_basename = .data$csv_basename,
      sheet_name = .data$sheet_name,
      csv_path = .data$csv_path_workbook
    )
  openxlsx::addWorksheet(wb, "TOC")
  openxlsx::writeData(wb, "TOC", toc, withFilter = TRUE)
  openxlsx::freezePane(wb, "TOC", firstActiveRow = 2, firstActiveCol = 1)

  for (i in seq_len(nrow(base_tbl))) {
    sh <- base_tbl$sheet_name[[i]]
    openxlsx::addWorksheet(wb, sh)
    df <- sieler2026_supp_read_csv_for_excel(base_tbl$csv_path[[i]])
    openxlsx::writeData(wb, sh, df, withFilter = TRUE)
    openxlsx::freezePane(wb, sh, firstActiveRow = 2, firstActiveCol = 1)
  }

  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

  manifest_tbl <- base_tbl
  manifest_tbl$csv_path <- manifest_tbl$csv_path_workbook
  manifest_tbl$csv_path_workbook <- NULL
  readr::write_csv(manifest_tbl, out_manifest)

  invisible(list(
    out_xlsx = out_xlsx,
    out_manifest = out_manifest,
    manifest_tbl = manifest_tbl
  ))
}

#' Slim INDEX workbook (META + master TOC only) from a combined manifest tibble.
sieler2026_supp_write_index_xlsx <- function(index_tbl, meta_tbl, out_xlsx) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Install openxlsx: install.packages('openxlsx')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install dplyr: install.packages('dplyr')")
  }
  toc <- index_tbl |>
    dplyr::transmute(
      supp_id = .data$supp_id,
      citation_table = .data$citation_table,
      module = .data$module_dir,
      group = .data$group,
      csv_basename = .data$csv_basename,
      sheet_name = .data$sheet_name,
      workbook_file = .data$workbook_file,
      csv_path = .data$csv_path
    )
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "META")
  openxlsx::writeData(wb, "META", meta_tbl, withFilter = TRUE)
  openxlsx::freezePane(wb, "META", firstActiveRow = 2, firstActiveCol = 1)
  openxlsx::addWorksheet(wb, "TOC")
  openxlsx::writeData(wb, "TOC", toc, withFilter = TRUE)
  openxlsx::freezePane(wb, "TOC", firstActiveRow = 2, firstActiveCol = 1)
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  invisible(out_xlsx)
}

#' Sorted `Results/NN__*/` module folder names (excludes `_Archive`).
sieler2026_supp_list_result_modules <- function(proj_path) {
  path.results <- file.path(proj_path, "Results")
  mods <- list.dirs(path.results, full.names = FALSE, recursive = FALSE)
  mods <- mods[grepl("^\\d{2}__", mods)]
  mods <- setdiff(mods, "_Archive")
  sort(mods)
}

#' Modules with at least one `Tables/*.csv` file.
sieler2026_supp_modules_with_tables <- function(proj_path) {
  mods <- sieler2026_supp_list_result_modules(proj_path)
  mods[vapply(mods, function(m) {
    pt <- file.path(proj_path, "Results", m, "Tables")
    dir.exists(pt) &&
      length(list.files(pt, pattern = "\\.csv$", ignore.case = TRUE)) >= 1L
  }, logical(1L))]
}

#' Git commit SHA for supplement META (best effort).
sieler2026_supp_git_sha <- function(proj_path) {
  tryCatch(
    {
      sha <- system2(
        "git",
        c("-C", proj_path, "rev-parse", "HEAD"),
        stdout = TRUE,
        stderr = FALSE
      )
      if (length(sha) == 1L && nzchar(sha[[1]])) sha[[1]] else NA_character_
    },
    warning = function(w) NA_character_,
    error = function(e) NA_character_
  )
}

#' Installed package version string for supplement META.
sieler2026_supp_pkg_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NA_character_)
  }
  as.character(utils::packageVersion(pkg))
}

#' Parse sessionInfo() package blocks into name/version rows (R 4.5+ list or legacy character).
sieler2026_session_info_parse_pkg_strings <- function(pkg_obj, group_label) {
  empty <- tibble::tibble(
    record_type = character(0),
    key = character(0),
    value = character(0),
    group = character(0),
    notes = character(0)
  )
  if (length(pkg_obj) < 1L) {
    return(empty)
  }

  if (is.list(pkg_obj)) {
    pkg_names <- vapply(
      pkg_obj,
      function(desc) {
        if (is.null(desc[["Package"]])) {
          return(NA_character_)
        }
        as.character(desc[["Package"]])
      },
      character(1)
    )
    pkg_versions <- vapply(
      pkg_obj,
      function(desc) {
        if (is.null(desc[["Version"]])) {
          return(NA_character_)
        }
        as.character(desc[["Version"]])
      },
      character(1)
    )
    if (any(is.na(pkg_names)) && !is.null(names(pkg_obj))) {
      pkg_names[is.na(pkg_names)] <- names(pkg_obj)[is.na(pkg_names)]
    }
    return(tibble::tibble(
      record_type = "package",
      key = pkg_names,
      value = pkg_versions,
      group = group_label,
      notes = NA_character_
    ))
  }

  if (!is.character(pkg_obj)) {
    return(empty)
  }

  # basePkgs: package names only (no embedded version string).
  pkg_versions <- vapply(
    pkg_obj,
    function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        return(NA_character_)
      }
      as.character(utils::packageVersion(pkg))
    },
    character(1)
  )
  tibble::tibble(
    record_type = "package",
    key = pkg_obj,
    value = pkg_versions,
    group = group_label,
    notes = NA_character_
  )
}

#' Canonical analysis-driver metadata for modules 01–08 (submission software provenance).
sieler2026_session_info_analysis_driver_registry <- function(proj_path) {
  proj_path <- normalizePath(proj_path, winslash = "/", mustWork = TRUE)
  rel <- function(p) sub(paste0("^", proj_path, "/?"), "", normalizePath(p, winslash = "/", mustWork = FALSE))

  drivers <- list(
    list(
      module = "01__Diversity",
      driver = "Code/01__Analysis/01__Diversity.R",
      bundle_glob = "*diversity*__bundle.rds",
      primary_packages = "glmmTMB; emmeans; car"
    ),
    list(
      module = "02__Composition",
      driver = "Code/01__Analysis/02__Composition.R",
      bundle_glob = "*composition*__bundle.rds",
      primary_packages = "vegan; microViz"
    ),
    list(
      module = "03__DiffAbund",
      driver = "Code/01__Analysis/03__DiffAbund.R",
      bundle_glob = "*diffabund*__bundle.rds",
      primary_packages = "maaslin3"
    ),
    list(
      module = "04__DiffGeneExp",
      driver = "Code/01__Analysis/04__DiffGeneExp.R",
      bundle_glob = "*diffgeneexp*__bundle.rds",
      primary_packages = "DESeq2"
    ),
    list(
      module = "05__Mort-Inf",
      driver = "Code/01__Analysis/05__Mort-Inf.R",
      bundle_glob = "*mortinf*__bundle.rds",
      primary_packages = "glmmTMB; emmeans"
    ),
    list(
      module = "06__Taxon-DEG-Mort",
      driver = "Code/01__Analysis/06__Taxon-DEG-Mort.R",
      bundle_glob = "*taxon_deg_mort*__bundle.rds",
      primary_packages = "nptest; DESeq2"
    ),
    list(
      module = "07__FunctionalAnno",
      driver = "Code/01__Analysis/07__FunctionalAnno.R",
      bundle_glob = "*functional_anno*__bundle.rds",
      primary_packages = "clusterProfiler; org.Dr.eg.db"
    ),
    list(
      module = "08__NeutralModel",
      driver = "Code/01__Analysis/08__NeutralModel.R",
      bundle_glob = "*neutral_model*__bundle.rds",
      primary_packages = "minpack.lm; Hmisc"
    )
  )

  dplyr::bind_rows(lapply(drivers, function(d) {
    stats_dir <- file.path(proj_path, "Results", d$module, "Stats")
    bundle_hits <- character(0)
    if (dir.exists(stats_dir)) {
      bundle_hits <- list.files(stats_dir, pattern = gsub("\\*", ".*", d$bundle_glob), full.names = TRUE)
    }
    bundle_rel <- if (length(bundle_hits) >= 1L) rel(bundle_hits[[1L]]) else NA_character_
    tibble::tibble(
      record_type = c("analysis_driver", "analysis_driver"),
      key = d$module,
      value = c(d$driver, bundle_rel),
      group = c("driver_script", "bundle_rds"),
      notes = c(d$primary_packages, NA_character_)
    )
  }))
}

#' Load project analysis libraries if not already attached (for sessionInfo capture).
sieler2026_session_info_ensure_analysis_env <- function(proj_path) {
  if ("package:phyloseq" %in% search()) {
    return(invisible(FALSE))
  }
  init_script <- file.path(proj_path, "Code", "00__Setup", "00__InitializeEnvironment.R")
  if (!file.exists(init_script)) {
    stop("Missing init script: ", init_script)
  }
  source(init_script, local = FALSE)
  invisible(TRUE)
}

#' Unified long-format session / software provenance tibble.
#'
#' @param proj_path Absolute project root.
#' @return Tibble with columns `record_type`, `key`, `value`, `group`, `notes`.
sieler2026_session_info_to_tibble <- function(proj_path) {
  proj_path <- normalizePath(proj_path, winslash = "/", mustWork = TRUE)
  sieler2026_session_info_ensure_analysis_env(proj_path)

  si <- utils::sessionInfo()
  captured_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  git_sha <- sieler2026_supp_git_sha(proj_path)
  locale_str <- paste(si$locale, collapse = "; ")

  env_tbl <- tibble::tibble(
    record_type = "environment",
    key = c(
      "r_version",
      "platform",
      "running_os",
      "locale",
      "timezone",
      "captured_at_utc",
      "proj_path",
      "git_commit"
    ),
    value = c(
      si$R.version$version.string,
      si$platform,
      si$running,
      locale_str,
      Sys.timezone(),
      captured_at,
      proj_path,
      git_sha
    ),
    group = NA_character_,
    notes = NA_character_
  )

  pkg_tbl <- dplyr::bind_rows(
    sieler2026_session_info_parse_pkg_strings(si$basePkgs, "base"),
    sieler2026_session_info_parse_pkg_strings(si$otherPkgs, "other_attached"),
    sieler2026_session_info_parse_pkg_strings(si$loadedOnly, "loaded_only")
  )

  driver_tbl <- sieler2026_session_info_analysis_driver_registry(proj_path)

  dplyr::bind_rows(env_tbl, pkg_tbl, driver_tbl)
}

#' Write Software_SessionInfo CSV to build dir and supplementary root.
#'
#' @param proj_path Absolute project root.
#' @param mode Build mode label (e.g. `submission`).
#' @param stamp Date stamp `YYYY-MM-DD`.
#' @return Absolute path to CSV in the build directory (for zip inclusion).
sieler2026_export_session_info_csv <- function(
    proj_path,
    mode,
    stamp = format(Sys.Date(), "%Y-%m-%d")
) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Install readr: install.packages('readr')")
  }
  paths <- sieler2026_supp_paths(proj_path, mode, stamp)
  base_name <- paste0("Software_SessionInfo__", mode, "__", stamp, ".csv")
  out_build <- file.path(paths$build, base_name)
  out_root <- file.path(paths$root, base_name)

  tbl <- sieler2026_session_info_to_tibble(proj_path)
  readr::write_csv(tbl, out_build)
  file.copy(out_build, out_root, overwrite = TRUE)

  message("Wrote software session info: ", out_build)
  invisible(normalizePath(out_build, winslash = "/", mustWork = TRUE))
}

#' Basename of session-info CSV for a build (INDEX meta and zip manifests).
sieler2026_session_info_csv_basename <- function(mode, stamp = format(Sys.Date(), "%Y-%m-%d")) {
  paste0("Software_SessionInfo__", mode, "__", stamp, ".csv")
}

#' Supplementary folder layout: root (docs + YAML), dated build dir, archive.
#'
#' @param proj_path Absolute project root.
#' @param mode Build mode label (e.g. `submission`).
#' @param stamp Date stamp `YYYY-MM-DD` (default: today).
#' @return Named list with `root`, `build`, `archive` absolute paths.
sieler2026_supp_paths <- function(proj_path, mode, stamp = format(Sys.Date(), "%Y-%m-%d")) {
  root <- file.path(proj_path, "Manuscript", "Supplementary")
  build <- file.path(root, "_build", paste0(mode, "__", stamp))
  archive <- file.path(root, "_build", "_archive")
  dir.create(build, recursive = TRUE, showWarnings = FALSE)
  dir.create(archive, recursive = TRUE, showWarnings = FALSE)
  list(root = root, build = build, archive = archive)
}

#' Write a zip archive with flat basenames (no directory paths inside the zip).
#'
#' @param files Character vector of absolute paths to existing files.
#' @param out_zip Absolute path for the output `.zip` file.
sieler2026_supp_write_zip <- function(files, out_zip) {
  files <- unique(normalizePath(files, winslash = "/", mustWork = FALSE))
  files <- files[file.exists(files)]
  if (length(files) < 1L) {
    stop("No files to zip for: ", out_zip)
  }
  out_zip <- normalizePath(out_zip, winslash = "/", mustWork = FALSE)
  dir.create(dirname(out_zip), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(out_zip)) {
    unlink(out_zip)
  }
  tmp_dir <- tempfile(pattern = "supp_zip_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  for (f in files) {
    file.copy(f, file.path(tmp_dir, basename(f)), overwrite = TRUE)
  }
  st <- utils::zip(
    zipfile = out_zip,
    files = list.files(tmp_dir, full.names = TRUE),
    flags = "-j -q"
  )
  if (!is.numeric(st) || length(st) != 1L || st != 0L) {
    stop("utils::zip failed (exit ", st, ") for: ", out_zip)
  }
  invisible(out_zip)
}

#' Move stray per-build xlsx/pdf/csv from supplementary root into `_build/_archive/`.
#'
#' @param paths Named list from [sieler2026_supp_paths()].
#' @return Invisibly, character vector of relocated basenames.
sieler2026_supp_relocate_stray_build_artifacts <- function(paths) {
  root <- paths$root
  archive <- paths$archive
  dir.create(archive, recursive = TRUE, showWarnings = FALSE)
  pat <- "^(Supplementary_Tables__|Supplementary_Figures__)"
  ext_ok <- "\\.(xlsx|pdf|csv)$"
  loose <- list.files(root, pattern = pat, full.names = TRUE)
  loose <- loose[grepl(ext_ok, loose, ignore.case = TRUE)]
  if (length(loose) < 1L) {
    return(invisible(character(0)))
  }
  moved <- character(0)
  for (f in loose) {
    dest <- file.path(archive, basename(f))
    if (file.exists(dest)) {
      dest <- file.path(archive, paste0(basename(f), "__", format(Sys.time(), "%H%M%S")))
    }
    ok <- file.rename(f, dest)
    if (isTRUE(ok)) {
      moved <- c(moved, basename(f))
    }
  }
  if (length(moved) > 0L) {
    message(
      "Relocated ",
      length(moved),
      " stray supplement artifact(s) from root -> ",
      archive
    )
  }
  invisible(moved)
}

#' Copy INDEX xlsx to root and create submission zip bundles for tables and figures.
#'
#' @param proj_path Absolute project root.
#' @param mode Build mode label.
#' @param stamp Date stamp `YYYY-MM-DD`.
#' @param combined If `TRUE`, also zip legacy `__ALL__` table/figure outputs when present.
#' @return Named list of output paths (`tables_zip`, `figures_zip`, `index_xlsx`, ...).
sieler2026_supp_finalize_packages <- function(
    proj_path,
    mode,
    stamp = format(Sys.Date(), "%Y-%m-%d"),
    combined = FALSE
) {
  paths <- sieler2026_supp_paths(proj_path, mode, stamp)
  build <- paths$build
  root <- paths$root

  index_xlsx_build <- file.path(
    build,
    paste0("Supplementary_Tables__INDEX__", mode, "__", stamp, ".xlsx")
  )
  index_manifest <- file.path(
    build,
    paste0("Supplementary_Tables__INDEX__", mode, "__", stamp, "__manifest.csv")
  )
  if (!file.exists(index_xlsx_build)) {
    stop("Missing INDEX workbook in build dir: ", index_xlsx_build)
  }

  session_csv <- sieler2026_export_session_info_csv(proj_path, mode, stamp)

  index_xlsx_root <- file.path(
    root,
    paste0("Supplementary_Tables__INDEX__", mode, "__", stamp, ".xlsx")
  )
  file.copy(index_xlsx_build, index_xlsx_root, overwrite = TRUE)

  mod_tbl_xlsx <- list.files(
    build,
    pattern = paste0(
      "^Supplementary_Tables__\\d{2}__.+__",
      gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", mode, perl = TRUE),
      "__",
      gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", stamp, perl = TRUE),
      "\\.xlsx$"
    ),
    full.names = TRUE
  )
  tables_zip_files <- c(mod_tbl_xlsx, index_xlsx_build, session_csv)
  if (file.exists(index_manifest)) {
    tables_zip_files <- c(tables_zip_files, index_manifest)
  }
  tables_zip_files <- unique(tables_zip_files[file.exists(tables_zip_files)])

  tables_zip <- file.path(
    root,
    paste0("Sieler2026_SupplementaryTables__", mode, "__", stamp, ".zip")
  )
  sieler2026_supp_write_zip(tables_zip_files, tables_zip)
  message("Wrote tables zip: ", tables_zip)

  fig_pdfs <- list.files(
    build,
    pattern = paste0(
      "^Supplementary_Figures__\\d{2}__.+__",
      gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", mode, perl = TRUE),
      "__",
      gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", stamp, perl = TRUE),
      "\\.pdf$"
    ),
    full.names = TRUE
  )
  fig_manifests <- sub("\\.pdf$", "__manifest.csv", fig_pdfs, ignore.case = TRUE)
  fig_zip_files <- unique(c(fig_pdfs, fig_manifests[file.exists(fig_manifests)]))

  figures_zip <- NULL
  if (length(fig_zip_files) >= 1L) {
    figures_zip <- file.path(
      root,
      paste0("Sieler2026_SupplementaryFigures__", mode, "__", stamp, ".zip")
    )
    sieler2026_supp_write_zip(fig_zip_files, figures_zip)
    message("Wrote figures zip: ", figures_zip)
  } else {
    message("No per-module figure PDFs found; skipped figures zip.")
  }

  tables_all_zip <- NULL
  figures_all_zip <- NULL
  if (isTRUE(combined)) {
    all_xlsx <- file.path(
      build,
      paste0("Supplementary_Tables__ALL__", mode, "__", stamp, ".xlsx")
    )
    all_xlsx_manifest <- file.path(
      build,
      paste0("Supplementary_Tables__ALL__", mode, "__", stamp, "__manifest.csv")
    )
    all_zip_files <- c(all_xlsx, all_xlsx_manifest, session_csv)
    all_zip_files <- all_zip_files[file.exists(all_zip_files)]
    if (length(all_zip_files) >= 1L) {
      tables_all_zip <- file.path(
        root,
        paste0("Sieler2026_SupplementaryTables__ALL__", mode, "__", stamp, ".zip")
      )
      sieler2026_supp_write_zip(all_zip_files, tables_all_zip)
      message("Wrote legacy ALL tables zip: ", tables_all_zip)
    }

    all_pdf <- file.path(
      build,
      paste0("Supplementary_Figures__ALL__", mode, "__", stamp, ".pdf")
    )
    all_pdf_manifest <- file.path(
      build,
      paste0("Supplementary_Figures__ALL__", mode, "__", stamp, "__manifest.csv")
    )
    all_fig_files <- c(all_pdf, all_pdf_manifest)
    all_fig_files <- all_fig_files[file.exists(all_fig_files)]
    if (length(all_fig_files) >= 1L) {
      figures_all_zip <- file.path(
        root,
        paste0("Sieler2026_SupplementaryFigures__ALL__", mode, "__", stamp, ".zip")
      )
      sieler2026_supp_write_zip(all_fig_files, figures_all_zip)
      message("Wrote legacy ALL figures zip: ", figures_all_zip)
    }
  }

  invisible(list(
    tables_zip = tables_zip,
    figures_zip = figures_zip,
    tables_all_zip = tables_all_zip,
    figures_all_zip = figures_all_zip,
    index_xlsx = index_xlsx_root
  ))
}

#' Table registry + `group_index` map (same rules as the supplement Excel builder).
#'
#' @param require_csv If `TRUE`, error when no `Tables/*.csv`. If `FALSE`, return empty
#'   registry when there are no CSVs (figure-only supplement assembly).
#'
#' @return Named list: `registry` (tibble through `supp_id` / `citation_table`),
#'   `grp_map` (tibble `group`, `group_index`), `module_index` (int).
sieler2026_supp_build_table_registry <- function(
    proj_path,
    module_dir,
    path_manuscript_supp,
    require_csv = TRUE
) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Install yaml: install.packages('yaml')")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Install tibble: install.packages('tibble')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install dplyr: install.packages('dplyr')")
  }

  mod_idx <- sieler2026_supp_module_index(module_dir)
  if (is.na(mod_idx)) {
    stop("Could not parse leading module number from: ", module_dir)
  }

  path_tables <- file.path(proj_path, "Results", module_dir, "Tables")
  csv_files <- if (dir.exists(path_tables)) {
    list.files(path_tables, pattern = "\\.csv$", full.names = TRUE)
  } else {
    character(0)
  }

  if (length(csv_files) == 0L) {
    if (isTRUE(require_csv)) {
      stop("No CSV files found in: ", path_tables)
    }
    return(list(
      registry = tibble::tibble(),
      grp_map = tibble::tibble(
        group = character(0),
        group_index = integer(0)
      ),
      module_index = mod_idx,
      path_tables = path_tables
    ))
  }

  overrides <- sieler2026_supp_read_table_overrides(path_manuscript_supp, module_dir)

  csv_basename <- basename(csv_files)
  base_tbl <- tibble::tibble(
    csv_path = csv_files,
    csv_basename = csv_basename,
    csv_stem = tools::file_path_sans_ext(csv_basename),
    group = vapply(csv_basename, function(bn) {
      sieler2026_supp_infer_group(tools::file_path_sans_ext(bn))
    }, character(1L)),
    module_index = mod_idx
  )

  # DiffAbund: include MaAsLin3 significant TSV outputs from Stats/ (not written under Tables/).
  # These are useful for full supplementary disclosure even if not used in the manuscript text.
  path.stats <- file.path(proj_path, "Results", module_dir, "Stats")
  maaslin_tsv <- if (dir.exists(path.stats)) {
    all_sig <- list.files(
      path.stats,
      pattern = "significant_results\\.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )
    all_sig[grepl("/maaslin_[^/]+/significant_results\\.tsv$", all_sig)]
  } else {
    character(0)
  }
  if (length(maaslin_tsv) > 0L) {
    maaslin_tsv <- sort(maaslin_tsv)
    maaslin_bn <- vapply(maaslin_tsv, function(p) {
      paste0(basename(dirname(p)), "__", basename(p))
    }, character(1L))
    base_tbl <- dplyr::bind_rows(
      base_tbl,
      tibble::tibble(
        csv_path = maaslin_tsv,
        csv_basename = maaslin_bn,
        csv_stem = tools::file_path_sans_ext(maaslin_bn),
        group = rep("zzz_maaslin3_significant_results", length(maaslin_tsv)),
        order_within_group = seq_along(maaslin_tsv),
        module_index = mod_idx
      )
    )
  }

  if (!is.null(overrides$tables) && is.list(overrides$tables)) {
    ov_list <- overrides$tables
    ov_tbl <- dplyr::bind_rows(lapply(ov_list, tibble::as_tibble))
    if (!"csv_basename" %in% colnames(ov_tbl)) {
      stop("supplement_overrides YAML: each `tables:` entry must include `csv_basename`.")
    }
    ov_cols <- intersect(
      colnames(ov_tbl),
      c("csv_basename", "group", "order_within_group", "supp_id_override", "exclude")
    )
    ov_tbl <- ov_tbl[, ov_cols, drop = FALSE]
    base_tbl <- dplyr::left_join(base_tbl, ov_tbl, by = "csv_basename", suffix = c("", "_ov"))
    if ("exclude" %in% colnames(base_tbl)) {
      base_tbl <- dplyr::filter(base_tbl, !isTRUE(.data$exclude))
      base_tbl <- dplyr::select(base_tbl, -dplyr::any_of("exclude"))
    }
    if ("group_ov" %in% colnames(base_tbl)) {
      base_tbl <- dplyr::mutate(
        base_tbl,
        group = dplyr::coalesce(.data$group_ov, .data$group)
      )
      base_tbl <- dplyr::select(base_tbl, -dplyr::any_of("group_ov"))
    }
    if ("order_within_group_ov" %in% colnames(base_tbl)) {
      base_tbl <- dplyr::mutate(
        base_tbl,
        order_within_group = dplyr::coalesce(.data$order_within_group_ov, .data$order_within_group)
      )
      base_tbl <- dplyr::select(base_tbl, -dplyr::any_of("order_within_group_ov"))
    }
  }

  grp_levels <- sort(unique(base_tbl$group))
  grp_map <- tibble::tibble(group = grp_levels, group_index = seq_along(grp_levels))

  base_tbl <- base_tbl |>
    dplyr::left_join(grp_map, by = "group")

  if (!"order_within_group" %in% colnames(base_tbl)) {
    base_tbl$order_within_group <- rep(NA_integer_, nrow(base_tbl))
  }

  base_tbl <- base_tbl |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(
      `_sort_key` = dplyr::if_else(
        is.na(.data$order_within_group),
        Inf,
        as.double(.data$order_within_group)
      )
    ) |>
    dplyr::arrange(.data$`_sort_key`, .data$csv_basename, .by_group = TRUE) |>
    dplyr::mutate(within_group_idx = dplyr::row_number()) |>
    dplyr::select(-dplyr::all_of("_sort_key")) |>
    dplyr::ungroup()

  base_tbl <- base_tbl |>
    dplyr::mutate(
      supp_id = paste(.data$module_index, .data$group_index, .data$within_group_idx, sep = "."),
      citation_table = paste0("Table S", .data$supp_id)
    )

  if ("supp_id_override" %in% colnames(base_tbl)) {
    base_tbl <- dplyr::mutate(
      base_tbl,
      supp_id = dplyr::coalesce(as.character(.data$supp_id_override), .data$supp_id),
      citation_table = paste0("Table S", .data$supp_id)
    )
  }

  base_tbl$module_dir <- module_dir

  list(
    registry = base_tbl,
    grp_map = grp_map,
    module_index = mod_idx,
    path_tables = path_tables
  )
}

#' Map a figure stem to a table `group` label (must align with `registry$group`).
#'
#' Uses optional YAML `group` when present; otherwise longest table group that
#' starts with the figure's `__`-based inferred group, then the same using
#' progressively shorter single-`_` prefixes of `fig_stem` (figures often use
#' one underscore where tables use `__`), then longest table group that the
#' figure stem starts with; otherwise the figure's inferred group (figure-only).
sieler2026_supp_match_figure_to_table_group <- function(fig_stem, table_groups, yaml_group) {
  yg <- yaml_group
  if (!is.null(yg) && !is.na(yg) && nzchar(as.character(yg))) {
    return(as.character(yg))
  }
  if (length(table_groups) == 0L) {
    return(sieler2026_supp_infer_group(fig_stem))
  }
  pick_longest_prefix_hit <- function(prefix) {
    hit <- table_groups[startsWith(table_groups, prefix)]
    if (length(hit) > 0L) {
      hit[which.max(nchar(hit))]
    } else {
      NA_character_
    }
  }
  fig_inferred <- sieler2026_supp_infer_group(fig_stem)
  hit_g <- pick_longest_prefix_hit(fig_inferred)
  if (!is.na(hit_g)) {
    return(hit_g)
  }
  parts_us <- strsplit(fig_stem, "_", fixed = TRUE)[[1]]
  if (length(parts_us) >= 2L) {
    for (k in seq(from = length(parts_us) - 1L, to = 1L, by = -1L)) {
      cand <- paste(parts_us[seq_len(k)], collapse = "_")
      hit_us <- pick_longest_prefix_hit(cand)
      if (!is.na(hit_us)) {
        return(hit_us)
      }
    }
  }
  hit2 <- table_groups[vapply(table_groups, function(g) {
    startsWith(fig_stem, g)
  }, logical(1L))]
  if (length(hit2) > 0L) {
    return(hit2[which.max(nchar(hit2))])
  }
  fig_inferred
}

#' Append figure-only groups to the table `grp_map` without renumbering table groups.
sieler2026_supp_merge_grp_map_with_figure_groups <- function(grp_map_tables, figure_groups) {
  if (length(figure_groups) == 0L) {
    return(grp_map_tables)
  }
  extra <- unique(figure_groups[!figure_groups %in% grp_map_tables$group])
  if (length(extra) == 0L) {
    return(grp_map_tables)
  }
  extra <- sort(extra)
  n0 <- nrow(grp_map_tables)
  dplyr::bind_rows(
    grp_map_tables,
    tibble::tibble(
      group = extra,
      group_index = seq_along(extra) + n0
    )
  )
}

# ----- Supplementary figure PDF cover (TOC) layout --------------------------------

#' Default layout for the text-only cover PDF (narrow margins, wrapped body lines).
#'
#' @return Named list used by `sieler2026_supp_cover_*` helpers.
sieler2026_supp_cover_layout_pars <- function() {
  list(
    mai = rep(0.12, 4L),
    x_left = 0.02,
    x_right = 0.995,
    title_y = 0.985,
    title_cex = 1.12,
    body_cex = 0.82,
    y_first = 0.905,
    y_next = 0.935,
    y_bottom = 0.05,
    step = 0.027
  )
}

#' Width of `txt` in user coordinates (active graphics device).
sieler2026_supp_cover_strwidth_u <- function(txt, cex) {
  suppressWarnings(graphics::strwidth(as.character(txt), cex = cex, units = "user"))
}

#' Break a single token (no spaces) into pieces each fitting `width_u` (measured).
sieler2026_supp_cover_break_long_token <- function(tok, width_u, cex) {
  tok <- as.character(tok)[[1L]]
  if (!nzchar(tok)) {
    return(character(0))
  }
  if (sieler2026_supp_cover_strwidth_u(tok, cex) <= width_u) {
    return(tok)
  }
  out <- character()
  n <- nchar(tok)
  i0 <- 1L
  while (i0 <= n) {
    lo <- i0
    hi <- n
    best <- i0
    while (lo <= hi) {
      mid <- lo + (hi - lo) %/% 2L
      ss <- substring(tok, i0, mid)
      if (sieler2026_supp_cover_strwidth_u(ss, cex) <= width_u) {
        best <- mid
        lo <- mid + 1L
      } else {
        hi <- mid - 1L
      }
    }
    if (best < i0) {
      best <- i0
    }
    out <- c(out, substring(tok, i0, best))
    i0 <- best + 1L
  }
  out
}

#' Split one logical TOC line into physical lines using measured widths (uses full `width_u`).
sieler2026_supp_cover_wrap_one <- function(s, width_u, cex) {
  sw <- function(txt) sieler2026_supp_cover_strwidth_u(txt, cex)
  s <- as.character(s)
  if (length(s) != 1L) {
    return(character(0))
  }
  s <- trimws(s[[1L]])
  if (!nzchar(s)) {
    return(character(0))
  }
  words <- strsplit(s, "\\s+", perl = TRUE)[[1L]]
  words <- words[nzchar(words)]
  if (length(words) == 0L) {
    return(character(0))
  }
  events <- list()
  for (w in words) {
    if (sw(w) <= width_u) {
      events[[length(events) + 1L]] <- list(t = w, glue = FALSE)
    } else {
      cs <- sieler2026_supp_cover_break_long_token(w, width_u = width_u, cex = cex)
      for (i in seq_along(cs)) {
        events[[length(events) + 1L]] <- list(t = cs[[i]], glue = i > 1L)
      }
    }
  }
  lines_out <- character()
  cur <- ""
  for (e in events) {
    trial <- if (!nzchar(cur)) {
      e$t
    } else if (isTRUE(e$glue)) {
      paste0(cur, e$t)
    } else {
      paste(cur, e$t, sep = " ")
    }
    if (sw(trial) <= width_u) {
      cur <- trial
    } else {
      if (nzchar(cur)) {
        lines_out <- c(lines_out, cur)
      }
      cur <- e$t
    }
  }
  if (nzchar(cur)) {
    lines_out <- c(lines_out, cur)
  }
  lines_out
}

#' Wrap every logical TOC line to one or more physical lines (requires active plot).
sieler2026_supp_cover_wrap_lines <- function(lines, width_u, cex) {
  out <- character()
  for (ln in lines) {
    out <- c(out, sieler2026_supp_cover_wrap_one(ln, width_u = width_u, cex = cex))
  }
  out
}

#' Number of physical lines after wrapping (opens a temporary PDF for `strwidth`).
sieler2026_supp_cover_n_physical_lines <- function(lines, layout = sieler2026_supp_cover_layout_pars()) {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf, width = 8.5, height = 11)
  tryCatch(
    {
      graphics::par(mai = layout$mai)
      graphics::plot.new()
      graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      wu <- layout$x_right - layout$x_left
      length(sieler2026_supp_cover_wrap_lines(lines, width_u = wu, cex = layout$body_cex))
    },
    finally = {
      grDevices::dev.off()
      if (file.exists(tf)) {
        unlink(tf)
      }
    }
  )
}

#' Pages required for the cover given how many physical (possibly wrapped) lines.
sieler2026_supp_cover_pages_needed <- function(
    n_lines,
    y_first,
    y_next,
    y_bottom,
    step
) {
  n_lines <- suppressWarnings(as.integer(n_lines))
  if (length(n_lines) != 1L || is.na(n_lines) || n_lines < 1L) {
    return(1L)
  }
  cap_first <- max(1L, floor((y_first - y_bottom) / step) + 1L)
  cap_next <- max(1L, floor((y_next - y_bottom) / step) + 1L)
  if (n_lines <= cap_first) {
    return(1L)
  }
  1L + as.integer(ceiling((n_lines - cap_first) / cap_next))
}

#' Iterate cover page count vs. wrapped TOC lines until `pdf_start_page` stabilizes.
#'
#' @param n_fig Integer number of figure rows (length of `page_counts_fig`).
#' @param page_counts_fig Integer vector of pages per figure PDF.
#' @param build_cover_lines Function accepting `pdf_start_page` (integer vector length
#'   `n_fig`) and returning a character vector of logical TOC lines (before wrapping).
#'
#' @return List: `n_cover_pages`, `pdf_start_page`, `cover_lines`.
sieler2026_supp_cover_resolve_pagination <- function(
    n_fig,
    page_counts_fig,
    build_cover_lines,
    layout = sieler2026_supp_cover_layout_pars(),
    max_iter = 15L
) {
  n_fig <- suppressWarnings(as.integer(n_fig))
  if (length(n_fig) != 1L || is.na(n_fig) || n_fig < 1L) {
    stop("n_fig must be a positive integer")
  }
  if (length(page_counts_fig) != n_fig) {
    stop("page_counts_fig length must equal n_fig")
  }
  pdf_start_page <- integer(n_fig)
  cover_lines <- character(0)
  n_cover <- 1L
  wu <- layout$x_right - layout$x_left
  for (iter in seq_len(max_iter)) {
    offset <- n_cover
    for (i in seq_len(n_fig)) {
      pdf_start_page[[i]] <- offset + 1L
      offset <- offset + as.integer(page_counts_fig[[i]])
    }
    cover_lines <- build_cover_lines(pdf_start_page)
    n_phys <- sieler2026_supp_cover_n_physical_lines(cover_lines, layout = layout)
    n_need <- sieler2026_supp_cover_pages_needed(
      n_phys,
      y_first = layout$y_first,
      y_next = layout$y_next,
      y_bottom = layout$y_bottom,
      step = layout$step
    )
    if (n_need == n_cover) {
      return(list(
        n_cover_pages = n_cover,
        pdf_start_page = pdf_start_page,
        cover_lines = cover_lines
      ))
    }
    n_cover <- n_need
  }
  stop(
    "Cover pagination did not converge after ",
    max_iter,
    " iterations (try shorter figure names or increase cover pages)."
  )
}

#' Write the supplementary-figures cover / TOC PDF (wrapped body text, tight margins).
sieler2026_supp_write_cover_pdf <- function(out_pdf, title, lines, layout = sieler2026_supp_cover_layout_pars()) {
  wu <- layout$x_right - layout$x_left
  grDevices::pdf(out_pdf, width = 8.5, height = 11)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mai = layout$mai)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  physical <- sieler2026_supp_cover_wrap_lines(lines, width_u = wu, cex = layout$body_cex)
  graphics::text(
    0.5,
    layout$title_y,
    title,
    adj = c(0.5, 1),
    cex = layout$title_cex
  )
  y <- layout$y_first
  for (ln in physical) {
    if (y < layout$y_bottom) {
      graphics::plot.new()
      graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      y <- layout$y_next
    }
    graphics::text(
      layout$x_left,
      y,
      ln,
      adj = c(0, 1),
      cex = layout$body_cex
    )
    y <- y - layout$step
  }
  invisible(out_pdf)
}
