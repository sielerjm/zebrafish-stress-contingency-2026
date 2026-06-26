# 99__BuildSupplementCombined.R
#
# Created by: Michael Sieler
# Date last updated: 2026-06-25
#
# Description:
#   Build **one** supplementary Excel workbook (all `Results/NN__*/Tables/*.csv`) and
#   **one** combined figure PDF (all `supplement_map__NN__*.yml` in module order).
#   Table `supp_id` values match per-module numbering (`module_index` from folder prefix).
#
# Expected input:
#   - CLI: `<mode>` (e.g. `submission`)
#
# Expected output:
#   - `Manuscript/Supplementary/_build/<mode>__<date>/Supplementary_Tables__ALL__<mode>__YYYY-MM-DD.xlsx` + manifest CSV
#   - `Manuscript/Supplementary/_build/<mode>__<date>/Supplementary_Figures__ALL__<mode>__YYYY-MM-DD.pdf` + manifest CSV (if any figure maps resolve)

options(stringsAsFactors = FALSE)

req <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install ", pkg, ": install.packages('", pkg, "')")
  }
}
req("here")
req("tibble")
req("dplyr")
req("readr")
req("yaml")
req("openxlsx")
req("pdftools")

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
  here::i_am("Code/01__Analysis/99__BuildSupplementCombined.R")
  proj.path <- as.character(here::here())
}

old_wd <- setwd(proj.path)
on.exit(setwd(old_wd), add = TRUE)

shared <- file.path(proj.path, "Code", "01__Analysis", "98__SupplementShared.R")
source(shared, local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
flag_tables_only <- "--tables-only" %in% args
flag_figures_only <- "--figures-only" %in% args
if (isTRUE(flag_tables_only) && isTRUE(flag_figures_only)) {
  stop("Use at most one of --tables-only and --figures-only.", call. = FALSE)
}
args_mode <- args[!grepl("^--", args)]
if (length(args_mode) < 1L) {
  stop(
    "Usage: Rscript Code/01__Analysis/99__BuildSupplementCombined.R <mode> [--tables-only] [--figures-only]\n",
    "Example: Rscript Code/01__Analysis/99__BuildSupplementCombined.R submission"
  )
}
mode <- args_mode[[1]]

run_tables <- !isTRUE(flag_figures_only)
run_figures <- !isTRUE(flag_tables_only)

path.results <- file.path(proj.path, "Results")
path.supp <- file.path(proj.path, "Manuscript", "Supplementary")
dir.create(path.supp, recursive = TRUE, showWarnings = FALSE)

mods <- sieler2026_supp_list_result_modules(proj.path)

git_sha <- sieler2026_supp_git_sha(proj.path)
pkg_v <- sieler2026_supp_pkg_version

stamp <- format(Sys.Date(), "%Y-%m-%d")
paths <- sieler2026_supp_paths(proj.path, mode, stamp)
path.build <- paths$build

# ----- Combined tables Excel -------------------------------------------------

if (isTRUE(run_tables)) {
regs <- list()
modules_with_tables <- character(0)
for (m in mods) {
  pt <- file.path(path.results, m, "Tables")
  n_csv <- if (dir.exists(pt)) {
    length(list.files(pt, pattern = "\\.csv$", ignore.case = TRUE))
  } else {
    0L
  }
  if (n_csv < 1L) {
    next
  }
  tr <- sieler2026_supp_build_table_registry(
    proj.path,
    m,
    path.supp,
    require_csv = TRUE
  )
  regs[[m]] <- tr$registry
  modules_with_tables <- c(modules_with_tables, m)
}

if (length(regs) == 0L) {
  stop("No module with Tables/*.csv found; nothing to build for combined Excel.")
}

base_tbl <- dplyr::bind_rows(regs) |>
  dplyr::arrange(.data$module_index, .data$group_index, .data$within_group_idx)

meta_tbl <- tibble::tibble(
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
    "tables_scope"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    git_sha,
    R.version.string,
    pkg_v("openxlsx"),
    pkg_v("readr"),
    pkg_v("yaml"),
    "ALL",
    mode,
    paste(modules_with_tables, collapse = ", "),
    "Results/*/Tables/*.csv + (03__DiffAbund) Stats/maaslin_*/significant_results.tsv"
  )
)

out_xlsx <- file.path(
  path.build,
  paste0("Supplementary_Tables__ALL__", mode, "__", stamp, ".xlsx")
)
out_manifest_tables <- file.path(
  path.build,
  paste0("Supplementary_Tables__ALL__", mode, "__", stamp, "__manifest.csv")
)

sieler2026_supp_write_excel_workbook(
  base_tbl,
  meta_tbl,
  out_xlsx,
  out_manifest_tables,
  proj.path
)

message("Wrote combined supplement Excel: ", out_xlsx)
message("Wrote tables manifest CSV:       ", out_manifest_tables)
}

# ----- Combined figures PDF ---------------------------------------------------

if (isTRUE(run_figures)) {

png_to_onepage_pdf <- function(png_path, out_pdf, dpi = 150) {
  if (requireNamespace("magick", quietly = TRUE)) {
    img <- magick::image_read(png_path)
    magick::image_write(img, path = out_pdf, format = "pdf")
    return(invisible(out_pdf))
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("For PNG inputs install `magick` or `png`: install.packages(c('magick','png'))")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Missing base recommended package `grid`.")
  }
  im <- png::readPNG(png_path)
  if (length(dim(im)) == 3L && dim(im)[[3L]] %in% c(3L, 4L)) {
    ht <- dim(im)[[1L]]
    wd <- dim(im)[[2L]]
  } else {
    stop("Unexpected PNG array shape for: ", png_path)
  }
  grDevices::pdf(out_pdf, width = wd / dpi, height = ht / dpi)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.raster(im, width = 1, height = 1)
  invisible(out_pdf)
}

supp_pdf_page_count <- function(path) {
  tryCatch(
    {
      np <- pdftools::pdf_info(path)$pages
      np <- suppressWarnings(as.integer(np))
      if (length(np) != 1L || is.na(np) || np < 1L) {
        1L
      } else {
        np
      }
    },
    error = function(e) 1L
  )
}

fig_blocks <- list()

for (m in mods) {
  map_path <- file.path(path.supp, paste0("supplement_map__", m, ".yml"))
  if (!file.exists(map_path)) {
    next
  }
  map_yaml <- yaml::read_yaml(map_path)
  if (is.null(map_yaml$figures) || length(map_yaml$figures) == 0L) {
    next
  }

  fig_dir <- file.path(path.results, m, "Figures")
  if (!is.null(map_yaml$figures_dir) && nzchar(map_yaml$figures_dir)) {
    p <- map_yaml$figures_dir
    if (!grepl("^/", p) && !grepl("^[A-Za-z]:", p)) {
      fig_dir <- file.path(proj.path, p)
    } else {
      fig_dir <- p
    }
  }
  if (!dir.exists(fig_dir)) {
    warning("Skipping figures for ", m, ": not a directory: ", fig_dir, call. = FALSE)
    next
  }

  fig_tbl <- sieler2026_supp_normalize_map_figures(map_yaml$figures)
  if (nrow(fig_tbl) == 0L) {
    next
  }
  fig_tbl$yaml_index <- seq_len(nrow(fig_tbl))
  fig_tbl$module_dir <- m
  fig_tbl$module_order <- match(m, mods)

  fig_tbl$resolved_path <- vapply(fig_tbl$basename, function(bn) {
    res <- sieler2026_supp_resolve_figure_path_pdf_first(fig_dir, bn)
    res$path
  }, character(1L))
  fig_tbl$resolved_ext <- vapply(fig_tbl$basename, function(bn) {
    res <- sieler2026_supp_resolve_figure_path_pdf_first(fig_dir, bn)
    res$ext
  }, character(1L))

  fig_ok <- dplyr::filter(fig_tbl, !is.na(.data$resolved_path))
  if (nrow(fig_ok) == 0L) {
    next
  }

  tr <- sieler2026_supp_build_table_registry(
    proj.path,
    m,
    path.supp,
    require_csv = FALSE
  )
  table_groups <- if (nrow(tr$registry) > 0L) {
    unique(tr$registry$group)
  } else {
    character(0)
  }

  fig_ok$matched_group <- vapply(seq_len(nrow(fig_ok)), function(i) {
    sieler2026_supp_match_figure_to_table_group(
      fig_ok$basename[[i]],
      table_groups,
      fig_ok$yaml_group[[i]]
    )
  }, character(1L))

  grp_full <- sieler2026_supp_merge_grp_map_with_figure_groups(tr$grp_map, fig_ok$matched_group)
  fig_ok <- fig_ok |>
    dplyr::left_join(
      dplyr::rename(grp_full, matched_group = group),
      by = "matched_group"
    ) |>
    dplyr::rename(group = matched_group)

  if (any(is.na(fig_ok$group_index))) {
    bad <- unique(fig_ok$group[is.na(fig_ok$group_index)])
    stop(
      "Missing group_index for module ", m, " figure group(s): ",
      paste(bad, collapse = ", ")
    )
  }

  if (!"order_within_group" %in% colnames(fig_ok)) {
    fig_ok$order_within_group <- rep(NA_integer_, nrow(fig_ok))
  }

  mod_idx <- tr$module_index
  if (is.na(mod_idx)) {
    stop("Could not parse module index for: ", m)
  }

  fig_ok <- fig_ok |>
    dplyr::mutate(module_index = mod_idx) |>
    dplyr::group_by(.data$module_dir, .data$group) |>
    dplyr::mutate(
      `_sort_key` = dplyr::if_else(
        is.na(.data$order_within_group),
        Inf,
        as.double(.data$order_within_group)
      )
    ) |>
    dplyr::arrange(.data$`_sort_key`, .data$yaml_index, .by_group = TRUE) |>
    dplyr::mutate(within_group_idx = dplyr::row_number()) |>
    dplyr::select(-dplyr::all_of("_sort_key")) |>
    dplyr::ungroup()

  fig_ok <- fig_ok |>
    dplyr::mutate(
      supp_id = paste(.data$module_index, .data$group_index, .data$within_group_idx, sep = "."),
      citation_figure = paste0("Figure S", .data$supp_id)
    )

  if (any(!is.na(fig_ok$supp_id_override))) {
    fig_ok <- dplyr::mutate(
      fig_ok,
      supp_id = dplyr::coalesce(as.character(.data$supp_id_override), .data$supp_id),
      citation_figure = paste0("Figure S", .data$supp_id)
    )
  }

  fig_ok$figure_path <- vapply(
    fig_ok$resolved_path,
    function(p) sieler2026_supp_path_for_workbook(p, proj.path),
    character(1L)
  )

  fig_blocks[[length(fig_blocks) + 1L]] <- fig_ok
}

if (length(fig_blocks) == 0L) {
  message("No supplementary figure maps resolved; skipping combined figure PDF.")
} else {
  fig_all <- dplyr::bind_rows(fig_blocks) |>
    dplyr::arrange(.data$module_order, .data$yaml_index)

  tmp_pdfs <- character(0)
  on.exit(
    {
      for (f in tmp_pdfs) {
        if (file.exists(f)) {
          unlink(f)
        }
      }
    },
    add = TRUE
  )

  parts_figures <- character(nrow(fig_all))
  for (i in seq_len(nrow(fig_all))) {
    res_path <- fig_all$resolved_path[[i]]
    res_ext <- fig_all$resolved_ext[[i]]
    mdir <- fig_all$module_dir[[i]]
    if (identical(res_ext, "pdf")) {
      parts_figures[[i]] <- res_path
      next
    }
    tmp_pdf <- tempfile(pattern = paste0("supp_", mdir, "_", i, "_"), fileext = ".pdf")
    tmp_pdfs <- c(tmp_pdfs, tmp_pdf)
    png_to_onepage_pdf(res_path, tmp_pdf)
    parts_figures[[i]] <- tmp_pdf
  }

  page_counts_fig <- vapply(parts_figures, supp_pdf_page_count, integer(1L))
  use_cover <- nrow(fig_all) > 0L
  cover_layout <- sieler2026_supp_cover_layout_pars()
  cover_res <- if (isTRUE(use_cover)) {
    sieler2026_supp_cover_resolve_pagination(
      n_fig = nrow(fig_all),
      page_counts_fig = page_counts_fig,
      build_cover_lines = function(pdf_start_page) {
        hdr <- c(
          paste0("Build (UTC): ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
          paste0("Git: ", git_sha),
          paste0("Bundle: ALL (", length(fig_blocks), " module map(s))"),
          paste0("Mode: ", mode),
          "Figures (combined PDF order; p. = first page of that figure in this PDF):"
        )
        cl2 <- paste0(
          "p. ",
          pdf_start_page,
          ": [",
          fig_all$module_dir,
          "] ",
          fig_all$citation_figure,
          " - ",
          fig_all$basename
        )
        c(hdr, cl2)
      },
      layout = cover_layout
    )
  } else {
    NULL
  }
  n_cover_pages <- if (isTRUE(use_cover)) {
    cover_res$n_cover_pages
  } else {
    0L
  }
  pdf_start_page <- if (isTRUE(use_cover)) {
    cover_res$pdf_start_page
  } else {
    integer(nrow(fig_all))
  }

  manifest_fig <- fig_all |>
    dplyr::transmute(
      figure_path = .data$figure_path,
      basename = .data$basename,
      fig_stem = .data$basename,
      module = .data$module_dir,
      group = .data$group,
      module_index = .data$module_index,
      group_index = .data$group_index,
      order_within_group = .data$order_within_group,
      within_group_idx = .data$within_group_idx,
      supp_id = .data$supp_id,
      citation_figure = .data$citation_figure,
      pdf_start_page = pdf_start_page,
      pdf_n_pages = page_counts_fig,
      resolved_ext = .data$resolved_ext,
      caption = .data$caption
    )

  out_pdf <- file.path(
    path.build,
    paste0("Supplementary_Figures__ALL__", mode, "__", stamp, ".pdf")
  )
  out_manifest_fig <- file.path(
    path.build,
    paste0("Supplementary_Figures__ALL__", mode, "__", stamp, "__manifest.csv")
  )

  parts <- character(0)
  if (use_cover) {
    cover_pdf <- tempfile(pattern = "supp_cover_ALL_", fileext = ".pdf")
    tmp_pdfs <- c(tmp_pdfs, cover_pdf)
    sieler2026_supp_write_cover_pdf(
      cover_pdf,
      title = "Supplementary figures (combined, all modules)",
      lines = cover_res$cover_lines,
      layout = cover_layout
    )
    parts <- c(parts, cover_pdf)
  }

  parts <- c(parts, parts_figures)
  invisible(pdftools::pdf_combine(parts, output = out_pdf))
  readr::write_csv(manifest_fig, out_manifest_fig)
  message("Wrote combined supplement figures PDF: ", out_pdf)
  message("Wrote figures manifest CSV:          ", out_manifest_fig)
}
}
