# 99__BuildSupplementFiguresPdf.R
#
# Created by: Michael Sieler
# Date last updated: 2026-06-25
#
# Description:
#   Combine supplementary figures into one PDF using a YAML map under
#   `Manuscript/Supplementary/`. For each logical basename, uses `.pdf` if present
#   in `Results/<module>/Figures/`, otherwise `.png` (rasterized to a one-page PDF).
#   The figure manifest mirrors the table supplement: same `supp_id` scheme
#   (`module_index.group_index.within_group_idx`), workbook-style paths, and
#   `group` aligned to table-registry groups (optional per-figure `group:` in YAML).
#
# Expected input:
#   - CLI: `<module_dir>` `<mode>` (e.g. `01__Diversity test`)
#   - Map file: `Manuscript/Supplementary/supplement_map__<module_dir>.yml`
#
# Expected output:
#   - `Manuscript/Supplementary/_build/<mode>__<date>/Supplementary_Figures__<module>__<mode>__YYYY-MM-DD.pdf`
#   - `Manuscript/Supplementary/_build/<mode>__<date>/Supplementary_Figures__<module>__<mode>__YYYY-MM-DD__manifest.csv`

options(stringsAsFactors = FALSE)

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install here: install.packages('here')")
}
if (!requireNamespace("tibble", quietly = TRUE)) {
  stop("Install tibble: install.packages('tibble')")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Install dplyr: install.packages('dplyr')")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Install readr: install.packages('readr')")
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Install yaml: install.packages('yaml')")
}
if (!requireNamespace("pdftools", quietly = TRUE)) {
  stop("Install pdftools: install.packages('pdftools')")
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
  here::i_am("Code/01__Analysis/99__BuildSupplementFiguresPdf.R")
  proj.path <- as.character(here::here())
}

old_wd <- setwd(proj.path)
on.exit(setwd(old_wd), add = TRUE)

shared <- file.path(proj.path, "Code", "01__Analysis", "98__SupplementShared.R")
source(shared, local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript Code/01__Analysis/99__BuildSupplementFiguresPdf.R <module_dir> <mode>\n",
    "Example: Rscript Code/01__Analysis/99__BuildSupplementFiguresPdf.R 01__Diversity test"
  )
}

module_dir <- args[[1]]
mode <- args[[2]]

path.results <- file.path(proj.path, "Results")
path.manuscript_supp <- file.path(proj.path, "Manuscript", "Supplementary")
dir.create(path.manuscript_supp, recursive = TRUE, showWarnings = FALSE)

map_path <- file.path(path.manuscript_supp, paste0("supplement_map__", module_dir, ".yml"))
if (!file.exists(map_path)) {
  stop("Missing supplement map YAML: ", map_path)
}

map_yaml <- yaml::read_yaml(map_path)
if (is.null(map_yaml$figures) || length(map_yaml$figures) == 0L) {
  stop("YAML must define non-empty `figures:` (list of basenames or `basename:` records).")
}

fig_dir_default <- file.path(path.results, module_dir, "Figures")
fig_dir <- if (!is.null(map_yaml$figures_dir) && nzchar(map_yaml$figures_dir)) {
  p <- map_yaml$figures_dir
  if (!grepl("^/", p) && !grepl("^[A-Za-z]:", p)) {
    file.path(proj.path, p)
  } else {
    p
  }
} else {
  fig_dir_default
}

if (!dir.exists(fig_dir)) {
  stop("Figures directory not found: ", fig_dir)
}

fig_tbl <- sieler2026_supp_normalize_map_figures(map_yaml$figures)
fig_tbl$yaml_index <- seq_len(nrow(fig_tbl))
fig_tbl$module_dir <- module_dir

tr <- sieler2026_supp_build_table_registry(
  proj.path,
  module_dir,
  path.manuscript_supp,
  require_csv = FALSE
)
table_groups <- if (nrow(tr$registry) > 0L) {
  unique(tr$registry$group)
} else {
  character(0)
}

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
  stop("No figure files resolved; nothing to combine.")
}

fig_ok$matched_group <- vapply(seq_len(nrow(fig_ok)), function(i) {
  sieler2026_supp_match_figure_to_table_group(
    fig_ok$basename[[i]],
    table_groups,
    fig_ok$yaml_group[[i]]
  )
}, character(1L))

if (any(!is.na(fig_ok$yaml_group))) {
  unk <- unique(fig_ok$yaml_group[!is.na(fig_ok$yaml_group) & !fig_ok$yaml_group %in% table_groups])
  unk <- unk[nzchar(unk)]
  if (length(unk) > 0L) {
    for (ug in unk) {
      warning(
        "YAML figure `group: ", ug, "` is not a table `group` in this module; ",
        "it will receive a new group_index after table groups."
      )
    }
  }
}

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
    "Missing group_index for figure group(s): ",
    paste(bad, collapse = ", "),
    " (check `group:` in supplement_map YAML matches a table registry group)."
  )
}

if (!"order_within_group" %in% colnames(fig_ok)) {
  fig_ok$order_within_group <- rep(NA_integer_, nrow(fig_ok))
}

mod_idx <- tr$module_index
if (is.na(mod_idx)) {
  stop("Could not parse module index for: ", module_dir)
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

git_sha <- tryCatch(
  {
    sha <- system2("git", c("-C", proj.path, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)
    if (length(sha) == 1L && nzchar(sha[[1]])) sha[[1]] else NA_character_
  },
  warning = function(w) NA_character_,
  error = function(e) NA_character_
)

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

stamp <- format(Sys.Date(), "%Y-%m-%d")
paths <- sieler2026_supp_paths(proj.path, mode, stamp)
out_pdf <- file.path(
  paths$build,
  paste0("Supplementary_Figures__", module_dir, "__", mode, "__", stamp, ".pdf")
)
out_manifest <- file.path(
  paths$build,
  paste0("Supplementary_Figures__", module_dir, "__", mode, "__", stamp, "__manifest.csv")
)

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

parts_figures <- character(nrow(fig_ok))
for (i in seq_len(nrow(fig_ok))) {
  res_path <- fig_ok$resolved_path[[i]]
  res_ext <- fig_ok$resolved_ext[[i]]

  if (identical(res_ext, "pdf")) {
    parts_figures[[i]] <- res_path
    next
  }

  tmp_pdf <- tempfile(pattern = paste0("supp_r_", i, "_"), fileext = ".pdf")
  tmp_pdfs <- c(tmp_pdfs, tmp_pdf)
  png_to_onepage_pdf(res_path, tmp_pdf)
  parts_figures[[i]] <- tmp_pdf
}

page_counts_fig <- vapply(parts_figures, supp_pdf_page_count, integer(1L))
include_cover <- sieler2026_supp_yaml_truthy(map_yaml$include_cover, default = FALSE)
cover_layout <- sieler2026_supp_cover_layout_pars()
cover_res <- if (isTRUE(include_cover)) {
  sieler2026_supp_cover_resolve_pagination(
    n_fig = nrow(fig_ok),
    page_counts_fig = page_counts_fig,
    build_cover_lines = function(pdf_start_page) {
      hdr <- c(
        paste0("Build (UTC): ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
        paste0("Git: ", git_sha),
        paste0("Module: ", module_dir),
        paste0("Figures dir: ", sieler2026_supp_path_for_workbook(fig_dir, proj.path)),
        "Figures (combined PDF order; p. = first page of that figure in this PDF):"
      )
      cl2 <- paste0(
        "p. ",
        pdf_start_page,
        ": ",
        fig_ok$citation_figure,
        " - ",
        fig_ok$basename
      )
      c(hdr, cl2)
    },
    layout = cover_layout
  )
} else {
  NULL
}
n_cover_pages <- if (isTRUE(include_cover)) {
  cover_res$n_cover_pages
} else {
  0L
}
pdf_start_page <- if (isTRUE(include_cover)) {
  cover_res$pdf_start_page
} else {
  off <- 0L
  pp <- integer(nrow(fig_ok))
  for (i in seq_len(nrow(fig_ok))) {
    pp[[i]] <- off + 1L
    off <- off + as.integer(page_counts_fig[[i]])
  }
  pp
}

manifest_out <- fig_ok |>
  dplyr::transmute(
    figure_path = .data$figure_path,
    basename = .data$basename,
    fig_stem = .data$basename,
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

parts <- character(0)

if (isTRUE(include_cover)) {
  cover_pdf <- tempfile(pattern = "supp_cover_", fileext = ".pdf")
  tmp_pdfs <- c(tmp_pdfs, cover_pdf)
  sieler2026_supp_write_cover_pdf(
    cover_pdf,
    title = "Supplementary figures (combined)",
    lines = cover_res$cover_lines,
    layout = cover_layout
  )
  parts <- c(parts, cover_pdf)
}

parts <- c(parts, parts_figures)

invisible(pdftools::pdf_combine(parts, output = out_pdf))

readr::write_csv(manifest_out, out_manifest)

message("Wrote combined supplementary PDF: ", out_pdf)
message("Wrote manifest CSV: ", out_manifest)
