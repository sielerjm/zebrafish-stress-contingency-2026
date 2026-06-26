# 03__HelperFunctions.R
# Created by: Michael Sieler
# Date last updated: 2026-04-25
#
# Description: Curated helper functions for microbiome preprocessing and downstream
#   analysis, organized by topic. Prefer sourcing `01__Libraries.R` and `02__PlotSettings.R`
#   first so packages and ggplot defaults are on the search path.
#
# Expected input:  none (defines functions in the global environment when sourced).
# Expected output: Function definitions listed in the INDEX below.
#
# INDEX — legacy file mapping (Code/99__Archive/Functions/HelperFunctions/)
#   General / tables ........ miscFunctions.R — p_val_format, SigStars, cutColNames,
#                             cutCellNames, set_GT, print_tables, save_tables_as_files, save_env
#   Phyloseq wrangling ...... miscFunctions.R (archived); active copies in sec. 2 — samdatAsDataframe,
#                             ps_rename, psObjToDfLong
#   DADA2 / paths ........... (this file, sec. 2a) — resolve_latest_rds() for dated *.rds in Data/DADA2
#   DADA2 post-processing ... (this file, sec. 2b) — screening + metadata for the RoL
#                             zebrafish experiment: apply_pre_analysis_filters,
#                             rename_asvs_sequential, ensure_parasite_column,
#                             augment_filtered_phyloseq_metadata
#   Alpha diversity ........ alphaFunctions.R — GLM / Tukey / Levene helpers; active: norm_scores,
#                             ps_calc_diversity.phy, populate_ps_list_alpha_diversity (sec. 3)
#   Beta diversity ........... betaFunctions.R — capscale, adonis2, …; active: build_beta_dist_matrices_for_ps_list (sec. 3)
#   Differential abundance .. diffAbundFunctions.R — MaAsLin2, ANCOM-BC stubs, run_ancombc2
#   Plotting ................ plotFunctions.R — gen_box_plot, gen_PCoA_*, ggplot helpers
#   Publication theme ....... theme_sieler2026_publication, theme_sieler2026_publication_with_grid,
#                             theme_sieler2026_trend_panel_grid_major_y, theme_sieler2026_composition_figure,
#                             SIELER2026_MIN_LINEWIDTH_MM (ISME / manuscript)
#   Export (Rmd / tables) ... export_results.R — get_doc_name, export helpers for Results/
#   GT tables (this file, sec. 5) ... patterns from Sieler_2026 Code/analysis/*.Rmd — see INDEX in sec. 5
#   Alpha beta-GLMM (sec. 6) ... build_alpha_data_model, make_model_coef_gt, trend plots for 01__Diversity.R,
#                             glmmtmb_history_parasite_interaction_coef_row / _interaction_p
#   Composition / beta (sec. 7) ... composition_inferential_gt, microviz_permanova_to_tidy,
#                             composition_betadisper_anova_to_tidy, composition_pcoa_parasite_faceted_by_history_plot,
#                             composition_pcoa_parasite_faceted_points_by_regime_plot, exposure_regime_colors,
#                             composition_betadisper_boxplot_history, composition_betadisper_boxplot_history_parasite,
#                             composition_betadisper_boxplot_atp, composition_add_significance_bars (02__Composition.R)
#   MaAsLin3 / diff. abund. (sec. 8) ... maaslin_read_all_results_tsv, maaslin_top_taxa_treatment_table (03__DiffAbund.R)
#   Host × taxa integration — DEG×DAT workflows in 04__TaxonGeneNetworkHelpers.R (module 06).
#
# Migration tip: keep preprocessing utilities in this file; copy statistical “stacks” from
#   the archive only when a pipeline needs them, then refactor globals into arguments.
#
# Dated files: this script does not embed calendar strings for RDS names. Use resolve_latest_rds()
#   or objects from 00__InitializeEnvironment.R (path.pseq.*, ps.list). Experimental day columns
#   (e.g. time_point levels 0, 14, …, 60) are sampling schedule, not file dates.

# --- Publication ggplot2 theme (manuscript / ISME J) ---------------------------------
# ggplot2 linewidth is in mm; ~0.35 mm ≈ 1 pt — ISME asks not to use lines thinner than 1 pt.
SIELER2026_MIN_LINEWIDTH_MM <- 0.35

# Shared theme for main-text figures: classic, consistent base size, bottom legend.
theme_sieler2026_publication <- function(base_size = 14, legend_position = "bottom") {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_position,
      plot.subtitle = ggplot2::element_text(lineheight = 1.15),
      plot.margin = ggplot2::margin(6, 8, 6, 8, "pt")
    )
}

# Extra layers for 02__Composition figures: light horizontal y-grid + bold legend titles.
theme_sieler2026_composition_layers <- function() {
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(
      colour = grDevices::grey(0.87),
      linewidth = 0.35
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(face = "bold")
  )
}

# Publication theme plus composition grid / legend-title styling (Michael Sieler; 2026-04-24).
theme_sieler2026_composition_figure <- function(base_size = 14, legend_position = "bottom") {
  theme_sieler2026_publication(base_size = base_size, legend_position = legend_position) +
    theme_sieler2026_composition_layers()
}

# Publication theme plus light major grid (replaces theme_minimal for scatter/dot figures).
theme_sieler2026_publication_with_grid <- function(base_size = 14, legend_position = "bottom") {
  theme_sieler2026_publication(base_size = base_size, legend_position = legend_position) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        colour = grDevices::grey(0.88),
        linewidth = max(0.25, SIELER2026_MIN_LINEWIDTH_MM * 0.65)
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Horizontal major grid at y-axis breaks only (discrete x unchanged). Matches
# theme_sieler2026_composition_layers major-y styling; for stressor-history trend
# and binomial tank trend plots on theme_classic (Michael Sieler; 2026-04-25).
theme_sieler2026_trend_panel_grid_major_y <- function() {
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(
      colour = grDevices::grey(0.87),
      linewidth = 0.35
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )
}

# =============================================================================
# 1. General — formatting, significance labels, column cleanup
# =============================================================================

# Format p-values for tables (values like "<0.001" via scales).
p_val_format <- function(x) {
  z <- scales::pvalue_format()(x)
  z[!is.finite(x)] <- ""
  z
}

# =============================================================================
# 1b. Results archiving — copy → verify → clear (used by Code/01__Analysis drivers)
# =============================================================================

# Archive existing Results outputs for a module directory, then (optionally) clear the module
# output folders so the next run produces a clean set of files.
#
# Archive layout:
#   Results/_Archive/<timestamp>/<module_name>/{Figures,Tables,Stats}/<copied outputs...>
#
# Safety:
# - Copies (never moves) existing outputs
# - Verifies that every top-level item copied exists in the archive before deleting originals
# - Skips when there are no outputs to archive (default)
sieler2026_archive_module_outputs <- function(
    path_res_module,
    module_name = basename(path_res_module),
    subdirs = c("Figures", "Tables", "Stats"),
    archive_root = file.path(path.results, "_Archive"),
    clear_original = TRUE,
    skip_if_empty = TRUE,
    timestamp = format(Sys.time(), "%Y-%m-%d__%H%M%S")
) {
  if (is.null(path_res_module) || !nzchar(path_res_module) || !dir.exists(path_res_module)) {
    warning("Archive skipped: module results directory missing: ", path_res_module)
    return(invisible(NULL))
  }

  # Gather top-level outputs per subdir (ignore any existing Archive folders).
  to_archive <- list()
  for (sd in subdirs) {
    d <- file.path(path_res_module, sd)
    if (!dir.exists(d)) {
      next
    }
    items <- list.files(d, full.names = TRUE, all.files = FALSE, no.. = TRUE)
    items <- items[basename(items) != "Archive"]
    if (length(items) > 0L) {
      to_archive[[sd]] <- items
    }
  }

  if (length(to_archive) == 0L) {
    if (isTRUE(skip_if_empty)) {
      return(invisible(NULL))
    }
    warning("No outputs found to archive for module: ", module_name)
    return(invisible(NULL))
  }

  archive_base <- file.path(archive_root, timestamp, module_name)
  dir.create(archive_base, recursive = TRUE, showWarnings = FALSE)

  copy_failures <- character()
  verify_failures <- character()

  copy_one <- function(src, dest_dir) {
    dest <- file.path(dest_dir, basename(src))
    if (dir.exists(src)) {
      # Use system cp for directory trees (handles large MaAsLin output folders reliably).
      # cp -R src dest_dir  -> creates dest_dir/<basename(src)>.
      ok <- tryCatch(
        {
          rc <- system2("cp", c("-R", src, dest_dir), stdout = FALSE, stderr = FALSE)
          isTRUE(rc == 0)
        },
        error = function(e) FALSE
      )
      return(list(ok = ok, dest = dest))
    }
    ok <- suppressWarnings(file.copy(src, dest, recursive = FALSE, copy.mode = TRUE, copy.date = TRUE))
    list(ok = isTRUE(ok), dest = dest)
  }

  for (sd in names(to_archive)) {
    src_items <- to_archive[[sd]]
    dest_dir <- file.path(archive_base, sd)
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

    for (src in src_items) {
      res <- copy_one(src, dest_dir)
      if (!isTRUE(res$ok)) {
        copy_failures <- c(copy_failures, src)
      }
    }

    # Verify: every top-level item now exists in the archive destination.
    for (src in src_items) {
      dest <- file.path(dest_dir, basename(src))
      if (!file.exists(dest) && !dir.exists(dest)) {
        verify_failures <- c(verify_failures, dest)
      }
    }
  }

  if (length(copy_failures) > 0L || length(verify_failures) > 0L) {
    stop(
      "Archiving failed for module ", module_name, ".\n",
      if (length(copy_failures) > 0L) paste0("Copy failures (first 10):\n- ", paste(utils::head(copy_failures, 10), collapse = "\n- "), "\n") else "",
      if (length(verify_failures) > 0L) paste0("Verification failures (first 10):\n- ", paste(utils::head(verify_failures, 10), collapse = "\n- "), "\n") else "",
      "No original outputs were deleted."
    )
  }

  if (isTRUE(clear_original)) {
    for (sd in names(to_archive)) {
      for (src in to_archive[[sd]]) {
        # Remove files/directories that were archived, leaving any existing Archive folder intact.
        unlink(src, recursive = TRUE, force = TRUE)
      }
    }
  }

  message("Archived previous outputs -> ", archive_base)
  invisible(list(timestamp = timestamp, archive_base = archive_base, module_name = module_name))
}

# =============================================================================
# 1c. Manuscript MainFigures — mirror PNG/PDF into Manuscript/MainFigures/Exports
# =============================================================================

# Portable path to the manuscript export folder (requires `proj.path` from init).
sieler2026_path_main_figure_exports <- function() {
  if (!exists("proj.path", inherits = TRUE) || !nzchar(as.character(proj.path))) {
    stop("proj.path not found: source Code/00__Setup/00__InitializeEnvironment.R from the repo root.")
  }
  file.path(proj.path, "Manuscript", "MainFigures", "Exports")
}

# Export filename under Manuscript/MainFigures/Exports/ (panel label + original basename).
# Example: panel_id "2.1" + "genus_relative_abundance_by_treatment.png" ->
#   "FIG_02-1__genus_relative_abundance_by_treatment.png"
sieler2026_main_figure_export_filename <- function(panel_id, src_basename) {
  pid <- trimws(as.character(panel_id))
  segs <- strsplit(pid, ".", fixed = TRUE)[[1L]]
  if (length(segs) < 2L) {
    stop("main_figures_manifest: panel_id must look like '2.1' for export naming, got: ", pid)
  }
  major <- suppressWarnings(as.integer(segs[[1L]]))
  if (!is.finite(major) || major < 0L) {
    stop("main_figures_manifest: invalid panel_id for export naming: ", pid)
  }
  minor <- segs[[2L]]
  sprintf("FIG_%02d-%s__%s", major, minor, src_basename)
}

# Copy main-text figure assets listed in Manuscript/MainFigures/main_figures_manifest.csv
# into Manuscript/MainFigures/Exports/ (names: FIG_XX-Y__<original_basename>.png|.pdf).
#
# @param driver_script Optional basename or path; only rows whose manifest `driver_script`
#   basename matches are copied (e.g. `"02__Composition.R"`).
# @param panel_ids Optional character vector of `panel_id` values (e.g. `"2.4"`); if NULL, all
#   matching-driver rows are synced.
# @param manifest_path Optional absolute path to the CSV; default uses `proj.path`.
sieler2026_sync_main_figures_from_manifest <- function(
    driver_script = NULL,
    panel_ids = NULL,
    manifest_path = NULL
) {
  if (!exists("proj.path", inherits = TRUE) || !nzchar(as.character(proj.path))) {
    stop("proj.path not found: source Code/00__Setup/00__InitializeEnvironment.R from the repo root.")
  }
  if (is.null(manifest_path)) {
    manifest_path <- file.path(proj.path, "Manuscript", "MainFigures", "main_figures_manifest.csv")
  }
  if (!file.exists(manifest_path)) {
    warning("Main-figures manifest not found; skip sync: ", manifest_path)
    return(invisible(NULL))
  }

  dest_dir <- sieler2026_path_main_figure_exports()
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  m <- utils::read.csv(
    manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c(panel_id = "character")
  )
  if (nrow(m) < 1L) {
    return(invisible(NULL))
  }

  if (!is.null(driver_script) && nzchar(as.character(driver_script))) {
    drv_b <- basename(as.character(driver_script))
    db <- vapply(as.character(m[["driver_script"]]), basename, character(1L))
    m <- m[db == drv_b, , drop = FALSE]
  }
  if (nrow(m) < 1L) {
    return(invisible(NULL))
  }

  if (!is.null(panel_ids) && length(panel_ids) > 0L) {
    pid <- as.character(panel_ids)
    m <- m[as.character(m[["panel_id"]]) %in% pid, , drop = FALSE]
  }
  if (nrow(m) < 1L) {
    return(invisible(NULL))
  }

  copy_one <- function(src, label, panel_id, dest_name_override = NULL) {
    src_chr <- as.character(src)
    if (length(src_chr) != 1L || is.na(src_chr) || !nzchar(src_chr)) {
      return(invisible(NULL))
    }
    full <- if (grepl("^/", src_chr) || grepl("^[A-Za-z]:", src_chr)) {
      src_chr
    } else {
      file.path(proj.path, src_chr)
    }
    if (!file.exists(full)) {
      warning("Main-figure sync: missing source (", label, "): ", full)
      return(invisible(NULL))
    }
    dest_name <- if (!is.null(dest_name_override) && nzchar(as.character(dest_name_override))) {
      as.character(dest_name_override)
    } else {
      sieler2026_main_figure_export_filename(panel_id, basename(full))
    }
    dest <- file.path(dest_dir, dest_name)
    ok <- suppressWarnings(
      file.copy(full, dest, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    )
    if (!isTRUE(ok)) {
      warning("Main-figure sync: file.copy failed for ", full)
    }
    invisible(dest)
  }

  for (i in seq_len(nrow(m))) {
    row <- m[i, , drop = FALSE]
    pid <- as.character(row[["panel_id"]])
    fig_id <- if ("figure_id" %in% names(row)) as.character(row[["figure_id"]]) else NA_character_
    drv <- if ("driver_script" %in% names(row)) tolower(trimws(as.character(row[["driver_script"]]))) else NA_character_

    # Special exceptions for external manuscript assets:
    # - Figure 01-1: ExperimentalDesignSchematic (keeps standard FIG_ naming via panel_id 1.1)
    # - Table 01-1: ExposureSchematic (panel_id 1.0, but exported as TABLE_01-1__... by override)
    is_external_asset <- isTRUE(drv == "external_design_asset")
    is_table_01_1 <- is_external_asset &&
      isTRUE(grepl("^table\\s*1\\b", fig_id, ignore.case = TRUE)) &&
      identical(pid, "1.0")

    if ("path_png" %in% names(row)) {
      png_src <- row[["path_png"]]
      png_dest_name_override <- NULL
      if (is_table_01_1) {
        png_dest_name_override <- "TABLE_01-1__ExposureSchematic.png"
      }
      png_dest <- copy_one(png_src, paste0("panel ", pid, " PNG"), pid, dest_name_override = png_dest_name_override)

      # For Table 01-1, additionally create a PDF for easy viewing/combining.
      if (is_table_01_1 && !is.null(png_dest) && file.exists(png_dest)) {
        pdf_dest <- file.path(dest_dir, "TABLE_01-1__ExposureSchematic.pdf")
        if (!requireNamespace("png", quietly = TRUE)) {
          stop(
            "Package 'png' is required to generate a PDF from the table PNG.\n",
            "Install with: install.packages('png')"
          )
        }
        img <- png::readPNG(png_dest)
        w_in <- max(1, ncol(img) / 100)
        h_in <- max(1, nrow(img) / 100)
        grDevices::pdf(pdf_dest, width = w_in, height = h_in, useDingbats = FALSE)
        grid::grid.newpage()
        grid::grid.draw(
          grid::rasterGrob(img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
        )
        grDevices::dev.off()
      }
    }
    if ("path_pdf" %in% names(row)) {
      # Table 01-1 has no source PDF in the manifest; it is generated from the PNG above.
      if (!is_table_01_1) {
        copy_one(row[["path_pdf"]], paste0("panel ", pid, " PDF"), pid)
      }
    }
  }

  message("Synced main-text figures for ", nrow(m), " manifest row(s) -> ", dest_dir)
  invisible(dest_dir)
}

# -----------------------------------------------------------------------------
# 1c.2 MaAsLin mortality — combined focal-taxa scatter helpers (03__DiffAbund.R)
# -----------------------------------------------------------------------------
#
# Created by: Michael Sieler
# Date last updated: 2026-04-27
#
# These helpers are shared by:
# - `Code/01__Analysis/03__DiffAbund.R` (`--figures-only` fast path)
# - `Code/01__Analysis/98__MainFiguresRefresh.R` (refresh panel_id `5.1` without a full 03 run)
#
# Expected globals: `treatment_color_scale` (optional) from `02__PlotSettings.R` via init.

get_treatment_colors_safe <- function(treatment_levels) {
  if (exists("treatment_color_scale", inherits = TRUE) &&
      is.atomic(treatment_color_scale) &&
      length(treatment_color_scale) >= 1L) {
    cols <- treatment_color_scale
    if (!is.null(names(cols)) && any(names(cols) %in% treatment_levels)) {
      cols <- cols[treatment_levels]
      cols[is.na(cols)] <- "grey60"
      return(cols)
    }
    if (length(cols) >= length(treatment_levels)) {
      names(cols) <- treatment_levels
      return(cols)
    }
  }
  cols <- grDevices::hcl.colors(length(treatment_levels), palette = "Dark 3")
  names(cols) <- treatment_levels
  cols
}

diffabund_build_mortality_combined_scatter_plot <- function(tank_df) {
  treatment_levels <- sort(unique(as.character(tank_df$Treatment)))
  cols_treat <- get_treatment_colors_safe(treatment_levels)
  taxon_ord <- c("Culicoidibacter", "Shewanella", "Flavobacterium", "Cetobacterium")
  taxon_present <- intersect(taxon_ord, unique(as.character(tank_df$Taxon)))
  tank_plot <- tank_df %>%
    dplyr::mutate(
      Taxon = factor(as.character(.data$Taxon), levels = taxon_present)
    )

  smooth_pred <- tank_plot %>%
    dplyr::group_by(.data$Taxon) %>%
    dplyr::group_modify(~ {
      ok <- stats::complete.cases(.x$log2_abund, .x$percent_mortality)
      if (sum(ok) < 3L) {
        return(tibble::tibble(log2_abund = numeric(), pred = numeric()))
      }
      d <- .x[ok, , drop = FALSE]
      fit <- stats::lm(percent_mortality ~ log2_abund, data = d)
      xr <- range(d$log2_abund, na.rm = TRUE)
      xg <- seq(xr[[1L]], xr[[2L]], length.out = 80L)
      tibble::tibble(
        log2_abund = xg,
        pred = as.numeric(stats::predict(fit, newdata = data.frame(log2_abund = xg)))
      )
    }) %>%
    dplyr::ungroup()

  lbl_df <- smooth_pred %>%
    dplyr::group_by(.data$Taxon) %>%
    dplyr::slice_max(order_by = .data$log2_abund, n = 1L, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      y_repel = pmin(.data$pred + 5, 96)
    )

  pal_set2 <- c(
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
  )
  set2_fill <- pal_set2[seq_along(treatment_levels)]
  names(set2_fill) <- treatment_levels

  taxon_line_colors <- c(
    Culicoidibacter = "#E41A1C",
    Shewanella = "#377EB8",
    Flavobacterium = "#4DAF4A",
    Cetobacterium = "#984EA3"
  )
  taxon_col_vec <- taxon_line_colors[taxon_present]
  taxon_col_vec[is.na(taxon_col_vec)] <- "grey35"
  names(taxon_col_vec) <- taxon_present

  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Install ggnewscale for the combined mortality scatter (see Code/00__Setup/01__Libraries.R).")
  }

  set.seed(42)
  ggplot2::ggplot(
    tank_plot,
    ggplot2::aes(
      x = .data$log2_abund,
      y = .data$percent_mortality
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.6) +
    ggplot2::geom_smooth(
      ggplot2::aes(group = .data$Taxon, fill = .data$Taxon),
      method = "lm",
      se = TRUE,
      linewidth = 0,
      color = NA,
      alpha = 0.32,
      show.legend = FALSE
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(group = .data$Taxon),
      method = "lm",
      se = FALSE,
      color = "white",
      linewidth = 2.08,
      show.legend = FALSE
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(group = .data$Taxon, color = .data$Taxon),
      method = "lm",
      se = FALSE,
      linewidth = 0.92,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = taxon_col_vec, drop = FALSE, guide = "none") +
    ggplot2::scale_colour_manual(
      values = taxon_col_vec,
      breaks = taxon_present,
      drop = FALSE,
      name = "Focal Taxa",
      guide = ggplot2::guide_legend(
        order = 2L,
        nrow = 1L,
        override.aes = list(
          shape = 21L,
          size = 3.2,
          alpha = 1,
          stroke = 0.45,
          fill = unname(taxon_col_vec[taxon_present]),
          colour = unname(taxon_col_vec[taxon_present])
        )
      )
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(
      shape = 21L,
      fill = "white",
      color = "white",
      stroke = 0.5,
      size = 3.45,
      alpha = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = .data$Treatment,
        color = .data$Taxon
      ),
      shape = 21L,
      stroke = 0.45,
      size = 2.8,
      alpha = 0.5
    ) +
    ggplot2::scale_fill_manual(
      values = set2_fill,
      drop = FALSE,
      name = "Exposure Regime",
      guide = ggplot2::guide_legend(
        order = 1L,
        nrow = 2L,
        byrow = TRUE,
        override.aes = list(
          alpha = 1,
          shape = 21L,
          size = 3.2,
          stroke = 0.45,
          colour = "grey30"
        )
      )
    ) +
    ggrepel::geom_label_repel(
      data = lbl_df,
      mapping = ggplot2::aes(
        x = .data$log2_abund,
        y = .data$y_repel,
        label = as.character(.data$Taxon),
        color = .data$Taxon
      ),
      inherit.aes = FALSE,
      fill = ggplot2::alpha("white", 0.94),
      size = 3.4,
      fontface = "bold",
      label.size = 0.35,
      label.padding = grid::unit(0.35, "lines"),
      min.segment.length = Inf,
      segment.size = 0.3,
      segment.color = ggplot2::alpha("grey40", 0.55),
      box.padding = ggplot2::unit(0.5, "lines"),
      point.padding = ggplot2::unit(0.35, "lines"),
      max.overlaps = Inf,
      seed = 42L,
      show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(limits = c(-25, 100), breaks = seq(0, 100, 25)) +
    ggplot2::labs(
      title = "Significantly associated genera with host mortality",
      subtitle = NULL,
      x = "Mean Log2 Taxon Abundance (per tank)",
      y = "Percent Mortality (per tank)"
    ) +
    theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.15, "cm"),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(t = 6, r = 8, b = 10, l = 8, unit = "pt")
    ) +
    ggplot2::coord_cartesian(clip = "off")
}

# Add significance stars from a p-value column (default name "p.value").
SigStars <- function(x, pval.var = "p.value") {
  x %>%
    dplyr::rename(p.value = !!rlang::sym(pval.var)) %>%
    dplyr::mutate(
      p.adj.sig = dplyr::case_when(
        p.value < 0.0001 ~ "****",
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value >= 0.05 ~ "ns"
      )
    )
}

# Strip suffix after a separator from selected column names (e.g. ASV__Genus -> ASV).
cutColNames <- function(df, cols, sep = "__") {
  df %>%
    dplyr::rename_with(~ sub(paste0(sep, ".*"), "", .x), .cols = dplyr::all_of(cols))
}

# Strip suffix after separator within selected columns (cell values).
cutCellNames <- function(df, col = c(), sep = "__") {
  df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(col), ~ sub(paste0(sep, ".*"), "", .x)))
}

# =============================================================================
# 2. Phyloseq — sample_data extraction, renames, long tables for modeling
#    (includes 2b: DADA2 post-processing filters + RoL metadata augmentation)
# =============================================================================

# Sample metadata as a base data.frame (avoids S4 quirks in some joins).
samdatAsDataframe <- function(ps) {
  samdat <- phyloseq::sample_data(ps)
  data.frame(samdat, check.names = FALSE, stringsAsFactors = FALSE)
}

# Rename sample_data columns via dplyr::rename semantics (legacy: miscFunctions.R).
# Prefer this over microViz::ps_rename for column renames so behavior matches archived pipelines.
ps_rename <- function(ps, ...) {
  ps <- microViz::ps_get(ps)
  df <- samdatAsDataframe(ps)
  df <- dplyr::rename(.data = df, ...)
  phyloseq::sample_data(ps) <- df
  ps
}

# Pivot normalized diversity columns long for stats (see microViz ps_mutate + norm_* patterns).
psObjToDfLong <- function(
    ps.obj,
    div.score,
    div.metric,
    pivot.long_col = "_norm"
) {
  ps.obj %>%
    microViz::samdat_tbl() %>%
    tidyr::pivot_longer(
      cols = dplyr::contains(pivot.long_col),
      names_to = div.metric,
      values_to = div.score
    ) %>%
    dplyr::ungroup()
}

# --- 2a. Newest dated RDS in a directory (same mtime rule as ps-list__*.rds in 00__InitializeEnvironment.R)
#
# Picks one path by latest file modification time — not by parsing the filename date — so re-touching
#   a file updates selection. For empty matches: returns NA when error_if_empty is FALSE; stops when TRUE.
resolve_latest_rds <- function(dir, pattern, error_if_empty = FALSE) {
  if (!dir.exists(dir)) {
    if (error_if_empty) {
      stop("Directory does not exist: ", dir)
    }
    return(NA_character_)
  }
  candidates <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(candidates) == 0L) {
    if (error_if_empty) {
      stop("No files matching ", encodeString(pattern, quote = "\""), " in ", dir)
    }
    return(NA_character_)
  }
  candidates[[which.max(file.info(candidates)$mtime)]]
}

# --- 2b. DADA2 post-processing (Sieler2026 / Rules-of-Life zebrafish gut 16S)
#
# Experiment (high level): adult zebrafish, factorial antibiotic (A) × temperature (T) × parasite (P)
#   exposure over 60 d; fecal 16S at several days. Design metadata columns expected in sample_data
#   include Antibiotics, Temperature, Parasite (0/1), Time or Timepoint (days), Tank.ID, fecal IDs, etc.
#   Functions below turn a raw merged phyloseq into analysis-ready objects used by 04__DataPreProcess.R
#   and downstream drivers. See Data/Context/ExperimentalDesignContext.md for full rationale.

# Drop sparse taxa, suspect organelle/host reads, placeholder phyla, then shallow samples.
# Order: (1) ASVs present in fewer than min_prevalence samples, (2) rows matching Mitochondria /
#   Chloroplast / Eukaryota anywhere in taxonomy, (3) unresolved bacterial phylum labels, (4) samples
#   with total reads below min_reads. Returns a smaller phyloseq; stops if fewer than two samples remain.
apply_pre_analysis_filters <- function(ps_input, min_prevalence = 2, min_reads = 1000) {
  otu <- as(phyloseq::otu_table(ps_input), "matrix")
  if (!phyloseq::taxa_are_rows(ps_input)) {
    otu <- t(otu)
  }
  taxa_prev <- rowSums(otu > 0)
  ps_out <- phyloseq::prune_taxa(taxa_prev >= min_prevalence, ps_input)

  # Taxonomy-driven removal only when tax_table exists (ASV-only pipelines skip this block).
  if (!is.null(phyloseq::tax_table(ps_out, errorIfNULL = FALSE))) {
    tax_df <- as.data.frame(as(phyloseq::tax_table(ps_out), "matrix"), stringsAsFactors = FALSE)
    rm_target <- apply(
      tax_df,
      1,
      function(x) any(grepl("Mitochondria|Chloroplast|Eukaryota", x, ignore.case = TRUE))
    )

    phylum_col <- if ("Phylum" %in% colnames(tax_df)) {
      "Phylum"
    } else if (ncol(tax_df) >= 2) {
      colnames(tax_df)[2]
    } else {
      NA_character_
    }
    unresolved <- rep(FALSE, nrow(tax_df))
    if (!is.na(phylum_col)) {
      phylum_vals <- tax_df[[phylum_col]]
      unresolved <- is.na(phylum_vals) |
        phylum_vals == "" |
        grepl("^Kingdom_Bacteria$", phylum_vals, ignore.case = TRUE) |
        grepl("^Bacteria Phylum$", phylum_vals, ignore.case = TRUE)
    }
    keep <- !(rm_target | unresolved)
    ps_out <- phyloseq::prune_taxa(rownames(tax_df)[keep], ps_out)
  }

  keep_samples <- phyloseq::sample_sums(ps_out) >= min_reads
  if (sum(keep_samples) < 2) {
    stop(
      "Low-read filter retained fewer than 2 samples. ",
      "Lower min_reads or inspect input phyloseq."
    )
  }
  ps_out <- phyloseq::prune_samples(keep_samples, ps_out)
  ps_out <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_out) > 0, ps_out)
  ps_out
}

# --- 2b. Contaminant removal (optional) ---------------------------------------
#
# Uses decontam::isContaminant(method="prevalence") to flag ASVs enriched in blanks
# (PCR blanks / kit blanks / other negative controls), then removes those taxa.
#
# Notes:
# - This does NOT require concentration data (unlike the "frequency" method).
# - Caller must supply a logical vector `neg` (length nsamples) marking blanks.
# - By default, it keeps blanks in the returned object unless remove_neg_samples is TRUE.
sieler2026_decontam_prevalence <- function(
    ps_in,
    neg,
    threshold = 0.5,
    remove_neg_samples = TRUE
) {
  if (!requireNamespace("decontam", quietly = TRUE)) {
    warning("Package `decontam` not installed; skipping contaminant removal.")
    return(list(ps = ps_in, contam_tbl = tibble::tibble()))
  }
  ps_in <- microViz::ps_get(ps_in)

  if (!is.logical(neg) || length(neg) != phyloseq::nsamples(ps_in)) {
    stop("`neg` must be a logical vector with length equal to nsamples(ps_in).")
  }
  if (sum(neg, na.rm = TRUE) < 1) {
    warning("No negative controls detected (sum(neg) < 1); skipping contaminant removal.")
    return(list(ps = ps_in, contam_tbl = tibble::tibble()))
  }

  contam <- decontam::isContaminant(
    ps_in,
    method = "prevalence",
    neg = neg,
    threshold = threshold
  )

  contam_tbl <- contam |>
    tibble::rownames_to_column("taxon") |>
    tibble::as_tibble() |>
    dplyr::arrange(dplyr::desc(.data$contaminant), dplyr::desc(.data$p))

  keep_taxa <- rownames(contam)[!contam$contaminant]
  ps_out <- phyloseq::prune_taxa(keep_taxa, ps_in)
  ps_out <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_out) > 0, ps_out)

  if (isTRUE(remove_neg_samples)) {
    ps_out <- phyloseq::prune_samples(!neg, ps_out)
  }

  list(ps = ps_out, contam_tbl = contam_tbl)
}

# Replace arbitrary ASV names (e.g. long hashes from DADA2) with stable ASV0001, ASV0002, …
#   so tables and plots stay compact. Skips if names already match ASV + four digits.
rename_asvs_sequential <- function(ps_in) {
  tn <- phyloseq::taxa_names(ps_in)
  if (all(grepl("^ASV[0-9]{4}$", tn))) {
    return(ps_in)
  }
  set.seed(42)
  new_names <- paste0("ASV", sprintf("%04d", seq_along(tn)))
  taxa_names(ps_in) <- new_names
  ps_in
}

# Older metadata spreadsheets used "Pathogen"; analysis code expects "Parasite". Rename only when
#   Pathogen is present and Parasite is not, so re-runs on already-fixed data are idempotent.
ensure_parasite_column <- function(ps_in) {
  cols <- colnames(microViz::samdat_tbl(ps_in))
  if ("Pathogen" %in% cols && !"Parasite" %in% cols) {
    return(ps_rename(ps_in, Parasite = Pathogen))
  }
  ps_in
}

# After read-depth / taxon filtering: add human-readable regime labels and model covariates.
#   Treatment: eight strings "A± T± P±" from binary Antibiotics, Temperature, Parasite.
#   HistoryLevel: count of *prior* stressors before the parasite phase (A and T only): 0, 1, or 2.
#     Explicit mapping by Treatment:
#       A- T- P- = 0, A- T- P+ = 0
#       A+ T- P- = 1, A+ T- P+ = 1
#       A- T+ P- = 1, A- T+ P+ = 1
#       A+ T+ P- = 2, A+ T+ P+ = 2
#   Parasite_Exposed / P in Treatment: parasite challenge vs not. Exp_Type groups regimes for plots.
#   Requires sample_data columns Parasite, Antibiotics, Temperature, fecal.sample.number, Sex, Time.
#   treatment_order_vec: factor levels for Treatment (must match project conventions).
augment_filtered_phyloseq_metadata <- function(ps_in, treatment_order_vec) {
  ps_step <- ps_in
  sd_cols <- colnames(microViz::samdat_tbl(ps_step))
  if ("Timepoint" %in% sd_cols && !"Time" %in% sd_cols) {
    ps_step <- ps_rename(ps_step, Time = Timepoint)
  }

  ps_step %>%
    microViz::ps_mutate(
      Treatment = dplyr::case_when(
        Antibiotics == 0 & Temperature == 0 & Parasite == 0 ~ "A- T- P-",
        Antibiotics == 0 & Temperature == 0 & Parasite == 1 ~ "A- T- P+",
        Antibiotics == 0 & Temperature == 1 & Parasite == 0 ~ "A- T+ P-",
        Antibiotics == 0 & Temperature == 1 & Parasite == 1 ~ "A- T+ P+",
        Antibiotics == 1 & Temperature == 0 & Parasite == 0 ~ "A+ T- P-",
        Antibiotics == 1 & Temperature == 0 & Parasite == 1 ~ "A+ T- P+",
        Antibiotics == 1 & Temperature == 1 & Parasite == 0 ~ "A+ T+ P-",
        Antibiotics == 1 & Temperature == 1 & Parasite == 1 ~ "A+ T+ P+",
        TRUE ~ "Unknown"
      ),
      .after = "Parasite"
    ) %>%
    microViz::ps_mutate(Sample = fecal.sample.number, .before = 1) %>%
    microViz::ps_mutate(Sample = gsub("^f", "", Sample)) %>%
    microViz::ps_filter(Treatment != "Unknown") %>%
    microViz::ps_mutate(
      History = dplyr::case_when(
        Treatment %in% c("A- T- P-", "A- T- P+") ~ 0L,
        Treatment %in% c("A+ T- P-", "A+ T- P+", "A- T+ P-", "A- T+ P+") ~ 1L,
        Treatment %in% c("A+ T+ P-", "A+ T+ P+") ~ 2L,
        TRUE ~ NA_integer_
      ),
      .after = "Treatment"
    ) %>%
    microViz::ps_mutate(
      treatment_code = dplyr::case_when(
        Antibiotics == 0 & Temperature == 0 & Parasite == 0 ~ "Aneg_Tneg_Pneg",
        Antibiotics == 0 & Temperature == 0 & Parasite == 1 ~ "Aneg_Tneg_Ppos",
        Antibiotics == 1 & Temperature == 0 & Parasite == 0 ~ "Apos_Tneg_Pneg",
        Antibiotics == 1 & Temperature == 0 & Parasite == 1 ~ "Apos_Tneg_Ppos",
        Antibiotics == 0 & Temperature == 1 & Parasite == 0 ~ "Aneg_Tpos_Pneg",
        Antibiotics == 0 & Temperature == 1 & Parasite == 1 ~ "Aneg_Tpos_Ppos",
        Antibiotics == 1 & Temperature == 1 & Parasite == 0 ~ "Apos_Tpos_Pneg",
        Antibiotics == 1 & Temperature == 1 & Parasite == 1 ~ "Apos_Tpos_Ppos"
      ),
      treatment_group = dplyr::case_when(
        Antibiotics == 0 & Temperature == 0 & Parasite == 1 ~ "Parasite",
        Antibiotics == 1 & Temperature == 0 & Parasite == 0 ~ "Antibiotics",
        Antibiotics == 1 & Temperature == 0 & Parasite == 1 ~ "Antibiotics_Parasite",
        Antibiotics == 0 & Temperature == 1 & Parasite == 0 ~ "Temperature",
        Antibiotics == 0 & Temperature == 1 & Parasite == 1 ~ "Temperature_Parasite",
        Antibiotics == 1 & Temperature == 1 & Parasite == 0 ~ "Antibiotics_Temperature",
        Antibiotics == 1 & Temperature == 1 & Parasite == 1 ~ "Antibiotics_Temperature_Parasite",
        TRUE ~ "Control"
      ),
      treatment_group = factor(
        treatment_group,
        levels = c(
          "Control", "Parasite",
          "Antibiotics", "Antibiotics_Parasite",
          "Temperature", "Temperature_Parasite",
          "Antibiotics_Temperature", "Antibiotics_Temperature_Parasite"
        )
      ),
      treatment_code = factor(treatment_code, levels = treatment_order_vec),
      # Experimental sampling days (days since study start), not RDS or calendar dates.
      time_point = factor(Time, levels = c(0, 14, 18, 25, 29, 60)),
      parasite_status = factor(
        ifelse(Parasite == 1, "Exposed", "Unexposed"),
        levels = c("Unexposed", "Exposed")
      ),
      sex = factor(Sex, levels = c("M", "F"))
    ) %>%
    microViz::ps_mutate(Treatment = factor(Treatment, levels = treatment_order_vec)) %>%
    microViz::ps_mutate(
      Exp_Type = dplyr::case_when(
        Treatment %in% c("A- T- P-", "A- T- P+") ~ "No prior stressor(s)",
        Treatment %in% c("A+ T- P-", "A+ T- P+") ~ "Antibiotics",
        Treatment %in% c("A- T+ P-", "A- T+ P+") ~ "Temperature",
        Treatment %in% c("A+ T+ P-", "A+ T+ P+") ~ "Combined"
      ),
      Exp_Type = factor(
        Exp_Type,
        levels = c("No prior stressor(s)", "Antibiotics", "Temperature", "Combined")
      )
    ) %>%
    microViz::ps_mutate(
      A = stringr::str_extract(Treatment, "A[+-]") %>% stringr::str_sub(2, 2),
      T = stringr::str_extract(Treatment, "T[+-]") %>% stringr::str_sub(2, 2),
      P = stringr::str_extract(Treatment, "P[+-]") %>% stringr::str_sub(2, 2)
    ) %>%
    microViz::ps_mutate(
      A = factor(A, levels = c("-", "+")),
      T = factor(T, levels = c("-", "+")),
      P = factor(P, levels = c("-", "+"))
    ) %>%
    microViz::ps_mutate(
      HistoryLevelNum = as.integer(.data$History),
      HistoryLevel = factor(HistoryLevelNum, levels = c(0, 1, 2)),
      Parasite_Exposed = dplyr::case_when(
        P == "-" ~ "Unexposed",
        P == "+" ~ "Exposed",
        TRUE ~ "Unknown"
      ),
      Parasite_Exposed = factor(Parasite_Exposed, levels = c("Unexposed", "Exposed")),
      History_Label = dplyr::case_when(
        HistoryLevelNum == 0 ~ "No prior stressors",
        HistoryLevelNum == 1 ~ "One prior stressor",
        HistoryLevelNum == 2 ~ "Two prior stressors",
        TRUE ~ "Unknown"
      ),
      History_Label = factor(
        History_Label,
        levels = c("No prior stressors", "One prior stressor", "Two prior stressors")
      ),
      Treatment. = as.numeric(factor(Treatment, levels = treatment_order_vec))
    ) %>%
    microViz::ps_mutate(
      History.Parasite = paste0(as.character(HistoryLevel), "_", Parasite)
    )
}

# OTU matrix for picante: samples as rows, taxa as columns (picante::pd expectation).
.ps_otu_mat_picante <- function(ps) {
  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) {
    t(mat)
  } else {
    mat
  }
}

# =============================================================================
# 3. Alpha diversity — normalization and phylogenetic diversity on phyloseq
# =============================================================================

# Normalize a numeric vector to [0, 1], using Tukey transform when Anderson–Darling rejects normality.
# Requires the `nortest` package (install.packages("nortest") if missing).
norm_scores <- function(x) {
  if (!requireNamespace("nortest", quietly = TRUE)) {
    stop("Package `nortest` is required for norm_scores(); install it with install.packages(\"nortest\")")
  }
  if (nortest::ad.test(x)$p.value <= 0.05) {
    x.trans <- rcompanion::transformTukey(x, plotit = FALSE, quiet = TRUE, statistic = 2)
    x.trans.norm <- (x.trans - min(x.trans, na.rm = TRUE)) /
      (max(x.trans, na.rm = TRUE) - min(x.trans, na.rm = TRUE))
    if (nortest::ad.test(x.trans.norm)$p.value < 0.05) {
      x.trans.norm <- (x.trans.norm - min(x.trans.norm, na.rm = TRUE)) /
        (max(x.trans.norm, na.rm = TRUE) - min(x.trans.norm, na.rm = TRUE))
      return(x.trans.norm)
    }
    return(x.trans.norm)
  }
  x.norm <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  x.norm
}

# Faith’s phylogenetic diversity (PD) via picante, joined into sample_data.
ps_calc_diversity.phy <- function(
    ps,
    index = "Phylogenetic",
    varname = index
) {
  df <- picante::pd(
    samp = .ps_otu_mat_picante(ps),
    tree = phyloseq::phy_tree(ps)
  ) %>%
    dplyr::rename(!!varname := PD) %>%
    dplyr::select(-SR)

  df[[".rownames."]] <- rownames(df)

  if (varname %in% phyloseq::sample_variables(ps)) {
    warning(varname, " is already in sample_data — overwriting.")
    phyloseq::sample_data(ps)[[varname]] <- NULL
  }

  microViz::ps_join(
    x = ps,
    y = df,
    type = "left",
    .keep_all_taxa = TRUE,
    match_sample_names = ".rownames.",
    keep_sample_name_col = FALSE
  )
}

# --- 3b. Batch alpha / beta for ps.list (replaces sourcing Code/99__Archive/.../AlphaDiversity.R and BetaDiversity.R in 04__DataPreProcess.R)

# Genus-level Shannon / Simpson / richness + norm_scores columns (legacy AlphaDiversity.R pipeline).
add_alpha_diversity_genus_norm_metrics <- function(ps, rank = "Genus") {
  sfx <- paste0("__", rank)
  ps %>%
    microViz::ps_calc_diversity(
      rank = rank,
      index = "shannon",
      varname = paste0("Shannon", sfx),
      exp = TRUE
    ) %>%
    microViz::ps_calc_diversity(
      rank = rank,
      index = "inverse_simpson",
      varname = paste0("Simpson", sfx)
    ) %>%
    microViz::ps_calc_richness(
      rank = rank,
      varname = paste0("Richness", sfx)
    ) %>%
    microViz::ps_mutate(
      dplyr::across(dplyr::contains(sfx), norm_scores, .names = "{.col}_norm")
    )
}

# Apply add_alpha_diversity_genus_norm_metrics to every element of ps_list (e.g. All, Unexposed, …).
populate_ps_list_alpha_diversity <- function(ps_list_in, rank = "Genus") {
  nm <- names(ps_list_in)
  out <- ps_list_in
  for (k in nm) {
    message("Calculating alpha diversity scores: ", k)
    out[[k]] <- add_alpha_diversity_genus_norm_metrics(ps_list_in[[k]], rank = rank)
  }
  out
}

# Parallel distance matrices per beta method (legacy BetaDiversity.R); returns a nested list named like ps_list.
build_beta_dist_matrices_for_ps_list <- function(
    ps_list_in,
    beta_methods,
    gunifrac_alpha = 0.5,
    rank_default = "Genus"
) {
  dc <- parallel::detectCores()
  n_workers <- max(1L, (if (is.na(dc)) 1L else as.integer(dc)) - 1L)
  old_plan <- future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(old_plan), add = TRUE)

  out <- list()
  for (nm in names(ps_list_in)) {
    message("Calculating beta diversity distance matrices: ", nm)
    tmp_ps <- ps_list_in[[nm]]
    out[[nm]] <- furrr::future_map(beta_methods, function(method) {
      rank_use <- if (method %in% c("gunifrac", "wunifrac", "unifrac")) {
        "unique"
      } else {
        rank_default
      }
      tmp_ps %>%
        microViz::tax_transform(trans = "identity", rank = rank_use) %>%
        microViz::dist_calc(method, gunifrac_alpha = gunifrac_alpha)
    })
    names(out[[nm]]) <- beta_methods
  }
  out
}

# =============================================================================
# 4. Optional — load legacy helper bundles from Archive (commented)
# =============================================================================
# When you need full beta-diversity or GLM helpers that still live in the archive,
# source the relevant file after `01__Libraries.R` (and plot settings if ggplot helpers), e.g.:
#
# legacy_hf <- here::here("Code", "99__Archive", "Functions", "HelperFunctions")
# source(file.path(legacy_hf, "betaFunctions.R"))
# source(file.path(legacy_hf, "alphaFunctions.R"))  # GLM / Tukey / Levene blocks
# source(file.path(legacy_hf, "plotFunctions.R"))
# source(file.path(legacy_hf, "diffAbundFunctions.R"))
# source(file.path(legacy_hf, "export_results.R"))
# DEG×DAT: Code/00__Setup/04__TaxonGeneNetworkHelpers.R (sourced by module 06).

# =============================================================================
# 5. GT tables — shared helpers for analysis Rmds (tidyverse + gt)
# =============================================================================
# Used across 01 Mortality, 02 Diversity, 03 Diff abund, 04 Taxon×mortality, 05 DEG, 06 networks,
# 07 Functional annotation: compose small helpers instead of duplicating pipelines.
#
# Typical pattern:
#   df %>% gt::gt() %>% gt::tab_header(...) %>% gt_fmt_numbers() %>% style_gt_significance(df)
#
# Analysis-specific layouts (row groups, data_color, HTML labels) stay in the Rmd; use the
# helpers below for p-value detection, significance fill, bold headers, and wide tables.

# Column names that look like p-values / FDR (for auto-styling or scientific notation).
gt_pvalue_column_names <- function(data_tbl) {
  cn <- names(data_tbl)
  cn[stringr::str_detect(
    stringr::str_to_lower(cn),
    "p[._]?value|^p$|^p[._]|q[._]?value|fdr|adj.*p|padj|pr\\(>|pval"
  )]
}

# Bold column labels (04 Taxon×mortality, 06 network summary tables).
gt_bold_column_labels <- function(gt_tbl) {
  gt_tbl %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    )
}

# Constrain table width for HTML overflow (07 Functional annotation uses lightblue + width).
gt_tab_width_px <- function(gt_tbl, width_px = 1200L) {
  gt_tbl %>%
    gt::tab_options(table.width = gt::px(width_px))
}

# Format all numeric columns to fixed decimals; optional scientific for tiny p-value columns.
gt_fmt_numbers <- function(
    gt_tbl,
    data_tbl,
    decimals = 3L,
    scientific_cols = NULL,
    scientific_decimals = 2L,
    pvalue_threshold = 1e-4
) {
  num_cols <- names(data_tbl)[vapply(data_tbl, is.numeric, logical(1L))]
  if (length(num_cols) == 0L) {
    return(gt_tbl)
  }
  if (is.null(scientific_cols)) {
    scientific_cols <- intersect(num_cols, gt_pvalue_column_names(data_tbl))
  }
  plain_cols <- setdiff(num_cols, scientific_cols)
  if (length(plain_cols) > 0L) {
    gt_tbl <- gt_tbl %>%
      gt::fmt_number(columns = dplyr::all_of(plain_cols), decimals = decimals)
  }
  if (length(scientific_cols) > 0L) {
    # For p-value-like columns, prefer "<0.0001" style instead of scientific notation.
    gt_tbl <- gt_tbl %>%
      gt::fmt(
        columns = dplyr::all_of(scientific_cols),
        fns = function(x) {
          x_num <- suppressWarnings(as.numeric(x))
          dplyr::case_when(
            is.na(x_num) ~ NA_character_,
            x_num < pvalue_threshold ~ paste0("<", sprintf("%.4f", pvalue_threshold)),
            TRUE ~ sprintf(paste0("%.", scientific_decimals + 2, "f"), x_num)
          )
        }
      )
  }
  gt_tbl
}

# Highlight significant p/q/FDR cells (rows where value < alpha) — from 01__Mortality_Infection.Rmd.
style_gt_significance <- function(
    gt_tbl,
    data_tbl,
    alpha = 0.05,
    sig_fill = "#e6f4ea",
    sig_text_color = "#1b7f3a",
    columns = NULL
) {
  sig_cols <- if (is.null(columns)) {
    gt_pvalue_column_names(data_tbl)
  } else {
    columns
  }

  if (length(sig_cols) == 0L) {
    return(gt_tbl)
  }

  for (sig_col in sig_cols) {
    sig_values <- suppressWarnings(as.numeric(data_tbl[[sig_col]]))
    sig_rows <- which(!is.na(sig_values) & sig_values < alpha)

    if (length(sig_rows) > 0L) {
      gt_tbl <- gt_tbl %>%
        gt::tab_style(
          style = list(
            gt::cell_fill(color = sig_fill),
            gt::cell_text(color = sig_text_color, weight = "bold")
          ),
          locations = gt::cells_body(
            columns = dplyr::all_of(sig_col),
            rows = sig_rows
          )
        )
    }
  }

  gt_tbl
}

# Data frame or model summary (e.g. emmeans::summary) → styled GT; same idea as 01__Mortality_Infection.
create_stat_gt <- function(
    result_object,
    table_title,
    table_subtitle = NULL,
    alpha = 0.05,
    sig_fill = "#e6f4ea",
    decimals = 3L,
    scientific_cols = NULL
) {
  result_df <- if (inherits(result_object, "data.frame")) {
    result_object
  } else {
    summary(result_object, infer = c(TRUE, TRUE)) %>% as.data.frame()
  }

  gt_tbl <- result_df %>%
    gt::gt() %>%
    gt::tab_header(title = table_title, subtitle = table_subtitle)

  gt_tbl <- gt_fmt_numbers(
    gt_tbl,
    data_tbl = result_df,
    decimals = decimals,
    scientific_cols = scientific_cols
  )

  style_gt_significance(
    gt_tbl = gt_tbl,
    data_tbl = result_df,
    alpha = alpha,
    sig_fill = sig_fill
  )
}

# One-call table for a ready-made tibble: header + number formats + significance highlighting.
gt_tidy_table <- function(
    data_tbl,
    title,
    subtitle = NULL,
    alpha = 0.05,
    sig_fill = "#e6f4ea",
    decimals = 3L,
    scientific_cols = NULL,
    bold_labels = FALSE,
    width_px = NULL
) {
  gt_tbl <- data_tbl %>%
    gt::gt() %>%
    gt::tab_header(title = title, subtitle = subtitle)

  gt_tbl <- gt_fmt_numbers(
    gt_tbl,
    data_tbl = data_tbl,
    decimals = decimals,
    scientific_cols = scientific_cols
  )

  gt_tbl <- style_gt_significance(
    gt_tbl = gt_tbl,
    data_tbl = data_tbl,
    alpha = alpha,
    sig_fill = sig_fill
  )

  if (isTRUE(bold_labels)) {
    gt_tbl <- gt_bold_column_labels(gt_tbl)
  }
  if (!is.null(width_px)) {
    gt_tbl <- gt_tab_width_px(gt_tbl, width_px = width_px)
  }
  gt_tbl
}

# =============================================================================
# 6. Alpha diversity — beta regression GLMMs (glmmTMB) for stress-history drivers
#    (used by Code/01__Analysis/01__Diversity.R; mirrors Sieler_2026 Code/analysis/02__Microbial_Diversity.Rmd)
# =============================================================================

# Model-ready sample table: clamp normalized Shannon / Simpson / richness to (eps, 1-eps) for beta_family().
# ps_list_obj: named list from preprocessing (e.g. ps.list); element_name: usually "TimeFinal" for Day 60 gut.
build_alpha_data_model <- function(ps_list_obj, element_name = "TimeFinal") {
  if (is.null(ps_list_obj[[element_name]])) {
    stop("ps_list_obj[[\"", element_name, "\"]] is missing. Check ps.list names after 04__DataPreProcess.R.")
  }
  ps_list_obj[[element_name]] %>%
    microViz::samdat_tbl() %>%
    dplyr::mutate(
      Shannon__Genus_norm = dplyr::case_when(
        Shannon__Genus_norm <= 0 ~ 0.0000001,
        Shannon__Genus_norm >= 1 ~ 0.999999,
        TRUE ~ Shannon__Genus_norm
      ),
      Simpson__Genus_norm = dplyr::case_when(
        Simpson__Genus_norm <= 0 ~ 0.0000001,
        Simpson__Genus_norm >= 1 ~ 0.999999,
        TRUE ~ Simpson__Genus_norm
      ),
      Richness__Genus_norm = dplyr::case_when(
        Richness__Genus_norm <= 0 ~ 0.0000001,
        Richness__Genus_norm >= 1 ~ 0.999999,
        TRUE ~ Richness__Genus_norm
      )
    )
}

# glmmTMB beta GLMM coefficient table (conditional component) as a styled gt table.
make_model_coef_gt <- function(fit_obj, title_text, subtitle_text, model_label) {
  coef_table <- summary(fit_obj)$coefficients$cond %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Term") %>%
    dplyr::rename(
      Estimate = "Estimate",
      SE = "Std. Error",
      Statistic = "z value",
      p_value = "Pr(>|z|)"
    ) %>%
    dplyr::mutate(
      Model = model_label,
      Estimate = round(Estimate, 4),
      SE = round(SE, 4),
      Statistic = round(Statistic, 3),
      p_value = round(p_value, 4)
    ) %>%
    dplyr::select(Model, Term, Estimate, SE, Statistic, p_value)

  coef_table %>%
    gt::gt() %>%
    gt::tab_header(
      title = title_text,
      subtitle = subtitle_text
    ) %>%
    gt::cols_label(
      Model = "Model",
      Term = "Term",
      Estimate = "Estimate",
      SE = "Std. Error",
      Statistic = "z",
      p_value = "p-value"
    ) %>%
    gt::tab_style(
      style = list(gt::cell_fill(color = "#e8f5e9")),
      locations = gt::cells_body(
        columns = p_value,
        rows = p_value < 0.05
      )
    )
}

# Highlight p-value / q-value columns in a gt table (legacy 02__Microbial_Diversity.Rmd pattern).
highlight_gt_significance <- function(gt_tbl, table_data, alpha = 0.05) {
  # Backward-compatible wrapper around the canonical `style_gt_significance()`.
  # Keep this name so older modules still work, but centralize behavior in one function.
  style_gt_significance(
    gt_tbl = gt_tbl,
    data_tbl = table_data,
    alpha = alpha,
    sig_fill = "#e8f5e9"
  )
}

# Predicted vs observed tank means across HistoryLevelNum for one normalized diversity column.
# Layer order (bottom → top): grey75 CI ribbon → outlined observed (tank means) → white-then-black
# marginal mean line → outlined predicted EMMs on top so predicted dots sit above the trend line.
# Quasirandom layout: ribbon/anchor → marginal lines → observed beeswarm → predicted EMMs on top.
# `observed_layout`: "point" (default stacked circles) or "quasirandom" (ggbeeswarm vertical spread).
# `subtitle` / `caption`: NULL or "" omits that `labs()` element.
alpha_diversity_stress_history_trend_plot <- function(
    alpha_data_model,
    fit_num,
    response_col,
    y_label,
    title,
    subtitle,
    caption = paste(
      "Colored points: observed tank-level means (normalized 0-1), fill by number of prior stressors.\n",
      "Black line: beta GLMM marginal means on the response scale (white underlay for contrast).\n",
      "Grey ribbon: 95% CI (grey75).",
      sep = ""
    ),
    observed_layout = c("point", "quasirandom")
) {
  observed_layout <- match.arg(observed_layout)

  # type = "response" back-transforms from logit so predicted means/CIs match observed y on (0, 1).
  emm_df <- emmeans::emmeans(
    fit_num,
    ~ HistoryLevelNum,
    at = list(HistoryLevelNum = c(0, 1, 2)),
    type = "response"
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(HistoryLevelNum = as.numeric(as.character(HistoryLevelNum))) %>%
    dplyr::rename(emmean = response)

  subtitle_wrapped <- NULL
  if (!is.null(subtitle) && nzchar(as.character(subtitle)[1])) {
    sub_lines <- strsplit(as.character(subtitle), "\n", fixed = TRUE)[[1]]
    sub_lines <- vapply(sub_lines, function(l) stringr::str_wrap(l, width = 72L), character(1))
    subtitle_wrapped <- paste(sub_lines, collapse = "\n")
  }

  caption_out <- caption
  if (is.null(caption_out) || !nzchar(trimws(as.character(caption_out)[1]))) {
    caption_out <- NULL
  }

  obs_df <- alpha_data_model %>%
    dplyr::group_by(HistoryLevelNum, Tank.ID) %>%
    dplyr::summarise(
      y_obs = mean(.data[[response_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      hist_f = factor(.data$HistoryLevelNum, levels = c(0, 1, 2))
    )

  # Same hex as history_color_scale / prior_stressor_history_colors_numeric (02__PlotSettings.R).
  hist_fill <- prior_stressor_history_colors_numeric

  # Invisible anchor so "Observed" appears in the Trend legend while real observed points use fill.
  legend_anchor_obs <- data.frame(
    HistoryLevelNum = 0.5,
    y_anchor = 0.5,
    dt = factor("Observed", levels = c("Observed", "Predicted"))
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = emm_df,
      ggplot2::aes(x = HistoryLevelNum, ymin = asymp.LCL, ymax = asymp.UCL),
      alpha = 0.42,
      fill = "grey75",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = legend_anchor_obs,
      ggplot2::aes(x = HistoryLevelNum, y = y_anchor, color = dt),
      alpha = 0,
      size = 0.001,
      inherit.aes = FALSE,
      show.legend = TRUE,
      na.rm = TRUE
    )

  if (!identical(observed_layout, "quasirandom")) {
    p <- p +
      ggplot2::geom_point(
        data = obs_df,
        ggplot2::aes(x = HistoryLevelNum, y = y_obs),
        shape = 21L,
        fill = "white",
        colour = "white",
        stroke = 0.5,
        size = 3.35,
        inherit.aes = FALSE,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = obs_df,
        ggplot2::aes(x = HistoryLevelNum, y = y_obs, fill = hist_f),
        shape = 21L,
        colour = "grey30",
        stroke = 0.35,
        size = 2.7,
        alpha = 0.88,
        inherit.aes = FALSE,
        show.legend = TRUE,
        na.rm = TRUE
      )
  }

  p <- p +
    ggplot2::geom_line(
      data = emm_df,
      ggplot2::aes(x = HistoryLevelNum, y = emmean, group = 1L),
      colour = "white",
      linewidth = 2.05,
      lineend = "round",
      inherit.aes = FALSE,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = emm_df,
      ggplot2::aes(x = HistoryLevelNum, y = emmean, group = 1L),
      colour = "black",
      linewidth = 0.95,
      lineend = "round",
      inherit.aes = FALSE,
      show.legend = FALSE,
      na.rm = TRUE
    )

  if (identical(observed_layout, "quasirandom")) {
    if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
      stop("Install package `ggbeeswarm` for observed_layout = \"quasirandom\".", call. = FALSE)
    }
    # width < default 0.4 keeps beeswarms narrow on discrete 0–2 x so clusters do not overlap.
    qr_width <- 0.12
    p <- p +
      ggbeeswarm::geom_quasirandom(
        data = obs_df,
        ggplot2::aes(x = HistoryLevelNum, y = y_obs),
        method = "quasirandom",
        width = qr_width,
        shape = 21L,
        fill = "white",
        colour = "white",
        stroke = 0.55,
        size = 3.75,
        groupOnX = TRUE,
        inherit.aes = FALSE,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      ggbeeswarm::geom_quasirandom(
        data = obs_df,
        ggplot2::aes(x = HistoryLevelNum, y = y_obs, fill = hist_f),
        method = "quasirandom",
        width = qr_width,
        shape = 21L,
        colour = "grey30",
        stroke = 0.35,
        size = 3.05,
        alpha = 0.88,
        groupOnX = TRUE,
        inherit.aes = FALSE,
        show.legend = TRUE,
        na.rm = TRUE
      )
  }

  p <- p +
    ggplot2::geom_point(
      data = emm_df,
      ggplot2::aes(x = HistoryLevelNum, y = emmean),
      shape = 21L,
      fill = "white",
      colour = "white",
      stroke = 0.55,
      size = 5.1,
      inherit.aes = FALSE,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = emm_df,
      ggplot2::aes(x = HistoryLevelNum, y = emmean, color = "Predicted"),
      shape = 19L,
      size = 3.5,
      inherit.aes = FALSE,
      show.legend = TRUE,
      na.rm = TRUE
    )

  p +
    ggplot2::scale_x_continuous(breaks = c(0, 1, 2)) +
    # Fixed (0, 1) y limits match normalized response scale and allow cross-metric comparison.
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_fill_manual(
      name = "Prior Stressor History",
      values = hist_fill,
      breaks = c("0", "1", "2"),
      labels = c("0", "1", "2"),
      guide = ggplot2::guide_legend(
        order = 2L,
        nrow = 1L,
        override.aes = list(
          shape = 21L,
          colour = "grey30",
          stroke = 0.35,
          size = 3.2,
          alpha = 1
        )
      )
    ) +
    ggplot2::scale_color_manual(
      name = "Trend",
      values = c(Observed = "grey45", Predicted = "black"),
      breaks = c("Observed", "Predicted"),
      guide = ggplot2::guide_legend(
        order = 1L,
        nrow = 1L,
        override.aes = list(
          shape = c(16L, 19L),
          size = c(2.8, 3.5),
          alpha = c(1, 1)
        )
      )
    ) +
    ggplot2::labs(
      x = "Number of prior stressors",
      y = y_label,
      title = title,
      subtitle = subtitle_wrapped,
      caption = caption_out
    ) +
    theme_sieler2026_publication(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold"),
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.18, "cm")
    )
}

# Tank-level binomial GLMM (HistoryLevelNum): observed rates vs emmeans on response scale.
# Theme: `theme_sieler2026_publication` + `theme_sieler2026_trend_panel_grid_major_y` (horizontal majors on y).
# Aesthetic match to `alpha_diversity_stress_history_trend_plot` (grey75 ribbon, history fill,
# white/black marginal line, predicted on top); y uses emmeans `prob` and 95% CI bounds.
# `observed_layout`: "point" (default) or "quasirandom" (ggbeeswarm; lines under beeswarm,
# predicted EMM dots on top; width 0.12).
# `y_lim_buffer`: padding below/above the 0-1 probability scale; limits use 100 * buffer when
# `y_as_percent` is TRUE (e.g. 0.05 -> +/-5 on a 0-100 axis).
# `y_as_percent`: multiply observed and emm response/CIs by 100 (axis 0-100%).
# `caption`: pass NULL (default) to omit plot caption.
# `subtitle`: pass NULL or "" to omit subtitle.
glmm_binomial_tank_history_numeric_trend_plot <- function(
    obs_df,
    obs_y_col,
    emm_df,
    y_label,
    title,
    subtitle,
    caption = NULL,
    y_breaks = NULL,
    y_lim_buffer = 0.05,
    y_as_percent = FALSE,
    observed_layout = c("point", "quasirandom")
) {
  observed_layout <- match.arg(observed_layout)
  if (!obs_y_col %in% names(obs_df)) {
    stop("obs_y_col '", obs_y_col, "' not found in obs_df.", call. = FALSE)
  }
  req_emm <- c("HistoryLevelNum", "prob", "asymp.LCL", "asymp.UCL")
  miss_emm <- setdiff(req_emm, names(emm_df))
  if (length(miss_emm) > 0L) {
    stop("emm_df missing columns: ", paste(miss_emm, collapse = ", "), call. = FALSE)
  }

  obs_plot <- obs_df %>%
    dplyr::mutate(
      y_obs = .data[[obs_y_col]],
      hist_f = factor(.data$HistoryLevelNum, levels = c(0, 1, 2))
    )

  emm_plot <- emm_df %>%
    dplyr::mutate(HistoryLevelNum = as.numeric(as.character(.data$HistoryLevelNum)))

  if (isTRUE(y_as_percent)) {
    obs_plot <- obs_plot %>%
      dplyr::mutate(y_obs = 100 * .data$y_obs)
    emm_plot <- emm_plot %>%
      dplyr::mutate(
        prob = 100 * .data$prob,
        asymp.LCL = 100 * .data$asymp.LCL,
        asymp.UCL = 100 * .data$asymp.UCL
      )
  }

  if (is.null(y_breaks)) {
    y_breaks <- if (isTRUE(y_as_percent)) {
      seq(0, 100, 20)
    } else {
      seq(0, 1, 0.2)
    }
  }

  subtitle_wrapped <- NULL
  if (!is.null(subtitle) && nzchar(as.character(subtitle))) {
    sub_lines <- strsplit(as.character(subtitle), "\n", fixed = TRUE)[[1]]
    sub_lines <- vapply(sub_lines, function(l) stringr::str_wrap(l, width = 72L), character(1))
    subtitle_wrapped <- paste(sub_lines, collapse = "\n")
  }

  y_lo <- if (isTRUE(y_as_percent)) {
    0 - 100 * y_lim_buffer
  } else {
    0 - y_lim_buffer
  }
  y_hi <- if (isTRUE(y_as_percent)) {
    100 + 100 * y_lim_buffer
  } else {
    1 + y_lim_buffer
  }

  hist_fill <- prior_stressor_history_colors_numeric

  legend_anchor_obs <- data.frame(
    HistoryLevelNum = 0.5,
    y_anchor = if (isTRUE(y_as_percent)) 50 else 0.5,
    dt = factor("Observed", levels = c("Observed", "Predicted"))
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = emm_plot,
      ggplot2::aes(x = HistoryLevelNum, ymin = asymp.LCL, ymax = asymp.UCL),
      alpha = 0.42,
      fill = "grey75",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = legend_anchor_obs,
      ggplot2::aes(x = HistoryLevelNum, y = y_anchor, color = dt),
      alpha = 0,
      size = 0.001,
      inherit.aes = FALSE,
      show.legend = TRUE,
      na.rm = TRUE
    )

  if (!identical(observed_layout, "quasirandom")) {
    p <- p +
      ggplot2::geom_point(
        data = obs_plot,
        ggplot2::aes(x = HistoryLevelNum, y = y_obs),
        shape = 21L,
        fill = "white",
        colour = "white",
        stroke = 0.5,
        size = 3.35,
        inherit.aes = FALSE,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = obs_plot,
        ggplot2::aes(x = HistoryLevelNum, y = y_obs, fill = hist_f),
        shape = 21L,
        colour = "grey30",
        stroke = 0.35,
        size = 2.7,
        alpha = 0.88,
        inherit.aes = FALSE,
        show.legend = TRUE,
        na.rm = TRUE
      )
  }

  p <- p +
    ggplot2::geom_line(
      data = emm_plot,
      ggplot2::aes(x = HistoryLevelNum, y = prob, group = 1L),
      colour = "white",
      linewidth = 2.05,
      lineend = "round",
      inherit.aes = FALSE,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = emm_plot,
      ggplot2::aes(x = HistoryLevelNum, y = prob, group = 1L),
      colour = "black",
      linewidth = 0.95,
      lineend = "round",
      inherit.aes = FALSE,
      show.legend = FALSE,
      na.rm = TRUE
    )

  if (identical(observed_layout, "quasirandom")) {
    if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
      stop("Install package `ggbeeswarm` for observed_layout = \"quasirandom\".", call. = FALSE)
    }
    qr_width <- 0.12
    p <- p +
      ggbeeswarm::geom_quasirandom(
        data = obs_plot,
        ggplot2::aes(x = HistoryLevelNum, y = y_obs),
        method = "quasirandom",
        width = qr_width,
        shape = 21L,
        fill = "white",
        colour = "white",
        stroke = 0.55,
        size = 3.75,
        groupOnX = TRUE,
        inherit.aes = FALSE,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      ggbeeswarm::geom_quasirandom(
        data = obs_plot,
        ggplot2::aes(x = HistoryLevelNum, y = y_obs, fill = hist_f),
        method = "quasirandom",
        width = qr_width,
        shape = 21L,
        colour = "grey30",
        stroke = 0.35,
        size = 3.05,
        alpha = 0.88,
        groupOnX = TRUE,
        inherit.aes = FALSE,
        show.legend = TRUE,
        na.rm = TRUE
      )
  }

  p <- p +
    ggplot2::geom_point(
      data = emm_plot,
      ggplot2::aes(x = HistoryLevelNum, y = prob),
      shape = 21L,
      fill = "white",
      colour = "white",
      stroke = 0.55,
      size = 5.1,
      inherit.aes = FALSE,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = emm_plot,
      ggplot2::aes(x = HistoryLevelNum, y = prob, color = "Predicted"),
      shape = 19L,
      size = 3.5,
      inherit.aes = FALSE,
      show.legend = TRUE,
      na.rm = TRUE
    )

  p +
    ggplot2::scale_x_continuous(breaks = c(0, 1, 2)) +
    ggplot2::scale_y_continuous(
      limits = c(y_lo, y_hi),
      breaks = y_breaks,
      expand = ggplot2::expansion(mult = c(0, 0)),
      oob = scales::squish
    ) +
    ggplot2::scale_fill_manual(
      name = "Prior Stressor History",
      values = hist_fill,
      breaks = c("0", "1", "2"),
      labels = c("0", "1", "2"),
      guide = ggplot2::guide_legend(
        order = 2L,
        nrow = 1L,
        override.aes = list(
          shape = 21L,
          colour = "grey30",
          stroke = 0.35,
          size = 3.2,
          alpha = 1
        )
      )
    ) +
    ggplot2::scale_color_manual(
      name = "Trend",
      values = c(Observed = "grey45", Predicted = "black"),
      breaks = c("Observed", "Predicted"),
      guide = ggplot2::guide_legend(
        order = 1L,
        nrow = 1L,
        override.aes = list(
          shape = c(16L, 19L),
          size = c(2.8, 3.5),
          alpha = c(1, 1)
        )
      )
    ) +
    ggplot2::labs(
      x = "Number of prior stressors",
      y = y_label,
      title = title,
      subtitle = subtitle_wrapped,
      caption = caption
    ) +
    theme_sieler2026_publication(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold"),
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.18, "cm")
    ) +
    theme_sieler2026_trend_panel_grid_major_y()
}

# Linear contrast (-1, 0, 1) on HistoryLevelNum EMMs; returns emmeans contrast object.
alpha_diversity_linear_trend_contrast <- function(fit_num) {
  emm_num <- emmeans::emmeans(
    fit_num,
    ~ HistoryLevelNum,
    at = list(HistoryLevelNum = c(0, 1, 2))
  )
  linear_contrast <- c(-1, 0, 1)
  emmeans::contrast(emm_num, list(linear = linear_contrast))
}

# One row: conditional (link-scale) fixed effect for HistoryLevelNum:Parasite (glmmTMB).
glmmtmb_history_parasite_interaction_coef_row <- function(fit_interaction) {
  co <- summary(fit_interaction)$coefficients$cond
  rn <- rownames(co)
  hit <- rn[
    stringr::str_detect(rn, "HistoryLevelNum") &
      stringr::str_detect(rn, "Parasite") &
      stringr::str_detect(rn, ":")
  ]
  if (length(hit) < 1L) {
    stop("No HistoryLevelNum:Parasite interaction row found in glmmTMB conditional coefficients.")
  }
  h <- hit[1]
  tibble::tibble(
    Term = h,
    Estimate = unname(co[h, "Estimate"]),
    SE = unname(co[h, "Std. Error"]),
    z = unname(co[h, "z value"]),
    p.value = unname(co[h, "Pr(>|z|)"])
  )
}

# Wald Anova (type II) p-value for the HistoryLevelNum:Parasite interaction in the interaction model.
glmmtmb_history_parasite_interaction_p <- function(fit_interaction) {
  av <- tryCatch(
    car::Anova(fit_interaction, type = "II"),
    error = function(e) NULL
  )
  if (!is.null(av)) {
    rn <- rownames(av)
    hit <- rn[stringr::str_detect(rn, "HistoryLevelNum") & stringr::str_detect(rn, "Parasite") & stringr::str_detect(rn, ":")]
    if (length(hit) == 1L) {
      pcol <- if ("Pr(>Chisq)" %in% colnames(av)) "Pr(>Chisq)" else if ("Pr(>|z|)" %in% colnames(av)) "Pr(>|z|)" else colnames(av)[ncol(av)]
      return(as.numeric(av[hit, pcol]))
    }
  }
  co <- summary(fit_interaction)$coefficients$cond
  rn2 <- rownames(co)
  hit2 <- rn2[stringr::str_detect(rn2, "HistoryLevelNum") & stringr::str_detect(rn2, "Parasite") & stringr::str_detect(rn2, ":")]
  if (length(hit2) >= 1L) {
    return(as.numeric(co[hit2[1], "Pr(>|z|)"]))
  }
  NA_real_
}

# Simple effects of Parasite (1 vs 0) at each HistoryLevelNum from the interaction beta-GLMM.
# Uses emmeans on the response scale (type = "response"); rows are pairwise Parasite contrasts within stratum.
# p_fdr_bh: Benjamini–Hochberg adjustment across the three strata (per metric / model).
alpha_diversity_parasite_simple_effects <- function(fit_interaction, metric_name = "Metric") {
  set.seed(42)
  rg <- emmeans::ref_grid(
    fit_interaction,
    at = list(HistoryLevelNum = c(0, 1, 2), Parasite = c(0, 1))
  )
  emm <- emmeans::emmeans(rg, ~ Parasite | HistoryLevelNum, type = "response")
  # Dispatch to pairs.emmGrid (emmeans); not exported as emmeans::pairs.
  pw <- pairs(emm, reverse = TRUE)
  summ <- summary(pw, infer = c(TRUE, TRUE))
  out <- summ %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      Metric = metric_name,
      HistoryLevelNum = as.numeric(as.character(.data$HistoryLevelNum)),
      p_fdr_bh = stats::p.adjust(.data$p.value, method = "BH")
    )
  # Beta GLMM pairs on response scale may label the effect `odds.ratio` (emmeans naming).
  if ("odds.ratio" %in% names(out) && !"estimate" %in% names(out)) {
    out <- out %>% dplyr::rename(estimate = odds.ratio)
  }
  ratio_col <- intersect(c("z.ratio", "t.ratio"), names(out))
  if (length(ratio_col) == 1L) {
    out <- out %>% dplyr::rename(z_or_t_ratio = dplyr::all_of(ratio_col))
  } else {
    out$z_or_t_ratio <- NA_real_
  }
  out %>%
    dplyr::select(
      "Metric",
      "HistoryLevelNum",
      "contrast",
      "estimate",
      "SE",
      "df",
      "z_or_t_ratio",
      "p.value",
      "p_fdr_bh",
      dplyr::any_of(c("lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL"))
    )
}

# Violin + jitter of per-sample normalized diversity by Parasite, faceted by prior stressor count.
alpha_diversity_parasite_within_history_plot <- function(
    alpha_data_model,
    response_col,
    y_label,
    title,
    subtitle,
    p_value_tbl = NULL
) {
  plot_df <- alpha_data_model %>%
    dplyr::mutate(
      Parasite_label = dplyr::if_else(.data$Parasite == 1L, "Exposed", "Unexposed"),
      Parasite_label = factor(.data$Parasite_label, levels = c("Unexposed", "Exposed")),
      History_facet = factor(
        .data$HistoryLevelNum,
        levels = c(0, 1, 2),
        labels = c("0 prior stressors", "1 prior stressor", "2 prior stressors")
      )
    )
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$Parasite_label,
      y = .data[[response_col]],
      fill = .data$Parasite_label
    )
  ) +
    ggplot2::geom_violin(alpha = 0.35, color = NA) +
    ggplot2::geom_jitter(
      width = 0.12,
      height = 0,
      alpha = 0.35,
      size = 1.4,
      color = "gray25"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(History_facet), nrow = 1L) +
    ggplot2::scale_fill_manual(values = c(Unexposed = "gray70", Exposed = "#E31A1C")) +
    ggplot2::labs(
      x = "Parasite Exposure",
      y = y_label,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    theme_sieler2026_publication(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_text(face = "bold")
    )

  # Optional: add per-facet comparison bars (Unexposed vs Exposed) using precomputed p-values.
  # p_value_tbl should have HistoryLevelNum and p.value columns (from emmeans pairs).
  if (!is.null(p_value_tbl) && nrow(p_value_tbl) > 0L) {
    y_by_facet <- plot_df %>%
      dplyr::group_by(.data$History_facet) %>%
      dplyr::summarise(y_max = max(.data[[response_col]], na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(y = pmin(0.98, .data$y_max + 0.06))

    p_annot <- p_value_tbl %>%
      dplyr::distinct(.data$HistoryLevelNum, .keep_all = TRUE) %>%
      dplyr::mutate(
        History_facet = factor(
          .data$HistoryLevelNum,
          levels = c(0, 1, 2),
          labels = c("0 prior stressors", "1 prior stressor", "2 prior stressors")
        ),
        label = dplyr::case_when(
          .data$p.value < 0.001 ~ "***",
          .data$p.value < 0.01 ~ "**",
          .data$p.value < 0.05 ~ "*",
          .data$p.value < 0.1 ~ ".",
          TRUE ~ "ns"
        ),
        x1 = 1L,
        x2 = 2L
      ) %>%
      dplyr::left_join(y_by_facet, by = "History_facet") %>%
      dplyr::filter(!is.na(.data$y)) %>%
      dplyr::select("History_facet", "label", "y", "x1", "x2")

    if (nrow(p_annot) > 0L) {
      p <- p +
        ggplot2::geom_segment(
          data = p_annot,
          ggplot2::aes(x = .data$x1, xend = .data$x2, y = .data$y, yend = .data$y),
          inherit.aes = FALSE,
          linewidth = SIELER2026_MIN_LINEWIDTH_MM
        ) +
        ggplot2::geom_segment(
          data = p_annot,
          ggplot2::aes(x = .data$x1, xend = .data$x1, y = .data$y, yend = .data$y - 0.02),
          inherit.aes = FALSE,
          linewidth = SIELER2026_MIN_LINEWIDTH_MM
        ) +
        ggplot2::geom_segment(
          data = p_annot,
          ggplot2::aes(x = .data$x2, xend = .data$x2, y = .data$y, yend = .data$y - 0.02),
          inherit.aes = FALSE,
          linewidth = SIELER2026_MIN_LINEWIDTH_MM
        ) +
        ggplot2::geom_text(
          data = p_annot,
          ggplot2::aes(x = 1.5, y = .data$y + 0.02, label = .data$label),
          inherit.aes = FALSE,
          size = 5
        )
    }
  }

  p
}


# =============================================================================
# 7) Composition / beta diversity (PERMANOVA, betadisper plots) — 02__Composition.R
# =============================================================================

# Eight exposure-regime labels (A± × T± × P±); must match 04__DataPreProcess.R / 02__PlotSettings.R.
exposure_regime_levels <- function() {
  c(
    "A- T- P-", "A- T- P+", "A+ T- P-", "A+ T- P+",
    "A- T+ P-", "A- T+ P+", "A+ T+ P-", "A+ T+ P+"
  )
}

# Colors for exposure regimes (must match `treatment_color_scale` in 02__PlotSettings.R).
exposure_regime_colors <- function() {
  lv <- exposure_regime_levels()
  cols <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666"
  )
  stats::setNames(cols, lv)
}

# microViz dist_permanova object -> tidy tibble (one row per model term).
microviz_permanova_to_tidy <- function(perma_obj, distance_label) {
  microViz::perm_get(perma_obj) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Term") %>%
    dplyr::filter(.data$Term != "Total") %>%
    dplyr::mutate(
      Distance = distance_label,
      R2 = round(.data$R2, 4),
      F = round(.data$F, 3),
      p_value = round(.data[["Pr(>F)"]], 4)
    ) %>%
    dplyr::select("Distance", "Term", "Df", "R2", "F", "p_value")
}

# Global PCoA (cmdscale on phyloseq distance) with facets by prior stressor history; points colored by parasite exposure.
composition_pcoa_parasite_faceted_by_history_plot <- function(
    ps_in,
    dist_method = "bray",
    title = "PCoA — parasite exposure within prior stressor history",
    subtitle = "Shared ordination across all samples; panels stratify prior stressor count."
) {
  ps_g <- ps_in %>%
    microViz::tax_transform("identity", rank = "Genus")
  dm <- phyloseq::distance(ps_g, method = dist_method)
  sn <- phyloseq::sample_names(ps_g)
  dmat <- as.matrix(dm)
  cds <- stats::cmdscale(dmat, k = 2L, eig = TRUE)
  pts <- cds$points
  if (is.null(pts) || ncol(pts) < 2L) {
    stop("PCoA: need at least two dimensions (check sample count).")
  }
  df_plot <- microViz::samdat_tbl(ps_g)
  if (!"Sample" %in% names(df_plot)) {
    df_plot$Sample <- sn
  }
  # Align rows to phyloseq sample order (cmdscale rows follow sample_names(ps_g)).
  m <- match(df_plot$Sample, sn)
  if (anyNA(m) || length(unique(stats::na.omit(m))) != nrow(df_plot)) {
    if (nrow(df_plot) == length(sn)) {
      df_plot$PC1 <- pts[, 1L]
      df_plot$PC2 <- pts[, 2L]
    } else {
      stop("composition_pcoa_parasite_faceted_by_history_plot: cannot align samples to PCoA scores.")
    }
  } else {
    pts_m <- pts[m, , drop = FALSE]
    df_plot$PC1 <- pts_m[, 1L]
    df_plot$PC2 <- pts_m[, 2L]
  }
  if (!("History_Label" %in% names(df_plot)) && "HistoryLevelNum" %in% names(df_plot)) {
    df_plot$History_Label <- factor(
      df_plot$HistoryLevelNum,
      levels = c(0, 1, 2),
      labels = c("No prior stressors", "One prior stressor", "Two prior stressors")
    )
  }
  if (!"Parasite_Exposed" %in% names(df_plot)) {
    df_plot$Parasite_Exposed <- dplyr::if_else(df_plot$Parasite == 1L, "Exposed", "Unexposed")
    df_plot$Parasite_Exposed <- factor(df_plot$Parasite_Exposed, levels = c("Unexposed", "Exposed"))
  }
  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = .data$PC1,
      y = .data$PC2,
      color = .data$Parasite_Exposed,
      fill = .data$Parasite_Exposed
    )
  ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = .data$Parasite_Exposed, color = .data$Parasite_Exposed),
      linewidth = 0.8,
      alpha = 0.4,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(shape = 21L, stroke = 0.4, size = 2.2, alpha = 0.85) +
    ggplot2::scale_color_manual(
      values = c(Unexposed = "gray35", Exposed = "#E31A1C"),
      name = "Parasite exposure"
    ) +
    ggplot2::scale_fill_manual(
      values = c(Unexposed = "white", Exposed = "#fccaca"),
      name = "Parasite exposure"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(History_Label), nrow = 1L) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "PC1",
      y = "PC2"
    ) +
    theme_sieler2026_composition_figure(base_size = 14) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
}

# Same global PCoA as `composition_pcoa_parasite_faceted_by_history_plot`, but points colored by
# eight-level exposure regime (`Treatment`) while ellipses stay on Parasite_Exposed (ggnewscale).
composition_pcoa_parasite_faceted_points_by_regime_plot <- function(
    ps_in,
    dist_method = "bray",
    title = "PCoA — exposure regime within prior stressor history",
    subtitle = "Points: exposure regime; ellipses: parasite exposure (shared cmdscale across facets)."
) {
  ps_g <- ps_in %>%
    microViz::tax_transform("identity", rank = "Genus")
  dm <- phyloseq::distance(ps_g, method = dist_method)
  sn <- phyloseq::sample_names(ps_g)
  dmat <- as.matrix(dm)
  cds <- stats::cmdscale(dmat, k = 2L, eig = TRUE)
  pts <- cds$points
  if (is.null(pts) || ncol(pts) < 2L) {
    stop("PCoA: need at least two dimensions (check sample count).")
  }
  df_plot <- microViz::samdat_tbl(ps_g)
  if (!"Sample" %in% names(df_plot)) {
    df_plot$Sample <- sn
  }
  m <- match(df_plot$Sample, sn)
  if (anyNA(m) || length(unique(stats::na.omit(m))) != nrow(df_plot)) {
    if (nrow(df_plot) == length(sn)) {
      df_plot$PC1 <- pts[, 1L]
      df_plot$PC2 <- pts[, 2L]
    } else {
      stop("composition_pcoa_parasite_faceted_points_by_regime_plot: cannot align samples to PCoA scores.")
    }
  } else {
    pts_m <- pts[m, , drop = FALSE]
    df_plot$PC1 <- pts_m[, 1L]
    df_plot$PC2 <- pts_m[, 2L]
  }
  if (!("History_Label" %in% names(df_plot)) && "HistoryLevelNum" %in% names(df_plot)) {
    df_plot$History_Label <- factor(
      df_plot$HistoryLevelNum,
      levels = c(0, 1, 2),
      labels = c("No prior stressors", "One prior stressor", "Two prior stressors")
    )
  }
  if (!"Parasite_Exposed" %in% names(df_plot)) {
    df_plot$Parasite_Exposed <- dplyr::if_else(df_plot$Parasite == 1L, "Exposed", "Unexposed")
    df_plot$Parasite_Exposed <- factor(df_plot$Parasite_Exposed, levels = c("Unexposed", "Exposed"))
  }
  if (!"Treatment" %in% names(df_plot)) {
    stop("composition_pcoa_parasite_faceted_points_by_regime_plot: sample_data needs a Treatment column.")
  }
  df_plot$Treatment <- factor(as.character(df_plot$Treatment), levels = exposure_regime_levels())
  regime_cols <- exposure_regime_colors()

  ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
    ggplot2::stat_ellipse(
      ggplot2::aes(color = .data$Parasite_Exposed, group = .data$Parasite_Exposed),
      linewidth = 0.8,
      alpha = 0.35,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c(Unexposed = "gray35", Exposed = "#E31A1C"),
      name = "Parasite exposure (ellipses)"
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$Treatment),
      shape = 19L,
      size = 2.2,
      alpha = 0.88
    ) +
    ggplot2::scale_color_manual(
      values = regime_cols[exposure_regime_levels()],
      breaks = exposure_regime_levels(),
      drop = FALSE,
      name = "Exposure regime"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(History_Label), nrow = 1L) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "PC1",
      y = "PC2"
    ) +
    theme_sieler2026_composition_figure(base_size = 14) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing = grid::unit(4, "pt")
    )
}

# GT table for inferential stats (PERMANOVA / ANOVA); styles small p-values.
composition_gt_format_p_chr <- function(x, threshold = 1e-4, digits = 4) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x < threshold ~ "<0.0001",
    TRUE ~ sprintf(paste0("%.", digits, "f"), x)
  )
}

composition_gt_style_significance_cells <- function(gt_tbl, alpha = 0.05, fill_color = "#d9f2d9") {
  # Prefer the canonical helper so styling stays consistent across modules.
  data_tbl <- gt_tbl[["_data"]]
  if (is.null(data_tbl) || nrow(data_tbl) == 0L) {
    return(gt_tbl)
  }
  style_gt_significance(
    gt_tbl = gt_tbl,
    data_tbl = data_tbl,
    alpha = alpha,
    sig_fill = fill_color
  )
}

composition_gt_format_pvalue_cols <- function(gt_tbl, threshold = 1e-4) {
  data_tbl <- gt_tbl[["_data"]]
  if (is.null(data_tbl) || nrow(data_tbl) == 0L) {
    return(gt_tbl)
  }
  sig_cols <- names(data_tbl)[purrr::map_lgl(names(data_tbl), function(col_name) {
    cl <- stringr::str_to_lower(col_name)
    stringr::str_detect(cl, "^p_value$|^p$|^p[._ -]?value$|^p-?adj$|^padj$|^q[._ -]?value$|fdr|^pr\\(>.*\\)$|adjusted[._ -]?p")
  })]
  for (col_name in sig_cols) {
    col_values <- suppressWarnings(as.numeric(data_tbl[[col_name]]))
    if (all(is.na(col_values))) {
      next
    }
    gt_tbl <- gt::fmt(
      gt_tbl,
      columns = dplyr::all_of(col_name),
      fns = function(x) composition_gt_format_p_chr(as.numeric(x), threshold = threshold)
    )
  }
  gt_tbl
}

composition_inferential_gt <- function(table_data, title, subtitle = NULL, alpha = 0.05) {
  gt_tbl <- table_data %>%
    gt::gt() %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    )
  composition_gt_style_significance_cells(gt_tbl, alpha = alpha) %>%
    composition_gt_format_pvalue_cols()
}

# Tukey table for betadisper post-hoc (significant rows only).
composition_build_tukey_sig_table <- function(tukey_obj, alpha = 0.05) {
  if (is.null(tukey_obj)) {
    return(tibble::tibble())
  }
  tukey_tbl <- as.data.frame(tukey_obj)
  if (!"p adj" %in% names(tukey_tbl)) {
    return(tibble::tibble())
  }
  tukey_tbl %>%
    tibble::rownames_to_column(var = "comparison") %>%
    dplyr::mutate(
      sig_label = dplyr::case_when(
        .data[["p adj"]] < 0.001 ~ "***",
        .data[["p adj"]] < 0.01 ~ "**",
        .data[["p adj"]] < 0.05 ~ "*",
        .data[["p adj"]] < 0.1 ~ ".",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::filter(.data[["p adj"]] < alpha)
}

composition_match_tukey_comparison <- function(comparison, level_values) {
  for (lhs in level_values) {
    for (rhs in level_values) {
      if (identical(lhs, rhs)) {
        next
      }
      if (comparison == paste0(lhs, "-", rhs)) {
        return(c(lhs, rhs))
      }
    }
  }
  NULL
}

# Add pairwise significance brackets above dispersion plots (betadisper violins / boxplots).
composition_add_significance_bars <- function(
    plot_obj,
    data_tbl,
    x_var,
    y_var,
    sig_tbl,
    level_values,
    level_map = NULL,
    bar_linewidth = 0.88,
    tip_length_factor = 0.22,
    stack_factor = 1.28,
    label_above_factor = 0.62,
    y_anchor_offset = 0.09,
    y_bar_per_row = NULL,
    y_label_per_row = NULL,
    y_hi_cap = NULL
) {
  if (nrow(sig_tbl) == 0L) {
    return(plot_obj)
  }
  y_max <- max(data_tbl[[y_var]], na.rm = TRUE)
  y_min <- min(data_tbl[[y_var]], na.rm = TRUE)
  y_range <- y_max - y_min
  y_step <- y_range * 0.05
  if (!is.finite(y_step) || y_step <= 0) {
    y_step <- 0.05
  }
  y_anchor <- y_max + y_range * y_anchor_offset
  tip_length <- y_step * tip_length_factor
  n_sig <- nrow(sig_tbl)
  use_fixed_y <- is.numeric(y_bar_per_row) &&
    length(y_bar_per_row) == n_sig &&
    all(is.finite(y_bar_per_row))
  use_fixed_lbl <- is.numeric(y_label_per_row) &&
    length(y_label_per_row) == n_sig &&
    all(is.finite(y_label_per_row))
  for (i in seq_len(n_sig)) {
    row <- sig_tbl[i, ]
    pair <- composition_match_tukey_comparison(row$comparison, level_values)
    if (is.null(pair)) {
      next
    }
    x1 <- pair[1]
    x2 <- pair[2]
    if (!is.null(level_map)) {
      x1 <- level_map[[x1]]
      x2 <- level_map[[x2]]
    }
    x_pos_1 <- match(x1, levels(data_tbl[[x_var]]))
    x_pos_2 <- match(x2, levels(data_tbl[[x_var]]))
    if (is.na(x_pos_1) || is.na(x_pos_2)) {
      next
    }
    y_pos <- if (isTRUE(use_fixed_y)) {
      y_bar_per_row[[i]]
    } else {
      y_anchor + (i * y_step * stack_factor)
    }
    y_txt <- if (isTRUE(use_fixed_lbl)) {
      y_label_per_row[[i]]
    } else {
      y_pos + (y_step * label_above_factor)
    }
    plot_obj <- plot_obj +
      ggplot2::annotate(
        "segment",
        x = x_pos_1,
        xend = x_pos_2,
        y = y_pos,
        yend = y_pos,
        color = "black",
        linewidth = bar_linewidth
      ) +
      ggplot2::annotate(
        "segment",
        x = x_pos_1,
        xend = x_pos_1,
        y = y_pos,
        yend = y_pos - tip_length,
        color = "black",
        linewidth = bar_linewidth
      ) +
      ggplot2::annotate(
        "segment",
        x = x_pos_2,
        xend = x_pos_2,
        y = y_pos,
        yend = y_pos - tip_length,
        color = "black",
        linewidth = bar_linewidth
      ) +
      ggplot2::annotate(
        "text",
        x = mean(c(x_pos_1, x_pos_2)),
        y = y_txt,
        label = row$sig_label,
        size = 4,
        color = "black"
      )
  }
  y_hi <- if (isTRUE(use_fixed_lbl)) {
    max(y_label_per_row, na.rm = TRUE) + y_step * 0.4
  } else {
    y_anchor + (n_sig * y_step * stack_factor) + (y_step * (label_above_factor + 0.35))
  }
  if (!is.null(y_hi_cap) && is.finite(y_hi_cap)) {
    y_hi <- min(y_hi, y_hi_cap)
  }
  plot_obj + ggplot2::expand_limits(y = y_hi)
}

# Boxplot of distance-to-centroid for HistoryLevel groups (betadisper); uses history_color_scale from 02__PlotSettings.R.
# Multi-layer violins, double quasirandom, white + history median via draw_quantiles; opaque history fills (Michael Sieler; 2026-04-24).
composition_betadisper_boxplot_history <- function(
    bdisp_obj,
    history_color_scale_vec,
    title_suffix
) {
  disp_data <- data.frame(
    Distance = bdisp_obj$HistoryLevel$model$distances,
    HistoryLevel = bdisp_obj$HistoryLevel$model$group
  ) %>%
    dplyr::mutate(
      History_Label = dplyr::case_when(
        .data$HistoryLevel == "0" ~ "No prior stressors",
        .data$HistoryLevel == "1" ~ "One prior stressor",
        .data$HistoryLevel == "2" ~ "Two prior stressors",
        TRUE ~ "Unknown"
      ),
      History_Label = factor(
        .data$History_Label,
        levels = c("No prior stressors", "One prior stressor", "Two prior stressors")
      )
    )
  tukey_sig <- composition_build_tukey_sig_table(bdisp_obj$HistoryLevel$tukeyHSD$group)
  p0 <- disp_data %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$History_Label, y = .data$Distance)) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = .data$History_Label, group = .data$History_Label),
      alpha = 0.18,
      colour = NA,
      linewidth = 0,
      trim = TRUE,
      width = 0.92,
      scale = "width",
      show.legend = FALSE
    ) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = .data$History_Label, group = .data$History_Label),
      colour = "white",
      linewidth = 1.18,
      trim = TRUE,
      width = 0.88,
      scale = "width",
      draw_quantiles = c(0.5),
      show.legend = FALSE
    ) +
    ggplot2::geom_violin(
      ggplot2::aes(
        colour = .data$History_Label,
        fill = .data$History_Label,
        group = .data$History_Label
      ),
      linewidth = 0.52,
      trim = TRUE,
      width = 0.76,
      scale = "width",
      draw_quantiles = c(0.5),
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(fill = .data$History_Label),
      colour = "white",
      stroke = 2.25,
      shape = 23L,
      groupOnX = TRUE,
      size = 1.26
    ) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(colour = .data$History_Label, fill = .data$History_Label),
      stroke = 1,
      shape = 23L,
      groupOnX = TRUE,
      size = 1.1
    ) +
    ggplot2::scale_fill_manual(
      values = history_color_scale_vec,
      name = "Stressor History",
      drop = FALSE
    ) +
    ggplot2::scale_colour_manual(values = history_color_scale_vec, drop = FALSE, guide = "none") +
    ggplot2::labs(
      x = "Stressor History",
      y = "Distance to centroid",
      title = paste0("Beta diversity dispersion by stressor history (", title_suffix, ")")
    ) +
    theme_sieler2026_composition_figure(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      legend.direction = "horizontal",
      legend.box = "horizontal",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  composition_add_significance_bars(
    plot_obj = p0,
    data_tbl = disp_data,
    x_var = "History_Label",
    y_var = "Distance",
    sig_tbl = tukey_sig,
    level_values = c("0", "1", "2"),
    level_map = c(
      "0" = "No prior stressors",
      "1" = "One prior stressor",
      "2" = "Two prior stressors"
    ),
    stack_factor = 1.52,
    label_above_factor = 0.9
  ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 36, r = 20, b = 88, l = 20, unit = "pt")
    )
}

# Betadisper boxplot for HistoryLevel × Parasite groups (six levels: 0_0 … 2_1); Used in 02__Composition.R for Bray-Curtis and Canberra dispersion figures.
# Multi-layer violins, double quasirandom, medians via draw_quantiles on outline layers; opaque white / grey75 fill; point fill RoL-style (Michael Sieler; 2026-04-24).
composition_betadisper_boxplot_history_parasite <- function(
    bdisp_obj,
    history_color_scale_vec,
    title_suffix
) {
  disp_data <- data.frame(
    Distance = bdisp_obj$History.Parasite$model$distances,
    Group = bdisp_obj$History.Parasite$model$group
  ) %>%
    dplyr::mutate(
      HistoryLevel = stringr::str_extract(.data$Group, "^[0-9]"),
      Parasite = stringr::str_extract(.data$Group, "[0-9]$"),
      History_Label = dplyr::case_when(
        .data$HistoryLevel == "0" ~ "No prior stressors",
        .data$HistoryLevel == "1" ~ "One prior stressor",
        .data$HistoryLevel == "2" ~ "Two prior stressors",
        TRUE ~ "Unknown"
      ),
      Parasite_Exposed = dplyr::case_when(
        .data$Parasite == "0" ~ "Unexposed",
        .data$Parasite == "1" ~ "Exposed",
        TRUE ~ "Unknown"
      ),
      History_Label = factor(
        .data$History_Label,
        levels = c("No prior stressors", "One prior stressor", "Two prior stressors")
      ),
      Parasite_Exposed = factor(.data$Parasite_Exposed, levels = c("Unexposed", "Exposed")),
      fill_cond = dplyr::if_else(
        .data$Parasite_Exposed == "Unexposed",
        "white",
        "grey75",
        missing = "white"
      ),
      # RoL analogue: Control=white fill, "Treatment" arm = hue fill — here Unexposed=white, Exposed=history hue.
      point_fill = dplyr::if_else(
        .data$Parasite_Exposed == "Unexposed",
        "white",
        as.character(.data$History_Label),
        missing = "white"
      )
    )

  tukey_within_history <- as.data.frame(bdisp_obj$History.Parasite$tukeyHSD$group) %>%
    tibble::rownames_to_column(var = "comparison") %>%
    dplyr::filter(
      stringr::str_detect(.data$comparison, "^0_1-0_0$|^0_0-0_1$|^1_1-1_0$|^1_0-1_1$|^2_1-2_0$|^2_0-2_1$")
    ) %>%
    dplyr::mutate(
      History_Label = dplyr::case_when(
        stringr::str_detect(.data$comparison, "^0_") ~ "No prior stressors",
        stringr::str_detect(.data$comparison, "^1_") ~ "One prior stressor",
        stringr::str_detect(.data$comparison, "^2_") ~ "Two prior stressors",
        TRUE ~ "Unknown"
      ),
      History_Label = factor(
        .data$History_Label,
        levels = c("No prior stressors", "One prior stressor", "Two prior stressors")
      ),
      sig_label = dplyr::case_when(
        .data[["p adj"]] < 0.001 ~ "***",
        .data[["p adj"]] < 0.01 ~ "**",
        .data[["p adj"]] < 0.05 ~ "*",
        .data[["p adj"]] < 0.1 ~ ".",
        TRUE ~ "ns"
      )
    )

  history_aov <- stats::aov(Distance ~ History_Label, data = disp_data)
  history_tukey_raw <- stats::TukeyHSD(history_aov, "History_Label")
  history_tukey <- composition_build_tukey_sig_table(history_tukey_raw[["History_Label"]]) %>%
    dplyr::arrange(
      match(
        .data$sig_label,
        c("***", "**", "*", ".", "ns"),
        nomatch = 99L
      ),
      .data$comparison
    )

  y_max <- max(disp_data$Distance, na.rm = TRUE)
  y_range <- max(disp_data$Distance, na.rm = TRUE) - min(disp_data$Distance, na.rm = TRUE)
  y_step <- y_range * 0.05
  if (!is.finite(y_step) || y_step <= 0) {
    y_step <- 0.05
  }

  dodge_w <- 0.98
  pos_d <- ggplot2::position_dodge(width = dodge_w)
  grey75_hex <- grDevices::grey(0.75)
  fill_parasite <- c(white = "#FFFFFF", grey75 = grey75_hex)

  p0 <- disp_data %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$History_Label,
        y = .data$Distance,
        colour = .data$History_Label,
        shape = .data$Parasite_Exposed,
        group = interaction(.data$History_Label, .data$Parasite_Exposed)
      )
    ) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = .data$History_Label),
      position = pos_d,
      alpha = 0.18,
      colour = NA,
      linewidth = 0,
      trim = TRUE,
      scale = "width",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = history_color_scale_vec, guide = "none") +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_violin(
      ggplot2::aes(fill = .data$fill_cond),
      position = pos_d,
      colour = NA,
      linewidth = 0,
      trim = TRUE,
      scale = "width",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fill_parasite, guide = "none") +
    ggplot2::geom_violin(
      ggplot2::aes(
        x = .data$History_Label,
        y = .data$Distance,
        group = interaction(.data$History_Label, .data$Parasite_Exposed)
      ),
      position = pos_d,
      inherit.aes = FALSE,
      fill = NA,
      colour = "white",
      linewidth = 2.24,
      trim = TRUE,
      scale = "width",
      draw_quantiles = c(0.5),
      show.legend = FALSE
    ) +
    ggplot2::geom_violin(
      ggplot2::aes(
        x = .data$History_Label,
        y = .data$Distance,
        colour = .data$History_Label,
        group = interaction(.data$History_Label, .data$Parasite_Exposed)
      ),
      position = pos_d,
      inherit.aes = FALSE,
      fill = NA,
      linewidth = 1.04,
      trim = TRUE,
      scale = "width",
      draw_quantiles = c(0.5),
      show.legend = FALSE
    ) +
    ggnewscale::new_scale_fill() +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(fill = .data$point_fill),
      colour = "white",
      stroke = 2.35,
      dodge.width = dodge_w,
      size = 0.97
    ) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(colour = .data$History_Label, fill = .data$point_fill),
      stroke = 1,
      dodge.width = dodge_w,
      size = 0.83
    ) +
    ggplot2::scale_fill_manual(
      values = c(history_color_scale_vec, white = "#FFFFFF"),
      guide = "none"
    ) +
    ggplot2::scale_colour_manual(
      values = history_color_scale_vec,
      name = "Prior Stressor History",
      labels = c(
        "No prior stressors" = "None",
        "One prior stressor" = "One",
        "Two prior stressors" = "Two"
      ),
      guide = ggplot2::guide_legend(
        order = 1L,
        # Filled keys (match alpha stress-history trend): fill = history hue per level; rim like Shannon legend.
        override.aes = list(
          fill = c(
            history_color_scale_vec[["No prior stressors"]],
            history_color_scale_vec[["One prior stressor"]],
            history_color_scale_vec[["Two prior stressors"]]
          ),
          shape = 21L,
          colour = "grey30",
          stroke = 0.35,
          size = 3.2
        )
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c("Unexposed" = 21L, "Exposed" = 23L),
      name = "Parasite Exposure",
      labels = c("Unexposed" = "Unexposed", "Exposed" = "Exposed"),
      guide = ggplot2::guide_legend(
        order = 2L,
        override.aes = list(
          colour = "black",
          stroke = 0.5,
          size = 3,
          fill = c("#FFFFFF", grey75_hex)
        )
      )
    ) +
    ggplot2::scale_x_discrete(
      labels = c(
        "No prior stressors" = "None",
        "One prior stressor" = "One",
        "Two prior stressors" = "Two"
      )
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.15),
      breaks = seq(0, 1, 0.25)
    ) +
    ggplot2::labs(
      x = "Stressor History",
      y = "Distance to centroid",
      title = stringr::str_wrap(
        paste0(
          "Beta diversity dispersion by stressor history and parasite exposure (",
          title_suffix,
          ")"
        ),
        width = 68
      )
    ) +
    theme_sieler2026_composition_figure(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(6, "pt"),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  tukey_ordered <- tukey_within_history %>%
    dplyr::mutate(.is_ns = .data$sig_label == "ns") %>%
    dplyr::arrange(.data$.is_ns, .data$History_Label) %>%
    dplyr::mutate(.is_ns = NULL)

  n_w <- nrow(tukey_ordered)
  y_bar_within <- 0.9
  y_lbl_within <- 0.95
  if (n_w > 0L) {
    tip_length <- y_range * 0.01
    for (i in seq_len(n_w)) {
      row <- tukey_ordered[i, ]
      x_pos <- as.numeric(row$History_Label)
      x_min <- x_pos - 0.2
      x_max <- x_pos + 0.2
      p0 <- p0 +
        ggplot2::annotate(
          "segment",
          x = x_min,
          xend = x_max,
          y = y_bar_within,
          yend = y_bar_within,
          color = "black",
          linewidth = 0.88
        ) +
        ggplot2::annotate(
          "segment",
          x = x_min,
          xend = x_min,
          y = y_bar_within,
          yend = y_bar_within - tip_length,
          color = "black",
          linewidth = 0.88
        ) +
        ggplot2::annotate(
          "segment",
          x = x_max,
          xend = x_max,
          y = y_bar_within,
          yend = y_bar_within - tip_length,
          color = "black",
          linewidth = 0.88
        ) +
        ggplot2::annotate(
          "text",
          x = (x_min + x_max) / 2,
          y = y_lbl_within,
          label = row$sig_label,
          size = 4,
          color = "black"
        )
    }
  }

  n_hist_sig <- nrow(history_tukey)
  y_bar_hist <- if (n_hist_sig == 1L) {
    1
  } else if (n_hist_sig == 2L) {
    c(1, 1.1)
  } else {
    NULL
  }
  y_lbl_hist <- if (n_hist_sig == 1L) {
    1.05
  } else if (n_hist_sig == 2L) {
    c(1.025, 1.125)
  } else {
    NULL
  }
  p0 <- composition_add_significance_bars(
    plot_obj = p0,
    data_tbl = disp_data %>% dplyr::mutate(Distance_for_sig = .data$Distance + (4 * y_step)),
    x_var = "History_Label",
    y_var = "Distance_for_sig",
    sig_tbl = history_tukey,
    level_values = c("No prior stressors", "One prior stressor", "Two prior stressors"),
    stack_factor = 1.1,
    label_above_factor = 1.02,
    y_bar_per_row = y_bar_hist,
    y_label_per_row = y_lbl_hist,
    y_hi_cap = 1.15
  )

  y_expand <- min(
    1.15,
    max(
      1.1,
      y_lbl_within + y_range * 0.05,
      if (n_hist_sig >= 1L) max(y_lbl_hist, na.rm = TRUE) + y_range * 0.06 else 0
    )
  )
  p0 <- p0 + ggplot2::expand_limits(y = y_expand) +
    ggplot2::theme(
      # Bottom margin: keep modest padding under stacked bottom legends. A very large b
      # (e.g. 110 pt) was tuned for taller exports and reads as excessive white space on square ggsave.
      plot.margin = ggplot2::margin(t = 40, r = 24, b = 44, l = 24, unit = "pt")
    )
  p0
}

# Betadisper violins for factorial A×T×P groups; P- = opaque white fill, P+ = grey75; outline = regime color; medians via draw_quantiles on outline violins.
# Uses `treatment_color_scale` and `treatment_order` from 02__PlotSettings.R (Michael Sieler; 2026-04-24).
composition_betadisper_boxplot_atp <- function(bdisp_obj, title_suffix) {
  disp_data <- data.frame(
    Distance = bdisp_obj$ATP_group$model$distances,
    Group = bdisp_obj$ATP_group$model$group
  ) %>%
    dplyr::mutate(
      Group = factor(as.character(.data$Group), levels = treatment_order),
      fill_cond = dplyr::if_else(
        stringr::str_detect(as.character(.data$Group), "P\\+$"),
        "grey75",
        "white",
        missing = "white"
      )
    )
  fill_parasite <- c(white = "#FFFFFF", grey75 = grDevices::grey(0.75))
  ggplot2::ggplot(
    disp_data,
    ggplot2::aes(
      x = .data$Group,
      y = .data$Distance,
      colour = .data$Group,
      group = .data$Group
    )
  ) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = .data$Group),
      alpha = 0.16,
      colour = NA,
      linewidth = 0,
      trim = TRUE,
      width = 0.92,
      scale = "width",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = treatment_color_scale, guide = "none") +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_violin(
      ggplot2::aes(fill = .data$fill_cond),
      colour = "white",
      linewidth = 1.08,
      trim = TRUE,
      width = 0.86,
      scale = "width",
      draw_quantiles = c(0.5),
      show.legend = FALSE
    ) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = .data$fill_cond),
      linewidth = 0.48,
      trim = TRUE,
      width = 0.74,
      scale = "width",
      draw_quantiles = c(0.5),
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fill_parasite, guide = "none") +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(fill = .data$fill_cond),
      colour = "white",
      stroke = 2.15,
      shape = 23L,
      groupOnX = TRUE,
      size = 1.18
    ) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(colour = .data$Group, fill = .data$fill_cond),
      stroke = 1,
      shape = 23L,
      groupOnX = TRUE,
      size = 1.05
    ) +
    ggplot2::scale_colour_manual(
      values = treatment_color_scale,
      drop = FALSE,
      name = "Exposure Regime",
      guide = ggplot2::guide_legend(
        nrow = 2L,
        byrow = TRUE,
        override.aes = list(
          alpha = 1,
          fill = "white",
          shape = 23L,
          size = 3
        )
      )
    ) +
    ggplot2::labs(
      x = "A × T × P group",
      y = "Distance to centroid",
      title = paste0("Beta diversity dispersion by factorial exposure regimes (", title_suffix, ")")
    ) +
    theme_sieler2026_composition_figure(base_size = 13, legend_position = "bottom") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# betadisper ANOVA table (Groups vs Residuals) -> tidy tibble with Distance label; used in 02__Composition.R.
composition_betadisper_anova_to_tidy <- function(anova_obj, distance_label) {
  anova_obj %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Term") %>%
    dplyr::mutate(
      Distance = distance_label,
      dplyr::across(where(is.numeric), ~ round(.x, 6))
    ) %>%
    dplyr::relocate("Distance", .before = 1) %>%
    dplyr::rename(p_value = `Pr(>F)`)
}

# =============================================================================
# 8. MaAsLin3 — read results, manuscript tables (03__DiffAbund.R)
# =============================================================================

# Read MaAsLin3 `all_results.tsv` (tab-separated) into a tibble.
maaslin_read_all_results_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("MaAsLin3 results not found: ", path)
  }
  readr::read_tsv(path, show_col_types = FALSE)
}

# Top taxa table for exposure-regime contrasts: q < q_max, non-reference regimes vs control,
# ranked by |coef|, then first n_top unique features (manuscript Table 1A-style).
# MaAsLin3 uses metadata == "Treatment" and `value` = regime label; reference is excluded via `value`.
maaslin_top_taxa_treatment_table <- function(
    results_tbl,
    q_col = "qval_individual",
    coef_col = "coef",
    feature_col = "feature",
    metadata_col = "metadata",
    value_col = "value",
    reference_treatment = "A- T- P-",
    n_top = 10L,
    q_max = 0.05
) {
  req <- c(q_col, coef_col, feature_col, metadata_col)
  if (!all(req %in% names(results_tbl))) {
    stop(
      "MaAsLin3 table missing expected columns. Need: ",
      paste(req, collapse = ", "),
      " — got: ",
      paste(names(results_tbl), collapse = ", ")
    )
  }
  sig <- results_tbl %>%
    dplyr::filter(
      !is.na(.data[[q_col]]),
      .data[[q_col]] < q_max,
      .data[[metadata_col]] == "Treatment",
      if (value_col %in% names(results_tbl)) {
        .data[[value_col]] != reference_treatment
      } else {
        TRUE
      }
    ) %>%
    dplyr::mutate(abs_coef = abs(.data[[coef_col]])) %>%
    dplyr::arrange(dplyr::desc(.data$abs_coef), .data[[q_col]])

  sig %>%
    dplyr::distinct(.data[[feature_col]], .keep_all = TRUE) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::select(-"abs_coef")
}

