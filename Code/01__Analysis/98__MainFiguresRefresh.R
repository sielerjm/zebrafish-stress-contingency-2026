# 98__MainFiguresRefresh.R
# Created by: Michael Sieler
# Date last updated: 2026-04-27
#
# Description: Re-`ggsave` main-text panels whose ggplot objects are stored in module bundle RDS
#   files, then mirror updated PNG/PDF into `Manuscript/MainFigures/Exports/` via the manifest.
#   Use after theme/size tweaks without re-running full analysis drivers.
#
# Limitation: This re-exports ggplots as last saved in each bundle. If you changed theme or layer
#   code but did not re-run the owning driver (so the bundle still holds old plot objects), run a
#   full driver once (or the relevant `--figures-only` driver for 03/04/05) before refresh.
#
# Expected input: Run from Sieler2026 repo root; existing Results/... bundle RDS and manifest CSV.
# Expected output: Overwrites listed figure files under Results/.../Figures/ and copies into Exports/.
#
# Examples:
#   Rscript Code/01__Analysis/98__MainFiguresRefresh.R --panels=2.2,4.1,4.2
#   Rscript Code/01__Analysis/98__MainFiguresRefresh.R --panels=2.2,2.3,2.4
#   Rscript Code/01__Analysis/98__MainFiguresRefresh.R --panels=all

init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
if (!file.exists(init_rel)) {
  stop(
    "Set working directory to the Sieler2026 repository root, or run:\n",
    "  Rscript Code/01__Analysis/98__MainFiguresRefresh.R\n",
    "from the project root so `", init_rel, "` resolves."
  )
}
source(init_rel)

manifest_path <- file.path(proj.path, "Manuscript", "MainFigures", "main_figures_manifest.csv")
if (!file.exists(manifest_path)) {
  stop("Manifest not found: ", manifest_path)
}
manifest <- utils::read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c(panel_id = "character")
)

args_tr <- commandArgs(trailingOnly = TRUE)
arg_panels <- args_tr[startsWith(args_tr, "--panels=")]
if (length(arg_panels) != 1L) {
  stop("Provide exactly one argument like --panels=2.2 or --panels=all")
}
panels_arg <- sub("^--panels=", "", arg_panels)
if (!nzchar(panels_arg)) {
  stop("Empty --panels= value.")
}

bundle_refreshable <- c("2.1", "2.2", "2.3", "2.4", "4.1", "4.2", "5.1", "5.2", "5.3", "6.1")
if (tolower(panels_arg) == "all") {
  panel_ids <- bundle_refreshable
} else {
  panel_ids <- trimws(strsplit(panels_arg, ",", fixed = TRUE)[[1L]])
  panel_ids <- panel_ids[nzchar(panel_ids)]
  bad <- setdiff(panel_ids, bundle_refreshable)
  if (length(bad) > 0L) {
    stop(
      "Unsupported panel_id(s): ", paste(bad, collapse = ", "),
      ". Supported: ", paste(bundle_refreshable, collapse = ", "), ", or all."
    )
  }
}

manifest_row <- function(pid) {
  m <- manifest[as.character(manifest$panel_id) == as.character(pid), , drop = FALSE]
  if (nrow(m) != 1L) {
    stop("Manifest must contain exactly one row for panel_id ", pid)
  }
  m[1L, , drop = FALSE]
}

save_panel_from_ggplot <- function(pid, plot_obj) {
  if (is.null(plot_obj)) {
    warning("Panel ", pid, ": ggplot object is NULL; skip.")
    return(invisible(NULL))
  }
  row <- manifest_row(pid)
  w <- as.numeric(row$width_in)
  h <- as.numeric(row$height_in)
  dpi <- as.numeric(row$dpi_png)
  png_rel <- row$path_png
  pdf_rel <- row$path_pdf
  png_abs <- file.path(proj.path, png_rel)
  pdf_abs <- file.path(proj.path, pdf_rel)
  dir.create(dirname(png_abs), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(pdf_abs), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(pdf_abs, plot = plot_obj, width = w, height = h, units = "in", device = "pdf")
  ggplot2::ggsave(png_abs, plot = plot_obj, width = w, height = h, units = "in", dpi = dpi)
  message("Refreshed panel ", pid, " -> ", png_rel)
  sieler2026_sync_main_figures_from_manifest(driver_script = basename(as.character(row$driver_script)), panel_ids = pid)
  invisible(NULL)
}

read_bundle <- function(rel) {
  p <- file.path(proj.path, rel)
  if (!file.exists(p)) {
    stop("Bundle RDS not found: ", p)
  }
  readRDS(p)
}

for (pid in panel_ids) {
  if (identical(pid, "2.1")) {
    b <- read_bundle("Results/02__Composition/Stats/composition__gut__bundle.rds")
    save_panel_from_ggplot(pid, b$modules$relative_abundance$treatment$figure)
  } else if (identical(pid, "2.2")) {
    b <- read_bundle("Results/01__Diversity/Stats/diversity__gut__bundle.rds")
    save_panel_from_ggplot(pid, b$modules$stress_history$Simpson$figure_quasirandom)
  } else if (identical(pid, "2.3")) {
    b <- read_bundle("Results/02__Composition/Stats/composition__gut__bundle.rds")
    save_panel_from_ggplot(pid, b$modules$beta_bray$pcoa_history_parasite$figure)
  } else if (identical(pid, "2.4")) {
    b <- read_bundle("Results/02__Composition/Stats/composition__gut__bundle.rds")
    save_panel_from_ggplot(pid, b$modules$beta_bray$betadisper_history_parasite$figure)
  } else if (identical(pid, "4.1") || identical(pid, "4.2")) {
    pl <- read_bundle("Results/05__Mort-Inf/Stats/mortinf_main_text_figure_ggplots.rds")
    if (is.null(pl[["p_mortality_prior_stressor_trend_qr"]])) {
      stop(
        "Panel 4.1/4.2 refresh requires `p_mortality_prior_stressor_trend_qr` in:\n",
        "  Results/05__Mort-Inf/Stats/mortinf_main_text_figure_ggplots.rds\n",
        "Re-run a full `05__Mort-Inf.R` once (checkpoint RDS written at end of driver)."
      )
    }
    if (identical(pid, "4.1")) {
      save_panel_from_ggplot(pid, pl$p_mortality_prior_stressor_trend_qr)
    } else {
      save_panel_from_ggplot(pid, pl$p_infection_prevalence_trend_qr)
    }
  } else if (identical(pid, "5.1")) {
    csv_mort <- file.path(proj.path, "Results", "03__DiffAbund", "Tables", "mortality_tank_taxon_log2abund_vs_percent__Tank.csv")
    if (!file.exists(csv_mort)) {
      stop("Panel 5.1 refresh requires:\n  ", csv_mort, "\nRun a full 03__DiffAbund.R once (or `03__DiffAbund.R --figures-only` after outputs exist).")
    }
    tank_df_fig <- readr::read_csv(csv_mort, show_col_types = FALSE)
    p51 <- diffabund_build_mortality_combined_scatter_plot(tank_df_fig)
    save_panel_from_ggplot(pid, p51)
  } else if (identical(pid, "5.2")) {
    b <- read_bundle("Results/06__Taxon-DEG-Mort/Stats/taxon_deg_mort__bundle.rds")
    save_panel_from_ggplot(pid, b$ggplots$all_sig_partial_cor_scatter)
  } else if (identical(pid, "5.3")) {
    b <- read_bundle("Results/07__FunctionalAnno/Stats/functional_anno__bundle.rds")
    save_panel_from_ggplot(pid, b$ggplots_focal_four_combined$kegg)
  } else if (identical(pid, "6.1")) {
    b <- read_bundle("Results/08__NeutralModel/Stats/neutral_model__bundle.rds")
    save_panel_from_ggplot(pid, b$ggplot_focal_four_genera)
  }
}

# Ensure external manuscript assets are mirrored into Exports/ with special naming.
# - Fig 01-1: ExperimentalDesignSchematic
# - Table 01-1: ExposureSchematic (also generates TABLE_01-1 PDF from PNG for combining)
sieler2026_sync_main_figures_from_manifest(panel_ids = c("1.1", "1.0"))

# Assemble a single multi-page PDF for easy viewing (one figure per page) in manuscript order:
#   1) Fig 01-1 (experimental design schematic)
#   2) Table 01-1 (exposure schematic)
#   3) Remaining manuscript figures (Fig 2 -> Fig 6)
export_dir <- sieler2026_path_main_figure_exports()
combined_pdf <- file.path(export_dir, "MainFigures__all_panels__one_per_page.pdf")

pdf_candidates <- c(
  file.path(export_dir, "FIG_01-1__ExperimentalDesignSchematic.pdf"),
  file.path(export_dir, "TABLE_01-1__ExposureSchematic.pdf")
)

manifest_remaining <- manifest[
  !is.na(manifest$path_pdf) &
    nzchar(as.character(manifest$path_pdf)) &
    grepl("^\\s*[2-9]\\.", as.character(manifest$panel_id)),
  ,
  drop = FALSE
]
for (i in seq_len(nrow(manifest_remaining))) {
  pid <- as.character(manifest_remaining$panel_id[[i]])
  src_basename <- basename(as.character(manifest_remaining$path_pdf[[i]]))
  pdf_candidates <- c(
    pdf_candidates,
    file.path(export_dir, sieler2026_main_figure_export_filename(pid, src_basename))
  )
}

pdf_existing <- pdf_candidates[file.exists(pdf_candidates)]
pdf_missing <- setdiff(pdf_candidates, pdf_existing)
if (length(pdf_missing) > 0L) {
  warning(
    "Combined main-figures PDF: missing ", length(pdf_missing), " expected export PDF(s). ",
    "They will be skipped.\n",
    paste0(" - ", pdf_missing, collapse = "\n")
  )
}

if (!requireNamespace("qpdf", quietly = TRUE)) {
  stop(
    "Package 'qpdf' is required to combine PDFs.\n",
    "Install with: install.packages('qpdf')"
  )
}
if (length(pdf_existing) < 1L) {
  stop("No export PDFs found to combine under: ", export_dir)
}
qpdf::pdf_combine(input = pdf_existing, output = combined_pdf)
message("Wrote combined main-figures PDF -> ", combined_pdf)

message("98__MainFiguresRefresh complete.")
