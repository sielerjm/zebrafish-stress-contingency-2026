# 03__DiffAbund.R
# Created by: Michael Sieler
# Date last updated: 2026-04-27
#
# Description: Genus-level differential abundance with MaAsLin3 (final timepoint): eight exposure
#   regimes vs reference, parasite exposure, stress-history linear trend, and stress×parasite
#   interaction. Writes MaAsLin output folders under Results/03__DiffAbund/Stats/, summary tables,
#   and diffabund__gut__bundle.rds for Code/02__Results/03__DiffAbund.Rmd.
#
# Expected input:  Run from repo root; ps.list[["TimeFinal"]] from 04__DataPreProcess.R / ps-list RDS.
# Expected output:  Results/03__DiffAbund/Stats/{maaslin_*,diffabund__gut__bundle.rds};
#   Tables/top10_significant_taxa_exposure_regimes__TankModel.csv (manuscript Table 1A-style).
#
# Note: Bipartite gene–taxon networks (partial correlations) live in 06__Taxon-DEG-Mort.R — not MaAsLin3.

init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
if (!file.exists(init_rel)) {
  stop(
    "Set working directory to the Sieler2026 repository root, or run:\n",
    "  Rscript Code/01__Analysis/03__DiffAbund.R\n",
    "from the project root so `", init_rel, "` resolves."
  )
}
source(init_rel)

path_res <- file.path(path.results, "03__DiffAbund")
path_fig <- file.path(path_res, "Figures")
path_tbl <- file.path(path_res, "Tables")
path_stats <- file.path(path_res, "Stats")
dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(path_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(path_stats, recursive = TRUE, showWarnings = FALSE)

args_03 <- commandArgs(trailingOnly = TRUE)
figures_only_03 <- isTRUE(any(args_03 == "--figures-only"))
if (isTRUE(figures_only_03)) {
  csv_mort <- file.path(path_tbl, "mortality_tank_taxon_log2abund_vs_percent__Tank.csv")
  if (!file.exists(csv_mort)) {
    stop("03__DiffAbund.R --figures-only requires:\n  ", csv_mort, "\nRun a full 03 driver first.")
  }
  tank_df_fig <- readr::read_csv(csv_mort, show_col_types = FALSE)
  p52 <- diffabund_build_mortality_combined_scatter_plot(tank_df_fig)
  ggplot2::ggsave(
    file.path(path_fig, "maaslin_mortality_top_taxa_scatter_combined__Tank.pdf"),
    p52,
    width = 8,
    height = 8,
    device = "pdf"
  )
  ggplot2::ggsave(
    file.path(path_fig, "maaslin_mortality_top_taxa_scatter_combined__Tank.png"),
    p52,
    width = 8,
    height = 8,
    dpi = 300L
  )
  sieler2026_sync_main_figures_from_manifest(driver_script = "03__DiffAbund.R", panel_ids = "5.1")
  message("03__DiffAbund.R --figures-only complete.")
  quit(save = "no", status = 0L)
}

# Archive previous outputs (copy → verify → clear) so each run is clean but reproducible.
# Note: Stats may include MaAsLin3 output directories; these are archived as folders.
sieler2026_archive_module_outputs(
  path_res_module = path_res,
  module_name = "03__DiffAbund",
  subdirs = c("Figures", "Tables", "Stats"),
  clear_original = TRUE
)

bundle_rds <- file.path(path_stats, "diffabund__gut__bundle.rds")

FORCE_RERUN_MAASLIN <- FALSE
RERUN_EXPOSURE_REGIME <- FALSE
RERUN_STRESS_HISTORY <- FALSE
RERUN_MORTALITY <- FALSE

check_maaslin_results <- function(output_path) {
  file.exists(file.path(output_path, "all_results.tsv"))
}

cleanup_maaslin_results <- function(output_path, force_rerun) {
  if (isTRUE(force_rerun) && dir.exists(output_path)) {
    message("Removing existing MaAsLin3 folder: ", output_path)
    unlink(output_path, recursive = TRUE)
  }
}

n_cores <- min(8L, max(1L, parallel::detectCores()))

if (!exists("ps.list", inherits = TRUE)) {
  stop("ps.list not found. Run 04__DataPreProcess.R and ensure ps-list__*.rds exists under Data/r_objects/.")
}

ps_list_element <- "TimeFinal"
ps_final <- ps.list[[ps_list_element]]

ps_genus <- microViz::tax_agg(ps_final, rank = "Genus")

feature_table <- microViz::otu_get(ps_genus) %>%
  as.data.frame()

metadata <- microViz::samdat_tbl(ps_final) %>%
  tibble::column_to_rownames(var = ".sample_name") %>%
  dplyr::filter(rownames(.) %in% rownames(feature_table)) %>%
  dplyr::mutate(
    Treatment = as.factor(.data$Treatment),
    Treatment = forcats::fct_relevel(.data$Treatment, "A- T- P-"),
    HistoryLevel_f = factor(
      .data$HistoryLevelNum,
      levels = c(0, 1, 2),
      labels = c("No prior stressors", "One prior stressor", "Two prior stressors")
    )
  )

feature_table <- feature_table[rownames(metadata), , drop = FALSE]

message(
  "Feature table ", paste(dim(feature_table), collapse = " x "),
  "; metadata ", paste(dim(metadata), collapse = " x ")
)

out_exp_tank <- file.path(path_stats, "maaslin_ExposureRegimes_Tank")
cleanup_maaslin_results(out_exp_tank, FORCE_RERUN_MAASLIN || RERUN_EXPOSURE_REGIME)
if (!check_maaslin_results(out_exp_tank) || FORCE_RERUN_MAASLIN || RERUN_EXPOSURE_REGIME) {
  message("Running MaAsLin3: Treatment + (1|Tank.ID) ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = ~ Treatment + (1 | Tank.ID),
    small_random_effects = TRUE,
    output = out_exp_tank,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 ExposureRegimes_Tank (results present).")
}

out_exp_notank <- file.path(path_stats, "maaslin_ExposureRegimes_noTank")
cleanup_maaslin_results(out_exp_notank, FORCE_RERUN_MAASLIN || RERUN_EXPOSURE_REGIME)
if (!check_maaslin_results(out_exp_notank) || FORCE_RERUN_MAASLIN || RERUN_EXPOSURE_REGIME) {
  message("Running MaAsLin3: ~ Treatment ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = ~ Treatment,
    output = out_exp_notank,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 ExposureRegimes_noTank (results present).")
}

out_par_tank <- file.path(path_stats, "maaslin_Parasite_Tank")
cleanup_maaslin_results(out_par_tank, FORCE_RERUN_MAASLIN || RERUN_EXPOSURE_REGIME)
if (!check_maaslin_results(out_par_tank) || FORCE_RERUN_MAASLIN || RERUN_EXPOSURE_REGIME) {
  message("Running MaAsLin3: Parasite_Exposed + (1|Tank.ID) ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = ~ Parasite_Exposed + (1 | Tank.ID),
    small_random_effects = TRUE,
    output = out_par_tank,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 Parasite_Tank (results present).")
}

out_par_notank <- file.path(path_stats, "maaslin_Parasite_noTank")
cleanup_maaslin_results(out_par_notank, FORCE_RERUN_MAASLIN || RERUN_EXPOSURE_REGIME)
if (!check_maaslin_results(out_par_notank) || FORCE_RERUN_MAASLIN || RERUN_EXPOSURE_REGIME) {
  message("Running MaAsLin3: ~ Parasite_Exposed ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = ~ Parasite_Exposed,
    output = out_par_notank,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 Parasite_noTank (results present).")
}

out_stress <- file.path(path_stats, "maaslin_StressHistory_Tank")
cleanup_maaslin_results(out_stress, FORCE_RERUN_MAASLIN || RERUN_STRESS_HISTORY)
if (!check_maaslin_results(out_stress) || FORCE_RERUN_MAASLIN || RERUN_STRESS_HISTORY) {
  message("Running MaAsLin3: HistoryLevelNum + (1|Tank.ID) ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = ~ HistoryLevelNum + (1 | Tank.ID),
    small_random_effects = TRUE,
    output = out_stress,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 StressHistory_Tank (results present).")
}

out_stress_factor <- file.path(path_stats, "maaslin_StressHistoryFactor_Tank")
cleanup_maaslin_results(out_stress_factor, FORCE_RERUN_MAASLIN || RERUN_STRESS_HISTORY)
if (!check_maaslin_results(out_stress_factor) || FORCE_RERUN_MAASLIN || RERUN_STRESS_HISTORY) {
  message("Running MaAsLin3: HistoryLevel_f + (1|Tank.ID) ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = ~ HistoryLevel_f + (1 | Tank.ID),
    small_random_effects = TRUE,
    output = out_stress_factor,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 StressHistoryFactor_Tank (results present).")
}

out_stresspath <- file.path(path_stats, "maaslin_StressPath_Tank")
cleanup_maaslin_results(out_stresspath, FORCE_RERUN_MAASLIN || RERUN_STRESS_HISTORY)
if (!check_maaslin_results(out_stresspath) || FORCE_RERUN_MAASLIN || RERUN_STRESS_HISTORY) {
  message("Running MaAsLin3: Parasite_Exposed * HistoryLevelNum + (1|Tank.ID) ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = ~ Parasite_Exposed * HistoryLevelNum + (1 | Tank.ID),
    small_random_effects = TRUE,
    output = out_stresspath,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 StressPath_Tank (results present).")
}

out_stresspath_factor <- file.path(path_stats, "maaslin_StressPathFactor_Tank")
cleanup_maaslin_results(out_stresspath_factor, FORCE_RERUN_MAASLIN || RERUN_STRESS_HISTORY)
if (!check_maaslin_results(out_stresspath_factor) || FORCE_RERUN_MAASLIN || RERUN_STRESS_HISTORY) {
  message("Running MaAsLin3: Parasite_Exposed * HistoryLevel_f + (1|Tank.ID) ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = ~ Parasite_Exposed * HistoryLevel_f + (1 | Tank.ID),
    small_random_effects = TRUE,
    output = out_stresspath_factor,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 StressPathFactor_Tank (results present).")
}

tbl_tank <- maaslin_read_all_results_tsv(file.path(out_exp_tank, "all_results.tsv"))
tbl_notank <- maaslin_read_all_results_tsv(file.path(out_exp_notank, "all_results.tsv"))

norm_names <- function(df) {
  if ("qval_individual" %in% names(df)) {
    return(df)
  }
  if ("qval" %in% names(df)) {
    return(dplyr::rename(df, qval_individual = "qval"))
  }
  df
}

tbl_tank <- norm_names(tbl_tank)
tbl_notank <- norm_names(tbl_notank)

top10_tbl <- tryCatch(
  maaslin_top_taxa_treatment_table(tbl_tank, n_top = 10L),
  error = function(e) {
    warning("Could not build top-10 table from Tank model: ", conditionMessage(e))
    tibble::tibble()
  }
)

readr::write_csv(
  top10_tbl,
  file.path(path_tbl, "top10_significant_taxa_exposure_regimes__TankModel.csv")
)

# ==============================================================================
# Mortality analysis (percent_mortality covariate; models on sample-level data)
# ==============================================================================

tank_mortality <- data.list[["Mortality"]] %>%
  dplyr::filter(Time == 60) %>%
  dplyr::group_by(Treatment, Tank.ID) %>%
  dplyr::summarise(
    alive = dplyr::n(),
    at_risk = 15, # 15 fish per tank at Day 60
    dead = at_risk - alive,
    percent_mortality = round((dead / at_risk) * 100, 1),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Treatment = as.character(Treatment),
    Tank.ID = as.numeric(as.character(Tank.ID))
  )

metadata_mort <- metadata %>%
  dplyr::mutate(
    Treatment = as.character(Treatment),
    Tank.ID = as.numeric(as.character(Tank.ID))
  )

# Add percent_mortality to sample-level metadata without disturbing rownames
sample_tank_key <- paste0(metadata_mort$Treatment, "_", metadata_mort$Tank.ID)
tank_key <- paste0(tank_mortality$Treatment, "_", tank_mortality$Tank.ID)
metadata_mort$percent_mortality <- tank_mortality$percent_mortality[match(sample_tank_key, tank_key)]

if (any(is.na(metadata_mort$percent_mortality))) {
  warning(
    "Some samples are missing percent_mortality after tank join; ",
    "check Treatment/Tank.ID matching. Missing: ",
    sum(is.na(metadata_mort$percent_mortality))
  )
}

out_mort_notank <- file.path(path_stats, "maaslin_Mortality_noTank")
cleanup_maaslin_results(out_mort_notank, FORCE_RERUN_MAASLIN || RERUN_MORTALITY)
if (!check_maaslin_results(out_mort_notank) || FORCE_RERUN_MAASLIN || RERUN_MORTALITY) {
  message("Running MaAsLin3: percent_mortality ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata_mort,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = "~ percent_mortality",
    output = out_mort_notank,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 Mortality_noTank (results present).")
}

out_mort_tank <- file.path(path_stats, "maaslin_Mortality_Tank")
cleanup_maaslin_results(out_mort_tank, FORCE_RERUN_MAASLIN || RERUN_MORTALITY)
if (!check_maaslin_results(out_mort_tank) || FORCE_RERUN_MAASLIN || RERUN_MORTALITY) {
  message("Running MaAsLin3: percent_mortality + (1|Tank.ID) ...")
  set.seed(42)
  maaslin3::maaslin3(
    input_data = feature_table,
    input_metadata = metadata_mort,
    min_abundance = 0.001,
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    correction = "BH",
    standardize = FALSE,
    formula = "~ percent_mortality + (1|Tank.ID)",
    small_random_effects = TRUE,
    output = out_mort_tank,
    cores = n_cores
  )
} else {
  message("Skipping MaAsLin3 Mortality_Tank (results present).")
}

# ==============================================================================
# Summary tables for Results Rmd (GT + CSV exports)
# ==============================================================================

# Build: (1) significant taxa counts by regime (noTank), (2) significant taxa counts by history (Tank),
# (3) significant taxa counts by history × parasite (Tank), (4) top mortality-associated taxa (Tank).
# These were previously computed in Code/02__Results/03__DiffAbund.Rmd; they now live in the driver so
# the Results Rmd is display-only.

maaslin_sig_counts_table <- function(all_results_tbl, metadata_name, value_levels, title, subtitle) {
  out <- all_results_tbl %>%
    dplyr::filter(.data$model == "abundance", .data$metadata == metadata_name) %>%
    dplyr::filter(!is.na(.data$qval_individual), .data$qval_individual < 0.05) %>%
    dplyr::group_by(.data$value) %>%
    dplyr::summarise(
      n_significant = dplyr::n(),
      n_positive = sum(.data$coef > 0, na.rm = TRUE),
      n_negative = sum(.data$coef < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::right_join(tibble::tibble(value = value_levels), by = "value") %>%
    dplyr::arrange(factor(.data$value, levels = value_levels)) %>%
    dplyr::mutate(
      n_significant = tidyr::replace_na(.data$n_significant, 0L),
      n_positive = tidyr::replace_na(.data$n_positive, 0L),
      n_negative = tidyr::replace_na(.data$n_negative, 0L)
    )

  out_total <- out %>%
    dplyr::summarise(
      value = "Total",
      n_significant = sum(.data$n_significant, na.rm = TRUE),
      n_positive = sum(.data$n_positive, na.rm = TRUE),
      n_negative = sum(.data$n_negative, na.rm = TRUE),
      .groups = "drop"
    )
  out2 <- dplyr::bind_rows(out, out_total)

  gt_tbl <- out2 %>%
    gt::gt() %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::cols_label(
      value = "Group",
      n_significant = "Total Significant",
      n_positive = "Positive Effect",
      n_negative = "Negative Effect"
    )

  list(data = out2, gt = gt_tbl)
}

maaslin_top_mortality_taxa_table <- function(all_results_tbl, n_top = 10L) {
  mort_sig <- all_results_tbl %>%
    dplyr::filter(
      .data$model == "abundance",
      .data$metadata == "percent_mortality",
      !is.na(.data$qval_individual),
      .data$qval_individual < 0.05
    )

  out <- mort_sig %>%
    dplyr::arrange(.data$qval_individual) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::mutate(
      Taxon = as.character(.data$feature),
      coef_val = .data$coef,
      q_value = .data$qval_individual,
      Association = dplyr::case_when(
        .data$coef_val > 0 ~ "Positive",
        .data$coef_val < 0 ~ "Negative",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select("Taxon", "coef_val", "q_value", "Association") %>%
    dplyr::mutate(
      coef_val = round(.data$coef_val, 3),
      q_value = as.numeric(.data$q_value)
    )

  gt_tbl <- out %>%
    gt::gt() %>%
    gt::tab_header(
      title = "Top significant taxa associated with mortality",
      subtitle = "MaAsLin3 (TSS + LOG): tank effect model; q < 0.05; top taxa by smallest q-value"
    ) %>%
    gt::cols_label(
      Taxon = "Taxon",
      coef_val = "Coefficient (log2 scale)",
      q_value = "q-value",
      Association = "Association"
    )

  gt_tbl <- style_gt_significance(gt_tbl, out, alpha = 0.05)

  list(data = out, gt = gt_tbl)
}

# Read MaAsLin3 results needed for summary tables
tbl_exp_notank <- maaslin_read_all_results_tsv(file.path(out_exp_notank, "all_results.tsv")) %>% norm_names()
tbl_hist_factor <- maaslin_read_all_results_tsv(file.path(out_stress_factor, "all_results.tsv")) %>% norm_names()
tbl_histpath_factor <- maaslin_read_all_results_tsv(file.path(out_stresspath_factor, "all_results.tsv")) %>% norm_names()
tbl_mort_tank <- maaslin_read_all_results_tsv(file.path(out_mort_tank, "all_results.tsv")) %>% norm_names()

treatment_order_no_ref <- c(
  "A- T- P+",
  "A- T+ P-",
  "A- T+ P+",
  "A+ T- P-",
  "A+ T- P+",
  "A+ T+ P-",
  "A+ T+ P+"
)
history_order_labels <- c("No prior stressors", "One prior stressor", "Two prior stressors")

tbl_counts_regime <- maaslin_sig_counts_table(
  all_results_tbl = tbl_exp_notank,
  metadata_name = "Treatment",
  value_levels = treatment_order_no_ref,
  title = "MaAsLin3: significant taxa counts by exposure regime",
  subtitle = "No-tank model; q < 0.05; counts of positive vs negative abundance coefficients (reference: A- T- P-)"
)
readr::write_csv(
  tbl_counts_regime$data,
  file.path(path_tbl, "significant_taxa_counts_by_exposure_regime__noTank.csv")
)
gt::gtsave(tbl_counts_regime$gt, file.path(path_tbl, "significant_taxa_counts_by_exposure_regime__noTank.html"))

tbl_counts_history <- maaslin_sig_counts_table(
  all_results_tbl = tbl_hist_factor,
  metadata_name = "HistoryLevel_f",
  value_levels = history_order_labels,
  title = "MaAsLin3: significant taxa counts by prior stressor history",
  subtitle = "Tank effect model; q < 0.05. Rows are stress-history levels (reference: No prior stressors)."
)
readr::write_csv(
  tbl_counts_history$data,
  file.path(path_tbl, "significant_taxa_counts_by_stress_history_levels__Tank.csv")
)
gt::gtsave(tbl_counts_history$gt, file.path(path_tbl, "significant_taxa_counts_by_stress_history_levels__Tank.html"))

tbl_counts_histpath <- maaslin_sig_counts_table(
  all_results_tbl = tbl_histpath_factor,
  metadata_name = "Parasite_Exposed:HistoryLevel_f",
  value_levels = history_order_labels,
  title = "MaAsLin3: significant taxa counts by stress history level × parasite exposure",
  subtitle = "Tank effect model; q < 0.05. Rows correspond to interaction terms by stress-history level."
)
readr::write_csv(
  tbl_counts_histpath$data,
  file.path(path_tbl, "significant_taxa_counts_by_stress_history_levels_x_parasite__Tank.csv")
)
gt::gtsave(
  tbl_counts_histpath$gt,
  file.path(path_tbl, "significant_taxa_counts_by_stress_history_levels_x_parasite__Tank.html")
)

tbl_top_mortality <- maaslin_top_mortality_taxa_table(tbl_mort_tank, n_top = 10L)
readr::write_csv(
  tbl_top_mortality$data,
  file.path(path_tbl, "top10_taxa_associated_with_mortality__Tank.csv")
)
gt::gtsave(
  tbl_top_mortality$gt,
  file.path(path_tbl, "top10_taxa_associated_with_mortality__Tank.html")
)

# ==============================================================================
# Mortality figures + backing tables (export to Results/03__DiffAbund/Figures + Tables)
# ==============================================================================

# These figures were previously generated only during knitting of
# Code/02__Results/03__DiffAbund.Rmd. To ensure they are included in the combined
# supplementary figures PDF, we export stable PDF/PNG artifacts from the driver.

save_plot_pair <- function(plot, stem_no_ext, w_in, h_in, dpi = 300L) {
  pdf_path <- file.path(path_fig, paste0(stem_no_ext, ".pdf"))
  png_path <- file.path(path_fig, paste0(stem_no_ext, ".png"))
  ggplot2::ggsave(pdf_path, plot = plot, width = w_in, height = h_in, units = "in", device = "pdf")
  ggplot2::ggsave(png_path, plot = plot, width = w_in, height = h_in, units = "in", dpi = dpi)
  list(pdf = pdf_path, png = png_path)
}

# Recompute the same top taxa list the GT table used (q < 0.05, smallest q)
mort_sig <- tbl_mort_tank %>%
  dplyr::filter(
    .data$model == "abundance",
    .data$metadata == "percent_mortality",
    !is.na(.data$qval_individual),
    .data$qval_individual < 0.05
  )

top10_mort <- mort_sig %>%
  dplyr::arrange(.data$qval_individual) %>%
  dplyr::slice_head(n = 10L) %>%
  dplyr::mutate(Taxon = as.character(.data$feature))

taxa_sig <- top10_mort %>%
  dplyr::arrange(.data$qval_individual) %>%
  dplyr::distinct(.data$Taxon) %>%
  dplyr::pull(.data$Taxon) %>%
  as.character()

# --- Figure 1: coefficient bar plot (tank random effect model) -----------------
bar_df <- top10_mort %>%
  dplyr::mutate(
    Taxon = as.character(.data$feature),
    Association = dplyr::case_when(
      .data$coef > 0 ~ "Positive",
      .data$coef < 0 ~ "Negative",
      TRUE ~ NA_character_
    ),
    Association = factor(.data$Association, levels = c("Positive", "Negative")),
    q_label = dplyr::if_else(
      .data$qval_individual < 0.0001,
      "q<0.0001",
      paste0("q=", sprintf("%.4f", .data$qval_individual))
    )
  )

taxon_order <- bar_df %>%
  dplyr::arrange(.data$qval_individual) %>%
  dplyr::pull(.data$Taxon) %>%
  unique()

bar_df <- bar_df %>%
  dplyr::filter(.data$Taxon %in% taxon_order) %>%
  dplyr::mutate(
    Taxon = factor(.data$Taxon, levels = rev(taxon_order)),
    label_offset = 0.02 * max(abs(.data$coef), na.rm = TRUE),
    label_offset_eff = pmin(.data$label_offset, abs(.data$coef) * 0.25),
    label_x = dplyr::case_when(
      .data$coef > 0 ~ .data$coef - .data$label_offset_eff,
      .data$coef < 0 ~ .data$coef + .data$label_offset_eff,
      TRUE ~ 0
    ),
    hjust_val = dplyr::case_when(
      .data$coef > 0 ~ 1,
      .data$coef < 0 ~ 0,
      TRUE ~ 0.5
    )
  )

readr::write_csv(
  bar_df %>%
    dplyr::select(
      Taxon,
      coef = .data$coef,
      qval_individual = .data$qval_individual,
      Association,
      q_label
    ),
  file.path(path_tbl, "mortality_maaslin_top_taxa_barplot_source__Tank.csv")
)

p_mort_bar <- ggplot2::ggplot(
  bar_df,
  ggplot2::aes(x = .data$coef, y = .data$Taxon, fill = .data$Association)
) +
  # Dummy layer ensures both legend categories display even if absent.
  ggplot2::geom_col(
    data = dplyr::tibble(
      Taxon = levels(bar_df$Taxon)[1],
      coef = 0,
      Association = factor("Positive", levels = c("Positive", "Negative"))
    ),
    ggplot2::aes(x = .data$coef, y = .data$Taxon, fill = .data$Association),
    width = 0.8,
    alpha = 0,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::scale_fill_manual(
    values = c("Positive" = "steelblue", "Negative" = "firebrick"),
    breaks = c("Positive", "Negative"),
    limits = c("Positive", "Negative"),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "MaAsLin3: Significant taxa vs percent mortality",
    subtitle = "Taxon Abundance ~ Percent Mortality + (1|Tank.ID)",
    x = "MaAsLin3 association coefficient (log2 scale)",
    y = NULL,
    fill = "Association"
  ) +
  theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$q_label, x = .data$label_x, hjust = .data$hjust_val),
    color = "white",
    fontface = "bold",
    size = 4,
    vjust = 0.5,
    inherit.aes = TRUE
  ) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(face = "bold", size = 14),
    axis.text.x = ggplot2::element_text(face = "bold", size = 12),
    axis.text.y = ggplot2::element_text(face = "bold", size = 12),
    plot.title = ggplot2::element_text(face = "bold", size = 14)
  )

fig_mort_bar_paths <- NULL
if (nrow(bar_df) > 0L) {
  fig_mort_bar_paths <- save_plot_pair(
    p_mort_bar,
    stem_no_ext = "maaslin_mortality_top_taxa_coef_bar__Tank",
    w_in = 8,
    h_in = 8
  )
}

# --- Figure 2: tank-level scatter (log2 abundance vs percent mortality) --------
ps_genus_scatter <- microViz::tax_agg(ps_final, rank = "Genus")
counts_scatter <- microViz::otu_get(ps_genus_scatter) %>%
  as.data.frame()

meta_scatter <- microViz::samdat_tbl(ps_final) %>%
  tibble::column_to_rownames(var = ".sample_name") %>%
  dplyr::filter(rownames(.) %in% rownames(counts_scatter)) %>%
  dplyr::mutate(
    Tank.ID = as.numeric(as.character(.data$Tank.ID)),
    Treatment = as.character(.data$Treatment)
  )

counts_scatter <- counts_scatter[rownames(meta_scatter), , drop = FALSE]

counts_mat <- as.matrix(counts_scatter)
rel_mat <- sweep(counts_mat, 1, rowSums(counts_mat), "/")

rel_sig <- rel_mat[, intersect(taxa_sig, colnames(rel_mat)), drop = FALSE]

tank_df <- tibble::tibble()
spearman_stats <- tibble::tibble()

if (ncol(rel_sig) > 0L) {
  min_nonzero <- suppressWarnings(min(rel_sig[rel_sig > 0], na.rm = TRUE))
  halfmin <- if (is.finite(min_nonzero) && min_nonzero > 0) min_nonzero / 2 else 1e-12
  log2_sig <- log2(rel_sig + halfmin)

  log2_df <- as.data.frame(log2_sig) %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of("sample_id"),
      names_to = "Taxon",
      values_to = "log2_abund"
    ) %>%
    dplyr::left_join(
      meta_scatter %>%
        tibble::rownames_to_column(var = "sample_id"),
      by = "sample_id"
    )

  tank_summary <- log2_df %>%
    dplyr::group_by(.data$Tank.ID, .data$Treatment, .data$Taxon) %>%
    dplyr::summarise(
      log2_abund = mean(.data$log2_abund, na.rm = TRUE),
      .groups = "drop"
    )

  rel_long <- as.data.frame(rel_sig) %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of("sample_id"),
      names_to = "Taxon",
      values_to = "rel_abund"
    ) %>%
    dplyr::left_join(
      meta_scatter %>%
        tibble::rownames_to_column(var = "sample_id"),
      by = "sample_id"
    )

  tank_presence <- rel_long %>%
    dplyr::group_by(.data$Tank.ID, .data$Treatment, .data$Taxon) %>%
    dplyr::summarise(
      detected = any(.data$rel_abund > 0),
      .groups = "drop"
    )

  tank_df <- tank_summary %>%
    dplyr::select(.data$Tank.ID, .data$Treatment, .data$Taxon, .data$log2_abund) %>%
    dplyr::left_join(
      tank_mortality %>% dplyr::select(.data$Tank.ID, .data$Treatment, .data$percent_mortality),
      by = c("Tank.ID", "Treatment")
    ) %>%
    dplyr::left_join(
      tank_presence %>%
        dplyr::group_by(.data$Taxon) %>%
        dplyr::summarise(prevalence = mean(.data$detected) * 100, .groups = "drop"),
      by = "Taxon"
    )

  spearman_stats <- tank_df %>%
    dplyr::group_by(.data$Taxon) %>%
    dplyr::summarise(
      rho = {
        d <- dplyr::cur_data()
        ok <- stats::complete.cases(d$log2_abund, d$percent_mortality)
        if (sum(ok) >= 3L) {
          stats::cor(d$log2_abund[ok], d$percent_mortality[ok], method = "spearman")
        } else {
          NA_real_
        }
      },
      p_value = {
        d <- dplyr::cur_data()
        ok <- stats::complete.cases(d$log2_abund, d$percent_mortality)
        if (sum(ok) >= 3L) {
          tryCatch(
            stats::cor.test(d$log2_abund[ok], d$percent_mortality[ok], method = "spearman", exact = FALSE)$p.value,
            error = function(e) NA_real_
          )
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(q_value = stats::p.adjust(.data$p_value, method = "BH"))

  tank_df <- tank_df %>%
    dplyr::left_join(spearman_stats, by = "Taxon") %>%
    dplyr::mutate(
      TaxonFacet = paste0(
        .data$Taxon,
        "\nSpearman rho = ", round(.data$rho, 2),
        ", q = ", signif(.data$q_value, 3),
        "\nPrevalence = ", round(.data$prevalence, 0), "%"
      )
    )

  readr::write_csv(tank_df, file.path(path_tbl, "mortality_tank_taxon_log2abund_vs_percent__Tank.csv"))
  readr::write_csv(spearman_stats, file.path(path_tbl, "mortality_tank_taxon_spearman_vs_percent__Tank.csv"))

  treatment_levels <- sort(unique(as.character(tank_df$Treatment)))
  cols_treat <- get_treatment_colors_safe(treatment_levels)

  p_mort_scatter <- ggplot2::ggplot(
    tank_df,
    ggplot2::aes(x = .data$log2_abund, y = .data$percent_mortality, color = .data$Treatment)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.6) +
    ggplot2::scale_color_manual(values = cols_treat, drop = FALSE) +
    ggplot2::labs(
      title = "Mortality scatter: significant taxa vs percent_mortality (tank-level means)",
      subtitle = "Points colored by exposure regime; dashed line + ribbon show overall trend (all treatments combined)",
      x = "Mean Log2 Taxon Abundance (per tank)",
      y = "Percent Mortality (per tank)"
    ) +
    ggplot2::geom_smooth(
      data = tank_df,
      mapping = ggplot2::aes(x = .data$log2_abund, y = .data$percent_mortality, group = 1),
      inherit.aes = FALSE,
      method = "lm",
      se = TRUE,
      color = "black",
      fill = "grey70",
      alpha = 0.35,
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_y_continuous(limits = c(-25, 100), breaks = seq(0, 100, 25)) +
    theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::facet_wrap(~TaxonFacet, scales = "free_x")

  fig_mort_scatter_paths <- save_plot_pair(
    p_mort_scatter,
    stem_no_ext = "maaslin_mortality_top_taxa_scatter_facets__Tank",
    w_in = 14,
    h_in = 10
  )

  # --- Figure 2b: same tank-level data, one panel — LM + SE per taxon, points above ribbons ----------
  p_mort_combined <- diffabund_build_mortality_combined_scatter_plot(tank_df)

  fig_mort_combined_paths <- save_plot_pair(
    p_mort_combined,
    stem_no_ext = "maaslin_mortality_top_taxa_scatter_combined__Tank",
    w_in = 8,
    h_in = 8
  )

  # --- Figure 3: Culicoidibacter-only panel ------------------------------------
  culico_taxon <- "Culicoidibacter"
  culico_df <- tank_df %>% dplyr::filter(.data$Taxon == culico_taxon)

  if (nrow(culico_df) > 0L) {
    culico_stats <- culico_df %>%
      dplyr::summarise(
        rho = dplyr::first(.data$rho),
        q_value = dplyr::first(.data$q_value),
        prevalence = dplyr::first(.data$prevalence),
        .groups = "drop"
      )

    p_culico_scatter <- ggplot2::ggplot(
      culico_df,
      ggplot2::aes(x = .data$log2_abund, y = .data$percent_mortality, color = .data$Treatment)
    ) +
      ggplot2::geom_hline(yintercept = 0, color = "grey40", linewidth = 0.7) +
      ggplot2::geom_smooth(
        data = culico_df,
        mapping = ggplot2::aes(x = .data$log2_abund, y = .data$percent_mortality, group = 1),
        inherit.aes = FALSE,
        method = "lm",
        se = TRUE,
        color = "black",
        fill = "grey70",
        alpha = 0.35,
        linetype = "dashed",
        linewidth = 1
      ) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::scale_color_manual(values = cols_treat, drop = FALSE) +
      ggplot2::labs(
        title = paste0(culico_taxon, " Abundance vs. Mortality"),
        subtitle = paste0(
          "Spearman rho = ", round(culico_stats$rho, 2),
          ", q = ", signif(culico_stats$q_value, 3),
          ", Prevalence = ", round(culico_stats$prevalence, 0), "%"
        ),
        x = "Mean Log2 Taxon Abundance",
        y = "Percent Mortality"
      ) +
      ggplot2::scale_y_continuous(limits = c(-25, 100), breaks = seq(0, 100, 25)) +
      theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 11),
        axis.title = ggplot2::element_text(face = "bold"),
        legend.title = ggplot2::element_text(face = "bold")
      )

    fig_culico_paths <- save_plot_pair(
      p_culico_scatter,
      stem_no_ext = "maaslin_mortality_culicoidibacter_scatter__Tank",
      w_in = 8,
      h_in = 6
    )
  }
}

bundle <- list(
  bundle_version = "1.0",
  description = "MaAsLin3 genus-level differential abundance (gut, TimeFinal)",
  paths = list(
    maaslin_ExposureRegimes_Tank = out_exp_tank,
    maaslin_ExposureRegimes_noTank = out_exp_notank,
    maaslin_Parasite_Tank = out_par_tank,
    maaslin_Parasite_noTank = out_par_notank,
    maaslin_StressHistory_Tank = out_stress,
    maaslin_StressHistoryFactor_Tank = out_stress_factor,
    maaslin_StressPath_Tank = out_stresspath,
    maaslin_StressPathFactor_Tank = out_stresspath_factor,
    maaslin_Mortality_noTank = out_mort_notank,
    maaslin_Mortality_Tank = out_mort_tank
  ),
  tables = list(
    top10_exposure_regimes = top10_tbl,
    all_results_exposure_tank_preview = head(tbl_tank, 20L),
    sig_counts_by_exposure_regime_noTank = tbl_counts_regime$data,
    sig_counts_by_stress_history_tank = tbl_counts_history$data,
    sig_counts_by_stress_history_x_parasite_tank = tbl_counts_histpath$data,
    top10_taxa_mortality_tank = tbl_top_mortality$data
  ),
  figures = list(
    mortality_top_taxa_coef_bar_tank = file.path("Results", "03__DiffAbund", "Figures", "maaslin_mortality_top_taxa_coef_bar__Tank.png"),
    mortality_top_taxa_scatter_facets_tank = file.path("Results", "03__DiffAbund", "Figures", "maaslin_mortality_top_taxa_scatter_facets__Tank.png"),
    mortality_top_taxa_scatter_combined_tank = file.path(
      "Results", "03__DiffAbund", "Figures", "maaslin_mortality_top_taxa_scatter_combined__Tank.png"
    ),
    mortality_culicoidibacter_scatter_tank = file.path("Results", "03__DiffAbund", "Figures", "maaslin_mortality_culicoidibacter_scatter__Tank.png")
  ),
  table_sig_counts_by_exposure_regime_noTank = tbl_counts_regime$gt,
  table_sig_counts_by_stress_history_tank = tbl_counts_history$gt,
  table_sig_counts_by_stress_history_x_parasite_tank = tbl_counts_histpath$gt,
  table_top10_taxa_associated_with_mortality_tank = tbl_top_mortality$gt,
  meta = list(
    reference_treatment = "A- T- P-",
    formula_exposure_tank = "~ Treatment + (1|Tank.ID)",
    formula_stress = "~ HistoryLevelNum + (1|Tank.ID)",
    formula_stress_factor = "~ HistoryLevel_f + (1|Tank.ID)",
    formula_stress_path = "~ Parasite_Exposed * HistoryLevelNum + (1|Tank.ID)",
    formula_stress_path_factor = "~ Parasite_Exposed * HistoryLevel_f + (1|Tank.ID)",
    mortality_covariate = "percent_mortality",
    formula_mortality_noTank = "~ percent_mortality",
    formula_mortality_tank = "~ percent_mortality + (1|Tank.ID)",
    ps_list_element = ps_list_element,
    min_abundance = 0.001,
    min_prevalence = 0.1
  )
)

saveRDS(bundle, bundle_rds)
message("Wrote bundle: ", bundle_rds)

sieler2026_sync_main_figures_from_manifest(driver_script = "03__DiffAbund.R")
