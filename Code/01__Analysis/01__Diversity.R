# 01__Diversity.R
# Created by: Michael Sieler
# Date last updated: 2026-04-27
#
# Description: Gut alpha-diversity analysis — beta regression GLMMs (glmmTMB) for stress-history
#   and stress-history × parasite effects (Shannon, inverse Simpson, genus richness norms).
#   Simple effects: emmeans pairwise Parasite (exposed vs unexposed) within each HistoryLevelNum; violin plots.
#   Writes figures, tables, and a serialized results bundle for manuscript Rmd (no re-fit).
#   gt table subtitles and bundle$meta carry human-readable model formulas (beta GLMM, logit link).
#
# Manuscript Rmd: load bundle without re-running this script:
#   div <- readRDS(here::here("Results", "01__Diversity", "Stats", "diversity__gut__bundle.rds"))
#   div$table_combined_interaction
#   div$table_combined_table2
#   div$modules$stress_history$Simpson$figure
#   knitr::include_graphics(div$modules$stress_history$Simpson$paths[["figure_pdf"]])
#
# Expected input:  Run from Sieler2026 project root; `Data/r_objects/ps-list__*.rds` from
#   `04__DataPreProcess.R` (sample metadata includes HistoryLevel, HistoryLevelNum, Parasite, Tank.ID).
# Expected output:  `Results/01__Diversity/Figures/*.pdf`, `Tables/*.csv`, `Stats/diversity__gut__bundle.rds`.

# --- Ensure project root (Sieler2026/) ----------------------------------------
init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
if (!file.exists(init_rel)) {
  stop(
    "Set working directory to the Sieler2026 repository root, or run:\n",
    "  Rscript Code/01__Analysis/01__Diversity.R\n",
    "from the project root so `", init_rel, "` resolves."
  )
}
source(init_rel)

# --- Paths --------------------------------------------------------------------
path_res_div <- file.path(path.results, "01__Diversity")
path_fig <- file.path(path_res_div, "Figures")
path_tbl <- file.path(path_res_div, "Tables")
path_stats <- file.path(path_res_div, "Stats")
dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(path_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(path_stats, recursive = TRUE, showWarnings = FALSE)

# Bundle paths should be portable across machines: store relative-to-repo paths.
rel_fig <- file.path("Results", "01__Diversity", "Figures")

# Archive previous outputs (copy → verify → clear) so each run is clean but reproducible.
sieler2026_archive_module_outputs(
  path_res_module = path_res_div,
  module_name = "01__Diversity",
  subdirs = c("Figures", "Tables", "Stats"),
  clear_original = TRUE
)

bundle_rds <- file.path(path_stats, "diversity__gut__bundle.rds")

# --- Constants ----------------------------------------------------------------
ps_list_element <- "TimeFinal"

# Human-readable model text (gt subtitles + bundle meta). y = Genus-level normalized index in (0,1).
model_subtitle_linear_trend <- paste(
  "Beta regression GLMM (glmmTMB; Beta family, logit link; random intercept per tank):",
  "y ~ HistoryLevelNum + (1 | Tank.ID).",
  "Rows: linear contrast on estimated means at 0, 1, and 2 prior stressors."
)
model_subtitle_table2 <- paste(
  "P-values only (no estimates). Stress-history trend p: linear contrast on EMMs from",
  "y ~ HistoryLevelNum + (1 | Tank.ID) (full table above).",
  "History × parasite p: HistoryLevelNum:Parasite in y ~ HistoryLevelNum * Parasite + (1 | Tank.ID)",
  "(Wald; type II Anova via car::Anova when available). See interaction table for coefficients."
)
model_subtitle_interaction <- paste(
  "Beta regression GLMM (glmmTMB; Beta family, logit link; random intercept per tank):",
  "y ~ HistoryLevelNum * Parasite + (1 | Tank.ID).",
  "Rows: Wald test for HistoryLevelNum:Parasite on the logit (link) scale (conditional fixed effects)."
)
model_subtitle_parasite_simple <- paste(
  "Simple effects from the interaction model (emmeans; response scale): pairwise Parasite contrast (1 − 0)",
  "within each prior stressor count (0, 1, 2). p_fdr_bh: Benjamini–Hochberg FDR across the three strata per metric."
)


if (!exists("ps.list", inherits = TRUE)) {
  stop("ps.list not found. Run 04__DataPreProcess.R and ensure ps-list__*.rds exists under Data/r_objects/.")
}

# --- Model-ready data ---------------------------------------------------------
alpha_data_model <- build_alpha_data_model(ps.list, element_name = ps_list_element)

message(
  "Beta-regression bounds check — Simpson__Genus_norm: ",
  paste(round(range(alpha_data_model$Simpson__Genus_norm, na.rm = TRUE), 6), collapse = " to ")
)

# --- Factorial exposure-regime covariates (A × T × P) --------------------------
# Antibiotics / Temperature / Parasite are coded 0/1 in sample_data; use factors for categorical contrasts.
alpha_data_model <- alpha_data_model %>%
  dplyr::mutate(
    Antibiotics_f = factor(.data$Antibiotics, levels = c(0, 1), labels = c("A-", "A+")),
    Temperature_f = factor(.data$Temperature, levels = c(0, 1), labels = c("T-", "T+")),
    Parasite_f = factor(.data$Parasite, levels = c(0, 1), labels = c("Unexposed", "Exposed"))
  )

# --- Fit models (reproducible) -------------------------------------------------
set.seed(42)

fit_shannon_history_factor <- glmmTMB::glmmTMB(
  Shannon__Genus_norm ~ HistoryLevel + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)
fit_shannon_history_num <- glmmTMB::glmmTMB(
  Shannon__Genus_norm ~ HistoryLevelNum + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)

fit_simpson_history_factor <- glmmTMB::glmmTMB(
  Simpson__Genus_norm ~ HistoryLevel + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)
fit_simpson_history_num <- glmmTMB::glmmTMB(
  Simpson__Genus_norm ~ HistoryLevelNum + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)

fit_richness_history_factor <- glmmTMB::glmmTMB(
  Richness__Genus_norm ~ HistoryLevel + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)
fit_richness_history_num <- glmmTMB::glmmTMB(
  Richness__Genus_norm ~ HistoryLevelNum + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)

fit_shannon_parasite <- glmmTMB::glmmTMB(
  Shannon__Genus_norm ~ HistoryLevelNum * Parasite + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)
fit_simpson_parasite <- glmmTMB::glmmTMB(
  Simpson__Genus_norm ~ HistoryLevelNum * Parasite + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)
fit_richness_parasite <- glmmTMB::glmmTMB(
  Richness__Genus_norm ~ HistoryLevelNum * Parasite + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)

# --- Factorial exposure-regime models (A × T × P) ------------------------------
# Use sum-to-zero contrasts for type-III Wald tests on factorial terms.
.old_contr <- options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))
on.exit(options(.old_contr), add = TRUE)

set.seed(42)
fit_shannon_atp <- glmmTMB::glmmTMB(
  Shannon__Genus_norm ~ Antibiotics_f * Temperature_f * Parasite_f + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)

set.seed(42)
fit_simpson_atp <- glmmTMB::glmmTMB(
  Simpson__Genus_norm ~ Antibiotics_f * Temperature_f * Parasite_f + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)

set.seed(42)
fit_richness_atp <- glmmTMB::glmmTMB(
  Richness__Genus_norm ~ Antibiotics_f * Temperature_f * Parasite_f + (1 | Tank.ID),
  data = alpha_data_model,
  family = glmmTMB::beta_family(link = "logit")
)

glmmtmb_type3_terms_tidy <- function(fit_obj, metric_label) {
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Missing package `car` required for type-III Wald tests (car::Anova).")
  }
  a3 <- car::Anova(fit_obj, type = 3)
  as.data.frame(a3) %>%
    tibble::rownames_to_column(var = "Term") %>%
    dplyr::rename(
      Chisq = "Chisq",
      Df = "Df",
      p_value = "Pr(>Chisq)"
    ) %>%
    dplyr::filter(.data$Term != "(Intercept)") %>%
    dplyr::mutate(
      Metric = metric_label,
      p_value = as.numeric(.data$p_value),
      p_value_round = round(.data$p_value, 4),
      p_signif = dplyr::case_when(
        .data$p_value < 0.001 ~ "***",
        .data$p_value < 0.01 ~ "**",
        .data$p_value < 0.05 ~ "*",
        .data$p_value < 0.1 ~ ".",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::select(Metric, Term, Chisq, Df, p_value = p_value_round, p_signif)
}

model_subtitle_atp_terms <- paste(
  "Type-III Wald tests from beta regression GLMMs (glmmTMB; Beta family, logit link; random intercept per tank):",
  "y ~ Antibiotics * Temperature * Parasite + (1 | Tank.ID).",
  "Antibiotics/Temperature/Parasite coded as 0/1 factors; sum-to-zero contrasts for type-III tests."
)

tbl_atp_terms <- dplyr::bind_rows(
  glmmtmb_type3_terms_tidy(fit_shannon_atp, "Shannon diversity"),
  glmmtmb_type3_terms_tidy(fit_simpson_atp, "Simpson diversity"),
  glmmtmb_type3_terms_tidy(fit_richness_atp, "Observed richness")
)

readr::write_csv(tbl_atp_terms, file.path(path_tbl, "alpha_factorial_ATP_terms__type3_wald.csv"))

gt_atp_terms <- tbl_atp_terms %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Exposure regimes (factorial): alpha diversity vs. Antibiotics × Temperature × Parasite",
    subtitle = model_subtitle_atp_terms
  ) %>%
  gt::cols_label(
    Metric = "Metric",
    Term = "Term",
    Chisq = "Chi-square",
    Df = "df",
    p_value = "p-value",
    p_signif = "Sig."
  ) %>%
  gt::tab_footnote(
    footnote = "Significance codes: *** <0.001, ** <0.01, * <0.05, . <0.1, ns otherwise.",
    locations = gt::cells_column_labels(columns = p_signif)
  )
gt_atp_terms <- highlight_gt_significance(gt_atp_terms, tbl_atp_terms)
gt::gtsave(gt_atp_terms, file.path(path_tbl, "alpha_factorial_ATP_terms__type3_wald.html"))

alpha_factorial_cell_plot <- function(fit_obj, y_lab, title_text) {
  emm <- emmeans::emmeans(
    fit_obj,
    ~ Antibiotics_f * Temperature_f * Parasite_f,
    type = "response"
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(
      AT = paste0(.data$Antibiotics_f, " ", .data$Temperature_f),
      AT = factor(.data$AT, levels = c("A- T-", "A+ T-", "A- T+", "A+ T+")),
      Parasite_f = factor(.data$Parasite_f, levels = c("Unexposed", "Exposed"))
    )

  ggplot2::ggplot(
    emm,
    ggplot2::aes(x = .data$AT, y = .data$response, color = .data$Parasite_f)
  ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.35), size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$asymp.LCL, ymax = .data$asymp.UCL),
      width = 0.15,
      position = ggplot2::position_dodge(width = 0.35),
      linewidth = SIELER2026_MIN_LINEWIDTH_MM
    ) +
    ggplot2::scale_color_manual(
      values = c("Unexposed" = "grey60", "Exposed" = "firebrick"),
      name = "Parasite Exposure"
    ) +
    ggplot2::labs(
      title = title_text,
      subtitle = "Estimated marginal means (response scale) ± 95% CI from GLMM: Antibiotics × Temperature × Parasite.",
      x = "Exposure regime (A × T)",
      y = y_lab,
      color = "Parasite Exposure"
    ) +
    theme_sieler2026_publication(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold")
    )
}

p_atp_shannon <- alpha_factorial_cell_plot(
  fit_shannon_atp,
  y_lab = "Shannon diversity (normalized)",
  title_text = "Alpha diversity across factorial exposure regimes (Shannon)"
)
p_atp_simpson <- alpha_factorial_cell_plot(
  fit_simpson_atp,
  y_lab = "Simpson diversity (normalized)",
  title_text = "Alpha diversity across factorial exposure regimes (Simpson)"
)
p_atp_richness <- alpha_factorial_cell_plot(
  fit_richness_atp,
  y_lab = "Observed richness (normalized)",
  title_text = "Alpha diversity across factorial exposure regimes (Richness)"
)

ggplot2::ggsave(file.path(path_fig, "alpha_factorial_ATP_shannon.pdf"), p_atp_shannon, width = 8, height = 6, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "alpha_factorial_ATP_shannon.png"), p_atp_shannon, width = 8, height = 6, dpi = 300)
ggplot2::ggsave(file.path(path_fig, "alpha_factorial_ATP_simpson.pdf"), p_atp_simpson, width = 8, height = 6, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "alpha_factorial_ATP_simpson.png"), p_atp_simpson, width = 8, height = 6, dpi = 300)
ggplot2::ggsave(file.path(path_fig, "alpha_factorial_ATP_richness.pdf"), p_atp_richness, width = 8, height = 6, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "alpha_factorial_ATP_richness.png"), p_atp_richness, width = 8, height = 6, dpi = 300)

# --- Linear trend contrasts (stress history) ----------------------------------
trend_shannon <- alpha_diversity_linear_trend_contrast(fit_shannon_history_num)
trend_simpson <- alpha_diversity_linear_trend_contrast(fit_simpson_history_num)
trend_richness <- alpha_diversity_linear_trend_contrast(fit_richness_history_num)

trend_shannon_summ <- summary(trend_shannon, infer = c(TRUE, TRUE)) %>% as.data.frame()
trend_simpson_summ <- summary(trend_simpson, infer = c(TRUE, TRUE)) %>% as.data.frame()
trend_richness_summ <- summary(trend_richness, infer = c(TRUE, TRUE)) %>% as.data.frame()

combined_trend_results_data <- dplyr::bind_rows(
  trend_shannon_summ %>% dplyr::filter(contrast == "linear") %>% dplyr::mutate(Metric = "Shannon diversity"),
  trend_simpson_summ %>% dplyr::filter(contrast == "linear") %>% dplyr::mutate(Metric = "Simpson diversity"),
  trend_richness_summ %>% dplyr::filter(contrast == "linear") %>% dplyr::mutate(Metric = "Observed richness")
) %>%
  dplyr::transmute(
    Metric = Metric,
    Contrast = "Linear trend (0 → 2)",
    Estimate = round(estimate, 3),
    SE = round(SE, 3),
    df = round(df, 1),
    z_ratio = round(z.ratio, 2),
    p_value = round(p.value, 3),
    p_signif = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

readr::write_csv(combined_trend_results_data, file.path(path_tbl, "combined_diversity_trends_stress_history.csv"))

combined_trend_gt <- combined_trend_results_data %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Linear trend: alpha diversity vs. prior stressor history",
    subtitle = model_subtitle_linear_trend
  ) %>%
  gt::cols_label(
    Metric = "Metric",
    Contrast = "Contrast",
    Estimate = "Estimate",
    SE = "SE",
    df = "df",
    z_ratio = "z",
    p_value = "p-value",
    p_signif = "Sig."
  ) %>%
  gt::tab_footnote(
    footnote = "Significance codes: *** <0.001, ** <0.01, * <0.05, . <0.1, ns otherwise.",
    locations = gt::cells_column_labels(columns = p_signif)
  )
combined_trend_gt <- highlight_gt_significance(combined_trend_gt, combined_trend_results_data)

gt::gtsave(combined_trend_gt, file.path(path_tbl, "combined_diversity_trends_stress_history.html"))

p_int_shannon <- glmmtmb_history_parasite_interaction_p(fit_shannon_parasite)
p_int_simpson <- glmmtmb_history_parasite_interaction_p(fit_simpson_parasite)
p_int_richness <- glmmtmb_history_parasite_interaction_p(fit_richness_parasite)

# --- Table B: History × parasite (interaction coefficient; link scale) ----------
interaction_coef_rows <- dplyr::bind_rows(
  glmmtmb_history_parasite_interaction_coef_row(fit_shannon_parasite) %>%
    dplyr::mutate(Metric = "Shannon diversity"),
  glmmtmb_history_parasite_interaction_coef_row(fit_simpson_parasite) %>%
    dplyr::mutate(Metric = "Simpson diversity"),
  glmmtmb_history_parasite_interaction_coef_row(fit_richness_parasite) %>%
    dplyr::mutate(Metric = "Observed richness")
)

combined_interaction_results_data <- interaction_coef_rows %>%
  dplyr::transmute(
    Metric = Metric,
    Term = Term,
    Estimate = round(Estimate, 4),
    SE = round(SE, 4),
    z_ratio = round(z, 3),
    p_value = round(p.value, 4),
    p_signif = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  )

readr::write_csv(combined_interaction_results_data, file.path(path_tbl, "combined_diversity_interaction_parasite.csv"))

combined_interaction_gt <- combined_interaction_results_data %>%
  gt::gt() %>%
  gt::tab_header(
    title = "History × parasite interaction: alpha diversity (link scale)",
    subtitle = model_subtitle_interaction
  ) %>%
  gt::cols_label(
    Metric = "Metric",
    Term = "Term",
    Estimate = "Estimate",
    SE = "SE",
    z_ratio = "z",
    p_value = "p-value",
    p_signif = "Sig."
  ) %>%
  gt::tab_footnote(
    footnote = "Significance codes: *** <0.001, ** <0.01, * <0.05, . <0.1, ns otherwise.",
    locations = gt::cells_column_labels(columns = p_signif)
  )
combined_interaction_gt <- highlight_gt_significance(combined_interaction_gt, combined_interaction_results_data)
gt::gtsave(combined_interaction_gt, file.path(path_tbl, "combined_diversity_interaction_parasite.html"))

# --- Table C: p-value summary only ----------------------------------------------
combined_table2_data <- dplyr::tibble(
  Metric = combined_trend_results_data$Metric,
  stress_history_trend_p = combined_trend_results_data$p_value,
  history_x_parasite_p = round(
    c(p_int_shannon, p_int_simpson, p_int_richness),
    4
  )
)

readr::write_csv(combined_table2_data, file.path(path_tbl, "combined_diversity_table2.csv"))

combined_table2_gt <- combined_table2_data %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Alpha diversity: p-value summary (stress history vs. history × parasite)",
    subtitle = model_subtitle_table2
  ) %>%
  gt::cols_label(
    Metric = "Metric",
    stress_history_trend_p = "Stress-history trend p",
    history_x_parasite_p = "History × parasite p"
  )
combined_table2_gt <- highlight_gt_significance(combined_table2_gt, combined_table2_data, alpha = 0.05)
gt::gtsave(combined_table2_gt, file.path(path_tbl, "combined_diversity_table2.html"))

# --- Table D: Parasite simple effects within each stressor history (emmeans) ---
parasite_simple_raw <- dplyr::bind_rows(
  alpha_diversity_parasite_simple_effects(fit_shannon_parasite, "Shannon diversity"),
  alpha_diversity_parasite_simple_effects(fit_simpson_parasite, "Simpson diversity"),
  alpha_diversity_parasite_simple_effects(fit_richness_parasite, "Observed richness")
)

readr::write_csv(parasite_simple_raw, file.path(path_tbl, "alpha_parasite_within_history__simple_effects.csv"))

parasite_simple_display <- parasite_simple_raw %>%
  dplyr::mutate(
    estimate = round(.data$estimate, 4),
    SE = round(.data$SE, 4),
    df = round(.data$df, 4),
    z_or_t_ratio = round(.data$z_or_t_ratio, 3),
    p_value = round(.data$p.value, 4),
    p_fdr_bh = round(.data$p_fdr_bh, 4),
    p_signif = dplyr::case_when(
      .data$p.value < 0.001 ~ "***",
      .data$p.value < 0.01 ~ "**",
      .data$p.value < 0.05 ~ "*",
      .data$p.value < 0.1 ~ ".",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::select(-dplyr::any_of(c("p.value")))

parasite_simple_gt <- parasite_simple_display %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Parasite exposure vs. unexposed within each prior stressor history (simple effects)",
    subtitle = model_subtitle_parasite_simple
  ) %>%
  gt::cols_label(
    Metric = "Metric",
    HistoryLevelNum = "Prior stressors",
    contrast = "Contrast",
    estimate = "Estimate",
    SE = "SE",
    df = "df",
    z_or_t_ratio = "z / t",
    p_value = "p-value",
    p_fdr_bh = "p (FDR-BH)",
    p_signif = "Sig."
  ) %>%
  gt::tab_footnote(
    footnote = "Contrast is Parasite 1 vs 0. emmeans may report the effect as odds.ratio on the response scale for beta models; z/t and p apply to that contrast. FDR-BH is across the three history strata within each metric.",
    locations = gt::cells_column_labels(columns = p_fdr_bh)
  ) %>%
  gt::tab_footnote(
    footnote = "Significance codes: *** <0.001, ** <0.01, * <0.05, . <0.1, ns otherwise.",
    locations = gt::cells_column_labels(columns = p_signif)
  )
parasite_simple_gt <- highlight_gt_significance(parasite_simple_gt, parasite_simple_display)
gt::gtsave(parasite_simple_gt, file.path(path_tbl, "alpha_parasite_within_history__simple_effects.html"))

# --- Figures: predicted vs observed trend lines --------------------------------
p_shannon <- alpha_diversity_stress_history_trend_plot(
  alpha_data_model,
  fit_shannon_history_num,
  "Shannon__Genus_norm",
  "Shannon diversity (normalized)",
  "Cumulative stress exposure effects on Shannon diversity",
  paste(
    "Beta GLMM: Shannon__Genus_norm ~ HistoryLevelNum + (1 | Tank.ID).",
    "Predicted means with 95% CI.",
    sep = "\n"
  )
)
p_simpson <- alpha_diversity_stress_history_trend_plot(
  alpha_data_model,
  fit_simpson_history_num,
  "Simpson__Genus_norm",
  "Simpson diversity (normalized)",
  "Cumulative stress exposure effects on Simpson diversity",
  paste(
    "Beta GLMM: Simpson__Genus_norm ~ HistoryLevelNum + (1 | Tank.ID).",
    "Predicted means with 95% CI.",
    sep = "\n"
  )
)
p_richness <- alpha_diversity_stress_history_trend_plot(
  alpha_data_model,
  fit_richness_history_num,
  "Richness__Genus_norm",
  "Observed richness (normalized)",
  "Cumulative stress exposure effects on observed richness",
  paste(
    "Beta GLMM: Richness__Genus_norm ~ HistoryLevelNum + (1 | Tank.ID).",
    "Predicted means with 95% CI.",
    sep = "\n"
  )
)

# Quasirandom observed layout: no subtitle/caption; parallel figure filenames.
p_shannon_qr <- alpha_diversity_stress_history_trend_plot(
  alpha_data_model,
  fit_shannon_history_num,
  "Shannon__Genus_norm",
  "Shannon diversity (normalized)",
  "Cumulative stress exposure effects on Shannon diversity",
  subtitle = NULL,
  caption = NULL,
  observed_layout = "quasirandom"
)
p_simpson_qr <- alpha_diversity_stress_history_trend_plot(
  alpha_data_model,
  fit_simpson_history_num,
  "Simpson__Genus_norm",
  "Simpson diversity (normalized)",
  "Cumulative stress exposure effects on Simpson diversity",
  subtitle = NULL,
  caption = NULL,
  observed_layout = "quasirandom"
) +
  theme_sieler2026_trend_panel_grid_major_y()
p_richness_qr <- alpha_diversity_stress_history_trend_plot(
  alpha_data_model,
  fit_richness_history_num,
  "Richness__Genus_norm",
  "Observed richness (normalized)",
  "Cumulative stress exposure effects on observed richness",
  subtitle = NULL,
  caption = NULL,
  observed_layout = "quasirandom"
)

fig_w <- 8
fig_h <- 8
# Main-text Simpson quasirandom (panel 2.2): same inches as Mort/Inf trend panels (Fig 4.1–4.2)
# so ggplot point sizes (mm) match visual scale across drivers (`main_figures_manifest.csv`).
fig_main_trend_in <- 6
ggplot2::ggsave(file.path(path_fig, "shannon_diversity_trend.pdf"), p_shannon, width = fig_w, height = fig_h, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "shannon_diversity_trend.png"), p_shannon, width = fig_w, height = fig_h, dpi = 300)
ggplot2::ggsave(file.path(path_fig, "simpson_diversity_trend.pdf"), p_simpson, width = fig_w, height = fig_h, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "simpson_diversity_trend.png"), p_simpson, width = fig_w, height = fig_h, dpi = 300)
ggplot2::ggsave(file.path(path_fig, "observed_diversity_trend.pdf"), p_richness, width = fig_w, height = fig_h, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "observed_diversity_trend.png"), p_richness, width = fig_w, height = fig_h, dpi = 300)
ggplot2::ggsave(
  file.path(path_fig, "shannon_diversity_trend_quasirandom.pdf"),
  p_shannon_qr,
  width = fig_w,
  height = fig_h,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "shannon_diversity_trend_quasirandom.png"),
  p_shannon_qr,
  width = fig_w,
  height = fig_h,
  dpi = 300
)
ggplot2::ggsave(
  file.path(path_fig, "simpson_diversity_trend_quasirandom.pdf"),
  p_simpson_qr,
  width = fig_main_trend_in,
  height = fig_main_trend_in,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "simpson_diversity_trend_quasirandom.png"),
  p_simpson_qr,
  width = fig_main_trend_in,
  height = fig_main_trend_in,
  dpi = 300
)
ggplot2::ggsave(
  file.path(path_fig, "observed_diversity_trend_quasirandom.pdf"),
  p_richness_qr,
  width = fig_w,
  height = fig_h,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "observed_diversity_trend_quasirandom.png"),
  p_richness_qr,
  width = fig_w,
  height = fig_h,
  dpi = 300
)

# --- Figures: parasite exposure within each prior stressor history -------------
fig_pw_w <- 12
fig_pw_h <- 5
p_parasite_shannon <- alpha_diversity_parasite_within_history_plot(
  alpha_data_model,
  "Shannon__Genus_norm",
  "Shannon diversity (normalized)",
  "Alpha diversity by parasite exposure within prior stressor history",
  paste(
    "Per-fish normalized values; facets = 0, 1, or 2 prior stressors.",
    "Same statistical model as interaction GLMM (Table D: simple effects).",
    sep = "\n"
  ),
  p_value_tbl = parasite_simple_raw %>%
    dplyr::filter(.data$Metric == "Shannon diversity") %>%
    dplyr::select("HistoryLevelNum", "p.value")
)
p_parasite_simpson <- alpha_diversity_parasite_within_history_plot(
  alpha_data_model,
  "Simpson__Genus_norm",
  "Simpson diversity (normalized)",
  "Alpha diversity by parasite exposure within prior stressor history",
  paste(
    "Per-fish normalized values; facets = 0, 1, or 2 prior stressors.",
    "Same statistical model as interaction GLMM (Table D: simple effects).",
    sep = "\n"
  ),
  p_value_tbl = parasite_simple_raw %>%
    dplyr::filter(.data$Metric == "Simpson diversity") %>%
    dplyr::select("HistoryLevelNum", "p.value")
)
p_parasite_richness <- alpha_diversity_parasite_within_history_plot(
  alpha_data_model,
  "Richness__Genus_norm",
  "Observed richness (normalized)",
  "Alpha diversity by parasite exposure within prior stressor history",
  paste(
    "Per-fish normalized values; facets = 0, 1, or 2 prior stressors.",
    "Same statistical model as interaction GLMM (Table D: simple effects).",
    sep = "\n"
  ),
  p_value_tbl = parasite_simple_raw %>%
    dplyr::filter(.data$Metric == "Observed richness") %>%
    dplyr::select("HistoryLevelNum", "p.value")
)
ggplot2::ggsave(
  file.path(path_fig, "shannon_parasite_within_history.pdf"),
  p_parasite_shannon,
  width = fig_pw_w,
  height = fig_pw_h,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "shannon_parasite_within_history.png"),
  p_parasite_shannon,
  width = fig_pw_w,
  height = fig_pw_h,
  dpi = 300
)
ggplot2::ggsave(
  file.path(path_fig, "simpson_parasite_within_history.pdf"),
  p_parasite_simpson,
  width = fig_pw_w,
  height = fig_pw_h,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "simpson_parasite_within_history.png"),
  p_parasite_simpson,
  width = fig_pw_w,
  height = fig_pw_h,
  dpi = 300
)
ggplot2::ggsave(
  file.path(path_fig, "richness_parasite_within_history.pdf"),
  p_parasite_richness,
  width = fig_pw_w,
  height = fig_pw_h,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "richness_parasite_within_history.png"),
  p_parasite_richness,
  width = fig_pw_w,
  height = fig_pw_h,
  dpi = 300
)

# --- Assemble bundle -----------------------------------------------------------
modules <- list(
  exposure_regime_factorial = list(
    Shannon = list(
      stats = list(fit = fit_shannon_atp),
      figure = p_atp_shannon,
      paths = c(
        figure_pdf = file.path(rel_fig, "alpha_factorial_ATP_shannon.pdf"),
        figure_png = file.path(rel_fig, "alpha_factorial_ATP_shannon.png")
      )
    ),
    Simpson = list(
      stats = list(fit = fit_simpson_atp),
      figure = p_atp_simpson,
      paths = c(
        figure_pdf = file.path(rel_fig, "alpha_factorial_ATP_simpson.pdf"),
        figure_png = file.path(rel_fig, "alpha_factorial_ATP_simpson.png")
      )
    ),
    Richness = list(
      stats = list(fit = fit_richness_atp),
      figure = p_atp_richness,
      paths = c(
        figure_pdf = file.path(rel_fig, "alpha_factorial_ATP_richness.pdf"),
        figure_png = file.path(rel_fig, "alpha_factorial_ATP_richness.png")
      )
    )
  ),
  stress_history = list(
    Shannon = list(
      stats = list(
        fit_factor = fit_shannon_history_factor,
        fit_num = fit_shannon_history_num,
        trend_contrast = trend_shannon
      ),
      figure = p_shannon,
      paths = c(
        figure_pdf = file.path(rel_fig, "shannon_diversity_trend.pdf"),
        figure_png = file.path(rel_fig, "shannon_diversity_trend.png")
      ),
      figure_quasirandom = p_shannon_qr,
      paths_quasirandom = c(
        figure_pdf = file.path(rel_fig, "shannon_diversity_trend_quasirandom.pdf"),
        figure_png = file.path(rel_fig, "shannon_diversity_trend_quasirandom.png")
      )
    ),
    Simpson = list(
      stats = list(
        fit_factor = fit_simpson_history_factor,
        fit_num = fit_simpson_history_num,
        trend_contrast = trend_simpson
      ),
      figure = p_simpson,
      paths = c(
        figure_pdf = file.path(rel_fig, "simpson_diversity_trend.pdf"),
        figure_png = file.path(rel_fig, "simpson_diversity_trend.png")
      ),
      figure_quasirandom = p_simpson_qr,
      paths_quasirandom = c(
        figure_pdf = file.path(rel_fig, "simpson_diversity_trend_quasirandom.pdf"),
        figure_png = file.path(rel_fig, "simpson_diversity_trend_quasirandom.png")
      )
    ),
    Richness = list(
      stats = list(
        fit_factor = fit_richness_history_factor,
        fit_num = fit_richness_history_num,
        trend_contrast = trend_richness
      ),
      figure = p_richness,
      paths = c(
        figure_pdf = file.path(rel_fig, "observed_diversity_trend.pdf"),
        figure_png = file.path(rel_fig, "observed_diversity_trend.png")
      ),
      figure_quasirandom = p_richness_qr,
      paths_quasirandom = c(
        figure_pdf = file.path(rel_fig, "observed_diversity_trend_quasirandom.pdf"),
        figure_png = file.path(rel_fig, "observed_diversity_trend_quasirandom.png")
      )
    )
  ),
  parasite_interaction = list(
    Shannon = list(
      stats = list(fit = fit_shannon_parasite, interaction_p = p_int_shannon)
    ),
    Simpson = list(
      stats = list(fit = fit_simpson_parasite, interaction_p = p_int_simpson)
    ),
    Richness = list(
      stats = list(fit = fit_richness_parasite, interaction_p = p_int_richness)
    )
  ),
  parasite_within_history = list(
    Shannon = list(
      stats = list(simple_effects_tidy = parasite_simple_raw %>% dplyr::filter(.data$Metric == "Shannon diversity")),
      figure = p_parasite_shannon,
      paths = c(
        figure_pdf = file.path(rel_fig, "shannon_parasite_within_history.pdf"),
        figure_png = file.path(rel_fig, "shannon_parasite_within_history.png")
      )
    ),
    Simpson = list(
      stats = list(simple_effects_tidy = parasite_simple_raw %>% dplyr::filter(.data$Metric == "Simpson diversity")),
      figure = p_parasite_simpson,
      paths = c(
        figure_pdf = file.path(rel_fig, "simpson_parasite_within_history.pdf"),
        figure_png = file.path(rel_fig, "simpson_parasite_within_history.png")
      )
    ),
    Richness = list(
      stats = list(simple_effects_tidy = parasite_simple_raw %>% dplyr::filter(.data$Metric == "Observed richness")),
      figure = p_parasite_richness,
      paths = c(
        figure_pdf = file.path(rel_fig, "richness_parasite_within_history.pdf"),
        figure_png = file.path(rel_fig, "richness_parasite_within_history.png")
      )
    )
  )
)

diversity_bundle <- list(
  meta = list(
    run_date = as.character(Sys.Date()),
    ps_list_element = ps_list_element,
    bundle_version = "1.5",
    script = "Code/01__Analysis/01__Diversity.R",
    model_factorial_atp_terms = model_subtitle_atp_terms,
    model_linear_trend = model_subtitle_linear_trend,
    model_table2 = model_subtitle_table2,
    model_interaction = model_subtitle_interaction,
    model_parasite_simple_effects = model_subtitle_parasite_simple,
    model_formulas = c(
      "Factorial exposure regimes (A×T×P)" = "y ~ Antibiotics * Temperature * Parasite + (1 | Tank.ID)",
      "Stress history (numeric)" = "y ~ HistoryLevelNum + (1 | Tank.ID)",
      "Stress history (factor)" = "y ~ HistoryLevel + (1 | Tank.ID)",
      "History × parasite" = "y ~ HistoryLevelNum * Parasite + (1 | Tank.ID)",
      "Parasite simple effects" = "emmeans pairwise Parasite | HistoryLevelNum (same interaction model)"
    )
  ),
  table_factorial_atp_terms = gt_atp_terms,
  table_factorial_atp_terms_tidy = tbl_atp_terms,
  table_combined_table2 = combined_table2_gt,
  table_combined_trends = combined_trend_gt,
  table_combined_interaction = combined_interaction_gt,
  table_parasite_within_history = parasite_simple_gt,
  tables_tidy = list(
    factorial_atp_terms = tbl_atp_terms,
    combined_table2 = combined_table2_data,
    combined_trends = combined_trend_results_data,
    combined_interaction = combined_interaction_results_data,
    parasite_within_history = parasite_simple_display
  ),
  modules = modules
)

saveRDS(diversity_bundle, bundle_rds)
message("Saved diversity bundle: ", bundle_rds)

sieler2026_sync_main_figures_from_manifest(driver_script = "01__Diversity.R")
