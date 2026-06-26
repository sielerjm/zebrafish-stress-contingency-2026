# 05__Mort-Inf.R
# Created by: Michael Sieler
# Date last updated: 2026-04-27
#
# Description: Day-60 mortality and infection (P+ survivors at Day 60) vs prior stressor history.
#   Mortality: binomial GLMM with HistoryLevelNum * Parasite + (1|Tank); infection: Day-60
#   survivor prevalence (unique Sample), tank-level binomial GLMMs, Treatment pairwise (FDR).
#   Exports trend figures (wide + square), bar charts (incl. burden by prior history), CSV tables,
#   Stats RDS, mortinf__host__bundle.rds.
#
# Expected input:  Run from Sieler2026 root; data.list from Data/r_objects/data-list__*.rds
#   (04__DataPreProcess.R). Uses Mortality and Infection_Tank elements.
# Expected output:  Results/05__Mort-Inf/{Figures,Tables,Stats} and mortinf__host__bundle.rds.
#   Optional: Rscript .../05__Mort-Inf.R --figures-only re-ggsaves manifest Fig 4 panels from
#   Stats/mortinf_main_text_figure_ggplots.rds (skips archive and the full modeling pipeline).

init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
if (!file.exists(init_rel)) {
  stop(
    "Set working directory to the Sieler2026 repository root, or run:\n",
    "  Rscript Code/01__Analysis/05__Mort-Inf.R\n",
    "from the project root so `", init_rel, "` resolves."
  )
}
source(init_rel)

args_05 <- commandArgs(trailingOnly = TRUE)
figures_only_05 <- isTRUE(any(args_05 == "--figures-only"))

if (!isTRUE(figures_only_05)) {
  if (!exists("data.list", inherits = TRUE) || is.null(data.list)) {
    stop(
      "data.list not found. Run Code/00__Setup/04__DataPreProcess.R and ensure ",
      "data-list__*.rds exists under Data/r_objects/."
    )
  }
}

path_res <- file.path(path.results, "05__Mort-Inf")
path_fig <- file.path(path_res, "Figures")
path_tbl <- file.path(path_res, "Tables")
path_stats <- file.path(path_res, "Stats")
dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(path_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(path_stats, recursive = TRUE, showWarnings = FALSE)

fig_in <- 6
fig_bar_sq <- 6

if (isTRUE(figures_only_05)) {
  chk_gg <- file.path(path_stats, "mortinf_main_text_figure_ggplots.rds")
  if (!file.exists(chk_gg)) {
    stop(
      "05__Mort-Inf.R --figures-only requires saved ggplot checkpoints:\n  ", chk_gg,
      "\nRun a full 05 driver once to create them."
    )
  }
  pl <- readRDS(chk_gg)
  req <- c(
    "p_mortality_exposure_by_prior_history",
    "p_infection_prevalence_exposure_by_prior_history",
    "p_mortality_prior_stressor_trend",
    "p_mortality_prior_stressor_trend_qr",
    "p_infection_prevalence_trend_qr"
  )
  miss <- req[vapply(req, function(nm) is.null(pl[[nm]]), logical(1L))]
  if (length(miss) > 0L) {
    stop("05 --figures-only: mortinf_main_text_figure_ggplots.rds is missing elements: ", paste(miss, collapse = ", "))
  }
  ggplot2::ggsave(
    file.path(path_fig, "mortality_percent_by_exposure_regime_prior_history_square.pdf"),
    pl$p_mortality_exposure_by_prior_history,
    width = fig_bar_sq,
    height = fig_bar_sq,
    device = "pdf"
  )
  ggplot2::ggsave(
    file.path(path_fig, "mortality_percent_by_exposure_regime_prior_history_square.png"),
    pl$p_mortality_exposure_by_prior_history,
    width = fig_bar_sq,
    height = fig_bar_sq,
    dpi = 300L
  )
  ggplot2::ggsave(
    file.path(path_fig, "infection_prevalence_by_exposure_regime_prior_history_square.pdf"),
    pl$p_infection_prevalence_exposure_by_prior_history,
    width = fig_bar_sq,
    height = fig_bar_sq,
    device = "pdf"
  )
  ggplot2::ggsave(
    file.path(path_fig, "infection_prevalence_by_exposure_regime_prior_history_square.png"),
    pl$p_infection_prevalence_exposure_by_prior_history,
    width = fig_bar_sq,
    height = fig_bar_sq,
    dpi = 300L
  )
  ggplot2::ggsave(
    file.path(path_fig, "mortality_prior_stressor_trend_square.pdf"),
    pl$p_mortality_prior_stressor_trend,
    width = fig_in,
    height = fig_in,
    device = "pdf"
  )
  ggplot2::ggsave(
    file.path(path_fig, "mortality_prior_stressor_trend_square.png"),
    pl$p_mortality_prior_stressor_trend,
    width = fig_in,
    height = fig_in,
    dpi = 300L
  )
  ggplot2::ggsave(
    file.path(path_fig, "mortality_prior_stressor_trend_square_quasirandom.pdf"),
    pl$p_mortality_prior_stressor_trend_qr,
    width = fig_in,
    height = fig_in,
    device = "pdf"
  )
  ggplot2::ggsave(
    file.path(path_fig, "mortality_prior_stressor_trend_square_quasirandom.png"),
    pl$p_mortality_prior_stressor_trend_qr,
    width = fig_in,
    height = fig_in,
    dpi = 300L
  )
  ggplot2::ggsave(
    file.path(path_fig, "infection_prevalence_trend_predicted_square_quasirandom.pdf"),
    pl$p_infection_prevalence_trend_qr,
    width = fig_in,
    height = fig_in,
    device = "pdf"
  )
  ggplot2::ggsave(
    file.path(path_fig, "infection_prevalence_trend_predicted_square_quasirandom.png"),
    pl$p_infection_prevalence_trend_qr,
    width = fig_in,
    height = fig_in,
    dpi = 300L
  )
  sieler2026_sync_main_figures_from_manifest(
    driver_script = "05__Mort-Inf.R",
    panel_ids = c("4.1", "4.2")
  )
  message("05__Mort-Inf.R --figures-only complete.")
  quit(save = "no", status = 0L)
}

# Archive previous outputs (copy → verify → clear) so each run is clean but reproducible.
sieler2026_archive_module_outputs(
  path_res_module = path_res,
  module_name = "05__Mort-Inf",
  subdirs = c("Figures", "Tables", "Stats"),
  clear_original = TRUE
)

bundle_rds <- file.path(path_stats, "mortinf__host__bundle.rds")

# --- Tank-level mortality (Day 60) -------------------------------------------------
tmp_mort_tank <- data.list[["Mortality"]] %>%
  dplyr::filter(.data$Time == 60L) %>%
  dplyr::group_by(.data$Treatment, .data$Tank.ID) %>%
  dplyr::summarise(
    n = dplyr::n(),
    HistoryLevel = dplyr::first(.data$HistoryLevel),
    HistoryLevelNum = dplyr::first(.data$HistoryLevelNum),
    Parasite = dplyr::first(.data$Parasite),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    HistoryLevel = dplyr::case_when(
      .data$HistoryLevelNum == 0 ~ "No prior stressors",
      .data$HistoryLevelNum == 1 ~ "One prior stressor",
      .data$HistoryLevelNum == 2 ~ "Two prior stressors",
      TRUE ~ as.character(.data$HistoryLevel)
    ),
    HistoryLevel = factor(.data$HistoryLevel, levels = history_order),
    alive = .data$n,
    at_risk = 15L,
    dead = .data$at_risk - .data$alive,
    percent_mortality = round((.data$dead / .data$at_risk) * 100, 1),
    Treatment = factor(.data$Treatment, levels = treatment_order)
  ) %>%
  dplyr::select(-"n") %>%
  dplyr::relocate("at_risk", .before = "alive") %>%
  dplyr::relocate("dead", .after = "alive")

set.seed(42)
fit_mort_history_num <- glmmTMB::glmmTMB(
  cbind(dead, at_risk - dead) ~ HistoryLevelNum + (1 | Tank.ID),
  family = stats::binomial(link = "logit"),
  data = tmp_mort_tank
)

set.seed(42)
fit_mort_history_factor <- glmmTMB::glmmTMB(
  cbind(dead, at_risk - dead) ~ HistoryLevel + (1 | Tank.ID),
  family = stats::binomial(link = "logit"),
  data = tmp_mort_tank
)

emm_mort_history_num <- emmeans::emmeans(
  fit_mort_history_num,
  ~ HistoryLevelNum,
  at = list(HistoryLevelNum = c(0, 1, 2)),
  type = "response"
)

emm_mort_history_factor <- emmeans::emmeans(
  fit_mort_history_factor,
  ~ HistoryLevel,
  type = "response"
)

set.seed(42)
fit_mort_history_parasite <- glmmTMB::glmmTMB(
  cbind(dead, at_risk - dead) ~ HistoryLevelNum * Parasite + (1 | Tank.ID),
  family = stats::binomial(link = "logit"),
  data = tmp_mort_tank
)

set.seed(42)
fit_mort_history_factor_parasite <- glmmTMB::glmmTMB(
  cbind(dead, at_risk - dead) ~ HistoryLevel * Parasite + (1 | Tank.ID),
  family = stats::binomial(link = "logit"),
  data = tmp_mort_tank
)

emm_mort_parasite_resp <- emmeans::emmeans(
  fit_mort_history_parasite,
  ~ HistoryLevelNum | Parasite,
  at = list(HistoryLevelNum = c(0, 1, 2)),
  type = "response"
)

emm_mort_factor_parasite_resp <- emmeans::emmeans(
  fit_mort_history_factor_parasite,
  ~ HistoryLevel | Parasite,
  type = "response"
)

joint_tests_mortality <- emmeans::joint_tests(fit_mort_history_parasite)
readr::write_csv(
  as.data.frame(joint_tests_mortality),
  file.path(path_tbl, "mortality_glmm_joint_tests.csv")
)

joint_tests_mortality_factor <- emmeans::joint_tests(fit_mort_history_factor_parasite)
readr::write_csv(
  as.data.frame(joint_tests_mortality_factor),
  file.path(path_tbl, "mortality_glmm_joint_tests_factor.csv")
)

trend_mortality_history <- emmeans::contrast(emm_mort_history_num, "poly")
pairs_mort_history <- emmeans::contrast(emm_mort_history_num, "pairwise", adjust = "fdr")
df_mort_hist_poly <- as.data.frame(trend_mortality_history)
df_mort_hist_pairs <- as.data.frame(pairs_mort_history)
df_mort_hist_emmeans <- as.data.frame(
  summary(emm_mort_history_num, type = "response", infer = c(TRUE, TRUE))
)
readr::write_csv(df_mort_hist_poly, file.path(path_tbl, "mortality_glmm_history_poly_contrast.csv"))
readr::write_csv(df_mort_hist_pairs, file.path(path_tbl, "mortality_glmm_history_pairs_fdr.csv"))
readr::write_csv(df_mort_hist_emmeans, file.path(path_tbl, "mortality_glmm_history_emmeans_response.csv"))

df_mort_hist_factor_emmeans <- as.data.frame(
  summary(emm_mort_history_factor, type = "response", infer = c(TRUE, TRUE))
)
df_mort_hist_factor_pairs <- as.data.frame(
  emmeans::contrast(emm_mort_history_factor, "pairwise", adjust = "fdr")
)
readr::write_csv(
  df_mort_hist_factor_emmeans,
  file.path(path_tbl, "mortality_glmm_history_factor_emmeans_response.csv")
)
readr::write_csv(
  df_mort_hist_factor_pairs,
  file.path(path_tbl, "mortality_glmm_history_factor_pairs_fdr.csv")
)
tmp_mort_tank_obs <- tmp_mort_tank %>%
  dplyr::mutate(
    mortality_obs = .data$dead / .data$at_risk,
    parasite_label = dplyr::if_else(
      .data$Parasite == 1L,
      "Parasite-exposed (P+)",
      "Unexposed (P-)"
    )
  )

emm_mort_parasite_df <- emm_mort_parasite_resp %>%
  as.data.frame() %>%
  dplyr::mutate(
    HistoryLevelNum = as.numeric(as.character(.data$HistoryLevelNum)),
    parasite_label = dplyr::if_else(
      .data$Parasite == 1L,
      "Parasite-exposed (P+)",
      "Unexposed (P-)"
    )
  )

emm_mort_hist_df <- emm_mort_history_num %>%
  as.data.frame() %>%
  dplyr::mutate(HistoryLevelNum = as.numeric(as.character(.data$HistoryLevelNum)))

p_mortality_prior_stressor_levels <- df_mort_hist_factor_emmeans %>%
  dplyr::mutate(HistoryLevel = factor(.data$HistoryLevel, levels = history_order)) %>%
  ggplot2::ggplot(
    ggplot2::aes(x = .data$HistoryLevel, y = .data$prob, color = .data$HistoryLevel)
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$asymp.LCL, ymax = .data$asymp.UCL),
    width = 0.18,
    linewidth = SIELER2026_MIN_LINEWIDTH_MM
  ) +
  ggplot2::scale_color_manual(values = history_color_scale, guide = "none") +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = ggplot2::expansion(mult = c(0, 0)),
    oob = scales::squish
  ) +
  ggplot2::labs(
    x = "Prior stressor history",
    y = "Mortality probability",
    title = "Mortality by prior stressor history (levels)",
    subtitle = "Tank-level binomial GLMM; emmeans ± 95% CI (response scale)",
    caption = "Model: cbind(dead, alive) ~ HistoryLevel + (1 | Tank.ID)."
  ) +
  theme_sieler2026_publication(base_size = 14)

p_mortality_prior_stressor_levels_x_parasite <- emm_mort_factor_parasite_resp %>%
  as.data.frame() %>%
  dplyr::mutate(
    HistoryLevel = factor(.data$HistoryLevel, levels = history_order),
    parasite_label = dplyr::if_else(
      .data$Parasite == 1L,
      "Parasite-exposed (P+)",
      "Unexposed (P-)"
    )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(x = .data$HistoryLevel, y = .data$prob, color = .data$HistoryLevel)
  ) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.35), size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$asymp.LCL, ymax = .data$asymp.UCL),
    position = ggplot2::position_dodge(width = 0.35),
    width = 0.18,
    linewidth = SIELER2026_MIN_LINEWIDTH_MM
  ) +
  ggplot2::facet_wrap(ggplot2::vars(parasite_label), nrow = 1) +
  ggplot2::scale_color_manual(values = history_color_scale, name = "Prior stressor history") +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = ggplot2::expansion(mult = c(0, 0)),
    oob = scales::squish
  ) +
  ggplot2::labs(
    x = "Prior stressor history",
    y = "Mortality probability",
    title = "Mortality by prior stressor history × parasite exposure",
    subtitle = "Tank-level binomial GLMM; emmeans ± 95% CI (response scale)",
    caption = "Model: cbind(dead, alive) ~ HistoryLevel * Parasite + (1 | Tank.ID)."
  ) +
  theme_sieler2026_publication(base_size = 14, legend_position = "bottom")

p_mortality_prior_stressor_trend <- glmm_binomial_tank_history_numeric_trend_plot(
  obs_df = tmp_mort_tank_obs,
  obs_y_col = "mortality_obs",
  emm_df = emm_mort_hist_df,
  y_label = "Mortality (%)",
  title = "Final Mortality by Prior Stressor History",
  subtitle = NULL,
  y_as_percent = TRUE
)
p_mortality_prior_stressor_trend_qr <- glmm_binomial_tank_history_numeric_trend_plot(
  obs_df = tmp_mort_tank_obs,
  obs_y_col = "mortality_obs",
  emm_df = emm_mort_hist_df,
  y_label = "Mortality (%)",
  title = "Final Mortality by Prior Stressor History",
  subtitle = NULL,
  y_as_percent = TRUE,
  observed_layout = "quasirandom"
)

fig_mort_hist_pdf <- file.path(path_fig, "mortality_prior_stressor_trend.pdf")
fig_mort_hist_png <- file.path(path_fig, "mortality_prior_stressor_trend.png")
fig_mort_hist_square_pdf <- file.path(path_fig, "mortality_prior_stressor_trend_square.pdf")
fig_mort_hist_square_png <- file.path(path_fig, "mortality_prior_stressor_trend_square.png")
ggplot2::ggsave(fig_mort_hist_pdf, p_mortality_prior_stressor_trend, width = fig_in, height = fig_in, device = "pdf")
ggplot2::ggsave(fig_mort_hist_png, p_mortality_prior_stressor_trend, width = fig_in, height = fig_in, dpi = 300)
ggplot2::ggsave(fig_mort_hist_square_pdf, p_mortality_prior_stressor_trend, width = fig_in, height = fig_in, device = "pdf")
ggplot2::ggsave(fig_mort_hist_square_png, p_mortality_prior_stressor_trend, width = fig_in, height = fig_in, dpi = 300)
fig_mort_hist_qr_pdf <- file.path(path_fig, "mortality_prior_stressor_trend_quasirandom.pdf")
fig_mort_hist_qr_png <- file.path(path_fig, "mortality_prior_stressor_trend_quasirandom.png")
fig_mort_hist_square_qr_pdf <- file.path(path_fig, "mortality_prior_stressor_trend_square_quasirandom.pdf")
fig_mort_hist_square_qr_png <- file.path(path_fig, "mortality_prior_stressor_trend_square_quasirandom.png")
ggplot2::ggsave(fig_mort_hist_qr_pdf, p_mortality_prior_stressor_trend_qr, width = fig_in, height = fig_in, device = "pdf")
ggplot2::ggsave(fig_mort_hist_qr_png, p_mortality_prior_stressor_trend_qr, width = fig_in, height = fig_in, dpi = 300)
ggplot2::ggsave(fig_mort_hist_square_qr_pdf, p_mortality_prior_stressor_trend_qr, width = fig_in, height = fig_in, device = "pdf")
ggplot2::ggsave(fig_mort_hist_square_qr_png, p_mortality_prior_stressor_trend_qr, width = fig_in, height = fig_in, dpi = 300)

fig_mort_hist_levels_pdf <- file.path(path_fig, "mortality_by_prior_stressor_history_levels.pdf")
fig_mort_hist_levels_png <- file.path(path_fig, "mortality_by_prior_stressor_history_levels.png")
ggplot2::ggsave(fig_mort_hist_levels_pdf, p_mortality_prior_stressor_levels, width = 10, height = 6, device = "pdf")
ggplot2::ggsave(fig_mort_hist_levels_png, p_mortality_prior_stressor_levels, width = 10, height = 6, dpi = 300)

fig_mort_hist_levels_parasite_pdf <- file.path(
  path_fig,
  "mortality_by_prior_stressor_history_levels_x_parasite.pdf"
)
fig_mort_hist_levels_parasite_png <- file.path(
  path_fig,
  "mortality_by_prior_stressor_history_levels_x_parasite.png"
)
ggplot2::ggsave(
  fig_mort_hist_levels_parasite_pdf,
  p_mortality_prior_stressor_levels_x_parasite,
  width = 12,
  height = 6,
  device = "pdf"
)
ggplot2::ggsave(
  fig_mort_hist_levels_parasite_png,
  p_mortality_prior_stressor_levels_x_parasite,
  width = 12,
  height = 6,
  dpi = 300
)

p_mortality_parasite_trend_comparison <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = tmp_mort_tank_obs,
    ggplot2::aes(
      x = .data$HistoryLevelNum,
      y = .data$mortality_obs,
      color = .data$parasite_label
    ),
    alpha = 0.75,
    size = 2.2,
    position = ggplot2::position_jitter(width = 0.04, height = 0)
  ) +
  ggplot2::geom_line(
    data = emm_mort_parasite_df,
    ggplot2::aes(
      x = .data$HistoryLevelNum,
      y = .data$prob,
      color = .data$parasite_label,
      group = .data$parasite_label
    ),
    linewidth = 1.2
  ) +
  ggplot2::geom_point(
    data = emm_mort_parasite_df,
    ggplot2::aes(
      x = .data$HistoryLevelNum,
      y = .data$prob,
      color = .data$parasite_label,
      group = .data$parasite_label
    ),
    size = 4,
    shape = 18L
  ) +
  ggplot2::geom_ribbon(
    data = emm_mort_parasite_df,
    ggplot2::aes(
      x = .data$HistoryLevelNum,
      ymin = .data$asymp.LCL,
      ymax = .data$asymp.UCL,
      fill = .data$parasite_label,
      group = .data$parasite_label
    ),
    alpha = 0.2,
    color = NA
  ) +
  ggplot2::scale_x_continuous(breaks = c(0, 1, 2)) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = ggplot2::expansion(mult = c(0, 0)),
    oob = scales::squish
  ) +
  ggplot2::scale_color_manual(
    name = "Parasite challenge",
    values = c("Parasite-exposed (P+)" = "#E41A1C", "Unexposed (P-)" = "#377EB8")
  ) +
  ggplot2::scale_fill_manual(
    name = "Parasite challenge",
    values = c("Parasite-exposed (P+)" = "#E41A1C", "Unexposed (P-)" = "#377EB8")
  ) +
  ggplot2::labs(
    x = "Number of prior stressors",
    y = "Mortality probability",
    title = "Mortality by prior stressors: parasite vs unexposed",
    subtitle = "Predicted probabilities with 95% confidence intervals",
    caption = paste(
      "Tank-level observed points (jittered x); diamonds: marginal predicted means.\n",
      "Model: cbind(dead, alive) ~ HistoryLevelNum * Parasite + (1 | Tank.ID).",
      sep = ""
    )
  ) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::theme(legend.position = "bottom")

fig_mort_para_pdf <- file.path(path_fig, "mortality_parasite_trend_comparison.pdf")
fig_mort_para_png <- file.path(path_fig, "mortality_parasite_trend_comparison.png")
fig_mort_para_square_pdf <- file.path(path_fig, "mortality_parasite_trend_comparison_square.pdf")
fig_mort_para_square_png <- file.path(path_fig, "mortality_parasite_trend_comparison_square.png")
ggplot2::ggsave(fig_mort_para_pdf, p_mortality_parasite_trend_comparison, width = fig_in + 0.5, height = fig_in + 0.35, device = "pdf")
ggplot2::ggsave(fig_mort_para_png, p_mortality_parasite_trend_comparison, width = fig_in + 0.5, height = fig_in + 0.35, dpi = 300)
ggplot2::ggsave(fig_mort_para_square_pdf, p_mortality_parasite_trend_comparison, width = fig_in, height = fig_in, device = "pdf")
ggplot2::ggsave(fig_mort_para_square_png, p_mortality_parasite_trend_comparison, width = fig_in, height = fig_in, dpi = 300)

# --- Mortality by exposure regime (descriptive) ------------------------------------
mortality_treatment_tbl <- data.list[["Mortality"]] %>%
  dplyr::filter(.data$Time == 60) %>%
  dplyr::group_by(.data$Treatment) %>%
  dplyr::count() %>%
  dplyr::mutate(
    alive = .data$n,
    at_risk = 45L,
    dead = .data$at_risk - .data$alive,
    percent_mortality = round((.data$dead / .data$at_risk) * 100, 1),
    Treatment = factor(.data$Treatment, levels = treatment_order)
  ) %>%
  dplyr::select(-"n") %>%
  dplyr::relocate("at_risk", .before = "alive") %>%
  dplyr::relocate("dead", .after = "alive") %>%
  dplyr::ungroup() %>%
  dplyr::arrange(.data$Treatment)

tbl_deaths <- mortality_treatment_tbl %>%
  dplyr::transmute(
    Treatment = as.character(.data$Treatment),
    deaths = .data$dead,
    at_risk = .data$at_risk,
    survived = .data$alive
  )

tbl_mort_pct <- mortality_treatment_tbl %>%
  dplyr::transmute(
    Treatment = as.character(.data$Treatment),
    mortality_percent = .data$percent_mortality
  )

readr::write_csv(tbl_deaths, file.path(path_tbl, "deaths_by_exposure_regime.csv"))
readr::write_csv(tbl_mort_pct, file.path(path_tbl, "mortality_percent_by_exposure_regime.csv"))
readr::write_csv(mortality_treatment_tbl, file.path(path_tbl, "mortality_summary_by_exposure_regime.csv"))

# Also export GT HTML versions for manuscript display (avoid re-rendering in Results Rmd).
gt_style_treatment_cells <- function(gt_tbl, df, col = "Treatment") {
  if (!exists("treatment_color_scale", inherits = TRUE)) {
    return(gt_tbl)
  }
  if (!(col %in% names(df))) {
    return(gt_tbl)
  }
  trts <- unique(as.character(df[[col]]))
  for (trt in trts) {
    if (trt %in% names(treatment_color_scale)) {
      trt_color <- treatment_color_scale[[trt]]
      gt_tbl <- gt_tbl %>%
        gt::tab_style(
          style = gt::cell_fill(color = trt_color),
          locations = gt::cells_body(columns = dplyr::all_of(col), rows = !!rlang::sym(col) == trt)
        ) %>%
        gt::tab_style(
          style = gt::cell_text(color = "white", weight = "bold"),
          locations = gt::cells_body(columns = dplyr::all_of(col), rows = !!rlang::sym(col) == trt)
        )
    }
  }
  gt_tbl
}

gt_deaths <- tbl_deaths %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Deaths by exposure regime (Day 60)",
    subtitle = "Tank-level summary; at_risk = 45 fish per regime"
  )
gt_deaths <- gt_style_treatment_cells(gt_deaths, tbl_deaths, col = "Treatment")

gt_mort_pct <- tbl_mort_pct %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Percent mortality by exposure regime (Day 60)",
    subtitle = "Tank-level summary; percent mortality = dead / at_risk × 100"
  )
gt_mort_pct <- gt_style_treatment_cells(gt_mort_pct, tbl_mort_pct, col = "Treatment")

gt_mort_summary <- mortality_treatment_tbl %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Mortality summary by exposure regime (Day 60)",
    subtitle = "Counts at risk, alive, dead, and percent mortality per exposure regime"
  )
gt_mort_summary <- gt_style_treatment_cells(gt_mort_summary, mortality_treatment_tbl, col = "Treatment")

gt::gtsave(gt_deaths, file.path(path_tbl, "deaths_by_exposure_regime.html"))
gt::gtsave(gt_mort_pct, file.path(path_tbl, "mortality_percent_by_exposure_regime.html"))
gt::gtsave(gt_mort_summary, file.path(path_tbl, "mortality_summary_by_exposure_regime.html"))

# --- Bar charts: mortality & infection by exposure regime, faceted by prior stressor history
#     (matches manuscript ggplot style: theme_sieler2026_publication, facet_grid + space = "free_x")
prior_stressor_history_from_treatment <- function(trt) {
  factor(
    dplyr::case_when(
      as.character(trt) %in% c("A- T- P-", "A- T- P+") ~ "Zero",
      as.character(trt) %in% c("A+ T- P-", "A+ T- P+", "A- T+ P-", "A- T+ P+") ~ "One",
      as.character(trt) %in% c("A+ T+ P-", "A+ T+ P+") ~ "Two",
      TRUE ~ NA_character_
    ),
    levels = c("Zero", "One", "Two")
  )
}

mort_bar_df <- mortality_treatment_tbl %>%
  dplyr::mutate(prior_stressor_history = prior_stressor_history_from_treatment(.data$Treatment)) %>%
  dplyr::group_by(.data$prior_stressor_history) %>%
  dplyr::mutate(
    Treatment = forcats::fct_reorder(.data$Treatment, .data$percent_mortality, .desc = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(`Prior Stressor History` = prior_stressor_history)

# --- NEW: Bar charts by A×T regime (pooling P) and by A×T regime × P -------------
# A×T regimes collapse parasite status (P-/P+) to show “regardless of parasite exposure”.
at_from_treatment <- function(trt) {
  stringr::str_replace(as.character(trt), "\\s+P[+-]$", "")
}

mort_at_hist_pooled <- tmp_mort_tank %>%
  dplyr::mutate(
    AT = at_from_treatment(.data$Treatment),
    AT = factor(.data$AT, levels = c("A- T-", "A+ T-", "A- T+", "A+ T+"))
  ) %>%
  dplyr::group_by(.data$HistoryLevel, .data$AT) %>%
  dplyr::summarise(
    dead = sum(.data$dead, na.rm = TRUE),
    at_risk = sum(.data$at_risk, na.rm = TRUE),
    percent_mortality = round((.data$dead / .data$at_risk) * 100, 1),
    .groups = "drop"
  )

p_mortality_at_by_history_pooled_parasite <- ggplot2::ggplot(
  mort_at_hist_pooled,
  ggplot2::aes(x = .data$AT, y = .data$percent_mortality, fill = .data$AT)
) +
  ggplot2::geom_col(color = "black", linewidth = SIELER2026_MIN_LINEWIDTH_MM) +
  ggplot2::geom_text(
    ggplot2::aes(
      y = .data$percent_mortality / 2,
      label = paste0(.data$percent_mortality, "%\n(", .data$dead, "/", .data$at_risk, ")")
    ),
    color = "white",
    vjust = 0.5,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(.data$HistoryLevel),
    scales = "free",
    space = "free_x"
  ) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::labs(
    title = "Final mortality by A×T regime (parasite pooled)",
    subtitle = "Day 60; parasite exposure pooled within each stress-history stratum",
    x = "Exposure regime (A × T; P pooled)",
    y = "Mortality (%)",
    caption = "Facets: prior stressor history. Bars pool P- and P+ tanks within each A×T regime."
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(hjust = 0, lineheight = 1.25)
  )

mort_at_hist_by_parasite <- tmp_mort_tank %>%
  dplyr::mutate(
    AT = at_from_treatment(.data$Treatment),
    AT = factor(.data$AT, levels = c("A- T-", "A+ T-", "A- T+", "A+ T+")),
    parasite_label = dplyr::if_else(.data$Parasite == 1L, "Exposed", "Unexposed")
  ) %>%
  dplyr::group_by(.data$HistoryLevel, .data$AT, .data$parasite_label) %>%
  dplyr::summarise(
    dead = sum(.data$dead, na.rm = TRUE),
    at_risk = sum(.data$at_risk, na.rm = TRUE),
    percent_mortality = round((.data$dead / .data$at_risk) * 100, 1),
    .groups = "drop"
  )

p_mortality_at_by_history_and_parasite <- ggplot2::ggplot(
  mort_at_hist_by_parasite,
  ggplot2::aes(
    x = .data$AT,
    y = .data$percent_mortality,
    fill = .data$parasite_label,
    color = .data$HistoryLevel
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.8),
    linewidth = SIELER2026_MIN_LINEWIDTH_MM
  ) +
  ggplot2::scale_fill_manual(
    values = c("Unexposed" = "grey60", "Exposed" = "firebrick"),
    name = "Parasite exposure"
  ) +
  ggplot2::scale_color_manual(
    values = history_color_scale,
    name = "Prior Stressor History"
  ) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(.data$HistoryLevel),
    scales = "free",
    space = "free_x"
  ) +
  theme_sieler2026_publication(base_size = 14, legend_position = "bottom") +
  ggplot2::labs(
    title = "Final mortality by A×T regime, stratified by parasite exposure",
    subtitle = "Day 60; within each prior stressor history level, compare P− vs P+ side-by-side",
    x = "Exposure regime (A × T)",
    y = "Mortality (%)"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold")
  )

p_mortality_exposure_by_prior_history <- ggplot2::ggplot(
  mort_bar_df,
  ggplot2::aes(x = .data$Treatment, y = .data$percent_mortality, fill = .data$Treatment)
) +
  ggplot2::geom_col(color = "black", linewidth = SIELER2026_MIN_LINEWIDTH_MM) +
  ggplot2::geom_text(
    ggplot2::aes(
      y = .data$percent_mortality / 2,
      label = paste0(
        .data$percent_mortality,
        "%\n(",
        .data$dead,
        "/",
        .data$at_risk,
        ")"
      )
    ),
    color = "white",
    vjust = 0.5,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = treatment_color_scale, breaks = treatment_order) +
  ggplot2::scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(`Prior Stressor History`),
    scales = "free",
    space = "free_x"
  ) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::labs(
    title = "Final Mortality by Prior Stressor History",
    x = "Exposure regime",
    y = "Mortality (%)"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold")
  )

fig_mort_exposure_hist_pdf <- file.path(path_fig, "mortality_percent_by_exposure_regime_prior_history.pdf")
fig_mort_exposure_hist_png <- file.path(path_fig, "mortality_percent_by_exposure_regime_prior_history.png")
ggplot2::ggsave(fig_mort_exposure_hist_pdf, p_mortality_exposure_by_prior_history, width = 14, height = 5.75, device = "pdf")
ggplot2::ggsave(fig_mort_exposure_hist_png, p_mortality_exposure_by_prior_history, width = 14, height = 5.75, dpi = 300)

fig_mort_exposure_hist_square_pdf <- file.path(path_fig, "mortality_percent_by_exposure_regime_prior_history_square.pdf")
fig_mort_exposure_hist_square_png <- file.path(path_fig, "mortality_percent_by_exposure_regime_prior_history_square.png")
ggplot2::ggsave(
  fig_mort_exposure_hist_square_pdf,
  p_mortality_exposure_by_prior_history,
  width = fig_bar_sq,
  height = fig_bar_sq,
  device = "pdf"
)
ggplot2::ggsave(
  fig_mort_exposure_hist_square_png,
  p_mortality_exposure_by_prior_history,
  width = fig_bar_sq,
  height = fig_bar_sq,
  dpi = 300
)

fig_mort_at_hist_pooled_pdf <- file.path(path_fig, "mortality_percent_by_AT_regime_prior_history_pooled_parasite.pdf")
fig_mort_at_hist_pooled_png <- file.path(path_fig, "mortality_percent_by_AT_regime_prior_history_pooled_parasite.png")
ggplot2::ggsave(
  fig_mort_at_hist_pooled_pdf,
  p_mortality_at_by_history_pooled_parasite,
  width = 14,
  height = 5.75,
  device = "pdf"
)
ggplot2::ggsave(
  fig_mort_at_hist_pooled_png,
  p_mortality_at_by_history_pooled_parasite,
  width = 14,
  height = 5.75,
  dpi = 300
)

fig_mort_at_hist_para_pdf <- file.path(path_fig, "mortality_percent_by_AT_regime_prior_history_by_parasite.pdf")
fig_mort_at_hist_para_png <- file.path(path_fig, "mortality_percent_by_AT_regime_prior_history_by_parasite.png")
ggplot2::ggsave(
  fig_mort_at_hist_para_pdf,
  p_mortality_at_by_history_and_parasite,
  width = 14,
  height = 5.75,
  device = "pdf"
)
ggplot2::ggsave(
  fig_mort_at_hist_para_png,
  p_mortality_at_by_history_and_parasite,
  width = 14,
  height = 5.75,
  dpi = 300
)

inf_bar_df <- data.list[["Infection"]] %>%
  dplyr::mutate(prior_stressor_history = prior_stressor_history_from_treatment(.data$Treatment)) %>%
  dplyr::group_by(.data$prior_stressor_history) %>%
  dplyr::mutate(
    Treatment = forcats::fct_reorder(.data$Treatment, .data$percent_infected, .desc = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(`Prior Stressor History` = prior_stressor_history)

p_infection_prevalence_exposure_by_prior_history <- ggplot2::ggplot(
  inf_bar_df,
  ggplot2::aes(x = .data$Treatment, y = .data$percent_infected, fill = .data$Treatment)
) +
  ggplot2::geom_col(color = "black", linewidth = SIELER2026_MIN_LINEWIDTH_MM) +
  ggplot2::geom_text(
    ggplot2::aes(
      y = .data$percent_infected / 2,
      label = paste0(
        .data$percent_infected,
        "%\n(",
        .data$n_infected,
        "/",
        .data$n_survivors_sampled,
        ")"
      )
    ),
    color = "white",
    vjust = 0.5,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = treatment_color_scale, breaks = treatment_order) +
  ggplot2::scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(`Prior Stressor History`),
    scales = "free",
    space = "free_x"
  ) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::labs(
    title = "Infection Prevalence by Prior Stressor History",
    x = "Exposure regime",
    y = "Infection prevalence (%)"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold")
  )

fig_inf_exposure_hist_pdf <- file.path(path_fig, "infection_prevalence_by_exposure_regime_prior_history.pdf")
fig_inf_exposure_hist_png <- file.path(path_fig, "infection_prevalence_by_exposure_regime_prior_history.png")
ggplot2::ggsave(fig_inf_exposure_hist_pdf, p_infection_prevalence_exposure_by_prior_history, width = 14, height = 5.75, device = "pdf")
ggplot2::ggsave(fig_inf_exposure_hist_png, p_infection_prevalence_exposure_by_prior_history, width = 14, height = 5.75, dpi = 300)

fig_inf_exposure_hist_square_pdf <- file.path(path_fig, "infection_prevalence_by_exposure_regime_prior_history_square.pdf")
fig_inf_exposure_hist_square_png <- file.path(path_fig, "infection_prevalence_by_exposure_regime_prior_history_square.png")
ggplot2::ggsave(
  fig_inf_exposure_hist_square_pdf,
  p_infection_prevalence_exposure_by_prior_history,
  width = fig_bar_sq,
  height = fig_bar_sq,
  device = "pdf"
)
ggplot2::ggsave(
  fig_inf_exposure_hist_square_png,
  p_infection_prevalence_exposure_by_prior_history,
  width = fig_bar_sq,
  height = fig_bar_sq,
  dpi = 300
)

readr::write_csv(
  mort_bar_df %>%
    dplyr::rename(Prior_Stressor_History = `Prior Stressor History`),
  file.path(path_tbl, "mortality_percent_by_exposure_regime_prior_history.csv")
)
readr::write_csv(
  inf_bar_df %>%
    dplyr::rename(Prior_Stressor_History = `Prior Stressor History`),
  file.path(path_tbl, "infection_prevalence_by_exposure_regime_prior_history.csv")
)

# Fish-level P+ Day 60 survivors (one row per Sample); used for burden summaries & bar plots
inf_fish_p1 <- data.list[["Mortality"]] %>%
  dplyr::filter(.data$Time == 60L, .data$Parasite == 1L) %>%
  dplyr::distinct(.data$Sample, .keep_all = TRUE)

history_stressor_colors <- c(
  "No prior stressors" = "#1B9E77",
  "One prior stressor" = "#D95F02",
  "Two prior stressors" = "#7570B3"
)

# --- Infection: tank-level prevalence GLMM + trend plot -----------------------------
tmp_infection_tank <- data.list[["Infection_Tank"]] %>%
  dplyr::filter(.data$n_survivors_sampled > 0L) %>%
  dplyr::mutate(Treatment = forcats::fct_drop(.data$Treatment))

set.seed(42)
fit_inf_history_num <- glmmTMB::glmmTMB(
  cbind(n_infected, n_survivors_sampled - n_infected) ~ HistoryLevelNum + (1 | Tank.ID),
  family = stats::binomial(link = "logit"),
  data = tmp_infection_tank
)

emm_inf_history_num <- emmeans::emmeans(
  fit_inf_history_num,
  ~ HistoryLevelNum,
  at = list(HistoryLevelNum = c(0, 1, 2)),
  type = "response"
)

trend_infection_history <- emmeans::contrast(emm_inf_history_num, "poly")

pairs_inf_history <- emmeans::contrast(emm_inf_history_num, "pairwise", adjust = "fdr")

df_inf_hist_poly <- as.data.frame(trend_infection_history)
df_inf_hist_pairs <- as.data.frame(pairs_inf_history)
df_inf_hist_emmeans <- as.data.frame(
  summary(emm_inf_history_num, type = "response", infer = c(TRUE, TRUE))
)

readr::write_csv(df_inf_hist_poly, file.path(path_tbl, "infection_glmm_history_poly_contrast.csv"))
readr::write_csv(df_inf_hist_pairs, file.path(path_tbl, "infection_glmm_history_pairs_fdr.csv"))
readr::write_csv(df_inf_hist_emmeans, file.path(path_tbl, "infection_glmm_history_emmeans_response.csv"))

set.seed(42)
fit_inf_treatment <- glmmTMB::glmmTMB(
  cbind(n_infected, n_survivors_sampled - n_infected) ~ Treatment + (1 | Tank.ID),
  family = stats::binomial(link = "logit"),
  data = tmp_infection_tank
)

emm_inf_treatment <- emmeans::emmeans(
  fit_inf_treatment,
  ~ Treatment,
  type = "response"
)

pairs_inf_treatment <- emmeans::contrast(emm_inf_treatment, "pairwise", adjust = "fdr")
readr::write_csv(
  as.data.frame(pairs_inf_treatment),
  file.path(path_tbl, "infection_glmm_treatment_pairs_fdr.csv")
)

emm_inf_plot <- emm_inf_history_num %>%
  as.data.frame() %>%
  dplyr::mutate(HistoryLevelNum = as.numeric(as.character(.data$HistoryLevelNum)))

infection_tank_obs <- tmp_infection_tank %>%
  dplyr::mutate(prevalence_obs = .data$n_infected / .data$n_survivors_sampled)

p_infection_prevalence_trend <- glmm_binomial_tank_history_numeric_trend_plot(
  obs_df = infection_tank_obs,
  obs_y_col = "prevalence_obs",
  emm_df = emm_inf_plot,
  y_label = "Infection prevalence (%)",
  title = "Infection Prevalence by Prior Stressor History",
  subtitle = NULL,
  y_as_percent = TRUE
)
p_infection_prevalence_trend_qr <- glmm_binomial_tank_history_numeric_trend_plot(
  obs_df = infection_tank_obs,
  obs_y_col = "prevalence_obs",
  emm_df = emm_inf_plot,
  y_label = "Infection prevalence (%)",
  title = "Infection Prevalence by Prior Stressor History",
  subtitle = NULL,
  y_as_percent = TRUE,
  observed_layout = "quasirandom"
)

fig_inf_pdf <- file.path(path_fig, "infection_prevalence_trend_predicted.pdf")
fig_inf_png <- file.path(path_fig, "infection_prevalence_trend_predicted.png")
fig_inf_square_pdf <- file.path(path_fig, "infection_prevalence_trend_predicted_square.pdf")
fig_inf_square_png <- file.path(path_fig, "infection_prevalence_trend_predicted_square.png")
ggplot2::ggsave(fig_inf_pdf, p_infection_prevalence_trend, width = fig_in, height = fig_in, device = "pdf")
ggplot2::ggsave(fig_inf_png, p_infection_prevalence_trend, width = fig_in, height = fig_in, dpi = 300)
ggplot2::ggsave(fig_inf_square_pdf, p_infection_prevalence_trend, width = fig_in, height = fig_in, device = "pdf")
ggplot2::ggsave(fig_inf_square_png, p_infection_prevalence_trend, width = fig_in, height = fig_in, dpi = 300)
fig_inf_qr_pdf <- file.path(path_fig, "infection_prevalence_trend_predicted_quasirandom.pdf")
fig_inf_qr_png <- file.path(path_fig, "infection_prevalence_trend_predicted_quasirandom.png")
fig_inf_square_qr_pdf <- file.path(path_fig, "infection_prevalence_trend_predicted_square_quasirandom.pdf")
fig_inf_square_qr_png <- file.path(path_fig, "infection_prevalence_trend_predicted_square_quasirandom.png")
ggplot2::ggsave(fig_inf_qr_pdf, p_infection_prevalence_trend_qr, width = fig_in, height = fig_in, device = "pdf")
ggplot2::ggsave(fig_inf_qr_png, p_infection_prevalence_trend_qr, width = fig_in, height = fig_in, dpi = 300)
ggplot2::ggsave(fig_inf_square_qr_pdf, p_infection_prevalence_trend_qr, width = fig_in, height = fig_in, device = "pdf")
ggplot2::ggsave(fig_inf_square_qr_png, p_infection_prevalence_trend_qr, width = fig_in, height = fig_in, dpi = 300)

# --- Infection prevalence & burden by prior stressor history (P+ Day 60 survivors) ---
infection_history <- inf_fish_p1 %>%
  dplyr::group_by(.data$History) %>%
  dplyr::summarise(
    n_survivors_sampled = dplyr::n(),
    n_infected = sum(.data$Total.Worm.Count > 0, na.rm = TRUE),
    mean_worm_burden = mean(.data$Total.Worm.Count, na.rm = TRUE),
    total_worm_count = sum(.data$Total.Worm.Count, na.rm = TRUE),
    mean_worm_burden_infected = mean(
      dplyr::if_else(.data$Total.Worm.Count > 0, .data$Total.Worm.Count, NA_real_),
      na.rm = TRUE
    ),
    percent_infected = round((.data$n_infected / .data$n_survivors_sampled) * 100, 1),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    History_Label = dplyr::case_when(
      .data$History == 0L ~ "No prior stressors",
      .data$History == 1L ~ "One prior stressor",
      .data$History == 2L ~ "Two prior stressors",
      TRUE ~ "Unknown"
    ),
    History_Label = factor(
      .data$History_Label,
      levels = c("No prior stressors", "One prior stressor", "Two prior stressors")
    )
  )

readr::write_csv(infection_history, file.path(path_tbl, "infection_prevalence_by_prior_stressor_history.csv"))

gt_infection_history <- infection_history %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Infection prevalence by prior stressor history (Day 60, P+ survivors)",
    subtitle = "One row per history level; prevalence = infected / sampled survivors"
  ) %>%
  gt::cols_label(
    History_Label = "Prior stressor history",
    n_survivors_sampled = "Survivors sampled",
    n_infected = "Infected",
    percent_infected = "Percent infected",
    mean_worm_burden = "Mean worm burden",
    total_worm_count = "Total worm count",
    mean_worm_burden_infected = "Mean burden (infected only)"
  )
gt::gtsave(gt_infection_history, file.path(path_tbl, "infection_prevalence_by_prior_stressor_history.html"))

# Burden by exposure regime (facets = prior stressor history), fish-level sums / conditional means
inf_burden_by_regime_prior_df <- inf_fish_p1 %>%
  dplyr::mutate(prior_stressor_history = prior_stressor_history_from_treatment(.data$Treatment)) %>%
  dplyr::group_by(.data$prior_stressor_history, .data$Treatment) %>%
  dplyr::summarise(
    n_survivors_sampled = dplyr::n(),
    n_infected = sum(.data$Total.Worm.Count > 0, na.rm = TRUE),
    percent_infected = round((.data$n_infected / .data$n_survivors_sampled) * 100, 1),
    total_worm_count = sum(.data$Total.Worm.Count, na.rm = TRUE),
    mean_worm_burden_infected = mean(
      dplyr::if_else(.data$Total.Worm.Count > 0, .data$Total.Worm.Count, NA_real_),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::group_by(.data$prior_stressor_history) %>%
  dplyr::mutate(
    Treatment = forcats::fct_reorder(.data$Treatment, .data$total_worm_count, .desc = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(`Prior Stressor History` = prior_stressor_history)

readr::write_csv(
  inf_burden_by_regime_prior_df %>%
    dplyr::rename(Prior_Stressor_History = `Prior Stressor History`),
  file.path(path_tbl, "infection_burden_by_exposure_regime_prior_stressor_history.csv")
)

gt_inf_burden <- inf_burden_by_regime_prior_df %>%
  dplyr::rename(Prior_Stressor_History = `Prior Stressor History`) %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Infection burden by exposure regime (P+; Day 60)",
    subtitle = "Faceted groupings by prior stressor history; totals and infected-only mean burden"
  ) %>%
  gt::cols_label(
    Prior_Stressor_History = "Prior stressor history",
    Treatment = "Exposure regime",
    n_survivors_sampled = "Survivors sampled",
    n_infected = "Infected",
    percent_infected = "Percent infected",
    total_worm_count = "Total worm count",
    mean_worm_burden_infected = "Mean burden (infected only)"
  )
gt_inf_burden <- gt_style_treatment_cells(gt_inf_burden, inf_burden_by_regime_prior_df, col = "Treatment")
gt::gtsave(gt_inf_burden, file.path(path_tbl, "infection_burden_by_exposure_regime_prior_stressor_history.html"))

# --- Bar charts: infection metrics by prior stressor history (single panel, 3 levels) ---
p_infection_prevalence_prior_history_bars <- ggplot2::ggplot(
  infection_history,
  ggplot2::aes(x = .data$History_Label, y = .data$percent_infected, fill = .data$History_Label)
) +
  ggplot2::geom_col(color = "black", linewidth = SIELER2026_MIN_LINEWIDTH_MM) +
  ggplot2::geom_text(
    ggplot2::aes(
      y = .data$percent_infected / 2,
      label = paste0(
        .data$percent_infected,
        "%\n(",
        .data$n_infected,
        "/",
        .data$n_survivors_sampled,
        ")"
      )
    ),
    color = "white",
    vjust = 0.5,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = history_stressor_colors) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::expand_limits(y = 0) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::labs(
    title = "Infection prevalence by prior stressor history",
    subtitle = "Day 60 P+ survivors; percent infected = infected / sampled fish per history level",
    x = "Prior stressor history",
    y = "Percent infected",
    caption = paste(
      "Percent infected = # fish with worms / total fish sampled (# infected / total sampled).",
      "History = count of prior stressors (antibiotic and/or temperature) before parasite phase.",
      sep = "\n"
    )
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.caption = ggplot2::element_text(hjust = 1, lineheight = 1.25)
  )

p_infection_total_worms_prior_history_bars <- ggplot2::ggplot(
  infection_history,
  ggplot2::aes(x = .data$History_Label, y = .data$total_worm_count, fill = .data$History_Label)
) +
  ggplot2::geom_col(color = "black", linewidth = SIELER2026_MIN_LINEWIDTH_MM) +
  ggplot2::geom_text(
    ggplot2::aes(
      y = .data$total_worm_count / 2,
      label = as.integer(.data$total_worm_count)
    ),
    color = "white",
    vjust = 0.5,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = history_stressor_colors) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::expand_limits(y = 0) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::labs(
    title = "Total worm counts by prior stressor history",
    subtitle = "Sum of all worms across all sampled P+ survivors per history level",
    x = "Prior stressor history",
    y = "Total worm count",
    caption = "Total worms = sum of worm counts across all fish sampled within each history level."
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.caption = ggplot2::element_text(hjust = 0, lineheight = 1.25)
  )

infection_history_mean_lab <- infection_history %>%
  dplyr::mutate(
    mean_worm_burden_infected_plot = dplyr::coalesce(.data$mean_worm_burden_infected, 0),
    label_mean_infected = dplyr::if_else(
      .data$n_infected > 0L,
      paste0(
        sprintf("%.1f", round(.data$mean_worm_burden_infected, 1)),
        " (n=",
        .data$n_infected,
        ")"
      ),
      "—"
    ),
    y_text = dplyr::if_else(
      .data$n_infected > 0L,
      .data$mean_worm_burden_infected_plot / 2,
      0
    )
  )

p_infection_mean_worms_infected_prior_history_bars <- ggplot2::ggplot(
  infection_history_mean_lab,
  ggplot2::aes(x = .data$History_Label, y = .data$mean_worm_burden_infected_plot, fill = .data$History_Label)
) +
  ggplot2::geom_col(color = "black", linewidth = SIELER2026_MIN_LINEWIDTH_MM) +
  ggplot2::geom_text(
    ggplot2::aes(y = .data$y_text, label = .data$label_mean_infected),
    color = "white",
    vjust = 0.5,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = history_stressor_colors) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::expand_limits(y = 0) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::labs(
    title = "Mean worm burden among infected fish",
    subtitle = "Average worm count per infected fish only (P+ survivors, Day 60)",
    x = "Prior stressor history",
    y = "Mean worm count (infected fish)",
    caption = paste(
      "Worm count = average worms per infected fish (excluding uninfected fish);",
      "n = number of infected fish in that history level.",
      sep = "\n"
    )
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.caption = ggplot2::element_text(hjust = 1, lineheight = 1.25)
  )

fig_inf_prev_hist_bars_pdf <- file.path(path_fig, "infection_prevalence_by_prior_stressor_history_bars.pdf")
fig_inf_prev_hist_bars_png <- file.path(path_fig, "infection_prevalence_by_prior_stressor_history_bars.png")
fig_inf_prev_hist_bars_square_pdf <- file.path(path_fig, "infection_prevalence_by_prior_stressor_history_bars_square.pdf")
fig_inf_prev_hist_bars_square_png <- file.path(path_fig, "infection_prevalence_by_prior_stressor_history_bars_square.png")
ggplot2::ggsave(fig_inf_prev_hist_bars_pdf, p_infection_prevalence_prior_history_bars, width = 7, height = 5.75, device = "pdf")
ggplot2::ggsave(fig_inf_prev_hist_bars_png, p_infection_prevalence_prior_history_bars, width = 7, height = 5.75, dpi = 300)
ggplot2::ggsave(fig_inf_prev_hist_bars_square_pdf, p_infection_prevalence_prior_history_bars, width = fig_bar_sq, height = fig_bar_sq, device = "pdf")
ggplot2::ggsave(fig_inf_prev_hist_bars_square_png, p_infection_prevalence_prior_history_bars, width = fig_bar_sq, height = fig_bar_sq, dpi = 300)

fig_inf_total_hist_bars_pdf <- file.path(path_fig, "infection_total_worm_count_by_prior_stressor_history_bars.pdf")
fig_inf_total_hist_bars_png <- file.path(path_fig, "infection_total_worm_count_by_prior_stressor_history_bars.png")
fig_inf_total_hist_bars_square_pdf <- file.path(path_fig, "infection_total_worm_count_by_prior_stressor_history_bars_square.pdf")
fig_inf_total_hist_bars_square_png <- file.path(path_fig, "infection_total_worm_count_by_prior_stressor_history_bars_square.png")
ggplot2::ggsave(fig_inf_total_hist_bars_pdf, p_infection_total_worms_prior_history_bars, width = 7, height = 5.75, device = "pdf")
ggplot2::ggsave(fig_inf_total_hist_bars_png, p_infection_total_worms_prior_history_bars, width = 7, height = 5.75, dpi = 300)
ggplot2::ggsave(fig_inf_total_hist_bars_square_pdf, p_infection_total_worms_prior_history_bars, width = fig_bar_sq, height = fig_bar_sq, device = "pdf")
ggplot2::ggsave(fig_inf_total_hist_bars_square_png, p_infection_total_worms_prior_history_bars, width = fig_bar_sq, height = fig_bar_sq, dpi = 300)

fig_inf_mean_hist_bars_pdf <- file.path(path_fig, "infection_mean_worm_burden_infected_by_prior_stressor_history_bars.pdf")
fig_inf_mean_hist_bars_png <- file.path(path_fig, "infection_mean_worm_burden_infected_by_prior_stressor_history_bars.png")
fig_inf_mean_hist_bars_square_pdf <- file.path(path_fig, "infection_mean_worm_burden_infected_by_prior_stressor_history_bars_square.pdf")
fig_inf_mean_hist_bars_square_png <- file.path(path_fig, "infection_mean_worm_burden_infected_by_prior_stressor_history_bars_square.png")
ggplot2::ggsave(fig_inf_mean_hist_bars_pdf, p_infection_mean_worms_infected_prior_history_bars, width = 7, height = 5.75, device = "pdf")
ggplot2::ggsave(fig_inf_mean_hist_bars_png, p_infection_mean_worms_infected_prior_history_bars, width = 7, height = 5.75, dpi = 300)
ggplot2::ggsave(fig_inf_mean_hist_bars_square_pdf, p_infection_mean_worms_infected_prior_history_bars, width = fig_bar_sq, height = fig_bar_sq, device = "pdf")
ggplot2::ggsave(fig_inf_mean_hist_bars_square_png, p_infection_mean_worms_infected_prior_history_bars, width = fig_bar_sq, height = fig_bar_sq, dpi = 300)

# --- Bar charts: total & mean burden by exposure regime (faceted by prior stressor history) ---
p_infection_total_worms_exposure_by_prior_history <- ggplot2::ggplot(
  inf_burden_by_regime_prior_df,
  ggplot2::aes(x = .data$Treatment, y = .data$total_worm_count, fill = .data$Treatment)
) +
  ggplot2::geom_col(color = "black", linewidth = SIELER2026_MIN_LINEWIDTH_MM) +
  ggplot2::geom_text(
    ggplot2::aes(
      y = .data$total_worm_count / 2,
      label = as.integer(.data$total_worm_count)
    ),
    color = "white",
    vjust = 0.5,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = treatment_color_scale, breaks = treatment_order) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(`Prior Stressor History`),
    scales = "free",
    space = "free_x"
  ) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::labs(
    title = "Total worm counts by exposure regime",
    subtitle = "Sum of worm counts across sampled P+ survivors; facets = prior stressor history",
    x = "Exposure regime",
    y = "Total worm count",
    caption = paste0(
      "Facets: prior stressor history. Colors match exposure regime.\n",
      "Total worms = sum of worm counts for all sampled fish in each regime × history cell."
    )
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(hjust = 0, lineheight = 1.25)
  )

inf_burden_mean_lab <- inf_burden_by_regime_prior_df %>%
  dplyr::mutate(
    mean_worm_burden_infected_plot = dplyr::coalesce(.data$mean_worm_burden_infected, 0),
    label_mean_infected = dplyr::if_else(
      .data$n_infected > 0L,
      paste0(
        sprintf("%.1f", round(.data$mean_worm_burden_infected, 1)),
        " (n=",
        .data$n_infected,
        ")"
      ),
      "—"
    ),
    y_text = dplyr::if_else(
      .data$n_infected > 0L,
      .data$mean_worm_burden_infected_plot / 2,
      0
    )
  ) %>%
  dplyr::group_by(.data$`Prior Stressor History`) %>%
  dplyr::mutate(
    Treatment = forcats::fct_reorder(.data$Treatment, .data$mean_worm_burden_infected, .desc = TRUE)
  ) %>%
  dplyr::ungroup()

p_infection_mean_worms_infected_exposure_by_prior_history <- ggplot2::ggplot(
  inf_burden_mean_lab,
  ggplot2::aes(x = .data$Treatment, y = .data$mean_worm_burden_infected_plot, fill = .data$Treatment)
) +
  ggplot2::geom_col(color = "black", linewidth = SIELER2026_MIN_LINEWIDTH_MM) +
  ggplot2::geom_text(
    ggplot2::aes(y = .data$y_text, label = .data$label_mean_infected),
    color = "white",
    vjust = 0.5,
    size = 4,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = treatment_color_scale, breaks = treatment_order) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(`Prior Stressor History`),
    scales = "free",
    space = "free_x"
  ) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::labs(
    title = "Mean worm burden among infected fish by exposure regime",
    subtitle = "Average worms per infected fish; facets = prior stressor history",
    x = "Exposure regime",
    y = "Mean worm count (infected fish)",
    caption = paste0(
      "Mean among fish with Total.Worm.Count > 0 only; n = infected fish per cell.\n",
      "Facets: prior stressor history; colors match exposure regime."
    )
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(hjust = 0, lineheight = 1.25)
  )

fig_inf_total_regime_hist_pdf <- file.path(path_fig, "infection_total_worm_count_by_exposure_regime_prior_history.pdf")
fig_inf_total_regime_hist_png <- file.path(path_fig, "infection_total_worm_count_by_exposure_regime_prior_history.png")
fig_inf_total_regime_hist_square_pdf <- file.path(path_fig, "infection_total_worm_count_by_exposure_regime_prior_history_square.pdf")
fig_inf_total_regime_hist_square_png <- file.path(path_fig, "infection_total_worm_count_by_exposure_regime_prior_history_square.png")
ggplot2::ggsave(fig_inf_total_regime_hist_pdf, p_infection_total_worms_exposure_by_prior_history, width = 14, height = 5.75, device = "pdf")
ggplot2::ggsave(fig_inf_total_regime_hist_png, p_infection_total_worms_exposure_by_prior_history, width = 14, height = 5.75, dpi = 300)
ggplot2::ggsave(fig_inf_total_regime_hist_square_pdf, p_infection_total_worms_exposure_by_prior_history, width = fig_bar_sq, height = fig_bar_sq, device = "pdf")
ggplot2::ggsave(fig_inf_total_regime_hist_square_png, p_infection_total_worms_exposure_by_prior_history, width = fig_bar_sq, height = fig_bar_sq, dpi = 300)

fig_inf_mean_regime_hist_pdf <- file.path(path_fig, "infection_mean_worm_burden_infected_by_exposure_regime_prior_history.pdf")
fig_inf_mean_regime_hist_png <- file.path(path_fig, "infection_mean_worm_burden_infected_by_exposure_regime_prior_history.png")
fig_inf_mean_regime_hist_square_pdf <- file.path(path_fig, "infection_mean_worm_burden_infected_by_exposure_regime_prior_history_square.pdf")
fig_inf_mean_regime_hist_square_png <- file.path(path_fig, "infection_mean_worm_burden_infected_by_exposure_regime_prior_history_square.png")
ggplot2::ggsave(fig_inf_mean_regime_hist_pdf, p_infection_mean_worms_infected_exposure_by_prior_history, width = 14, height = 5.75, device = "pdf")
ggplot2::ggsave(fig_inf_mean_regime_hist_png, p_infection_mean_worms_infected_exposure_by_prior_history, width = 14, height = 5.75, dpi = 300)
ggplot2::ggsave(fig_inf_mean_regime_hist_square_pdf, p_infection_mean_worms_infected_exposure_by_prior_history, width = fig_bar_sq, height = fig_bar_sq, device = "pdf")
ggplot2::ggsave(fig_inf_mean_regime_hist_square_png, p_infection_mean_worms_infected_exposure_by_prior_history, width = fig_bar_sq, height = fig_bar_sq, dpi = 300)

saveRDS(
  list(
    p_mortality_exposure_by_prior_history = p_mortality_exposure_by_prior_history,
    p_infection_prevalence_exposure_by_prior_history = p_infection_prevalence_exposure_by_prior_history,
    p_mortality_prior_stressor_trend = p_mortality_prior_stressor_trend,
    p_mortality_prior_stressor_trend_qr = p_mortality_prior_stressor_trend_qr,
    p_infection_prevalence_trend_qr = p_infection_prevalence_trend_qr
  ),
  file.path(path_stats, "mortinf_main_text_figure_ggplots.rds")
)

# --- Save RDS ----------------------------------------------------------------------
saveRDS(tmp_mort_tank, file.path(path_stats, "tmp.Mort_Tank.rds"))
saveRDS(tmp_infection_tank, file.path(path_stats, "tmp.Infection_Tank.rds"))
saveRDS(fit_mort_history_num, file.path(path_stats, "fit_mort_history_num.rds"))
saveRDS(fit_mort_history_parasite, file.path(path_stats, "fit_mort_history_parasite.rds"))
saveRDS(fit_mort_history_factor, file.path(path_stats, "fit_mort_history_factor.rds"))
saveRDS(fit_mort_history_factor_parasite, file.path(path_stats, "fit_mort_history_factor_parasite.rds"))
saveRDS(fit_inf_history_num, file.path(path_stats, "fit_inf_history_num.rds"))
saveRDS(fit_inf_treatment, file.path(path_stats, "fit_inf_treatment.rds"))
saveRDS(emm_mort_history_num, file.path(path_stats, "emm_mort_history_num.rds"))
saveRDS(emm_mort_history_factor, file.path(path_stats, "emm_mort_history_factor.rds"))
saveRDS(emm_inf_history_num, file.path(path_stats, "emm_inf_history_num.rds"))
saveRDS(emm_inf_treatment, file.path(path_stats, "emm_inf_treatment.rds"))
saveRDS(trend_infection_history, file.path(path_stats, "trend_infection_history.rds"))

# Project-relative paths (portable; work with knitr root.dir = project root)
path_res_rel <- file.path("Results", "05__Mort-Inf")
path_fig_rel <- file.path(path_res_rel, "Figures")
path_tbl_rel <- file.path(path_res_rel, "Tables")
path_stats_rel <- file.path(path_res_rel, "Stats")

bundle <- list(
  meta = list(
    run_date = as.character(Sys.Date()),
    script = "Code/01__Analysis/05__Mort-Inf.R",
    models = c(
      mortality_main = "cbind(dead, at_risk - dead) ~ HistoryLevelNum * Parasite + (1 | Tank.ID), binomial glmmTMB",
      mortality_null_history = "cbind(dead, at_risk - dead) ~ HistoryLevelNum + (1 | Tank.ID), binomial glmmTMB",
      infection_history_linear = "cbind(n_infected, n_survivors_sampled - n_infected) ~ HistoryLevelNum + (1 | Tank.ID), P+ Day 60 survivors, binomial glmmTMB",
      infection_treatment = "cbind(n_infected, n_survivors_sampled - n_infected) ~ Treatment + (1 | Tank.ID), P+ Day 60 survivors, binomial glmmTMB; pairs(..., adjust = FDR)"
    )
  ),
  paths = list(
    figures = list(
      mortality_prior_stressor_trend_pdf = file.path(path_fig_rel, "mortality_prior_stressor_trend.pdf"),
      mortality_prior_stressor_trend_png = file.path(path_fig_rel, "mortality_prior_stressor_trend.png"),
      mortality_prior_stressor_trend_square_pdf = file.path(path_fig_rel, "mortality_prior_stressor_trend_square.pdf"),
      mortality_prior_stressor_trend_square_png = file.path(path_fig_rel, "mortality_prior_stressor_trend_square.png"),
      mortality_prior_stressor_trend_quasirandom_pdf = file.path(
        path_fig_rel,
        "mortality_prior_stressor_trend_quasirandom.pdf"
      ),
      mortality_prior_stressor_trend_quasirandom_png = file.path(
        path_fig_rel,
        "mortality_prior_stressor_trend_quasirandom.png"
      ),
      mortality_prior_stressor_trend_square_quasirandom_pdf = file.path(
        path_fig_rel,
        "mortality_prior_stressor_trend_square_quasirandom.pdf"
      ),
      mortality_prior_stressor_trend_square_quasirandom_png = file.path(
        path_fig_rel,
        "mortality_prior_stressor_trend_square_quasirandom.png"
      ),
      mortality_by_prior_stressor_history_levels_pdf = file.path(
        path_fig_rel,
        "mortality_by_prior_stressor_history_levels.pdf"
      ),
      mortality_by_prior_stressor_history_levels_png = file.path(
        path_fig_rel,
        "mortality_by_prior_stressor_history_levels.png"
      ),
      mortality_by_prior_stressor_history_levels_x_parasite_pdf = file.path(
        path_fig_rel,
        "mortality_by_prior_stressor_history_levels_x_parasite.pdf"
      ),
      mortality_by_prior_stressor_history_levels_x_parasite_png = file.path(
        path_fig_rel,
        "mortality_by_prior_stressor_history_levels_x_parasite.png"
      ),
      mortality_parasite_trend_comparison_pdf = file.path(path_fig_rel, "mortality_parasite_trend_comparison.pdf"),
      mortality_parasite_trend_comparison_png = file.path(path_fig_rel, "mortality_parasite_trend_comparison.png"),
      mortality_parasite_trend_comparison_square_pdf = file.path(
        path_fig_rel,
        "mortality_parasite_trend_comparison_square.pdf"
      ),
      mortality_parasite_trend_comparison_square_png = file.path(
        path_fig_rel,
        "mortality_parasite_trend_comparison_square.png"
      ),
      infection_prevalence_trend_predicted_pdf = file.path(
        path_fig_rel,
        "infection_prevalence_trend_predicted.pdf"
      ),
      infection_prevalence_trend_predicted_png = file.path(
        path_fig_rel,
        "infection_prevalence_trend_predicted.png"
      ),
      infection_prevalence_trend_predicted_square_pdf = file.path(
        path_fig_rel,
        "infection_prevalence_trend_predicted_square.pdf"
      ),
      infection_prevalence_trend_predicted_square_png = file.path(
        path_fig_rel,
        "infection_prevalence_trend_predicted_square.png"
      ),
      infection_prevalence_trend_predicted_quasirandom_pdf = file.path(
        path_fig_rel,
        "infection_prevalence_trend_predicted_quasirandom.pdf"
      ),
      infection_prevalence_trend_predicted_quasirandom_png = file.path(
        path_fig_rel,
        "infection_prevalence_trend_predicted_quasirandom.png"
      ),
      infection_prevalence_trend_predicted_square_quasirandom_pdf = file.path(
        path_fig_rel,
        "infection_prevalence_trend_predicted_square_quasirandom.pdf"
      ),
      infection_prevalence_trend_predicted_square_quasirandom_png = file.path(
        path_fig_rel,
        "infection_prevalence_trend_predicted_square_quasirandom.png"
      ),
      mortality_percent_by_exposure_regime_prior_history_pdf = file.path(
        path_fig_rel,
        "mortality_percent_by_exposure_regime_prior_history.pdf"
      ),
      mortality_percent_by_exposure_regime_prior_history_png = file.path(
        path_fig_rel,
        "mortality_percent_by_exposure_regime_prior_history.png"
      ),
      infection_prevalence_by_exposure_regime_prior_history_pdf = file.path(
        path_fig_rel,
        "infection_prevalence_by_exposure_regime_prior_history.pdf"
      ),
      infection_prevalence_by_exposure_regime_prior_history_png = file.path(
        path_fig_rel,
        "infection_prevalence_by_exposure_regime_prior_history.png"
      ),
      mortality_percent_by_exposure_regime_prior_history_square_pdf = file.path(
        path_fig_rel,
        "mortality_percent_by_exposure_regime_prior_history_square.pdf"
      ),
      mortality_percent_by_exposure_regime_prior_history_square_png = file.path(
        path_fig_rel,
        "mortality_percent_by_exposure_regime_prior_history_square.png"
      ),
      mortality_percent_by_AT_regime_prior_history_pooled_parasite_pdf = file.path(
        path_fig_rel,
        "mortality_percent_by_AT_regime_prior_history_pooled_parasite.pdf"
      ),
      mortality_percent_by_AT_regime_prior_history_pooled_parasite_png = file.path(
        path_fig_rel,
        "mortality_percent_by_AT_regime_prior_history_pooled_parasite.png"
      ),
      mortality_percent_by_AT_regime_prior_history_by_parasite_pdf = file.path(
        path_fig_rel,
        "mortality_percent_by_AT_regime_prior_history_by_parasite.pdf"
      ),
      mortality_percent_by_AT_regime_prior_history_by_parasite_png = file.path(
        path_fig_rel,
        "mortality_percent_by_AT_regime_prior_history_by_parasite.png"
      ),
      infection_prevalence_by_exposure_regime_prior_history_square_pdf = file.path(
        path_fig_rel,
        "infection_prevalence_by_exposure_regime_prior_history_square.pdf"
      ),
      infection_prevalence_by_exposure_regime_prior_history_square_png = file.path(
        path_fig_rel,
        "infection_prevalence_by_exposure_regime_prior_history_square.png"
      ),
      infection_prevalence_by_prior_stressor_history_bars_pdf = file.path(
        path_fig_rel,
        "infection_prevalence_by_prior_stressor_history_bars.pdf"
      ),
      infection_prevalence_by_prior_stressor_history_bars_png = file.path(
        path_fig_rel,
        "infection_prevalence_by_prior_stressor_history_bars.png"
      ),
      infection_prevalence_by_prior_stressor_history_bars_square_pdf = file.path(
        path_fig_rel,
        "infection_prevalence_by_prior_stressor_history_bars_square.pdf"
      ),
      infection_prevalence_by_prior_stressor_history_bars_square_png = file.path(
        path_fig_rel,
        "infection_prevalence_by_prior_stressor_history_bars_square.png"
      ),
      infection_total_worm_count_by_prior_stressor_history_bars_pdf = file.path(
        path_fig_rel,
        "infection_total_worm_count_by_prior_stressor_history_bars.pdf"
      ),
      infection_total_worm_count_by_prior_stressor_history_bars_png = file.path(
        path_fig_rel,
        "infection_total_worm_count_by_prior_stressor_history_bars.png"
      ),
      infection_total_worm_count_by_prior_stressor_history_bars_square_pdf = file.path(
        path_fig_rel,
        "infection_total_worm_count_by_prior_stressor_history_bars_square.pdf"
      ),
      infection_total_worm_count_by_prior_stressor_history_bars_square_png = file.path(
        path_fig_rel,
        "infection_total_worm_count_by_prior_stressor_history_bars_square.png"
      ),
      infection_mean_worm_burden_infected_by_prior_stressor_history_bars_pdf = file.path(
        path_fig_rel,
        "infection_mean_worm_burden_infected_by_prior_stressor_history_bars.pdf"
      ),
      infection_mean_worm_burden_infected_by_prior_stressor_history_bars_png = file.path(
        path_fig_rel,
        "infection_mean_worm_burden_infected_by_prior_stressor_history_bars.png"
      ),
      infection_mean_worm_burden_infected_by_prior_stressor_history_bars_square_pdf = file.path(
        path_fig_rel,
        "infection_mean_worm_burden_infected_by_prior_stressor_history_bars_square.pdf"
      ),
      infection_mean_worm_burden_infected_by_prior_stressor_history_bars_square_png = file.path(
        path_fig_rel,
        "infection_mean_worm_burden_infected_by_prior_stressor_history_bars_square.png"
      ),
      infection_total_worm_count_by_exposure_regime_prior_history_pdf = file.path(
        path_fig_rel,
        "infection_total_worm_count_by_exposure_regime_prior_history.pdf"
      ),
      infection_total_worm_count_by_exposure_regime_prior_history_png = file.path(
        path_fig_rel,
        "infection_total_worm_count_by_exposure_regime_prior_history.png"
      ),
      infection_total_worm_count_by_exposure_regime_prior_history_square_pdf = file.path(
        path_fig_rel,
        "infection_total_worm_count_by_exposure_regime_prior_history_square.pdf"
      ),
      infection_total_worm_count_by_exposure_regime_prior_history_square_png = file.path(
        path_fig_rel,
        "infection_total_worm_count_by_exposure_regime_prior_history_square.png"
      ),
      infection_mean_worm_burden_infected_by_exposure_regime_prior_history_pdf = file.path(
        path_fig_rel,
        "infection_mean_worm_burden_infected_by_exposure_regime_prior_history.pdf"
      ),
      infection_mean_worm_burden_infected_by_exposure_regime_prior_history_png = file.path(
        path_fig_rel,
        "infection_mean_worm_burden_infected_by_exposure_regime_prior_history.png"
      ),
      infection_mean_worm_burden_infected_by_exposure_regime_prior_history_square_pdf = file.path(
        path_fig_rel,
        "infection_mean_worm_burden_infected_by_exposure_regime_prior_history_square.pdf"
      ),
      infection_mean_worm_burden_infected_by_exposure_regime_prior_history_square_png = file.path(
        path_fig_rel,
        "infection_mean_worm_burden_infected_by_exposure_regime_prior_history_square.png"
      )
    ),
    tables = list(
      deaths_by_exposure_regime = file.path(path_tbl_rel, "deaths_by_exposure_regime.csv"),
      deaths_by_exposure_regime_html = file.path(path_tbl_rel, "deaths_by_exposure_regime.html"),
      mortality_percent_by_exposure_regime = file.path(
        path_tbl_rel,
        "mortality_percent_by_exposure_regime.csv"
      ),
      mortality_percent_by_exposure_regime_html = file.path(
        path_tbl_rel,
        "mortality_percent_by_exposure_regime.html"
      ),
      mortality_summary_by_exposure_regime = file.path(
        path_tbl_rel,
        "mortality_summary_by_exposure_regime.csv"
      ),
      mortality_summary_by_exposure_regime_html = file.path(
        path_tbl_rel,
        "mortality_summary_by_exposure_regime.html"
      ),
      mortality_percent_by_exposure_regime_prior_history = file.path(
        path_tbl_rel,
        "mortality_percent_by_exposure_regime_prior_history.csv"
      ),
      infection_prevalence_by_exposure_regime_prior_history = file.path(
        path_tbl_rel,
        "infection_prevalence_by_exposure_regime_prior_history.csv"
      ),
      infection_prevalence_by_prior_stressor_history = file.path(
        path_tbl_rel,
        "infection_prevalence_by_prior_stressor_history.csv"
      ),
      infection_prevalence_by_prior_stressor_history_html = file.path(
        path_tbl_rel,
        "infection_prevalence_by_prior_stressor_history.html"
      ),
      infection_burden_by_exposure_regime_prior_stressor_history = file.path(
        path_tbl_rel,
        "infection_burden_by_exposure_regime_prior_stressor_history.csv"
      ),
      infection_burden_by_exposure_regime_prior_stressor_history_html = file.path(
        path_tbl_rel,
        "infection_burden_by_exposure_regime_prior_stressor_history.html"
      ),
      mortality_glmm_joint_tests = file.path(path_tbl_rel, "mortality_glmm_joint_tests.csv"),
      mortality_glmm_joint_tests_factor = file.path(path_tbl_rel, "mortality_glmm_joint_tests_factor.csv"),
      mortality_glmm_history_poly_contrast = file.path(
        path_tbl_rel,
        "mortality_glmm_history_poly_contrast.csv"
      ),
      mortality_glmm_history_pairs_fdr = file.path(path_tbl_rel, "mortality_glmm_history_pairs_fdr.csv"),
      mortality_glmm_history_emmeans_response = file.path(
        path_tbl_rel,
        "mortality_glmm_history_emmeans_response.csv"
      ),
      mortality_glmm_history_factor_emmeans_response = file.path(
        path_tbl_rel,
        "mortality_glmm_history_factor_emmeans_response.csv"
      ),
      mortality_glmm_history_factor_pairs_fdr = file.path(
        path_tbl_rel,
        "mortality_glmm_history_factor_pairs_fdr.csv"
      ),
      infection_glmm_treatment_pairs_fdr = file.path(
        path_tbl_rel,
        "infection_glmm_treatment_pairs_fdr.csv"
      ),
      infection_glmm_history_poly_contrast = file.path(
        path_tbl_rel,
        "infection_glmm_history_poly_contrast.csv"
      ),
      infection_glmm_history_pairs_fdr = file.path(path_tbl_rel, "infection_glmm_history_pairs_fdr.csv"),
      infection_glmm_history_emmeans_response = file.path(
        path_tbl_rel,
        "infection_glmm_history_emmeans_response.csv"
      )
    ),
    stats_dir = path_stats_rel
  ),
  tables = list(
    deaths_by_exposure_regime = tbl_deaths,
    mortality_percent_by_exposure_regime = tbl_mort_pct,
    mortality_summary_by_exposure_regime = mortality_treatment_tbl,
    mortality_percent_by_exposure_regime_prior_history = mort_bar_df %>%
      dplyr::rename(Prior_Stressor_History = `Prior Stressor History`),
    infection_prevalence_by_exposure_regime_prior_history = inf_bar_df %>%
      dplyr::rename(Prior_Stressor_History = `Prior Stressor History`),
    infection_prevalence_by_prior_stressor_history = infection_history,
    infection_burden_by_exposure_regime_prior_stressor_history = inf_burden_by_regime_prior_df %>%
      dplyr::rename(Prior_Stressor_History = `Prior Stressor History`),
    mortality_glmm_joint_tests = as.data.frame(joint_tests_mortality),
    mortality_glmm_joint_tests_factor = as.data.frame(joint_tests_mortality_factor),
    mortality_glmm_history_poly_contrast = df_mort_hist_poly,
    mortality_glmm_history_pairs_fdr = df_mort_hist_pairs,
    mortality_glmm_history_emmeans_response = df_mort_hist_emmeans,
    mortality_glmm_history_factor_emmeans_response = df_mort_hist_factor_emmeans,
    mortality_glmm_history_factor_pairs_fdr = df_mort_hist_factor_pairs,
    infection_glmm_treatment_pairs_fdr = as.data.frame(pairs_inf_treatment),
    infection_glmm_history_poly_contrast = df_inf_hist_poly,
    infection_glmm_history_pairs_fdr = df_inf_hist_pairs,
    infection_glmm_history_emmeans_response = df_inf_hist_emmeans
  ),
  table_deaths_by_exposure_regime = gt_deaths,
  table_mortality_percent_by_exposure_regime = gt_mort_pct,
  table_mortality_summary_by_exposure_regime = gt_mort_summary,
  table_infection_prevalence_by_prior_stressor_history = gt_infection_history,
  table_infection_burden_by_exposure_regime_prior_stressor_history = gt_inf_burden,
  models = list(
    fit_mort_history_num = fit_mort_history_num,
    fit_mort_history_parasite = fit_mort_history_parasite,
    fit_mort_history_factor = fit_mort_history_factor,
    fit_mort_history_factor_parasite = fit_mort_history_factor_parasite,
    fit_inf_history_num = fit_inf_history_num,
    fit_inf_treatment = fit_inf_treatment,
    emm_mort_history_num = emm_mort_history_num,
    emm_mort_history_factor = emm_mort_history_factor,
    trend_mortality_history = trend_mortality_history,
    pairs_mort_history = pairs_mort_history,
    emm_mort_factor_parasite_resp = emm_mort_factor_parasite_resp,
    emm_inf_history_num = emm_inf_history_num,
    emm_inf_treatment = emm_inf_treatment,
    trend_infection_history = trend_infection_history,
    pairs_inf_history = pairs_inf_history,
    pairs_inf_treatment = pairs_inf_treatment
  )
)

saveRDS(bundle, bundle_rds)
message("Saved bundle: ", bundle_rds)

sieler2026_sync_main_figures_from_manifest(driver_script = "05__Mort-Inf.R")
