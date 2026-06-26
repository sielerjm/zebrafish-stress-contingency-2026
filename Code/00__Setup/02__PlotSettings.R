# 02__PlotSettings.R
# Created by: Michael Sieler
# Date last updated: 2026-04-24
#
# Description: Defines ggplot2 default theme, RColorBrewer-based palettes, and named
#   color scales for experiment factors (temperature/DPE, treatments, stress history).
#
# Expected input:  Packages available — ggplot2, RColorBrewer, purrr, grid (grid is
#   typically attached with ggplot2; purrr comes with tidyverse).
# Expected output:  Objects in the global environment: palette vectors, factor orders,
#   named color scales, `col.list`, and side effect `ggplot2::theme_update()` for
#   session-wide plot defaults.

# --- Optional: ensure dependencies when this file is sourced standalone
suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
  library(purrr)
})

# Color palettes ---------------------------------------------------------------
# Reference: RColorBrewer — https://www.datanovia.com/en/blog/the-a-z-of-rcolorbrewer-palette/

# Sequential (ordinal): low → high
pal.BuPu <- RColorBrewer::brewer.pal(9, "BuPu")
pal.Greys <- RColorBrewer::brewer.pal(9, "Greys")
pal.Blues <- RColorBrewer::brewer.pal(9, "Blues")
pal.Greens <- RColorBrewer::brewer.pal(9, "Greens")
pal.Reds <- RColorBrewer::brewer.pal(9, "Reds")
pal.Purples <- RColorBrewer::brewer.pal(9, "Purples")
pal.Oranges <- RColorBrewer::brewer.pal(9, "Oranges")

# Qualitative (categorical)
pal.Set1 <- RColorBrewer::brewer.pal(9, "Set1") # not colorblind-friendly
pal.Set2 <- RColorBrewer::brewer.pal(8, "Set2")
pal.Dark2 <- RColorBrewer::brewer.pal(8, "Dark2")
pal.Paired <- RColorBrewer::brewer.pal(12, "Paired")

# Diverging
pal.RdYlGn <- RColorBrewer::brewer.pal(11, "RdYlGn")
pal.Spectral <- RColorBrewer::brewer.pal(9, "Spectral")
pal.BrBg <- RColorBrewer::brewer.pal(11, "BrBG")

# Sloan neutral-model partitions (above / neutral / below Wilson band) -----------------------
# Used by Code/01__Analysis/08__NeutralModel.R; keep in sync with neutral-model figures.
partition_colors_neutral_model <- c(
  above = "firebrick",
  below = "steelblue",
  neutral = "grey35"
)

# Plot-specific palettes -------------------------------------------------------

col.Temp <- c("#3B65DB", "#7AB84C", "#A03022")
col.DPE <- RColorBrewer::brewer.pal(9, "YlOrRd")
col.Treat <- pal.Dark2[c(3, 5)]
col.Worm <- RColorBrewer::brewer.pal(9, "YlOrBr")

# Five evenly spaced steps from Blues, Greens, Reds brewer scales (for Temp × DPE)
pal.TempDPE <- c(
  RColorBrewer::brewer.pal(9, "Blues")[c(1, 3, 5, 7, 9)],
  RColorBrewer::brewer.pal(9, "Greens")[c(1, 3, 5, 7, 9)],
  RColorBrewer::brewer.pal(9, "Reds")[c(1, 3, 5, 7, 9)]
)

TempDPE.breaks <- c(
  "28°C_0DPE", "32°C_0DPE", "35°C_0DPE",
  "28°C_14DPE", "32°C_14DPE", "35°C_14DPE",
  "28°C_21DPE", "32°C_21DPE", "35°C_21DPE",
  "28°C_28DPE", "32°C_28DPE", "35°C_28DPE",
  "28°C_42DPE", "32°C_42DPE", "35°C_42DPE"
)

col.TempDPE <- c(
  pal.TempDPE[1], pal.TempDPE[6], pal.TempDPE[11],
  pal.TempDPE[2], pal.TempDPE[7], pal.TempDPE[12],
  pal.TempDPE[3], pal.TempDPE[8], pal.TempDPE[13],
  pal.TempDPE[4], pal.TempDPE[9], pal.TempDPE[14],
  pal.TempDPE[5], pal.TempDPE[10], pal.TempDPE[15]
)

col.TempDPE_v2 <- c(
  pal.TempDPE[1:5],
  pal.TempDPE[6:10],
  pal.TempDPE[11:15]
)

TreatTempDPE.breaks <- c(
  "28°C_0DPE_Control", "32°C_0DPE_Control", "35°C_0DPE_Control",
  "28°C_14DPE_Control", "32°C_14DPE_Control", "35°C_14DPE_Control",
  "28°C_21DPE_Control", "32°C_21DPE_Control", "35°C_21DPE_Control",
  "28°C_28DPE_Control", "32°C_28DPE_Control", "35°C_28DPE_Control",
  "28°C_42DPE_Control", "32°C_42DPE_Control", "35°C_42DPE_Control",
  "28°C_0DPE_Exposed", "32°C_0DPE_Exposed", "35°C_0DPE_Exposed",
  "28°C_14DPE_Exposed", "32°C_14DPE_Exposed", "35°C_14DPE_Exposed",
  "28°C_21DPE_Exposed", "32°C_21DPE_Exposed", "35°C_21DPE_Exposed",
  "28°C_28DPE_Exposed", "32°C_28DPE_Exposed", "35°C_28DPE_Exposed",
  "28°C_42DPE_Exposed", "32°C_42DPE_Exposed", "35°C_42DPE_Exposed"
)

TempDPE.breaks_v2 <- c(
  "28°C_0DPE", "28°C_14DPE", "28°C_21DPE", "28°C_28DPE", "28°C_42DPE",
  "32°C_0DPE", "32°C_14DPE", "32°C_21DPE", "32°C_28DPE", "32°C_42DPE",
  "35°C_0DPE", "35°C_14DPE", "35°C_21DPE", "35°C_28DPE", "35°C_42DPE"
)

col.TreatTempDPE <- c(
  pal.TempDPE[1:5], pal.TempDPE[6:10], pal.TempDPE[11:15],
  pal.TempDPE[1:5], pal.TempDPE[6:10], pal.TempDPE[11:15]
)

# Named list of default palettes for looping over treatment dimensions
col.list <- list(
  Temperature = col.Temp,
  DPE = col.DPE,
  Treatment = col.Treat
)

# ggplot2 theme (session default) --------------------------------------------
# Adapted from Keaton Stagaman; updates the active ggplot2 theme for the session.

# Session defaults sized for manuscript figures (~14 pt base) so partial themes are not dominated
# by oversized inherited text. Main-text drivers should still add theme_sieler2026_* explicitly.
my_theme <- ggplot2::theme_update(
  legend.position = "bottom",
  legend.box = "vertical",
  legend.box.just = "center",
  legend.title = ggplot2::element_text(size = 12, face = "bold"),
  legend.text = ggplot2::element_text(size = 11),
  legend.key = ggplot2::element_rect(fill = "white"),
  legend.key.size = grid::unit(1, "line"),
  legend.spacing.y = grid::unit(0, "cm"),
  strip.text = ggplot2::element_text(size = 12, face = "bold"),
  plot.caption = ggplot2::element_text(hjust = 0, size = 10),
  axis.text = ggplot2::element_text(size = 11),
  axis.title = ggplot2::element_text(size = 12, face = "bold"),
  panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
  panel.background = ggplot2::element_rect(fill = "white"),
  panel.grid.major = ggplot2::element_line(colour = pal.Greys[4], linewidth = 0.3),
  panel.grid.minor = ggplot2::element_line(colour = pal.Greys[4], linewidth = 0.2),
  strip.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.4)
)

# Treatment & stress-history scales --------------------------------------------

treatment_order <- c(
  "A- T- P-", # Control
  "A- T- P+", # Parasite
  "A+ T- P-", # Antibiotics
  "A+ T- P+", # Antibiotics_Parasite
  "A- T+ P-", # Temperature
  "A- T+ P+", # Temperature_Parasite
  "A+ T+ P-", # Antibiotics_Temperature
  "A+ T+ P+" # Antibiotics_Temperature_Parasite
)

treatment_colors <- c(
  "#1B9E77", # A- T- P- (Control)
  "#D95F02", # A- T- P+ (Parasite)
  "#7570B3", # A+ T- P- (Antibiotics)
  "#E7298A", # A+ T- P+ (Antibiotics_Parasite)
  "#66A61E", # A- T+ P- (Temperature)
  "#E6AB02", # A- T+ P+ (Temperature_Parasite)
  "#A6761D", # A+ T+ P- (Antibiotics_Temperature)
  "#666666" # A+ T+ P+ (Antibiotics_Temperature_Parasite)
)

treatment_color_scale <- purrr::set_names(treatment_colors, treatment_order)

# Prior stressor count (HistoryLevelNum 0/1/2): single source for trend plots, PCoA, betadisper, etc.
prior_stressor_history_colors_numeric <- c(
  `0` = "#1B9E77",
  `1` = "#D95F02",
  `2` = "#7570B3"
)

history_order <- c(
  "No prior stressors",
  "One prior stressor",
  "Two prior stressors"
)

history_colors <- c(
  "No prior stressors" = prior_stressor_history_colors_numeric[["0"]],
  "One prior stressor" = prior_stressor_history_colors_numeric[["1"]],
  "Two prior stressors" = prior_stressor_history_colors_numeric[["2"]]
)

history_color_scale <- purrr::set_names(history_colors[history_order], history_order)
