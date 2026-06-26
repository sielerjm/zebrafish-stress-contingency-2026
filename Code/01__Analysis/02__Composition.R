# 02__Composition.R
# Created by: Michael Sieler
# Date last updated: 2026-04-25 (composition theme: y-grid + bold legend titles; betadisper layout in helpers;
#   betadisper history×parasite ggsave 9×9 in; PCoA coord_fixed(1, expand=FALSE); history-only legends right;
#   history×parasite PCoA legends bottom (stacked); compact caption)
#
# Description: Gut microbial **composition** - relative abundance by exposure regime (Fig. A-style);
#   PERMANOVA (Tables C-D) and betadisper ANOVA tables on Bray-Curtis and Canberra distances;
#   PCoA and betadisper violin figures for both metrics (stress history; history × parasite).
#   Stratified: PERMANOVA + betadisper for Parasite within each HistoryLevelNum; faceted PCoA (shared axes).
#   Table D PERMANOVA: distance ~ HistoryLevelNum * Parasite (parallels alpha GLMM fixed effects; no tank RE).
#   Writes figures, CSV tables, gt HTML, and `composition__gut__bundle.rds` for `Code/02__Results/02__Composition.Rmd`.
#
# Expected input:  Run from Sieler2026 project root; `ps.list` from `04__DataPreProcess.R` (element `TimeFinal`).
# Expected output:  `Results/02__Composition/Figures/`, `Tables/` (e.g. `genus_mean_relative_abundance_global.csv`),
#   `Stats/composition__gut__bundle.rds`.
#
# Manuscript Rmd: load bundle without re-running this script:
#   comp <- readRDS(here::here("Results", "02__Composition", "Stats", "composition__gut__bundle.rds"))

# --- Ensure project root (Sieler2026/) ----------------------------------------
init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
if (!file.exists(init_rel)) {
  stop(
    "Set working directory to the Sieler2026 repository root, or run:\n",
    "  Rscript Code/01__Analysis/02__Composition.R\n",
    "from the project root so `", init_rel, "` resolves."
  )
}
source(init_rel)

# --- Paths --------------------------------------------------------------------
path_res_comp <- file.path(path.results, "02__Composition")
path_fig <- file.path(path_res_comp, "Figures")
path_tbl <- file.path(path_res_comp, "Tables")
path_stats <- file.path(path_res_comp, "Stats")
dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(path_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(path_stats, recursive = TRUE, showWarnings = FALSE)

# Archive previous outputs (copy → verify → clear) so each run is clean but reproducible.
sieler2026_archive_module_outputs(
  path_res_module = path_res_comp,
  module_name = "02__Composition",
  subdirs = c("Figures", "Tables", "Stats"),
  clear_original = TRUE
)

bundle_rds <- file.path(path_stats, "composition__gut__bundle.rds")

# --- Constants ----------------------------------------------------------------
ps_list_element <- "TimeFinal"
n_perm <- 999L
# Cap workers for microViz PERMANOVA permutations (avoid oversubscription on laptops).
n_proc <- min(8L, max(1L, parallel::detectCores()))

if (!exists("ps.list", inherits = TRUE)) {
  stop("ps.list not found. Run 04__DataPreProcess.R and ensure ps-list__*.rds exists under Data/r_objects/.")
}

ps_final <- ps.list[[ps_list_element]]

# Human-readable model text for gt subtitles and bundle meta.
model_subtitle_permanova_history <- paste(
  "PERMANOVA on genus-level Bray-Curtis and Canberra distances (microViz::dist_permanova):",
  "distance ~ HistoryLevelNum.",
  paste0(n_perm, " permutations; seed 42.")
)
model_subtitle_permanova_history_parasite <- paste(
  "PERMANOVA on genus-level Bray-Curtis and Canberra distances (microViz::dist_permanova / vegan::adonis2):",
  "distance ~ HistoryLevelNum * Parasite (main effects + interaction; multivariate analogue of the alpha-diversity",
  "GLMM y ~ HistoryLevelNum * Parasite + (1 | Tank.ID), without random effects).",
  "adonis2 by = terms (sequential sums of squares; terms added in order: HistoryLevelNum, Parasite, interaction).",
  paste0(n_perm, " permutations; seed 42.")
)
model_subtitle_permanova_parasite_within_history <- paste(
  "Stratified PERMANOVA: within each prior stressor count (HistoryLevelNum = 0, 1, 2), distance ~ Parasite.",
  paste0(n_perm, " permutations; seed 42. p_fdr_bh: Benjamini–Hochberg across the three strata per distance metric.")
)
model_subtitle_betadisper_parasite_within_history <- paste(
  "Stratified betadisper ANOVA: grouping by Parasite within each HistoryLevelNum stratum.",
  "Rows: Groups vs Residuals; p_fdr_bh is across the three strata per distance metric."
)

# --- Factorial exposure regimes (A × T × P; inferential models) ----------------
# Antibiotics / Temperature / Parasite are 0/1 in sample_data; use factors for categorical contrasts.
ps_final_fac <- ps_final %>%
  microViz::ps_mutate(
    Antibiotics_f = factor(Antibiotics, levels = c(0, 1), labels = c("A-", "A+")),
    Temperature_f = factor(Temperature, levels = c(0, 1), labels = c("T-", "T+")),
    Parasite_f = factor(Parasite, levels = c(0, 1), labels = c("Unexposed", "Exposed")),
    # Use the canonical 8-level exposure-regime factor for plotting/dispersion palettes.
    ATP_group = factor(as.character(Treatment), levels = treatment_order)
  )

model_subtitle_permanova_factorial_atp <- paste(
  "PERMANOVA on genus-level Bray-Curtis and Canberra distances (microViz::dist_permanova / vegan::adonis2):",
  "distance ~ Antibiotics * Temperature * Parasite (A×T×P; term-level tests via adonis2 by = terms).",
  paste0(n_perm, " permutations; seed 42.")
)

model_subtitle_betadisper_factorial_atp <- paste(
  "ANOVA on distance-to-centroid (betadisper); grouping: Antibiotics × Temperature × Parasite (8 levels).",
  "Rows: Groups vs Residuals."
)

# --- Fig. A: relative abundance by exposure regime (merged samples) ----------
# Aggregate to compositional abundances per treatment (mean of merged samples), then top genera barplot.
ps_genus_comp <- ps_final %>%
  microViz::tax_agg(rank = "Genus") %>%
  microViz::tax_transform("compositional")

ps_treatment_agg <- phyloseq::merge_samples(ps_genus_comp, "Treatment", fun = mean)
sample_data_treatment <- methods::as(
  phyloseq::sample_data(ps_treatment_agg),
  "data.frame"
) %>%
  tibble::rownames_to_column(var = "SampleID") %>%
  dplyr::mutate(
    Treatment = as.character(.data$SampleID),
    PriorStressorHistory = dplyr::case_when(
      .data$Treatment %in% c("A- T- P-", "A- T- P+") ~ "None",
      .data$Treatment %in% c("A+ T- P-", "A+ T- P+", "A- T+ P-", "A- T+ P+") ~ "One",
      .data$Treatment %in% c("A+ T+ P-", "A+ T+ P+") ~ "Two",
      TRUE ~ "Unknown"
    ),
    PriorStressorHistory = factor(
      .data$PriorStressorHistory,
      levels = c("None", "One", "Two", "Unknown")
    ),
    Treatment = factor(.data$Treatment, levels = treatment_order)
  ) %>%
  dplyr::filter(.data$PriorStressorHistory != "Unknown")
rownames(sample_data_treatment) <- sample_data_treatment$SampleID
sample_data_treatment$SampleID <- NULL
phyloseq::sample_data(ps_treatment_agg) <- phyloseq::sample_data(sample_data_treatment)

# --- Fig. A helpers: global genus means (CSV + palette ordering) ----------------
# Bar plot taxa: top 8 genera (by summed compositional abundance on treatment means) plus
# Flavobacterium if it is not already in that top 8; all remaining genera -> other_genus_label.
n_genus_core_bar <- 8L
forced_genus_bar <- "Flavobacterium"
other_genus_label <- "Other"

otu_genus_comp <- phyloseq::otu_table(ps_genus_comp) %>%
  methods::as("matrix")
if (!phyloseq::taxa_are_rows(ps_genus_comp)) {
  otu_genus_comp <- t(otu_genus_comp)
}
genus_labels_comp <- as.character(phyloseq::tax_table(ps_genus_comp)[, "Genus"])
# Sum compositional rows that share a genus label, then mean across samples (true genus totals).
genus_mat_global <- rowsum(otu_genus_comp, group = genus_labels_comp)
tbl_genus_global_mean <- tibble::tibble(
  Genus = rownames(genus_mat_global),
  mean_relative_abundance = rowMeans(genus_mat_global, na.rm = TRUE)
) %>%
  dplyr::arrange(dplyr::desc(.data$mean_relative_abundance)) %>%
  dplyr::mutate(rank = dplyr::row_number())

readr::write_csv(
  tbl_genus_global_mean,
  file.path(path_tbl, "genus_mean_relative_abundance_global.csv")
)

# Focal genera colors (match 03__DiffAbund.R); each hex must appear in palette_pool_genus_bar below.
focal_genus_colors <- c(
  Culicoidibacter = "#E41A1C",
  Shewanella = "#377EB8",
  Flavobacterium = "#4DAF4A",
  Cetobacterium = "#984EA3"
)
# Curated palette order (non-focal taxa draw from this list after focal + Other are reserved).
palette_pool_genus_bar <- c(
  "#A6CEE3",
  "#FC8D62",
  "#E41A1C",
  "#984EA3",
  "#377EB8",
  "#FFFF99",
  "#E5C494",
  "#FB9A99",
  "#4DAF4A",
  "grey85"
)

other_hex_genus_bar <- "grey85"
stopifnot(
  identical(other_hex_genus_bar, palette_pool_genus_bar[length(palette_pool_genus_bar)]),
  all(unname(focal_genus_colors) %in% palette_pool_genus_bar)
)

# Optional: pool swatches withheld from non-focal when a focal is present (hue clashes). Empty with the
# curated 10-color palette so five non-focal slots can still use FB9A99 / A6CEE3 / etc. as needed.
focal_genus_bar_pool_conflicts <- list()

# If any two taxa still share a hex after pool assignment, non-focal taxa yield to focal/Other first.
spare_hex_genus_bar_dedupe <- c(
  "#1B9E77",
  "#D95F02",
  "#A6761D",
  "#7570B3",
  "#E6AB02",
  "#A6D854",
  "#E7298A",
  "#CE1256",
  "#6A3D9A"
)

# Pre-merge to top 8 + Flavobacterium (+ Other); comp_barplot n_taxa must equal ntaxa() so nothing is re-cut.
ps_bar_g <- ps_treatment_agg %>%
  microViz::ps_get() %>%
  microViz::tax_agg(rank = "Genus", add_unique = TRUE) %>%
  microViz::ps_get() %>%
  microViz::tax_sort(by = sum, use_counts = TRUE)

ps_bar_taxa <- phyloseq::taxa_names(ps_bar_g)
n_core_take <- min(n_genus_core_bar, length(ps_bar_taxa))
top_core <- ps_bar_taxa[seq_len(n_core_take)]
keep_genera_bar <- unique(c(top_core, forced_genus_bar))
keep_genera_bar <- keep_genera_bar[keep_genera_bar %in% ps_bar_taxa]

otu_bar <- phyloseq::otu_table(ps_bar_g) %>%
  methods::as("matrix")
if (!phyloseq::taxa_are_rows(ps_bar_g)) {
  otu_bar <- t(otu_bar)
}
is_keep_bar <- rownames(otu_bar) %in% keep_genera_bar
keep_otu_bar <- otu_bar[is_keep_bar, , drop = FALSE]
other_vec_bar <- colSums(otu_bar[!is_keep_bar, , drop = FALSE])
other_row_bar <- matrix(
  other_vec_bar,
  nrow = 1L,
  dimnames = list(other_genus_label, colnames(keep_otu_bar))
)
otu_bar_merged <- rbind(keep_otu_bar, other_row_bar)
genus_order_bar <- order(-rowSums(otu_bar_merged), rownames(otu_bar_merged))
otu_bar_merged <- otu_bar_merged[genus_order_bar, , drop = FALSE]

tt_bar_src <- phyloseq::tax_table(ps_bar_g) %>%
  methods::as("matrix")
rank_names_bar <- colnames(tt_bar_src)
tax_bar_merged <- matrix(
  NA_character_,
  nrow = nrow(otu_bar_merged),
  ncol = ncol(tt_bar_src),
  dimnames = list(rownames(otu_bar_merged), rank_names_bar)
)
for (gi in rownames(otu_bar_merged)) {
  if (identical(gi, other_genus_label)) {
    tax_bar_merged[gi, ] <- tt_bar_src[1L, ]
    tax_bar_merged[gi, "Genus"] <- other_genus_label
  } else {
    src_i <- match(gi, rownames(tt_bar_src))
    tax_bar_merged[gi, ] <- tt_bar_src[src_i, ]
  }
}

ps_bar_for_plot <- phyloseq::phyloseq(
  phyloseq::otu_table(otu_bar_merged, taxa_are_rows = TRUE),
  phyloseq::tax_table(tax_bar_merged),
  phyloseq::sample_data(ps_bar_g)
)

# Mirror microViz::comp_barplot post-agg level names so palette matches the rendered fill.
ps_bar_top_work <- ps_bar_for_plot %>%
  microViz::ps_get() %>%
  microViz::tax_agg(rank = "Genus", add_unique = TRUE) %>%
  microViz::ps_get() %>%
  microViz::tax_sort(by = sum, use_counts = TRUE)
n_bar_taxa <- phyloseq::ntaxa(ps_bar_top_work)
phyloseq::tax_table(ps_bar_top_work) <- microViz:::tt_add_topN_var(
  phyloseq::tax_table(ps_bar_top_work),
  N = n_bar_taxa,
  other = other_genus_label,
  varname = ".top"
)
ps_bar_top_agg <- ps_bar_top_work %>%
  microViz::tax_agg(rank = ".top", force = TRUE, add_unique = TRUE) %>%
  microViz::ps_get()
genus_bar_levels <- unique(as.character(phyloseq::tax_table(ps_bar_top_agg)[, ".top"]))

focal_in_bar <- intersect(names(focal_genus_colors), genus_bar_levels)
focal_hex_used <- unique(unname(focal_genus_colors[focal_in_bar]))
pool_hex_blocked_near_focal <- unique(unlist(
  focal_genus_bar_pool_conflicts[
    intersect(focal_in_bar, names(focal_genus_bar_pool_conflicts))
  ],
  use.names = FALSE
))
remaining_pool_genus_bar <- palette_pool_genus_bar[
  !palette_pool_genus_bar %in% unique(c(
    focal_hex_used,
    other_hex_genus_bar,
    pool_hex_blocked_near_focal
  ))
]
nonfocal_in_bar <- setdiff(
  genus_bar_levels,
  c(other_genus_label, focal_in_bar)
)
mean_lookup <- stats::setNames(
  tbl_genus_global_mean$mean_relative_abundance,
  tbl_genus_global_mean$Genus
)
nf_order <- order(
  -unname(mean_lookup[nonfocal_in_bar]),
  nonfocal_in_bar
)
nonfocal_sorted <- nonfocal_in_bar[nf_order]
if (length(remaining_pool_genus_bar) < length(nonfocal_sorted)) {
  remaining_pool_genus_bar <- c(remaining_pool_genus_bar, spare_hex_genus_bar_dedupe)
}
if (length(remaining_pool_genus_bar) < length(nonfocal_sorted)) {
  stop(
    "Not enough distinct palette colors for non-focal genera in the bar plot. ",
    "Extend palette_pool_genus_bar or reduce n_genus_core_bar."
  )
}
comp_palette_genus_bar <- stats::setNames(
  rep(NA_character_, length(genus_bar_levels)),
  genus_bar_levels
)
comp_palette_genus_bar[[other_genus_label]] <- other_hex_genus_bar
comp_palette_genus_bar[focal_in_bar] <- unname(focal_genus_colors[focal_in_bar])
comp_palette_genus_bar[nonfocal_sorted] <- remaining_pool_genus_bar[seq_along(nonfocal_sorted)]
if (any(is.na(comp_palette_genus_bar))) {
  stop("Internal error: incomplete comp_palette_genus_bar for genus bar plot.")
}

# Enforce unique hex per taxon (focal + Other fixed; reassign duplicate non-focal from spares).
dedupe_genus_bar_palette <- function(
    pal,
    genus_order,
    focal_nms,
    other_nm,
    spare
) {
  pal_out <- pal
  for (.iter in seq_len(30L)) {
    vals <- unname(pal_out)
    dup_mask <- duplicated(vals) | duplicated(vals, fromLast = TRUE)
    if (!any(dup_mask)) {
      break
    }
    dup_vals <- unique(vals[dup_mask])
    for (hx in dup_vals) {
      nm_h <- names(pal_out)[unname(pal_out) == hx]
      keep_nm <- NA_character_
      if (any(nm_h %in% focal_nms)) {
        keep_nm <- nm_h[nm_h %in% focal_nms][[1L]]
      } else if (other_nm %in% nm_h) {
        keep_nm <- other_nm
      } else {
        ord <- match(nm_h, genus_order)
        keep_nm <- nm_h[which.min(replace(ord, is.na(ord), Inf))]
      }
      for (nm in setdiff(nm_h, keep_nm)) {
        if (nm %in% focal_nms) {
          next
        }
        if (identical(nm, other_nm)) {
          next
        }
        cand <- spare[!spare %in% unname(pal_out)]
        if (length(cand) < 1L) {
          stop("Ran out of spare hex colors for genus bar palette deduplication.")
        }
        pal_out[[nm]] <- cand[[1L]]
      }
    }
  }
  pal_out
}

comp_palette_genus_bar <- dedupe_genus_bar_palette(
  comp_palette_genus_bar,
  genus_bar_levels,
  focal_in_bar,
  other_genus_label,
  spare_hex_genus_bar_dedupe
)
if (length(unique(unname(comp_palette_genus_bar))) != length(comp_palette_genus_bar)) {
  stop("Internal error: genus bar palette still has duplicate hex values after deduplication.")
}

p_rel_abund <- ps_bar_for_plot %>%
  microViz::ps_get() %>%
  microViz::comp_barplot(
    tax_level = "Genus",
    n_taxa = n_bar_taxa,
    sample_order = treatment_order,
    bar_outline_colour = "black",
    merge_other = TRUE,
    other_name = other_genus_label,
    palette = comp_palette_genus_bar
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(PriorStressorHistory),
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Relative Abundance by exposure regime (genus level)",
    x = "Exposure regime",
    y = "Mean relative abundance",
    tag = "Prior Stressor History"
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3L, byrow = TRUE)) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    axis.title.x = ggplot2::element_text(size = 25),
    axis.title.y = ggplot2::element_text(size = 25),
    legend.title = ggplot2::element_text(size = 20, face = "bold"),
    legend.text = ggplot2::element_text(size = 18),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.width = grid::unit(0.9, "lines"),
    legend.key.height = grid::unit(0.9, "lines"),
    strip.placement = "outside",
    strip.text.y.right = ggplot2::element_text(angle = 270),
    plot.tag = ggplot2::element_text(
      angle = 270,
      size = 25,
      face = "bold",
      margin = ggplot2::margin(l = 6)
    ),
    plot.tag.position = c(1.025, 0.62),
    plot.margin = ggplot2::margin(8, 52, 24, 8)
  ) +
  theme_sieler2026_composition_layers()

fig_w_bar <- 10
fig_h_bar <- 9
ggplot2::ggsave(
  file.path(path_fig, "genus_relative_abundance_by_treatment.pdf"),
  p_rel_abund,
  width = fig_w_bar,
  height = fig_h_bar,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "genus_relative_abundance_by_treatment.png"),
  p_rel_abund,
  width = fig_w_bar,
  height = fig_h_bar,
  dpi = 300
)

# --- PERMANOVA: stress history and stress history + parasite (Bray + Canberra) -
set.seed(42)
permanova_stress_bray <- ps_final %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("bray") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("HistoryLevelNum"),
    n_perms = n_perm,
    n_processes = n_proc
  )

set.seed(42)
permanova_stress_canberra <- ps_final %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("canberra") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("HistoryLevelNum"),
    n_perms = n_perm,
    n_processes = n_proc
  )

# Table D: full factorial on distances — same model as variables = "HistoryLevelNum * Parasite"
# (see microViz dist_permanova: https://david-barnett.github.io/microViz/reference/dist_permanova.html).
# Using variables + interactions matches the package example (PERM3); identical to a single RHS string.
# by = "terms" -> vegan::adonis2 sequential SS so HistoryLevelNum, Parasite, and interaction rows appear.
set.seed(42)
permanova_stress_parasite_bray <- ps_final %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("bray") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("HistoryLevelNum", "Parasite"),
    interactions = "HistoryLevelNum * Parasite",
    n_perms = n_perm,
    n_processes = n_proc,
    by = "terms"
  )

set.seed(42)
permanova_stress_parasite_canberra <- ps_final %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("canberra") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("HistoryLevelNum", "Parasite"),
    interactions = "HistoryLevelNum * Parasite",
    n_perms = n_perm,
    n_processes = n_proc,
    by = "terms"
  )

tbl_perma_hist <- dplyr::bind_rows(
  microviz_permanova_to_tidy(permanova_stress_bray, "Bray-Curtis"),
  microviz_permanova_to_tidy(permanova_stress_canberra, "Canberra")
)
readr::write_csv(tbl_perma_hist, file.path(path_tbl, "permanova_stress_history__bray_canberra.csv"))

tbl_perma_hp <- dplyr::bind_rows(
  microviz_permanova_to_tidy(permanova_stress_parasite_bray, "Bray-Curtis"),
  microviz_permanova_to_tidy(permanova_stress_parasite_canberra, "Canberra")
)
readr::write_csv(tbl_perma_hp, file.path(path_tbl, "permanova_stress_history_parasite__bray_canberra.csv"))

gt_perma_hist <- tbl_perma_hist %>%
  composition_inferential_gt(
    title = "Table C - PERMANOVA: stress history",
    subtitle = model_subtitle_permanova_history,
    alpha = 0.05
  ) %>%
  gt::tab_footnote(
    footnote = "Rows list model terms; p_value is permutation p-value (Pr(>F)).",
    locations = gt::cells_column_labels(columns = p_value)
  )

gt_perma_hp <- tbl_perma_hp %>%
  composition_inferential_gt(
    title = "Table D - PERMANOVA: stress history × parasite (interaction model)",
    subtitle = model_subtitle_permanova_history_parasite,
    alpha = 0.05
  ) %>%
  gt::tab_footnote(
    footnote = "Rows list model terms; p_value is permutation p-value (Pr(>F)).",
    locations = gt::cells_column_labels(columns = p_value)
  )

gt::gtsave(gt_perma_hist, file.path(path_tbl, "permanova_stress_history.html"))
gt::gtsave(gt_perma_hp, file.path(path_tbl, "permanova_stress_history_parasite.html"))

# --- PERMANOVA: factorial exposure regimes (A × T × P; Bray + Canberra) --------
set.seed(42)
permanova_factorial_bray <- ps_final_fac %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("bray") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("Antibiotics_f", "Temperature_f", "Parasite_f"),
    interactions = "Antibiotics_f * Temperature_f * Parasite_f",
    n_perms = n_perm,
    n_processes = n_proc,
    by = "terms"
  )

set.seed(42)
permanova_factorial_canberra <- ps_final_fac %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("canberra") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("Antibiotics_f", "Temperature_f", "Parasite_f"),
    interactions = "Antibiotics_f * Temperature_f * Parasite_f",
    n_perms = n_perm,
    n_processes = n_proc,
    by = "terms"
  )

tbl_perma_atp <- dplyr::bind_rows(
  microviz_permanova_to_tidy(permanova_factorial_bray, "Bray-Curtis"),
  microviz_permanova_to_tidy(permanova_factorial_canberra, "Canberra")
) %>%
  dplyr::mutate(
    Term = dplyr::recode(
      .data$Term,
      "Antibiotics_f" = "Antibiotics",
      "Temperature_f" = "Temperature",
      "Parasite_f" = "Parasite",
      "Antibiotics_f:Temperature_f" = "Antibiotics:Temperature",
      "Antibiotics_f:Parasite_f" = "Antibiotics:Parasite",
      "Temperature_f:Parasite_f" = "Temperature:Parasite",
      "Antibiotics_f:Temperature_f:Parasite_f" = "Antibiotics:Temperature:Parasite"
    )
  )

readr::write_csv(tbl_perma_atp, file.path(path_tbl, "permanova_factorial_ATP__bray_canberra.csv"))

gt_perma_atp <- tbl_perma_atp %>%
  composition_inferential_gt(
    title = "PERMANOVA: factorial exposure regimes (Antibiotics × Temperature × Parasite)",
    subtitle = model_subtitle_permanova_factorial_atp,
    alpha = 0.05
  ) %>%
  gt::tab_footnote(
    footnote = "Rows list model terms; p_value is permutation p-value (Pr(>F)) from adonis2 (by = terms).",
    locations = gt::cells_column_labels(columns = p_value)
  )
gt::gtsave(gt_perma_atp, file.path(path_tbl, "permanova_factorial_ATP.html"))

# --- Bray distances: PCoA + betadisper (legacy microViz pipelines) ----------------
# PCoA: stress history (color by History_Label, ellipses).
set.seed(42)
p_pcoa_bray_history <- ps_final %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("bray") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("HistoryLevel"),
    n_perms = n_perm,
    n_processes = n_proc
  ) %>%
  microViz::ord_calc(method = "PCoA") %>%
  microViz::ord_plot(color = "History_Label") +
  ggplot2::scale_color_manual(
    values = history_color_scale,
    name = "Prior Stressor History",
    labels = c(
      "No prior stressors" = "None",
      "One prior stressor" = "One",
      "Two prior stressors" = "Two"
    )
  ) +
  ggnewscale::new_scale_color() +
  ggplot2::stat_ellipse(ggplot2::aes(color = History_Label), linewidth = 1, alpha = 0.7) +
  ggplot2::scale_color_manual(
    values = history_color_scale,
    name = "Prior Stressor History",
    labels = c(
      "No prior stressors" = "None",
      "One prior stressor" = "One",
      "Two prior stressors" = "Two"
    ),
    guide = "none"
  ) +
  ggplot2::coord_fixed(ratio = 1, expand = FALSE, clip = "off") +
  ggplot2::labs(title = "PCoA - stressor history (Bray-Curtis)") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.margin = ggplot2::margin(0, 0, 0, 6, "pt"),
    legend.spacing.y = grid::unit(4, "pt"),
    legend.key.height = grid::unit(0.5, "cm"),
    plot.margin = ggplot2::margin(t = 10, r = 8, b = 6, l = 8, unit = "pt"),
    plot.caption = ggplot2::element_text(
      hjust = 1,
      size = rel(0.68),
      colour = "grey35",
      margin = ggplot2::margin(t = 2, r = 2, b = 0, l = 2, unit = "pt")
    )
  ) +
  theme_sieler2026_composition_layers()

# PCoA: stress history with parasite stratification (point fill + ellipse linetype).
set.seed(42)
p_pcoa_bray_hp <- ps_final %>%
  microViz::ps_mutate(
    Fill_Group = dplyr::if_else(Parasite_Exposed == "Unexposed", "white", as.character(History_Label))
  ) %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("bray") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("HistoryLevel"),
    n_perms = n_perm,
    n_processes = n_proc
  ) %>%
  microViz::ord_calc(method = "PCoA") %>%
  microViz::ord_plot(
    color = "History_Label",
    fill = "Fill_Group",
    shape = "Parasite_Exposed",
    stroke = 1,
    size = 2
  ) +
  ggplot2::scale_color_manual(
    values = history_color_scale,
    name = "Prior Stressor History",
    labels = c(
      "No prior stressors" = "None",
      "One prior stressor" = "One",
      "Two prior stressors" = "Two"
    ),
    guide = ggplot2::guide_legend(order = 1L)
  ) +
  ggplot2::scale_fill_manual(
    values = c("white" = "#FFFFFF", history_color_scale),
    guide = "none"
  ) +
  ggplot2::scale_shape_manual(
    values = c("Unexposed" = 21L, "Exposed" = 23L),
    name = "Parasite Exposure",
    # Merged with linetype scale (same title): combined keys show point + line (legend fills: white circle, black diamond).
    guide = ggplot2::guide_legend(
      order = 2L,
      override.aes = list(
        fill = c("#FFFFFF", "#000000"),
        colour = c("black", "black"),
        stroke = c(0.9, 0.9),
        size = 3.2
      )
    )
  ) +
  ggnewscale::new_scale_color() +
  ggplot2::stat_ellipse(
    ggplot2::aes(color = History_Label, linetype = Parasite_Exposed),
    linewidth = 1,
    alpha = 0.7
  ) +
  ggplot2::scale_color_manual(
    values = history_color_scale,
    name = "Prior Stressor History",
    labels = c(
      "No prior stressors" = "None",
      "One prior stressor" = "One",
      "Two prior stressors" = "Two"
    ),
    guide = "none"
  ) +
  ggplot2::scale_linetype_manual(
    values = c("Unexposed" = "solid", "Exposed" = "dashed"),
    name = "Parasite Exposure",
    guide = ggplot2::guide_legend(
      order = 2L,
      override.aes = list(linewidth = 0.32, linetype = c("solid", "dashed"))
    )
  ) +
  ggplot2::coord_fixed(ratio = 1, expand = FALSE, clip = "off") +
  ggplot2::labs(title = "PCoA - stressor history and parasite exposure (Bray-Curtis)") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(6, "pt"),
    legend.key.width = grid::unit(1.6, "cm"),
    legend.key.height = grid::unit(0.55, "cm"),
    plot.margin = ggplot2::margin(t = 10, r = 12, b = 48, l = 10, unit = "pt"),
    plot.caption = ggplot2::element_text(
      hjust = 1,
      size = rel(0.68),
      colour = "grey35",
      margin = ggplot2::margin(t = 2, r = 2, b = 0, l = 2, unit = "pt")
    )
  ) +
  theme_sieler2026_composition_layers()

fig_w <- 8
fig_h <- 8
ggplot2::ggsave(file.path(path_fig, "pcoa_bray_stress_history.pdf"), p_pcoa_bray_history, width = fig_w, height = fig_h, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "pcoa_bray_stress_history.png"), p_pcoa_bray_history, width = fig_w, height = fig_h, dpi = 300)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_bray_stress_history_parasite.pdf"),
  p_pcoa_bray_hp,
  width = 9,
  height = 9,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_bray_stress_history_parasite.png"),
  p_pcoa_bray_hp,
  width = 9,
  height = 9,
  dpi = 300
)

# Betadisper on Bray-Curtis distances.
disp_stress_bray <- ps_final %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("bray") %>%
  microViz::dist_bdisp(variables = c("HistoryLevel")) %>%
  microViz::bdisp_get()

disp_stress_parasite_bray <- ps_final %>%
  microViz::ps_mutate(History.Parasite = paste0(HistoryLevel, "_", Parasite)) %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("bray") %>%
  microViz::dist_bdisp(variables = c("History.Parasite")) %>%
  microViz::bdisp_get()

p_bdisp_hist <- composition_betadisper_boxplot_history(
  disp_stress_bray,
  history_color_scale_vec = history_color_scale,
  title_suffix = "Bray-Curtis"
)

p_bdisp_hp <- composition_betadisper_boxplot_history_parasite(
  disp_stress_parasite_bray,
  history_color_scale_vec = history_color_scale,
  title_suffix = "Bray-Curtis"
)

ggplot2::ggsave(file.path(path_fig, "betadisper_bray_stress_history.pdf"), p_bdisp_hist, width = fig_w, height = fig_h, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "betadisper_bray_stress_history.png"), p_bdisp_hist, width = fig_w, height = fig_h, dpi = 300)
ggplot2::ggsave(
  file.path(path_fig, "betadisper_bray_stress_history_parasite.pdf"),
  p_bdisp_hp,
  width = 9,
  height = 9,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "betadisper_bray_stress_history_parasite.png"),
  p_bdisp_hp,
  width = 9,
  height = 9,
  dpi = 300
)

# --- Factorial exposure regimes: PCoA + betadisper (Bray + Canberra) ------------
set.seed(42)
p_pcoa_bray_atp <- ps_final_fac %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("bray") %>%
  microViz::ord_calc(method = "PCoA") %>%
  microViz::ord_plot(color = "Parasite_f", shape = "Antibiotics_f") +
  ggplot2::scale_color_manual(
    values = c("Unexposed" = "grey60", "Exposed" = "firebrick"),
    name = "Parasite exposure"
  ) +
  ggplot2::facet_grid(rows = ggplot2::vars(Temperature_f), cols = ggplot2::vars(Antibiotics_f)) +
  ggplot2::labs(
    title = "PCoA (Bray-Curtis): factorial exposure regimes (A×T×P)",
    subtitle = "Points colored by parasite exposure; facets show Antibiotics (columns) and Temperature (rows)."
  ) +
  theme_sieler2026_composition_figure(base_size = 14, legend_position = "bottom") +
  ggplot2::theme(
    legend.box = "vertical",
    legend.spacing.y = grid::unit(6, "pt")
  )

set.seed(42)
p_pcoa_canberra_atp <- ps_final_fac %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("canberra") %>%
  microViz::ord_calc(method = "PCoA") %>%
  microViz::ord_plot(color = "Parasite_f", shape = "Antibiotics_f") +
  ggplot2::scale_color_manual(
    values = c("Unexposed" = "grey60", "Exposed" = "firebrick"),
    name = "Parasite exposure"
  ) +
  ggplot2::facet_grid(rows = ggplot2::vars(Temperature_f), cols = ggplot2::vars(Antibiotics_f)) +
  ggplot2::labs(
    title = "PCoA (Canberra): factorial exposure regimes (A×T×P)",
    subtitle = "Points colored by parasite exposure; facets show Antibiotics (columns) and Temperature (rows)."
  ) +
  theme_sieler2026_composition_figure(base_size = 14, legend_position = "bottom") +
  ggplot2::theme(
    legend.box = "vertical",
    legend.spacing.y = grid::unit(6, "pt")
  )

ggplot2::ggsave(file.path(path_fig, "pcoa_bray_factorial_ATP.pdf"), p_pcoa_bray_atp, width = 10, height = 7, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "pcoa_bray_factorial_ATP.png"), p_pcoa_bray_atp, width = 10, height = 7, dpi = 300)
ggplot2::ggsave(file.path(path_fig, "pcoa_canberra_factorial_ATP.pdf"), p_pcoa_canberra_atp, width = 10, height = 7, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "pcoa_canberra_factorial_ATP.png"), p_pcoa_canberra_atp, width = 10, height = 7, dpi = 300)

disp_atp_bray <- ps_final_fac %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("bray") %>%
  microViz::dist_bdisp(variables = c("ATP_group")) %>%
  microViz::bdisp_get()

disp_atp_canberra <- ps_final_fac %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("canberra") %>%
  microViz::dist_bdisp(variables = c("ATP_group")) %>%
  microViz::bdisp_get()

tbl_bdisp_atp <- dplyr::bind_rows(
  composition_betadisper_anova_to_tidy(disp_atp_bray$ATP_group$anova, "Bray-Curtis"),
  composition_betadisper_anova_to_tidy(disp_atp_canberra$ATP_group$anova, "Canberra")
)
readr::write_csv(tbl_bdisp_atp, file.path(path_tbl, "betadisper_anova_factorial_ATP__bray_canberra.csv"))

gt_bdisp_atp <- tbl_bdisp_atp %>%
  composition_inferential_gt(
    title = "Betadisper ANOVA: factorial exposure regimes (A×T×P)",
    subtitle = model_subtitle_betadisper_factorial_atp,
    alpha = 0.05
  ) %>%
  gt::tab_footnote(
    footnote = "p_value is Pr(>F) from the anova() table on betadisper distances.",
    locations = gt::cells_column_labels(columns = p_value)
  )
gt::gtsave(gt_bdisp_atp, file.path(path_tbl, "betadisper_anova_factorial_ATP.html"))

p_bdisp_bray_atp <- composition_betadisper_boxplot_atp(disp_atp_bray, "Bray-Curtis")
p_bdisp_canberra_atp <- composition_betadisper_boxplot_atp(disp_atp_canberra, "Canberra")

ggplot2::ggsave(file.path(path_fig, "betadisper_bray_factorial_ATP.pdf"), p_bdisp_bray_atp, width = 10, height = 7, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "betadisper_bray_factorial_ATP.png"), p_bdisp_bray_atp, width = 10, height = 7, dpi = 300)
ggplot2::ggsave(file.path(path_fig, "betadisper_canberra_factorial_ATP.pdf"), p_bdisp_canberra_atp, width = 10, height = 7, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "betadisper_canberra_factorial_ATP.png"), p_bdisp_canberra_atp, width = 10, height = 7, dpi = 300)

# --- Canberra distances: PCoA + betadisper (mirrors Bray) -----------------------
set.seed(42)
p_pcoa_canberra_history <- ps_final %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("canberra") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("HistoryLevel"),
    n_perms = n_perm,
    n_processes = n_proc
  ) %>%
  microViz::ord_calc(method = "PCoA") %>%
  microViz::ord_plot(color = "History_Label") +
  ggplot2::scale_color_manual(
    values = history_color_scale,
    name = "Prior Stressor History",
    labels = c(
      "No prior stressors" = "None",
      "One prior stressor" = "One",
      "Two prior stressors" = "Two"
    )
  ) +
  ggnewscale::new_scale_color() +
  ggplot2::stat_ellipse(ggplot2::aes(color = History_Label), linewidth = 1, alpha = 0.7) +
  ggplot2::scale_color_manual(
    values = history_color_scale,
    name = "Prior Stressor History",
    labels = c(
      "No prior stressors" = "None",
      "One prior stressor" = "One",
      "Two prior stressors" = "Two"
    ),
    guide = "none"
  ) +
  ggplot2::coord_fixed(ratio = 1, expand = FALSE, clip = "off") +
  ggplot2::labs(title = "PCoA - stressor history (Canberra)") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.margin = ggplot2::margin(0, 0, 0, 6, "pt"),
    legend.spacing.y = grid::unit(4, "pt"),
    legend.key.height = grid::unit(0.5, "cm"),
    plot.margin = ggplot2::margin(t = 10, r = 8, b = 6, l = 8, unit = "pt"),
    plot.caption = ggplot2::element_text(
      hjust = 1,
      size = rel(0.68),
      colour = "grey35",
      margin = ggplot2::margin(t = 2, r = 2, b = 0, l = 2, unit = "pt")
    )
  ) +
  theme_sieler2026_composition_layers()

set.seed(42)
p_pcoa_canberra_hp <- ps_final %>%
  microViz::ps_mutate(
    Fill_Group = dplyr::if_else(Parasite_Exposed == "Unexposed", "white", as.character(History_Label))
  ) %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("canberra") %>%
  microViz::dist_permanova(
    seed = 42,
    variables = c("HistoryLevel"),
    n_perms = n_perm,
    n_processes = n_proc
  ) %>%
  microViz::ord_calc(method = "PCoA") %>%
  microViz::ord_plot(
    color = "History_Label",
    fill = "Fill_Group",
    shape = "Parasite_Exposed",
    stroke = 1,
    size = 2
  ) +
  ggplot2::scale_color_manual(
    values = history_color_scale,
    name = "Prior Stressor History",
    labels = c(
      "No prior stressors" = "None",
      "One prior stressor" = "One",
      "Two prior stressors" = "Two"
    ),
    guide = ggplot2::guide_legend(order = 1L)
  ) +
  ggplot2::scale_fill_manual(
    values = c("white" = "#FFFFFF", history_color_scale),
    guide = "none"
  ) +
  ggplot2::scale_shape_manual(
    values = c("Unexposed" = 21L, "Exposed" = 23L),
    name = "Parasite Exposure",
    guide = ggplot2::guide_legend(
      order = 2L,
      override.aes = list(
        fill = c("#FFFFFF", "#000000"),
        colour = c("black", "black"),
        stroke = c(0.9, 0.9),
        size = 3.2
      )
    )
  ) +
  ggnewscale::new_scale_color() +
  ggplot2::stat_ellipse(
    ggplot2::aes(color = History_Label, linetype = Parasite_Exposed),
    linewidth = 1,
    alpha = 0.7
  ) +
  ggplot2::scale_color_manual(
    values = history_color_scale,
    name = "Prior Stressor History",
    labels = c(
      "No prior stressors" = "None",
      "One prior stressor" = "One",
      "Two prior stressors" = "Two"
    ),
    guide = "none"
  ) +
  ggplot2::scale_linetype_manual(
    values = c("Unexposed" = "solid", "Exposed" = "dashed"),
    name = "Parasite Exposure",
    guide = ggplot2::guide_legend(
      order = 2L,
      override.aes = list(linewidth = 0.32, linetype = c("solid", "dashed"))
    )
  ) +
  ggplot2::coord_fixed(ratio = 1, expand = FALSE, clip = "off") +
  ggplot2::labs(title = "PCoA - stressor history and parasite exposure (Canberra)") +
  ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(6, "pt"),
    legend.key.width = grid::unit(1.6, "cm"),
    legend.key.height = grid::unit(0.55, "cm"),
    plot.margin = ggplot2::margin(t = 10, r = 12, b = 48, l = 10, unit = "pt"),
    plot.caption = ggplot2::element_text(
      hjust = 1,
      size = rel(0.68),
      colour = "grey35",
      margin = ggplot2::margin(t = 2, r = 2, b = 0, l = 2, unit = "pt")
    )
  ) +
  theme_sieler2026_composition_layers() +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(
      colour = grDevices::grey(0.87),
      linewidth = 0.35
    ),
    panel.grid.minor.x = ggplot2::element_line(
      colour = grDevices::grey(0.93),
      linewidth = 0.25
    ),
    panel.grid.minor.y = ggplot2::element_line(
      colour = grDevices::grey(0.93),
      linewidth = 0.25
    )
  )

ggplot2::ggsave(file.path(path_fig, "pcoa_canberra_stress_history.pdf"), p_pcoa_canberra_history, width = fig_w, height = fig_h, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "pcoa_canberra_stress_history.png"), p_pcoa_canberra_history, width = fig_w, height = fig_h, dpi = 300)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_canberra_stress_history_parasite.pdf"),
  p_pcoa_canberra_hp,
  width = 9,
  height = 9,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_canberra_stress_history_parasite.png"),
  p_pcoa_canberra_hp,
  width = 9,
  height = 9,
  dpi = 300
)

disp_stress_canberra <- ps_final %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("canberra") %>%
  microViz::dist_bdisp(variables = c("HistoryLevel")) %>%
  microViz::bdisp_get()

disp_stress_parasite_canberra <- ps_final %>%
  microViz::ps_mutate(History.Parasite = paste0(HistoryLevel, "_", Parasite)) %>%
  microViz::tax_transform("identity", rank = "Genus") %>%
  microViz::dist_calc("canberra") %>%
  microViz::dist_bdisp(variables = c("History.Parasite")) %>%
  microViz::bdisp_get()

p_bdisp_canberra_hist <- composition_betadisper_boxplot_history(
  disp_stress_canberra,
  history_color_scale_vec = history_color_scale,
  title_suffix = "Canberra"
)

p_bdisp_canberra_hp <- composition_betadisper_boxplot_history_parasite(
  disp_stress_parasite_canberra,
  history_color_scale_vec = history_color_scale,
  title_suffix = "Canberra"
)

ggplot2::ggsave(file.path(path_fig, "betadisper_canberra_stress_history.pdf"), p_bdisp_canberra_hist, width = fig_w, height = fig_h, device = "pdf")
ggplot2::ggsave(file.path(path_fig, "betadisper_canberra_stress_history.png"), p_bdisp_canberra_hist, width = fig_w, height = fig_h, dpi = 300)
ggplot2::ggsave(
  file.path(path_fig, "betadisper_canberra_stress_history_parasite.pdf"),
  p_bdisp_canberra_hp,
  width = 9,
  height = 9,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "betadisper_canberra_stress_history_parasite.png"),
  p_bdisp_canberra_hp,
  width = 9,
  height = 9,
  dpi = 300
)

# --- Betadisper ANOVA tables (Bray + Canberra) ---------------------------------
model_subtitle_betadisper_history <- paste(
  "ANOVA on distance-to-centroid (betadisper / vegan); grouping: HistoryLevel.",
  "Rows: Groups vs Residuals."
)
model_subtitle_betadisper_parasite <- paste(
  "ANOVA on distance-to-centroid (betadisper); grouping: HistoryLevel × Parasite (six levels)."
)

tbl_bdisp_hist <- dplyr::bind_rows(
  composition_betadisper_anova_to_tidy(disp_stress_bray$HistoryLevel$anova, "Bray-Curtis"),
  composition_betadisper_anova_to_tidy(disp_stress_canberra$HistoryLevel$anova, "Canberra")
)
readr::write_csv(tbl_bdisp_hist, file.path(path_tbl, "betadisper_anova_stress_history__bray_canberra.csv"))

tbl_bdisp_hp <- dplyr::bind_rows(
  composition_betadisper_anova_to_tidy(disp_stress_parasite_bray$History.Parasite$anova, "Bray-Curtis"),
  composition_betadisper_anova_to_tidy(disp_stress_parasite_canberra$History.Parasite$anova, "Canberra")
)
readr::write_csv(tbl_bdisp_hp, file.path(path_tbl, "betadisper_anova_stress_history_parasite__bray_canberra.csv"))

gt_bdisp_hist <- tbl_bdisp_hist %>%
  composition_inferential_gt(
    title = "Betadisper ANOVA: stress history (homogeneity of dispersions)",
    subtitle = model_subtitle_betadisper_history,
    alpha = 0.05
  ) %>%
  gt::tab_footnote(
    footnote = "p_value is Pr(>F) from the anova() table on betadisper distances.",
    locations = gt::cells_column_labels(columns = p_value)
  )

gt_bdisp_hp <- tbl_bdisp_hp %>%
  composition_inferential_gt(
    title = "Betadisper ANOVA: stress history × parasite",
    subtitle = model_subtitle_betadisper_parasite,
    alpha = 0.05
  ) %>%
  gt::tab_footnote(
    footnote = "p_value is Pr(>F) from the anova() table on betadisper distances.",
    locations = gt::cells_column_labels(columns = p_value)
  )

gt::gtsave(gt_bdisp_hist, file.path(path_tbl, "betadisper_anova_stress_history.html"))
gt::gtsave(gt_bdisp_hp, file.path(path_tbl, "betadisper_anova_stress_history_parasite.html"))

# --- Stratified: Parasite effect within each prior stressor history (0, 1, 2) ---
# PERMANOVA / betadisper on subsets; faceted PCoA uses one global ordination (cmdscale) for comparability.

stratified_permanova_parasite_row <- function(ps_obj, dist_method, distance_label, h) {
  ps_h <- ps_obj %>%
    microViz::ps_filter(HistoryLevelNum == h)
  sdt <- microViz::samdat_tbl(ps_h)
  n0 <- sum(sdt$Parasite == 0L, na.rm = TRUE)
  n1 <- sum(sdt$Parasite == 1L, na.rm = TRUE)
  n_all <- nrow(sdt)
  if (n_all < 4L || n0 < 1L || n1 < 1L) {
    warning(
      "Stratified PERMANOVA skipped for HistoryLevelNum = ", h,
      " (n = ", n_all, ", P0 = ", n0, ", P1 = ", n1, ")."
    )
    return(NULL)
  }
  set.seed(42)
  perma_h <- ps_h %>%
    microViz::tax_transform("identity", rank = "Genus") %>%
    microViz::dist_calc(dist_method) %>%
    microViz::dist_permanova(
      seed = 42L,
      variables = c("Parasite"),
      n_perms = n_perm,
      n_processes = n_proc
    )
  p_df <- microViz::perm_get(perma_h) %>%
    as.data.frame()
  if (!"Parasite" %in% rownames(p_df)) {
    warning("Stratified PERMANOVA: no Parasite term in adonis output for stratum ", h, ".")
    return(NULL)
  }
  f_col <- intersect(c("F", "F.Model"), colnames(p_df))
  f_col <- if (length(f_col) >= 1L) f_col[[1]] else NA_character_
  f_val <- if (!is.na(f_col)) as.numeric(p_df["Parasite", f_col]) else NA_real_
  tibble::tibble(
    stratum = h,
    Distance = distance_label,
    Term = "Parasite",
    Df = as.numeric(p_df["Parasite", "Df"]),
    R2 = round(as.numeric(p_df["Parasite", "R2"]), 4),
    F = round(f_val, 4),
    p_value = as.numeric(p_df["Parasite", "Pr(>F)"]),
    n_samples = n_all,
    n_parasite_0 = n0,
    n_parasite_1 = n1
  )
}

stratified_betadisper_parasite_row <- function(ps_obj, dist_method, distance_label, h) {
  ps_h <- ps_obj %>%
    microViz::ps_filter(HistoryLevelNum == h)
  sdt <- microViz::samdat_tbl(ps_h)
  n0 <- sum(sdt$Parasite == 0L, na.rm = TRUE)
  n1 <- sum(sdt$Parasite == 1L, na.rm = TRUE)
  n_all <- nrow(sdt)
  if (n_all < 4L || n0 < 1L || n1 < 1L) {
    warning(
      "Stratified betadisper skipped for HistoryLevelNum = ", h,
      " (n = ", n_all, ", P0 = ", n0, ", P1 = ", n1, ")."
    )
    return(NULL)
  }
  set.seed(42)
  # betadisper requires a grouping factor; Parasite is 0/1 numeric in sample_data.
  ps_h <- ps_h %>%
    microViz::ps_mutate(
      Parasite_f = factor(Parasite, levels = c(0L, 1L), labels = c("Unexposed", "Exposed"))
    )
  disp <- ps_h %>%
    microViz::tax_transform("identity", rank = "Genus") %>%
    microViz::dist_calc(dist_method) %>%
    microViz::dist_bdisp(variables = c("Parasite_f")) %>%
    microViz::bdisp_get()
  bd <- disp$Parasite_f
  aov_df <- bd$anova %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Term")
  row_g <- aov_df %>% dplyr::filter(.data$Term == "Groups")
  if (nrow(row_g) != 1L) {
    warning("Stratified betadisper: unexpected ANOVA table for stratum ", h, ".")
    return(NULL)
  }
  p_g <- suppressWarnings(as.numeric(row_g[["Pr(>F)"]]))
  f_col <- intersect(c("F value", "F"), names(row_g))
  f_val <- if (length(f_col) >= 1L) suppressWarnings(as.numeric(row_g[[f_col[1]]])) else NA_real_
  tibble::tibble(
    stratum = h,
    Distance = distance_label,
    Term = "Groups",
    Df = suppressWarnings(as.numeric(row_g[["Df"]])),
    sum_sq = suppressWarnings(as.numeric(row_g[["Sum Sq"]])),
    mean_sq = suppressWarnings(as.numeric(row_g[["Mean Sq"]])),
    F = f_val,
    p_value = p_g,
    n_samples = n_all,
    n_parasite_0 = n0,
    n_parasite_1 = n1
  )
}

stratified_rows_permanova <- function(ps_obj, dist_method, distance_label) {
  out <- list()
  for (h in c(0L, 1L, 2L)) {
    row_h <- stratified_permanova_parasite_row(ps_obj, dist_method, distance_label, h)
    if (!is.null(row_h)) {
      out[[length(out) + 1L]] <- row_h
    }
  }
  if (length(out) == 0L) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(out) %>%
    dplyr::group_by(.data$Distance) %>%
    dplyr::mutate(p_fdr_bh = stats::p.adjust(.data$p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p_value = round(.data$p_value, 4),
      p_fdr_bh = round(.data$p_fdr_bh, 4)
    )
}

stratified_rows_betadisper <- function(ps_obj, dist_method, distance_label) {
  out <- list()
  for (h in c(0L, 1L, 2L)) {
    row_h <- stratified_betadisper_parasite_row(ps_obj, dist_method, distance_label, h)
    if (!is.null(row_h)) {
      out[[length(out) + 1L]] <- row_h
    }
  }
  if (length(out) == 0L) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(out) %>%
    dplyr::group_by(.data$Distance) %>%
    dplyr::mutate(p_fdr_bh = stats::p.adjust(.data$p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p_value = round(.data$p_value, 4),
      p_fdr_bh = round(.data$p_fdr_bh, 4)
    )
}

tbl_perm_parasite_within <- dplyr::bind_rows(
  stratified_rows_permanova(ps_final, "bray", "Bray-Curtis"),
  stratified_rows_permanova(ps_final, "canberra", "Canberra")
)
readr::write_csv(tbl_perm_parasite_within, file.path(path_tbl, "permanova_parasite_within_history__bray_canberra.csv"))

tbl_bdisp_parasite_within <- dplyr::bind_rows(
  stratified_rows_betadisper(ps_final, "bray", "Bray-Curtis"),
  stratified_rows_betadisper(ps_final, "canberra", "Canberra")
)
readr::write_csv(tbl_bdisp_parasite_within, file.path(path_tbl, "betadisper_anova_parasite_within_history__bray_canberra.csv"))

if (nrow(tbl_perm_parasite_within) > 0L) {
  gt_perm_parasite_within <- tbl_perm_parasite_within %>%
    composition_inferential_gt(
      title = "Table E — PERMANOVA: Parasite within each prior stressor history (stratified)",
      subtitle = model_subtitle_permanova_parasite_within_history,
      alpha = 0.05
    ) %>%
    gt::tab_footnote(
      footnote = "p_value is permutation p (Pr(>F)); p_fdr_bh is Benjamini–Hochberg across strata 0–2 within each distance metric.",
      locations = gt::cells_column_labels(columns = p_fdr_bh)
    )
  gt::gtsave(gt_perm_parasite_within, file.path(path_tbl, "permanova_parasite_within_history.html"))
} else {
  gt_perm_parasite_within <- NULL
  warning("Stratified PERMANOVA table is empty (no valid strata).")
}

if (nrow(tbl_bdisp_parasite_within) > 0L) {
  gt_bdisp_parasite_within <- tbl_bdisp_parasite_within %>%
    composition_inferential_gt(
      title = "Table F — Betadisper ANOVA: Parasite within each prior stressor history (stratified)",
      subtitle = model_subtitle_betadisper_parasite_within_history,
      alpha = 0.05
    ) %>%
    gt::tab_footnote(
      footnote = "p_value is Pr(>F) on Groups (dispersion); p_fdr_bh is Benjamini–Hochberg across strata 0–2 within each distance metric.",
      locations = gt::cells_column_labels(columns = p_fdr_bh)
    )
  gt::gtsave(gt_bdisp_parasite_within, file.path(path_tbl, "betadisper_anova_parasite_within_history.html"))
} else {
  gt_bdisp_parasite_within <- NULL
  warning("Stratified betadisper table is empty (no valid strata).")
}

p_pcoa_bray_parasite_faceted <- composition_pcoa_parasite_faceted_by_history_plot(
  ps_final,
  dist_method = "bray",
  title = "PCoA (Bray-Curtis): parasite exposure, faceted by prior stressor history",
  subtitle = "Single ordination for all samples (cmdscale); compare P+ vs P- within each panel."
)
p_pcoa_canberra_parasite_faceted <- composition_pcoa_parasite_faceted_by_history_plot(
  ps_final,
  dist_method = "canberra",
  title = "PCoA (Canberra): parasite exposure, faceted by prior stressor history",
  subtitle = "Single ordination for all samples (cmdscale); compare P+ vs P- within each panel."
)

fig_fac_w <- 12
fig_fac_h <- 5
ggplot2::ggsave(
  file.path(path_fig, "pcoa_bray_parasite_faceted_by_history.pdf"),
  p_pcoa_bray_parasite_faceted,
  width = fig_fac_w,
  height = fig_fac_h,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_bray_parasite_faceted_by_history.png"),
  p_pcoa_bray_parasite_faceted,
  width = fig_fac_w,
  height = fig_fac_h,
  dpi = 300
)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_canberra_parasite_faceted_by_history.pdf"),
  p_pcoa_canberra_parasite_faceted,
  width = fig_fac_w,
  height = fig_fac_h,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_canberra_parasite_faceted_by_history.png"),
  p_pcoa_canberra_parasite_faceted,
  width = fig_fac_w,
  height = fig_fac_h,
  dpi = 300
)

# Faceted PCoA: points by eight exposure regimes (Treatment); ellipses by parasite exposure.
p_pcoa_bray_regime_pts <- composition_pcoa_parasite_faceted_points_by_regime_plot(
  ps_final,
  dist_method = "bray",
  title = "PCoA (Bray-Curtis): exposure regime, faceted by prior stressor history",
  subtitle = "Points: eight A/T/P regimes; ellipses: parasite exposure (same cmdscale as parasite-colored figure)."
)
p_pcoa_canberra_regime_pts <- composition_pcoa_parasite_faceted_points_by_regime_plot(
  ps_final,
  dist_method = "canberra",
  title = "PCoA (Canberra): exposure regime, faceted by prior stressor history",
  subtitle = "Points: eight A/T/P regimes; ellipses: parasite exposure (same cmdscale as parasite-colored figure)."
)

fig_fac_reg_w <- 12
fig_fac_reg_h <- 6.5
ggplot2::ggsave(
  file.path(path_fig, "pcoa_bray_parasite_faceted_by_history_regime_points.pdf"),
  p_pcoa_bray_regime_pts,
  width = fig_fac_reg_w,
  height = fig_fac_reg_h,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_bray_parasite_faceted_by_history_regime_points.png"),
  p_pcoa_bray_regime_pts,
  width = fig_fac_reg_w,
  height = fig_fac_reg_h,
  dpi = 300
)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_canberra_parasite_faceted_by_history_regime_points.pdf"),
  p_pcoa_canberra_regime_pts,
  width = fig_fac_reg_w,
  height = fig_fac_reg_h,
  device = "pdf"
)
ggplot2::ggsave(
  file.path(path_fig, "pcoa_canberra_parasite_faceted_by_history_regime_points.png"),
  p_pcoa_canberra_regime_pts,
  width = fig_fac_reg_w,
  height = fig_fac_reg_h,
  dpi = 300
)

# --- Assemble bundle -----------------------------------------------------------
modules <- list(
  relative_abundance = list(
    treatment = list(
      figure = p_rel_abund,
      paths = c(
        figure_pdf = file.path(path_fig, "genus_relative_abundance_by_treatment.pdf"),
        figure_png = file.path(path_fig, "genus_relative_abundance_by_treatment.png")
      )
    )
  ),
  beta_bray = list(
    factorial_atp = list(
      pcoa = list(
        figure = p_pcoa_bray_atp,
        paths = c(
          figure_pdf = file.path(path_fig, "pcoa_bray_factorial_ATP.pdf"),
          figure_png = file.path(path_fig, "pcoa_bray_factorial_ATP.png")
        )
      ),
      betadisper = list(
        figure = p_bdisp_bray_atp,
        paths = c(
          figure_pdf = file.path(path_fig, "betadisper_bray_factorial_ATP.pdf"),
          figure_png = file.path(path_fig, "betadisper_bray_factorial_ATP.png")
        )
      )
    ),
    pcoa_history = list(
      figure = p_pcoa_bray_history,
      paths = c(
        figure_pdf = file.path(path_fig, "pcoa_bray_stress_history.pdf"),
        figure_png = file.path(path_fig, "pcoa_bray_stress_history.png")
      )
    ),
    pcoa_history_parasite = list(
      figure = p_pcoa_bray_hp,
      paths = c(
        figure_pdf = file.path(path_fig, "pcoa_bray_stress_history_parasite.pdf"),
        figure_png = file.path(path_fig, "pcoa_bray_stress_history_parasite.png")
      )
    ),
    betadisper_history = list(
      figure = p_bdisp_hist,
      paths = c(
        figure_pdf = file.path(path_fig, "betadisper_bray_stress_history.pdf"),
        figure_png = file.path(path_fig, "betadisper_bray_stress_history.png")
      )
    ),
    betadisper_history_parasite = list(
      figure = p_bdisp_hp,
      paths = c(
        figure_pdf = file.path(path_fig, "betadisper_bray_stress_history_parasite.pdf"),
        figure_png = file.path(path_fig, "betadisper_bray_stress_history_parasite.png")
      )
    ),
    pcoa_parasite_faceted_history = list(
      figure = p_pcoa_bray_parasite_faceted,
      paths = c(
        figure_pdf = file.path(path_fig, "pcoa_bray_parasite_faceted_by_history.pdf"),
        figure_png = file.path(path_fig, "pcoa_bray_parasite_faceted_by_history.png")
      )
    ),
    pcoa_parasite_faceted_history_regime_points = list(
      figure = p_pcoa_bray_regime_pts,
      paths = c(
        figure_pdf = file.path(path_fig, "pcoa_bray_parasite_faceted_by_history_regime_points.pdf"),
        figure_png = file.path(path_fig, "pcoa_bray_parasite_faceted_by_history_regime_points.png")
      )
    )
  ),
  beta_canberra = list(
    factorial_atp = list(
      pcoa = list(
        figure = p_pcoa_canberra_atp,
        paths = c(
          figure_pdf = file.path(path_fig, "pcoa_canberra_factorial_ATP.pdf"),
          figure_png = file.path(path_fig, "pcoa_canberra_factorial_ATP.png")
        )
      ),
      betadisper = list(
        figure = p_bdisp_canberra_atp,
        paths = c(
          figure_pdf = file.path(path_fig, "betadisper_canberra_factorial_ATP.pdf"),
          figure_png = file.path(path_fig, "betadisper_canberra_factorial_ATP.png")
        )
      )
    ),
    pcoa_history = list(
      figure = p_pcoa_canberra_history,
      paths = c(
        figure_pdf = file.path(path_fig, "pcoa_canberra_stress_history.pdf"),
        figure_png = file.path(path_fig, "pcoa_canberra_stress_history.png")
      )
    ),
    pcoa_history_parasite = list(
      figure = p_pcoa_canberra_hp,
      paths = c(
        figure_pdf = file.path(path_fig, "pcoa_canberra_stress_history_parasite.pdf"),
        figure_png = file.path(path_fig, "pcoa_canberra_stress_history_parasite.png")
      )
    ),
    betadisper_history = list(
      figure = p_bdisp_canberra_hist,
      paths = c(
        figure_pdf = file.path(path_fig, "betadisper_canberra_stress_history.pdf"),
        figure_png = file.path(path_fig, "betadisper_canberra_stress_history.png")
      )
    ),
    betadisper_history_parasite = list(
      figure = p_bdisp_canberra_hp,
      paths = c(
        figure_pdf = file.path(path_fig, "betadisper_canberra_stress_history_parasite.pdf"),
        figure_png = file.path(path_fig, "betadisper_canberra_stress_history_parasite.png")
      )
    ),
    pcoa_parasite_faceted_history = list(
      figure = p_pcoa_canberra_parasite_faceted,
      paths = c(
        figure_pdf = file.path(path_fig, "pcoa_canberra_parasite_faceted_by_history.pdf"),
        figure_png = file.path(path_fig, "pcoa_canberra_parasite_faceted_by_history.png")
      )
    ),
    pcoa_parasite_faceted_history_regime_points = list(
      figure = p_pcoa_canberra_regime_pts,
      paths = c(
        figure_pdf = file.path(path_fig, "pcoa_canberra_parasite_faceted_by_history_regime_points.pdf"),
        figure_png = file.path(path_fig, "pcoa_canberra_parasite_faceted_by_history_regime_points.png")
      )
    )
  )
)

composition_bundle <- list(
  meta = list(
    run_date = as.character(Sys.Date()),
    bundle_version = "1.4",
    script = "Code/01__Analysis/02__Composition.R",
    ps_list_element = ps_list_element,
    description = paste(
      "Composition bundle: genus relative abundance by eight exposure regimes (merged samples);",
      "PERMANOVA (factorial A×T×P; history; history × parasite) and betadisper ANOVA tables on Bray-Curtis and Canberra;",
      "PERMANOVA (history + history × parasite) and betadisper ANOVA tables on Bray-Curtis and Canberra;",
      "PCoA and betadisper boxplots (stress history; history × parasite);",
      "stratified PERMANOVA/betadisper (Parasite within each HistoryLevelNum) and faceted PCoA (shared axes);",
      "optional faceted PCoA with points by eight exposure regimes and ellipses by parasite exposure."
    ),
    model_table_factorial_atp = model_subtitle_permanova_factorial_atp,
    model_table_c = model_subtitle_permanova_history,
    model_table_d = model_subtitle_permanova_history_parasite,
    model_table_dispersion_factorial_atp = model_subtitle_betadisper_factorial_atp,
    model_table_dispersion_history = model_subtitle_betadisper_history,
    model_table_dispersion_parasite = model_subtitle_betadisper_parasite,
    model_table_parasite_within_permanova = model_subtitle_permanova_parasite_within_history,
    model_table_parasite_within_betadisper = model_subtitle_betadisper_parasite_within_history,
    model_formulas = c(
      "PERMANOVA (factorial A×T×P)" = "distance ~ Antibiotics * Temperature * Parasite",
      "PERMANOVA (history)" = "distance ~ HistoryLevelNum",
      "PERMANOVA (history × parasite)" = "distance ~ HistoryLevelNum * Parasite",
      "PERMANOVA (stratified)" = "distance ~ Parasite within each HistoryLevelNum"
    ),
    n_perm = n_perm,
    n_processes = n_proc
  ),
  table_combined_trends = NULL,
  table_combined_interaction = NULL,
  table_combined_table2 = NULL,
  table_permanova_factorial_atp = gt_perma_atp,
  table_betadisper_anova_factorial_atp = gt_bdisp_atp,
  table_permanova_stress_history = gt_perma_hist,
  table_permanova_stress_history_parasite = gt_perma_hp,
  table_betadisper_anova_stress_history = gt_bdisp_hist,
  table_betadisper_anova_stress_history_parasite = gt_bdisp_hp,
  table_permanova_parasite_within_history = gt_perm_parasite_within,
  table_betadisper_parasite_within_history = gt_bdisp_parasite_within,
  tables_tidy = list(
    permanova_factorial_atp = tbl_perma_atp,
    betadisper_anova_factorial_atp = tbl_bdisp_atp,
    permanova_stress_history = tbl_perma_hist,
    permanova_stress_history_parasite = tbl_perma_hp,
    betadisper_anova_stress_history = tbl_bdisp_hist,
    betadisper_anova_stress_history_parasite = tbl_bdisp_hp,
    permanova_parasite_within_history = tbl_perm_parasite_within,
    betadisper_parasite_within_history = tbl_bdisp_parasite_within
  ),
  modules = modules
)

saveRDS(composition_bundle, bundle_rds)
message("Saved composition bundle: ", bundle_rds)

sieler2026_sync_main_figures_from_manifest(driver_script = "02__Composition.R")
