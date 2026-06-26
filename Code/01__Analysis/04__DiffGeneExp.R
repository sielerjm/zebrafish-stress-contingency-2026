# 04__DiffGeneExp.R
# Created by: Michael Sieler
# Date last updated: 2026-04-24
#
# Description: Host intestinal differential gene expression (DESeq2) at the final timepoint.
#   Runs:
#   1) Exposure-regime contrasts vs control (A- T- P-)
#   2) Prior stressor history pairwise + linear-trend models
#   3) Parasite exposure effects within each history level + interaction model
#   Exports manuscript-priority figures (top-50 DEG heatmap across parasite-by-history and
#   significant DEG bar chart across exposure regimes), summary tables, and a serialized bundle.
#
# Expected input:
#   - Data/DEG/salmon.merged.gene_counts_length_scaled.tsv (or corrected variant)
#   - Data/Metadata/metadata.tsv
#
# Expected output:
#   - Results/04__DiffGeneExp/Figures/*
#   - Results/04__DiffGeneExp/Tables/*
#   - Results/04__DiffGeneExp/Stats/diffgeneexp__host__bundle.rds + core DESeq2 objects

init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
if (!file.exists(init_rel)) {
  stop(
    "Set working directory to the Sieler2026 repository root, or run:\n",
    "  Rscript Code/01__Analysis/04__DiffGeneExp.R\n",
    "from the project root so `", init_rel, "` resolves."
  )
}
source(init_rel)

path_res <- file.path(path.results, "04__DiffGeneExp")
path_fig <- file.path(path_res, "Figures")
path_tbl <- file.path(path_res, "Tables")
path_stats <- file.path(path_res, "Stats")
dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(path_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(path_stats, recursive = TRUE, showWarnings = FALSE)

args_04 <- commandArgs(trailingOnly = TRUE)
figures_only_04 <- isTRUE(any(args_04 == "--figures-only"))
if (isTRUE(figures_only_04)) {
  rds_fig31 <- file.path(path_stats, "fig31_significant_genes_by_treatment_bar.rds")
  if (!file.exists(rds_fig31)) {
    stop("04__DiffGeneExp.R --figures-only requires:\n  ", rds_fig31, "\nRun a full 04 driver first.")
  }
  p31 <- readRDS(rds_fig31)
  fig_treatment_deg_pdf <- file.path(path_fig, "significant_genes_by_treatment_bar.pdf")
  fig_treatment_deg_png <- file.path(path_fig, "significant_genes_by_treatment_bar.png")
  ggplot2::ggsave(fig_treatment_deg_pdf, p31, width = 14, height = 5.75, device = "pdf")
  ggplot2::ggsave(fig_treatment_deg_png, p31, width = 14, height = 5.75, dpi = 300)
  sieler2026_sync_main_figures_from_manifest(driver_script = "04__DiffGeneExp.R", panel_ids = "3.1")
  message("04__DiffGeneExp.R --figures-only complete.")
  quit(save = "no", status = 0L)
}

# Archive previous outputs (copy → verify → clear) so each run is clean but reproducible.
sieler2026_archive_module_outputs(
  path_res_module = path_res,
  module_name = "04__DiffGeneExp",
  subdirs = c("Figures", "Tables", "Stats"),
  clear_original = TRUE
)

bundle_rds <- file.path(path_stats, "diffgeneexp__host__bundle.rds")

pick_counts_path <- function(path_deg_dir) {
  preferred <- file.path(path_deg_dir, "salmon.merged.gene_counts_length_scale__Corrected_f136-f138.tsv")
  fallback <- file.path(path_deg_dir, "salmon.merged.gene_counts_length_scaled.tsv")
  if (file.exists(preferred)) {
    return(preferred)
  }
  if (file.exists(fallback)) {
    return(fallback)
  }
  stop(
    "Could not find DEG count TSV in ", path_deg_dir,
    ". Expected one of:\n - ", preferred, "\n - ", fallback
  )
}

counts_path <- pick_counts_path(path.deg)
metadata_path <- file.path(path.data, "Metadata", "metadata.tsv")
if (!file.exists(metadata_path)) {
  stop("Metadata not found: ", metadata_path)
}

set.seed(42)

counts_tbl <- readr::read_tsv(counts_path, show_col_types = FALSE)
count_cols <- names(counts_tbl)[stringr::str_detect(names(counts_tbl), "^TS047_RoL_RNA_")]
if (length(count_cols) == 0L) {
  stop("No transcriptomics sample columns found in counts TSV: ", counts_path)
}

gene_ids <- counts_tbl$gene_id
gene_names <- counts_tbl$gene_name
counts_mat <- counts_tbl %>%
  dplyr::select(dplyr::all_of(count_cols)) %>%
  as.matrix()
mode(counts_mat) <- "numeric"
rownames(counts_mat) <- gene_ids

sample_id_num <- stringr::str_remove(colnames(counts_mat), "^TS047_RoL_RNA_")

metadata_raw <- readr::read_tsv(metadata_path, show_col_types = FALSE) %>%
  dplyr::rename(Time = Timepoint, Parasite = Pathogen) %>%
  dplyr::mutate(
    gut.sample.number = as.character(.data$gut.sample.number),
    gut_id_num = stringr::str_remove(.data$gut.sample.number, "^g"),
    Sample_RNA = paste0("TS047_RoL_RNA_", .data$gut_id_num)
  )

# Final-timepoint host transcriptomics are represented by RNA sample IDs present in counts.
metadata_final <- metadata_raw %>%
  dplyr::filter(.data$Time == 60, .data$gut_id_num %in% sample_id_num) %>%
  dplyr::mutate(
    Treatment = dplyr::case_when(
      .data$Antibiotics == 0 & .data$Temperature == 0 & .data$Parasite == 0 ~ "A- T- P-",
      .data$Antibiotics == 0 & .data$Temperature == 0 & .data$Parasite == 1 ~ "A- T- P+",
      .data$Antibiotics == 1 & .data$Temperature == 0 & .data$Parasite == 0 ~ "A+ T- P-",
      .data$Antibiotics == 1 & .data$Temperature == 0 & .data$Parasite == 1 ~ "A+ T- P+",
      .data$Antibiotics == 0 & .data$Temperature == 1 & .data$Parasite == 0 ~ "A- T+ P-",
      .data$Antibiotics == 0 & .data$Temperature == 1 & .data$Parasite == 1 ~ "A- T+ P+",
      .data$Antibiotics == 1 & .data$Temperature == 1 & .data$Parasite == 0 ~ "A+ T+ P-",
      .data$Antibiotics == 1 & .data$Temperature == 1 & .data$Parasite == 1 ~ "A+ T+ P+",
      TRUE ~ "Unknown"
    ),
    Treatment = factor(.data$Treatment, levels = treatment_order),
    Treatment_DESeq = dplyr::case_when(
      .data$Treatment == "A- T- P-" ~ "Aneg_Tneg_Pneg",
      .data$Treatment == "A- T- P+" ~ "Aneg_Tneg_Ppos",
      .data$Treatment == "A+ T- P-" ~ "Apos_Tneg_Pneg",
      .data$Treatment == "A+ T- P+" ~ "Apos_Tneg_Ppos",
      .data$Treatment == "A- T+ P-" ~ "Aneg_Tpos_Pneg",
      .data$Treatment == "A- T+ P+" ~ "Aneg_Tpos_Ppos",
      .data$Treatment == "A+ T+ P-" ~ "Apos_Tpos_Pneg",
      .data$Treatment == "A+ T+ P+" ~ "Apos_Tpos_Ppos",
      TRUE ~ "Unknown"
    ),
    Treatment_DESeq = factor(
      .data$Treatment_DESeq,
      levels = c(
        "Aneg_Tneg_Pneg", "Aneg_Tneg_Ppos", "Apos_Tneg_Pneg", "Apos_Tneg_Ppos",
        "Aneg_Tpos_Pneg", "Aneg_Tpos_Ppos", "Apos_Tpos_Pneg", "Apos_Tpos_Ppos"
      )
    ),
    HistoryLevelNum = dplyr::case_when(
      .data$Treatment %in% c("A- T- P-", "A- T- P+") ~ 0L,
      .data$Treatment %in% c("A+ T- P-", "A+ T- P+", "A- T+ P-", "A- T+ P+") ~ 1L,
      .data$Treatment %in% c("A+ T+ P-", "A+ T+ P+") ~ 2L,
      TRUE ~ NA_integer_
    ),
    HistoryLevel = factor(.data$HistoryLevelNum, levels = c(0, 1, 2)),
    Parasite_Exposed = dplyr::if_else(.data$Parasite == 1, "Exposed", "Unexposed"),
    Parasite_Exposed = factor(.data$Parasite_Exposed, levels = c("Unexposed", "Exposed"))
  ) %>%
  dplyr::filter(.data$Treatment != "Unknown") %>%
  dplyr::distinct(.data$Sample_RNA, .keep_all = TRUE) %>%
  as.data.frame()

if (nrow(metadata_final) == 0L) {
  stop("No final-timepoint samples matched between metadata and DEG counts.")
}

common_samples <- intersect(colnames(counts_mat), metadata_final$Sample_RNA)
counts_final <- counts_mat[, common_samples, drop = FALSE]
metadata_final <- metadata_final %>%
  dplyr::mutate(Sample_RNA = factor(.data$Sample_RNA, levels = colnames(counts_final))) %>%
  dplyr::arrange(.data$Sample_RNA) %>%
  dplyr::mutate(Sample_RNA = as.character(.data$Sample_RNA))

if (!all(colnames(counts_final) == metadata_final$Sample_RNA)) {
  stop("Sample alignment failed between counts and metadata.")
}

counts_int <- round(counts_final)
rownames(metadata_final) <- metadata_final$Sample_RNA

# Pre-filter genes: >=2 CPM in >=3 samples.
cpm <- sweep(counts_int, 2, colSums(counts_int), "/") * 1e6
keep <- rowSums(cpm >= 2) >= 3
counts_int_filt <- counts_int[keep, , drop = FALSE]
gene_map <- tibble::tibble(gene_id = gene_ids, gene_name = gene_names)

treatment_map <- metadata_final %>%
  dplyr::select("Treatment_DESeq", "Treatment") %>%
  dplyr::distinct()

# 1) Exposure regime model
dds_treatment <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_int_filt,
  colData = metadata_final,
  design = ~ Treatment_DESeq
)
dds_treatment$Treatment_DESeq <- stats::relevel(dds_treatment$Treatment_DESeq, ref = "Aneg_Tneg_Pneg")
set.seed(42)
dds_treatment <- DESeq2::DESeq(dds_treatment)

treatment_levels <- levels(dds_treatment$Treatment_DESeq)
treatment_contrasts <- treatment_levels[treatment_levels != "Aneg_Tneg_Pneg"]
treatment_results <- purrr::map(
  treatment_contrasts,
  function(lvl) {
    DESeq2::results(dds_treatment, contrast = c("Treatment_DESeq", lvl, "Aneg_Tneg_Pneg"), alpha = 0.05)
  }
) %>%
  purrr::set_names(treatment_contrasts)

all_treatment_results <- purrr::imap_dfr(
  treatment_results,
  function(res, lvl) {
    as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gene_map, by = "gene_id") %>%
      dplyr::mutate(
        Treatment_DESeq = lvl,
        Treatment = treatment_map$Treatment[match(lvl, treatment_map$Treatment_DESeq)],
        Treatment = dplyr::coalesce(.data$Treatment, lvl),
        Comparison = paste0(.data$Treatment, " vs A- T- P-")
      )
  }
)

treatment_deg_counts <- all_treatment_results %>%
  dplyr::filter(!is.na(.data$padj)) %>%
  dplyr::group_by(.data$Treatment) %>%
  dplyr::summarise(
    n_significant = sum(.data$padj < 0.05, na.rm = TRUE),
    n_up = sum(.data$padj < 0.05 & .data$log2FoldChange > 0, na.rm = TRUE),
    n_down = sum(.data$padj < 0.05 & .data$log2FoldChange < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Treatment = factor(.data$Treatment, levels = treatment_order)) %>%
  dplyr::arrange(.data$Treatment) %>%
  dplyr::mutate(
    total = .data$n_up + .data$n_down,
    prior_stressor_history = factor(
      dplyr::case_when(
        as.character(.data$Treatment) %in% c("A- T- P-", "A- T- P+") ~ "Zero",
        as.character(.data$Treatment) %in% c("A+ T- P-", "A+ T- P+", "A- T+ P-", "A- T+ P+") ~ "One",
        as.character(.data$Treatment) %in% c("A+ T+ P-", "A+ T+ P+") ~ "Two",
        TRUE ~ NA_character_
      ),
      levels = c("Zero", "One", "Two")
    )
  )

treatment_order_by_total <- treatment_deg_counts %>%
  dplyr::arrange(dplyr::desc(.data$total)) %>%
  dplyr::pull(.data$Treatment) %>%
  as.character()

treatment_deg_long <- treatment_deg_counts %>%
  dplyr::select("Treatment", "prior_stressor_history", "total", "n_up", "n_down") %>%
  tidyr::pivot_longer(cols = c("n_up", "n_down"), names_to = "direction", values_to = "n") %>%
  dplyr::mutate(
    direction = dplyr::recode(.data$direction, n_up = "Upregulated", n_down = "Downregulated")
  ) %>%
  dplyr::group_by(.data$prior_stressor_history) %>%
  dplyr::mutate(
    Treatment = forcats::fct_reorder(.data$Treatment, .data$total, .desc = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(`Prior Stressor History` = prior_stressor_history)

# Fixed y for count labels (~200) so they align across facets; white when bar extends past that height.
treatment_deg_labels <- treatment_deg_long %>%
  dplyr::mutate(
    label_text = as.character(.data$n),
    label_y = 200,
    label_color = dplyr::if_else(.data$n >= 200, "white", "black"),
    label_vjust = 0.5
  )

p_treatment_deg <- ggplot2::ggplot(
  treatment_deg_long,
  ggplot2::aes(x = .data$Treatment, y = .data$n, fill = .data$direction)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.9),
    color = "black",
    linewidth = SIELER2026_MIN_LINEWIDTH_MM
  ) +
  ggplot2::geom_text(
    data = treatment_deg_labels,
    ggplot2::aes(
      y = .data$label_y,
      label = .data$label_text,
      color = .data$label_color,
      vjust = .data$label_vjust
    ),
    position = ggplot2::position_dodge(width = 0.9),
    size = 4.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = c("Downregulated" = "steelblue", "Upregulated" = "firebrick"),
    breaks = c("Downregulated", "Upregulated")
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::expand_limits(y = 0) +
  # space = "free_x" makes panel width proportional to # of x categories so bar width matches across facets.
  ggplot2::facet_grid(
    cols = ggplot2::vars(`Prior Stressor History`),
    scales = "free",
    space = "free_x"
  ) +
  theme_sieler2026_publication(base_size = 14) +
  ggplot2::labs(
    title = "Number of Differentially Expressed Genes by Exposure Regime",
    subtitle = "Upregulated and downregulated genes (padj < 0.05), ordered by total significant genes within each prior stressor history",
    x = "Exposure Regimes",
    y = "Number of Genes",
    fill = "Association",
    caption = "All treatments compared to reference: A- T- P-"
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    strip.text = ggplot2::element_text(face = "bold")
  )

fig_treatment_deg_pdf <- file.path(path_fig, "significant_genes_by_treatment_bar.pdf")
fig_treatment_deg_png <- file.path(path_fig, "significant_genes_by_treatment_bar.png")
ggplot2::ggsave(fig_treatment_deg_pdf, p_treatment_deg, width = 14, height = 5.75, device = "pdf")
ggplot2::ggsave(fig_treatment_deg_png, p_treatment_deg, width = 14, height = 5.75, dpi = 300)
saveRDS(p_treatment_deg, file.path(path_stats, "fig31_significant_genes_by_treatment_bar.rds"))

# 2) Stress-history models
dds_history <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_int_filt,
  colData = metadata_final,
  design = ~ HistoryLevel
)
dds_history$HistoryLevel <- stats::relevel(dds_history$HistoryLevel, ref = "0")
set.seed(42)
dds_history <- DESeq2::DESeq(dds_history)

history_pairs <- list(
  "1_vs_0" = c("HistoryLevel", "1", "0"),
  "2_vs_0" = c("HistoryLevel", "2", "0"),
  "2_vs_1" = c("HistoryLevel", "2", "1")
)
history_results <- purrr::map(history_pairs, ~ DESeq2::results(dds_history, contrast = .x, alpha = 0.05))

dds_history_num <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_int_filt,
  colData = metadata_final %>% dplyr::mutate(HistoryLevelNum = as.numeric(.data$HistoryLevelNum)),
  design = ~ HistoryLevelNum
)
set.seed(42)
dds_history_num <- DESeq2::DESeq(dds_history_num)
history_trend_results <- DESeq2::results(dds_history_num, name = "HistoryLevelNum", alpha = 0.05)

history_deg_counts <- purrr::imap_dfr(
  history_results,
  function(res, comp) {
    df <- as.data.frame(res)
    tibble::tibble(
      comparison = comp,
      n_significant = sum(df$padj < 0.05, na.rm = TRUE),
      n_up = sum(df$padj < 0.05 & df$log2FoldChange > 0, na.rm = TRUE),
      n_down = sum(df$padj < 0.05 & df$log2FoldChange < 0, na.rm = TRUE)
    )
  }
)

# 3) Parasite-by-history models
history_levels <- c(0, 1, 2)
parasite_history_results <- purrr::map(
  history_levels,
  function(h) {
    meta_h <- metadata_final %>%
      dplyr::filter(.data$HistoryLevelNum == h) %>%
      as.data.frame()
    if (nrow(meta_h) == 0L || dplyr::n_distinct(meta_h$Parasite_Exposed) < 2L) {
      return(NULL)
    }
    counts_h <- counts_int_filt[, meta_h$Sample_RNA, drop = FALSE]
    rownames(meta_h) <- meta_h$Sample_RNA
    dds_h <- DESeq2::DESeqDataSetFromMatrix(countData = counts_h, colData = meta_h, design = ~ Parasite_Exposed)
    dds_h$Parasite_Exposed <- stats::relevel(dds_h$Parasite_Exposed, ref = "Unexposed")
    set.seed(42)
    dds_h <- DESeq2::DESeq(dds_h)
    DESeq2::results(dds_h, contrast = c("Parasite_Exposed", "Exposed", "Unexposed"), alpha = 0.05)
  }
) %>%
  purrr::set_names(paste0("History_", history_levels)) %>%
  purrr::discard(is.null)

dds_parasite_history <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_int_filt,
  colData = metadata_final,
  design = ~ Parasite_Exposed * HistoryLevelNum
)
dds_parasite_history$Parasite_Exposed <- stats::relevel(dds_parasite_history$Parasite_Exposed, ref = "Unexposed")
set.seed(42)
dds_parasite_history <- DESeq2::DESeq(dds_parasite_history)
interaction_name <- grep("Parasite_Exposed.*HistoryLevelNum", DESeq2::resultsNames(dds_parasite_history), value = TRUE)[1]
interaction_results <- DESeq2::results(dds_parasite_history, name = interaction_name, alpha = 0.05)

parasite_history_deg_counts <- purrr::imap_dfr(
  parasite_history_results,
  function(res, nm) {
    df <- as.data.frame(res)
    tibble::tibble(
      history_level = nm,
      n_significant = sum(df$padj < 0.05, na.rm = TRUE),
      n_up = sum(df$padj < 0.05 & df$log2FoldChange > 0, na.rm = TRUE),
      n_down = sum(df$padj < 0.05 & df$log2FoldChange < 0, na.rm = TRUE)
    )
  }
)

# Manuscript figure: top-50 DEG heatmap (parasite exposure across history levels)
top_genes_parasite <- purrr::map_dfr(
  parasite_history_results,
  function(res) {
    as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::filter(!is.na(.data$padj), .data$padj < 0.05) %>%
      dplyr::arrange(.data$padj) %>%
      dplyr::slice_head(n = 50) %>%
      dplyr::select("gene_id")
  }
) %>%
  dplyr::pull("gene_id") %>%
  unique() %>%
  head(50)

fig_heatmap_pdf <- file.path(path_fig, "heatmap_top_50_genes_parasite_exposure.pdf")
fig_heatmap_png <- file.path(path_fig, "heatmap_top_50_genes_parasite_exposure.png")

if (length(top_genes_parasite) > 0L) {
  vsd_parasite <- DESeq2::vst(dds_parasite_history, blind = FALSE)
  heatmap_mat <- SummarizedExperiment::assay(vsd_parasite)[top_genes_parasite, , drop = FALSE]
  annotation_col <- data.frame(
    `Stressor History` = factor(
      metadata_final$HistoryLevelNum,
      levels = c(0, 1, 2),
      labels = c("No prior stressors", "One prior stressor", "Two prior stressors")
    ),
    Parasite = factor(
      metadata_final$Parasite_Exposed,
      levels = c("Unexposed", "Exposed"),
      labels = c("Unexposed", "Exposed")
    ),
    row.names = colnames(heatmap_mat),
    check.names = FALSE
  )
  dark2_cols <- RColorBrewer::brewer.pal(3, "Dark2")
  stressor_history_colors <- c(
    "No prior stressors" = dark2_cols[1],
    "One prior stressor" = dark2_cols[2],
    "Two prior stressors" = dark2_cols[3]
  )
  ann_colors <- list(
    `Stressor History` = stressor_history_colors,
    Parasite = c("Unexposed" = "white", "Exposed" = "#E31A1C")
  )

  pheatmap::pheatmap(
    heatmap_mat,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = "Top 50 Differentially Expressed Genes: Parasite Exposure by Stressor History",
    fontsize = 10,
    annotation_legend = TRUE,
    border_color = "black",
    treeheight_row = 0,
    treeheight_col = 40,
    filename = fig_heatmap_pdf,
    width = 12,
    height = 9
  )
  pheatmap::pheatmap(
    heatmap_mat,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = "Top 50 Differentially Expressed Genes: Parasite Exposure by Stressor History",
    fontsize = 10,
    annotation_legend = TRUE,
    border_color = "black",
    treeheight_row = 0,
    treeheight_col = 40,
    filename = fig_heatmap_png,
    width = 12,
    height = 9
  )
}

# Export tables
readr::write_csv(treatment_deg_counts, file.path(path_tbl, "significant_genes_by_treatment.csv"))
readr::write_csv(history_deg_counts, file.path(path_tbl, "significant_genes_by_history.csv"))
readr::write_csv(parasite_history_deg_counts, file.path(path_tbl, "significant_genes_by_parasite_history.csv"))

# Also export GT HTML versions for the Results Rmd (avoid `knitr::kable()`).
gt_deg_counts_table <- function(df, title, subtitle) {
  gt_tbl <- df %>%
    gt::gt() %>%
    gt::tab_header(title = title, subtitle = subtitle)
  gt_tbl <- style_gt_significance(gt_tbl, df, alpha = 0.05)
  gt_tbl
}

gt_treatment <- gt_deg_counts_table(
  df = treatment_deg_counts,
  title = "Significant DEGs by exposure regime",
  subtitle = "DESeq2 contrasts vs A- T- P-; padj < 0.05"
)
gt_history <- gt_deg_counts_table(
  df = history_deg_counts,
  title = "Significant DEGs by prior stressor history",
  subtitle = "DESeq2: history comparisons; padj < 0.05"
)
gt_parasite_history <- gt_deg_counts_table(
  df = parasite_history_deg_counts,
  title = "Significant DEGs by parasite exposure within history level",
  subtitle = "DESeq2: parasite within history; padj < 0.05"
)

gt::gtsave(gt_treatment, file.path(path_tbl, "significant_genes_by_treatment.html"))
gt::gtsave(gt_history, file.path(path_tbl, "significant_genes_by_history.html"))
gt::gtsave(gt_parasite_history, file.path(path_tbl, "significant_genes_by_parasite_history.html"))

# Save stats objects
saveRDS(metadata_final, file.path(path_stats, "metadata_final.rds"))
saveRDS(counts_int_filt, file.path(path_stats, "counts_int_filtered.rds"))
saveRDS(keep, file.path(path_stats, "gene_filter_keep.rds"))
saveRDS(dds_treatment, file.path(path_stats, "dds_treatment.rds"))
saveRDS(treatment_results, file.path(path_stats, "treatment_results.rds"))
saveRDS(all_treatment_results, file.path(path_stats, "all_treatment_results.rds"))
saveRDS(dds_history, file.path(path_stats, "dds_history.rds"))
saveRDS(dds_history_num, file.path(path_stats, "dds_history_num.rds"))
saveRDS(history_results, file.path(path_stats, "history_results.rds"))
saveRDS(history_trend_results, file.path(path_stats, "history_trend_results.rds"))
saveRDS(dds_parasite_history, file.path(path_stats, "dds_parasite_history.rds"))
saveRDS(parasite_history_results, file.path(path_stats, "parasite_history_results.rds"))
saveRDS(interaction_results, file.path(path_stats, "interaction_results.rds"))
saveRDS(top_genes_parasite, file.path(path_stats, "top_genes_parasite_heatmap.rds"))

bundle <- list(
  meta = list(
    run_date = as.character(Sys.Date()),
    script = "Code/01__Analysis/04__DiffGeneExp.R",
    counts_path = counts_path,
    metadata_path = metadata_path,
    n_samples = ncol(counts_int_filt),
    n_genes_after_filter = nrow(counts_int_filt),
    model_formulas = c(
      treatment = "~ Treatment_DESeq",
      history_factor = "~ HistoryLevel",
      history_numeric = "~ HistoryLevelNum",
      parasite_by_history = "~ Parasite_Exposed * HistoryLevelNum"
    )
  ),
  paths = list(
    figures = list(
      significant_genes_by_treatment_bar_pdf = fig_treatment_deg_pdf,
      significant_genes_by_treatment_bar_png = fig_treatment_deg_png,
      heatmap_top_50_genes_parasite_exposure_pdf = fig_heatmap_pdf,
      heatmap_top_50_genes_parasite_exposure_png = fig_heatmap_png
    ),
    tables = list(
      significant_genes_by_treatment = file.path(path_tbl, "significant_genes_by_treatment.csv"),
      significant_genes_by_history = file.path(path_tbl, "significant_genes_by_history.csv"),
      significant_genes_by_parasite_history = file.path(path_tbl, "significant_genes_by_parasite_history.csv"),
      significant_genes_by_treatment_html = file.path(path_tbl, "significant_genes_by_treatment.html"),
      significant_genes_by_history_html = file.path(path_tbl, "significant_genes_by_history.html"),
      significant_genes_by_parasite_history_html = file.path(path_tbl, "significant_genes_by_parasite_history.html")
    ),
    stats_dir = path_stats
  ),
  tables = list(
    treatment_deg_counts = treatment_deg_counts,
    history_deg_counts = history_deg_counts,
    parasite_history_deg_counts = parasite_history_deg_counts
  ),
  table_treatment_deg_counts = gt_treatment,
  table_history_deg_counts = gt_history,
  table_parasite_history_deg_counts = gt_parasite_history
)
saveRDS(bundle, bundle_rds)
message("Saved bundle: ", bundle_rds)

sieler2026_sync_main_figures_from_manifest(driver_script = "04__DiffGeneExp.R")
