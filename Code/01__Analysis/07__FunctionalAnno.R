# 07__FunctionalAnno.R
# Created by: Michael Sieler
# Date last updated: 2026-06-26 (focal dotplots: cap tick = true value; y-axis rel 0.85)
#
# Description: GO and KEGG enrichment (clusterProfiler) for host genes with FDR-significant
#   partial correlations to genera from module 06: (1) union of genes linked to the top 10
#   genera by association count and (2) genes linked specifically to Culicoidibacter.
#   Gene symbols are mapped to Entrez IDs via org.Dr.eg.db (Danio rerio). KEGG uses organism
#   dre and may require network access on first run.
#
# Expected input: Run from Sieler2026 root; Results/06__Taxon-DEG-Mort/Tables/
#   combined_sig_partial_correlations.csv from Code/01__Analysis/06__Taxon-DEG-Mort.R.
# Expected output: Results/07__FunctionalAnno/{Figures,Tables,Stats}/ and
#   Stats/functional_anno__bundle.rds (focal-four combined dotplots + term-by-taxon wide tables;
#   focal-four gene-set Venns live in module 06).

init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
if (!file.exists(init_rel)) {
  stop(
    "Set working directory to the Sieler2026 repository root, or run:\n",
    "  Rscript Code/01__Analysis/07__FunctionalAnno.R\n",
    "from the project root so `", init_rel, "` resolves."
  )
}
source(init_rel)

path_res <- file.path(path.results, "07__FunctionalAnno")
path_fig <- file.path(path_res, "Figures")
path_tbl <- file.path(path_res, "Tables")
path_stats <- file.path(path_res, "Stats")
dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(path_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(path_stats, recursive = TRUE, showWarnings = FALSE)

# Archive previous outputs (copy → verify → clear) so each run is clean but reproducible.
sieler2026_archive_module_outputs(
  path_res_module = path_res,
  module_name = "07__FunctionalAnno",
  subdirs = c("Figures", "Tables", "Stats"),
  clear_original = TRUE
)

bundle_rds <- file.path(path_stats, "functional_anno__bundle.rds")
path_csv_06 <- file.path(path.results, "06__Taxon-DEG-Mort", "Tables", "combined_sig_partial_correlations.csv")

top_n_taxa <- 10L
top_n_plot_terms <- 15L

# Focal taxa for reviewer-ready enrichment replication
focal_genera <- c("Shewanella", "Culicoidibacter", "Flavobacterium", "Cetobacterium")
# Display / legend order (matches mortality combined scatter palette)
focal_genera_display <- c("Culicoidibacter", "Shewanella", "Flavobacterium", "Cetobacterium")
taxon_outline_colors <- c(
  Culicoidibacter = "#E41A1C",
  Shewanella = "#377EB8",
  Flavobacterium = "#4DAF4A",
  Cetobacterium = "#984EA3"
)

focal_slug <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
  x <- gsub("^_+|_+$", "", x, perl = TRUE)
  tolower(x)
}

# --- Helpers: plot labels -------------------------------------------------------
# KEGG y-axis: pathway name only (strip dre#####: if present)
kegg_pathway_label <- function(description) {
  sub("^dre[0-9]+:\\s*", "", description, perl = TRUE)
}

# GO MF: collapse long oxidoreductase activity terms for readability
shorten_go_mf_description <- function(description) {
  dplyr::if_else(
    grepl("(?i)^oxidoreductase activity", description, perl = TRUE),
    "Oxidoreductase activity",
    description
  )
}

# --- Helpers: enrichment tables and dotplots ------------------------------------

empty_go_df <- function() {
  tibble::tibble(
    Description = character(),
    p.adjust = numeric(),
    Count = integer(),
    geneID = character(),
    log_padj = numeric()
  )
}

empty_kegg_df <- function() {
  tibble::tibble(
    ID = character(),
    Description = character(),
    p.adjust = numeric(),
    Count = integer(),
    geneID = character(),
    log_padj = numeric()
  )
}

go_result_to_df <- function(res) {
  if (is.null(res)) {
    return(empty_go_df())
  }
  df <- as.data.frame(res)
  if (nrow(df) == 0L) {
    return(empty_go_df())
  }
  df %>%
    dplyr::select("ID", "Description", "p.adjust", "Count", "geneID") %>%
    dplyr::mutate(log_padj = -log10(.data$p.adjust)) %>%
    dplyr::arrange(.data$p.adjust)
}

kegg_result_to_df <- function(res) {
  if (is.null(res)) {
    return(empty_kegg_df())
  }
  df <- as.data.frame(res)
  if (nrow(df) == 0L) {
    return(empty_kegg_df())
  }
  df %>%
    dplyr::select("ID", "Description", "p.adjust", "Count", "geneID") %>%
    dplyr::mutate(log_padj = -log10(.data$p.adjust)) %>%
    dplyr::arrange(.data$p.adjust)
}

run_go_kegg <- function(entrez_ids, block_label) {
  out <- list(
    go_bp = NULL,
    go_mf = NULL,
    go_cc = NULL,
    kegg = NULL,
    go_bp_df = empty_go_df(),
    go_mf_df = empty_go_df(),
    go_cc_df = empty_go_df(),
    kegg_df = empty_kegg_df()
  )
  if (length(entrez_ids) < 1L) {
    return(out)
  }

  # Biological Process
  cat("\n  [", block_label, "] GO BP ...\n", sep = "")
  set.seed(42)
  out$go_bp <- tryCatch(
    clusterProfiler::enrichGO(
      gene = entrez_ids,
      OrgDb = org.Dr.eg.db::org.Dr.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) {
      message("  enrichGO BP failed: ", conditionMessage(e))
      NULL
    }
  )
  out$go_bp_df <- go_result_to_df(out$go_bp)

  cat("  [", block_label, "] GO MF ...\n", sep = "")
  set.seed(42)
  out$go_mf <- tryCatch(
    clusterProfiler::enrichGO(
      gene = entrez_ids,
      OrgDb = org.Dr.eg.db::org.Dr.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) {
      message("  enrichGO MF failed: ", conditionMessage(e))
      NULL
    }
  )
  out$go_mf_df <- go_result_to_df(out$go_mf)

  cat("  [", block_label, "] GO CC ...\n", sep = "")
  set.seed(42)
  out$go_cc <- tryCatch(
    clusterProfiler::enrichGO(
      gene = entrez_ids,
      OrgDb = org.Dr.eg.db::org.Dr.eg.db,
      keyType = "ENTREZID",
      ont = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) {
      message("  enrichGO CC failed: ", conditionMessage(e))
      NULL
    }
  )
  out$go_cc_df <- go_result_to_df(out$go_cc)

  cat("  [", block_label, "] KEGG ...\n", sep = "")
  set.seed(42)
  out$kegg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = entrez_ids,
      organism = "dre",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    ),
    error = function(e) {
      message("  enrichKEGG failed (network may be required): ", conditionMessage(e))
      NULL
    }
  )
  out$kegg_df <- kegg_result_to_df(out$kegg)

  out
}

# Dotplot: top terms by significance (data.frame from enrich, already sorted by p.adjust)
save_enrichment_dotplot <- function(
    df,
    title,
    subtitle,
    path_base_no_ext,
    w_in = 12,
    h_in = 8,
    shorten_mf_oxidoreductase = FALSE
) {
  if (nrow(df) < 1L) {
    return(invisible(NULL))
  }
  n_show <- min(top_n_plot_terms, nrow(df))
  plot_data <- df %>%
    dplyr::mutate(Description_plot = .data$Description)
  if (isTRUE(shorten_mf_oxidoreductase)) {
    plot_data <- plot_data %>%
      dplyr::mutate(
        Description_short = shorten_go_mf_description(.data$Description),
        # If shortening collapses distinct GO IDs to the same label, append the GO ID
        # so each dot corresponds to a unique functional annotation row.
        Description_plot = dplyr::if_else(
          duplicated(.data$Description_short) | duplicated(.data$Description_short, fromLast = TRUE),
          paste0(.data$Description_short, " (", .data$ID, ")"),
          .data$Description_short
        )
      )
  }
  plot_data <- plot_data %>%
    dplyr::slice_head(n = n_show) %>%
    dplyr::mutate(
      Description_wrapped = stringr::str_wrap(.data$Description_plot, width = 50),
      Description_wrapped = forcats::fct_reorder(.data$Description_wrapped, .data$log_padj)
    )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$log_padj, y = .data$Description_wrapped)) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$Count, fill = .data$log_padj),
      shape = 21L,
      colour = "black",
      stroke = max(0.45, SIELER2026_MIN_LINEWIDTH_MM)
    ) +
    ggplot2::scale_fill_viridis_c(option = "viridis", name = "-log10(adj. p)") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "-log10(Adjusted P-value)",
      y = NULL,
      size = "Gene count"
    ) +
    theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.85)),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.65), hjust = 0),
      plot.margin = ggplot2::margin(l = 20, r = 10, t = 10, b = 10)
    )

  p_png <- paste0(path_base_no_ext, ".png")
  p_pdf <- paste0(path_base_no_ext, ".pdf")
  ggplot2::ggsave(p_png, p, width = w_in, height = h_in, units = "in", dpi = 300L)
  ggplot2::ggsave(p_pdf, p, width = w_in, height = h_in, units = "in", device = "pdf")
  invisible(p)
}

save_gg_list <- function(enrich_df, prefix, title_prefix, subtitle) {
  if (nrow(enrich_df$go_bp_df) > 0L) {
    save_enrichment_dotplot(
      enrich_df$go_bp_df,
      title = paste0(title_prefix, " - GO Biological Process"),
      subtitle = subtitle,
      path_base_no_ext = file.path(path_fig, paste0("go_bp_", prefix, "_dotplot"))
    )
  }
  if (nrow(enrich_df$go_mf_df) > 0L) {
    save_enrichment_dotplot(
      enrich_df$go_mf_df,
      title = paste0(title_prefix, " - GO Molecular Function"),
      subtitle = subtitle,
      path_base_no_ext = file.path(path_fig, paste0("go_mf_", prefix, "_dotplot")),
      shorten_mf_oxidoreductase = TRUE
    )
  }
  if (nrow(enrich_df$go_cc_df) > 0L) {
    save_enrichment_dotplot(
      enrich_df$go_cc_df,
      title = paste0(title_prefix, " - GO Cellular Component"),
      subtitle = subtitle,
      path_base_no_ext = file.path(path_fig, paste0("go_cc_", prefix, "_dotplot"))
    )
  }
  if (nrow(enrich_df$kegg_df) > 0L) {
    # KEGG uses ID + Description; build a label column for plotting
    n_k <- min(top_n_plot_terms, nrow(enrich_df$kegg_df))
    kdf <- enrich_df$kegg_df %>%
      dplyr::slice_head(n = n_k) %>%
      dplyr::mutate(
        lab = kegg_pathway_label(paste0(.data$ID, ": ", .data$Description)),
        Description_wrapped = stringr::str_wrap(.data$lab, width = 50),
        Description_wrapped = forcats::fct_reorder(.data$Description_wrapped, .data$log_padj)
      )
    p <- ggplot2::ggplot(kdf, ggplot2::aes(x = .data$log_padj, y = .data$Description_wrapped)) +
      ggplot2::geom_point(
        ggplot2::aes(size = .data$Count, fill = .data$log_padj),
        shape = 21L,
        colour = "black",
        stroke = max(0.45, SIELER2026_MIN_LINEWIDTH_MM)
      ) +
      ggplot2::scale_fill_viridis_c(option = "viridis", name = "-log10(adj. p)") +
      ggplot2::labs(
        title = paste0(title_prefix, " - KEGG pathways"),
        subtitle = subtitle,
        x = "-log10(Adjusted P-value)",
        y = NULL,
        size = "Gene count"
      ) +
      theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.85)),
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.65), hjust = 0),
        plot.margin = ggplot2::margin(l = 20, r = 10, t = 10, b = 10)
      )
    ggplot2::ggsave(
      file.path(path_fig, paste0("kegg_", prefix, "_dotplot.png")),
      p,
      width = 12,
      height = 8,
      units = "in",
      dpi = 300L
    )
    ggplot2::ggsave(
      file.path(path_fig, paste0("kegg_", prefix, "_dotplot.pdf")),
      p,
      width = 12,
      height = 8,
      units = "in",
      device = "pdf"
    )
  }
}

# --- Focal four: combined enrichment dotplots + gene-set Venns -----------------
# Stack per-taxon enrichment tables; outline colour = taxon (same palette as mortality scatter).

bind_focal_enrichment <- function(focal_enrich_lst, display_order, extract_df) {
  pieces <- list()
  for (g in display_order) {
    if (!g %in% names(focal_enrich_lst)) {
      next
    }
    d <- extract_df(focal_enrich_lst[[g]])
    if (is.null(d) || nrow(d) < 1L) {
      next
    }
    pieces[[length(pieces) + 1L]] <- d %>%
      dplyr::mutate(Taxon = factor(g, levels = display_order))
  }
  if (length(pieces) < 1L) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(pieces)
}

prepare_combined_dotplot_data <- function(
    combined_df,
    shorten_mf_oxidoreductase = FALSE,
    is_kegg = FALSE
) {
  if (nrow(combined_df) < 1L) {
    return(combined_df)
  }
  plot_data <- combined_df %>%
    dplyr::mutate(Description_plot = .data$Description)
  if (isTRUE(is_kegg)) {
    plot_data <- plot_data %>%
      dplyr::mutate(
        Description_plot = kegg_pathway_label(paste0(.data$ID, ": ", .data$Description))
      )
  }
  if (isTRUE(shorten_mf_oxidoreductase)) {
    plot_data <- plot_data %>%
      dplyr::mutate(
        Description_short = shorten_go_mf_description(.data$Description),
        Description_plot = dplyr::if_else(
          duplicated(.data$Description_short) | duplicated(.data$Description_short, fromLast = TRUE),
          paste0(.data$Description_short, " (", .data$ID, ")"),
          .data$Description_short
        )
      )
  }
  top_ids <- plot_data %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::summarise(mx = max(.data$log_padj, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$mx)) %>%
    dplyr::slice_head(n = top_n_plot_terms) %>%
    dplyr::pull(.data$ID)
  plot_data %>%
    dplyr::filter(.data$ID %in% top_ids) %>%
    dplyr::mutate(
      Description_wrapped = stringr::str_wrap(.data$Description_plot, width = 50L),
      Description_wrapped = forcats::fct_reorder(.data$Description_wrapped, .data$log_padj)
    )
}

save_focal_combined_dotplot <- function(
    plot_data,
    title,
    subtitle,
    stem_no_ext,
    w_in = 12,
    h_in = 9,
    taxon_factor_levels = NULL,
    x_display_max = NULL
) {
  if (nrow(plot_data) < 1L) {
    return(invisible(NULL))
  }
  lvl <- taxon_factor_levels
  if (is.null(lvl)) {
    lvl <- names(taxon_outline_colors)
  }
  lvl <- lvl[lvl %in% names(taxon_outline_colors)]
  if (length(lvl) < 1L) {
    stop("save_focal_combined_dotplot: no valid taxon names in taxon_factor_levels.")
  }

  plot_data <- plot_data %>%
    dplyr::mutate(Taxon = factor(as.character(.data$Taxon), levels = lvl))

  # Optional x-axis display cap: squish outliers; right tick shows true -log10(adj. p).
  if (!is.null(x_display_max)) {
    x_cap <- as.numeric(x_display_max)[1L]
    if (!is.finite(x_cap) || x_cap <= 0) {
      stop("save_focal_combined_dotplot: x_display_max must be a positive finite number.")
    }
    plot_data <- plot_data %>%
      dplyr::mutate(log_padj_plot = pmin(.data$log_padj, x_cap))
  } else {
    plot_data <- plot_data %>%
      dplyr::mutate(log_padj_plot = .data$log_padj)
  }

  # Gene counts are integers; size legend breaks must be whole numbers (not 2.5, 5, …).
  cnt <- plot_data$Count
  cnt <- cnt[is.finite(cnt)]
  cmin <- max(1L, as.integer(floor(min(cnt))))
  cmax <- as.integer(ceiling(max(cnt)))
  if (!length(cmax) || !is.finite(cmax)) {
    cmax <- cmin
  }
  if (cmax < cmin) {
    cmax <- cmin
  }
  if (cmax - cmin <= 6L) {
    size_breaks <- seq.int(cmin, cmax, by = 1L)
  } else {
    size_breaks <- unique(as.integer(pretty(c(cmin, cmax), n = 6L)))
    size_breaks <- size_breaks[size_breaks >= cmin & size_breaks <= cmax]
  }
  if (length(size_breaks) < 1L) {
    size_breaks <- cmin
  }

  # Black outline on shape-21 points; stroke matches Gene count legend.
  point_stroke <- max(0.45, SIELER2026_MIN_LINEWIDTH_MM)
  guide_linewidth <- max(0.25, SIELER2026_MIN_LINEWIDTH_MM * 0.75)
  fill_vals <- unname(taxon_outline_colors[lvl])
  names(fill_vals) <- lvl

  term_levels <- levels(plot_data$Description_wrapped)
  n_terms <- length(term_levels)
  plot_data <- plot_data %>%
    dplyr::mutate(y_num = as.numeric(.data$Description_wrapped))

  y_lim_lo <- if (!is.null(x_display_max)) 0.35 else 0.5
  y_lim_hi <- n_terms + 0.5
  axis_y <- y_lim_lo

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$log_padj_plot, y = .data$y_num))

  # Capped x-axis: vertical guide at squished position + discontinuity marker on axis.
  if (!is.null(x_display_max)) {
    capped_rows <- plot_data %>%
      dplyr::filter(.data$log_padj > x_cap)
    if (nrow(capped_rows) > 0L) {
      p <- p +
        ggplot2::geom_segment(
          data = capped_rows,
          ggplot2::aes(
            x = .data$log_padj_plot,
            xend = .data$log_padj_plot,
            y = .data$y_num,
            yend = axis_y
          ),
          inherit.aes = FALSE,
          linewidth = guide_linewidth,
          colour = "grey55",
          linetype = "dashed"
        )
    }
    p <- p +
      ggplot2::annotate(
        "text",
        x = (x_cap - 1) + 0.5,
        y = axis_y,
        label = "...",
        size = 3.2,
        colour = "grey40",
        vjust = -0.6
      )
  }

  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$Count,
        fill = .data$Taxon
      ),
      shape = 21L,
      colour = "black",
      stroke = point_stroke,
      alpha = 0.92,
      position = ggplot2::position_jitter(width = 0, height = 0.14, seed = 42L)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(n_terms),
      labels = term_levels,
      limits = c(y_lim_lo, y_lim_hi),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::scale_fill_manual(
      name = "Focal taxa",
      values = fill_vals,
      limits = lvl,
      breaks = lvl,
      drop = FALSE
    ) +
    ggplot2::scale_size_continuous(
      range = c(2.2, 9),
      breaks = size_breaks,
      labels = as.character(size_breaks)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1L,
        ncol = 2L,
        byrow = TRUE,
        override.aes = list(
          size = 5.5,
          colour = "black",
          stroke = point_stroke,
          alpha = 1
        )
      ),
      size = ggplot2::guide_legend(
        order = 2L,
        override.aes = list(
          fill = "white",
          colour = "black",
          stroke = point_stroke,
          alpha = 1
        )
      )
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "-log10(Adjusted P-value)",
      y = NULL,
      size = "Gene count"
    ) +
    theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.85)),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.85), hjust = 0),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.box = "horizontal",
      legend.spacing.x = grid::unit(2, "pt"),
      legend.key.spacing.x = grid::unit(0.5, "pt"),
      legend.box.margin = ggplot2::margin(t = 2, r = 0, b = 2, l = 0),
      plot.margin = ggplot2::margin(
        l = 24,
        r = 18,
        t = 10,
        b = if (!is.null(x_display_max)) 14 else 10
      )
    ) +
    ggplot2::coord_cartesian(clip = if (!is.null(x_display_max)) "off" else "on")

  if (!is.null(x_display_max)) {
    cap_tick_label <- plot_data %>%
      dplyr::filter(.data$log_padj > x_cap) %>%
      dplyr::summarise(
        lbl = sprintf("%.1f", max(.data$log_padj, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::pull(.data$lbl)
    if (!length(cap_tick_label) || !nzchar(cap_tick_label)) {
      cap_tick_label <- as.character(x_cap)
    }
    x_breaks <- c(0, 1, 2, 3, 4, x_cap)
    x_labels <- c("0", "1", "2", "3", "4", cap_tick_label)
    p <- p +
      ggplot2::scale_x_continuous(
        limits = c(0, x_cap),
        breaks = x_breaks,
        labels = x_labels,
        oob = scales::squish,
        expand = ggplot2::expansion(mult = c(0.02, 0.08))
      )
  }

  ggplot2::ggsave(paste0(stem_no_ext, ".png"), p, width = w_in, height = h_in, units = "in", dpi = 300L)
  ggplot2::ggsave(paste0(stem_no_ext, ".pdf"), p, width = w_in, height = h_in, units = "in", device = "pdf")
  invisible(p)
}

# --- Focal four: term-by-taxon wide tables (GO/KEGG enrichment outputs) ---------

build_focal_term_by_taxon_wide <- function(plot_data, display_order) {
  if (is.null(plot_data) || nrow(plot_data) < 1L) {
    return(tibble::tibble())
  }
  if (!"Description_plot" %in% names(plot_data)) {
    stop("build_focal_term_by_taxon_wide: expected column Description_plot.")
  }
  desc_tbl <- plot_data %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::slice_max(order_by = .data$log_padj, n = 1L, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select("ID", Description = "Description_plot")

  cnt_wide <- plot_data %>%
    dplyr::mutate(Taxon = as.character(.data$Taxon)) %>%
    dplyr::select("ID", "Taxon", "Count") %>%
    tidyr::pivot_wider(
      id_cols = "ID",
      names_from = "Taxon",
      values_from = "Count",
      values_fill = NA_integer_,
      names_prefix = "Count_"
    )

  lp_wide <- plot_data %>%
    dplyr::mutate(Taxon = as.character(.data$Taxon)) %>%
    dplyr::select("ID", "Taxon", "log_padj") %>%
    tidyr::pivot_wider(
      id_cols = "ID",
      names_from = "Taxon",
      values_from = "log_padj",
      values_fill = NA_real_,
      names_prefix = "log_padj_"
    )

  out <- desc_tbl %>%
    dplyr::left_join(cnt_wide, by = "ID") %>%
    dplyr::left_join(lp_wide, by = "ID")

  for (t in display_order) {
    cc <- paste0("Count_", t)
    lp <- paste0("log_padj_", t)
    if (!cc %in% names(out)) {
      out[[cc]] <- NA_integer_
    }
    if (!lp %in% names(out)) {
      out[[lp]] <- NA_real_
    }
  }

  pref_counts <- paste0("Count_", display_order)
  pref_lp <- paste0("log_padj_", display_order)
  out %>%
    dplyr::select(
      "ID",
      "Description",
      dplyr::all_of(pref_counts),
      dplyr::all_of(pref_lp)
    )
}

gt_focal_term_by_taxon_wide <- function(wide_tbl, title, subtitle) {
  if (nrow(wide_tbl) < 1L) {
    return(NULL)
  }
  count_cols <- grep("^Count_", names(wide_tbl), value = TRUE)
  lp_cols <- grep("^log_padj_", names(wide_tbl), value = TRUE)
  g <- gt::gt(wide_tbl) %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::tab_options(table.font.size = gt::px(11))
  if (length(count_cols) > 0L) {
    g <- g %>%
      gt::fmt_integer(columns = dplyr::all_of(count_cols)) %>%
      gt::tab_spanner(label = "Gene count (enriched term)", columns = dplyr::all_of(count_cols))
  }
  if (length(lp_cols) > 0L) {
    g <- g %>%
      gt::fmt_number(columns = dplyr::all_of(lp_cols), decimals = 2) %>%
      gt::tab_spanner(label = "-log10(adj. p)", columns = dplyr::all_of(lp_cols))
  }
  g
}

run_focal_four_combined_outputs <- function(
    focal_enrich_lst,
    display_order
) {
  out <- list(go_bp = NULL, go_mf = NULL, go_cc = NULL, kegg = NULL)
  subt_dot <- paste0(
    "Focal genera (combined): genes with significant partial correlations (BH FDR < 0.1); ",
    "dot fill = focal taxon; x-axis = -log10(adj. p)"
  )

  # --- Combined dotplots (GO BP / MF / CC + KEGG) ---
  comb_bp <- bind_focal_enrichment(focal_enrich_lst, display_order, function(e) e$go_bp_df)
  pd_bp <- prepare_combined_dotplot_data(comb_bp, shorten_mf_oxidoreductase = FALSE, is_kegg = FALSE)
  if (nrow(pd_bp) > 0L) {
    readr::write_csv(pd_bp, file.path(path_tbl, "go_bp_focal_four_combined_source.csv"))
    wide_bp <- build_focal_term_by_taxon_wide(pd_bp, display_order)
    readr::write_csv(wide_bp, file.path(path_tbl, "go_bp_focal_four_term_by_taxon_wide.csv"))
    gt_bp_w <- gt_focal_term_by_taxon_wide(
      wide_bp,
      title = "GO Biological Process (focal four): terms by taxon",
      subtitle = "Rows = terms shown in combined dotplot; Count / -log10(adj. p) per taxon (NA = term not enriched for that taxon)."
    )
    if (!is.null(gt_bp_w)) {
      gt::gtsave(gt_bp_w, file.path(path_tbl, "go_bp_focal_four_term_by_taxon_wide.html"))
    }
    out$go_bp <- save_focal_combined_dotplot(
      pd_bp,
      title = "Functional enrichment - GO Biological Process",
      subtitle = subt_dot,
      stem_no_ext = file.path(path_fig, "go_bp_focal_four_dotplot"),
      taxon_factor_levels = display_order
    )
  }

  comb_mf <- bind_focal_enrichment(focal_enrich_lst, display_order, function(e) e$go_mf_df)
  pd_mf <- prepare_combined_dotplot_data(comb_mf, shorten_mf_oxidoreductase = TRUE, is_kegg = FALSE)
  if (nrow(pd_mf) > 0L) {
    readr::write_csv(pd_mf, file.path(path_tbl, "go_mf_focal_four_combined_source.csv"))
    wide_mf <- build_focal_term_by_taxon_wide(pd_mf, display_order)
    readr::write_csv(wide_mf, file.path(path_tbl, "go_mf_focal_four_term_by_taxon_wide.csv"))
    gt_mf_w <- gt_focal_term_by_taxon_wide(
      wide_mf,
      title = "GO Molecular Function (focal four): terms by taxon",
      subtitle = "Rows = terms shown in combined dotplot; MF oxidoreductase labels shortened as in the figure."
    )
    if (!is.null(gt_mf_w)) {
      gt::gtsave(gt_mf_w, file.path(path_tbl, "go_mf_focal_four_term_by_taxon_wide.html"))
    }
    out$go_mf <- save_focal_combined_dotplot(
      pd_mf,
      title = "Functional enrichment - GO Molecular Function",
      subtitle = subt_dot,
      stem_no_ext = file.path(path_fig, "go_mf_focal_four_dotplot"),
      taxon_factor_levels = display_order
    )
  }

  comb_cc <- bind_focal_enrichment(focal_enrich_lst, display_order, function(e) e$go_cc_df)
  pd_cc <- prepare_combined_dotplot_data(comb_cc, shorten_mf_oxidoreductase = FALSE, is_kegg = FALSE)
  if (nrow(pd_cc) > 0L) {
    readr::write_csv(pd_cc, file.path(path_tbl, "go_cc_focal_four_combined_source.csv"))
    wide_cc <- build_focal_term_by_taxon_wide(pd_cc, display_order)
    readr::write_csv(wide_cc, file.path(path_tbl, "go_cc_focal_four_term_by_taxon_wide.csv"))
    gt_cc_w <- gt_focal_term_by_taxon_wide(
      wide_cc,
      title = "GO Cellular Component (focal four): terms by taxon",
      subtitle = "Rows = terms shown in combined dotplot; Count / -log10(adj. p) per taxon (NA = not enriched)."
    )
    if (!is.null(gt_cc_w)) {
      gt::gtsave(gt_cc_w, file.path(path_tbl, "go_cc_focal_four_term_by_taxon_wide.html"))
    }
    out$go_cc <- save_focal_combined_dotplot(
      pd_cc,
      title = "Functional enrichment - GO Cellular Component",
      subtitle = subt_dot,
      stem_no_ext = file.path(path_fig, "go_cc_focal_four_dotplot"),
      taxon_factor_levels = display_order
    )
  }

  comb_k <- bind_focal_enrichment(focal_enrich_lst, display_order, function(e) e$kegg_df)
  pd_k <- prepare_combined_dotplot_data(comb_k, shorten_mf_oxidoreductase = FALSE, is_kegg = TRUE)
  if (nrow(pd_k) > 0L) {
    readr::write_csv(pd_k, file.path(path_tbl, "kegg_focal_four_combined_source.csv"))
    wide_k <- build_focal_term_by_taxon_wide(pd_k, display_order)
    readr::write_csv(wide_k, file.path(path_tbl, "kegg_focal_four_term_by_taxon_wide.csv"))
    gt_k_w <- gt_focal_term_by_taxon_wide(
      wide_k,
      title = "KEGG pathways (focal four): pathways by taxon",
      subtitle = "Rows = pathways shown in combined dotplot; pathway labels match the figure (dre prefix stripped in Description)."
    )
    if (!is.null(gt_k_w)) {
      gt::gtsave(gt_k_w, file.path(path_tbl, "kegg_focal_four_term_by_taxon_wide.html"))
    }
    out$kegg <- save_focal_combined_dotplot(
      pd_k,
      title = "Functional enrichment - KEGG pathways",
      subtitle = subt_dot,
      stem_no_ext = file.path(path_fig, "kegg_focal_four_dotplot"),
      taxon_factor_levels = display_order,
      x_display_max = 5
    )
  }

  invisible(out)
}

# --- Load partial-correlation edges ---------------------------------------------
if (!file.exists(path_csv_06)) {
  stop("Required file missing: ", path_csv_06, "\nRun Code/01__Analysis/06__Taxon-DEG-Mort.R first.")
}

correlation_data <- readr::read_csv(path_csv_06, show_col_types = FALSE)

if (!"gene_id" %in% names(correlation_data)) {
  if ("gene" %in% names(correlation_data)) {
    correlation_data <- correlation_data %>% dplyr::mutate(gene_id = .data$gene)
  } else {
    stop("Cannot find gene_id (or gene) column in combined_sig_partial_correlations.csv")
  }
}

if (!"TaxaID" %in% names(correlation_data)) {
  stop("Cannot find TaxaID column in combined_sig_partial_correlations.csv")
}

if (!"correlation" %in% names(correlation_data)) {
  stop("Cannot find correlation column in combined_sig_partial_correlations.csv")
}

set.seed(42)

if (nrow(correlation_data) < 1L) {
  message("No significant partial-correlation edges in input; writing minimal bundle only.")
  bundle <- list(
    meta = list(
      run_date = as.character(Sys.Date()),
      script = "Code/01__Analysis/07__FunctionalAnno.R",
      skipped = TRUE,
      reason = "zero rows in combined_sig_partial_correlations.csv",
      n_sig_edges = 0L,
      n_genes_top10_unique = 0L,
      n_genes_culicoidibacter_unique = 0L,
      mapping_rate_top10 = NA_real_,
      mapping_rate_culicoidibacter = NA_real_
    ),
    tables = list(),
    paths = list(
      figures = list(),
      tables = list(),
      stats_dir = path_stats
    ),
    enrich = list(top10 = NULL, culicoidibacter = NULL),
    ggplots_focal_four_combined = NULL
  )
  saveRDS(bundle, bundle_rds)
  message("Saved minimal bundle: ", bundle_rds)
  message("07__FunctionalAnno complete (skipped enrichment).")
} else {

# --- Top 10 genera (same ranking as module 06 network/scatter) ----------------
top_10_taxa_summary <- correlation_data %>%
  dplyr::group_by(.data$TaxaID) %>%
  dplyr::summarise(
    n_correlations = dplyr::n(),
    n_positive = sum(.data$correlation > 0, na.rm = TRUE),
    n_negative = sum(.data$correlation < 0, na.rm = TRUE),
    mean_correlation = mean(.data$correlation, na.rm = TRUE),
    mean_abs_correlation = mean(abs(.data$correlation), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$n_correlations), dplyr::desc(.data$mean_abs_correlation)) %>%
  dplyr::slice_head(n = top_n_taxa)

top_10_genus_names <- top_10_taxa_summary$TaxaID

readr::write_csv(top_10_taxa_summary, file.path(path_tbl, "top_10_taxa_summary_for_anno.csv"))

gene_counts_by_taxa <- correlation_data %>%
  dplyr::filter(.data$TaxaID %in% top_10_genus_names) %>%
  dplyr::group_by(.data$TaxaID) %>%
  dplyr::summarise(
    n_genes = length(unique(.data$gene_id)),
    n_positive = sum(.data$correlation > 0, na.rm = TRUE),
    n_negative = sum(.data$correlation < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$n_genes))

readr::write_csv(gene_counts_by_taxa, file.path(path_tbl, "gene_counts_by_taxa_top10.csv"))

genes_for_top_taxa <- correlation_data %>%
  dplyr::filter(.data$TaxaID %in% top_10_genus_names)

readr::write_csv(genes_for_top_taxa, file.path(path_tbl, "genes_for_top10_taxa_edges.csv"))

unique_genes_top10 <- unique(genes_for_top_taxa$gene_id)

# --- Focal genus edges (per-genus replication of prior Culicoidibacter outputs) -----------------
focal_edges <- list()
focal_genes <- list()
focal_by_comparison <- list()

for (g in focal_genera) {
  slug <- focal_slug(g)
  df_g <- correlation_data %>% dplyr::filter(.data$TaxaID == g)
  focal_edges[[g]] <- df_g
  focal_genes[[g]] <- unique(df_g$gene_id)

  # Keep prior naming scheme (culicoidibacter_...) by using slug in filename.
  readr::write_csv(df_g, file.path(path_tbl, paste0(slug, "_partial_correlation_edges.csv")))

  if ("comparison_name" %in% names(df_g) && nrow(df_g) > 0L) {
    by_comp <- df_g %>%
      dplyr::group_by(.data$comparison_name) %>%
      dplyr::summarise(
        n_genes = length(unique(.data$gene_id)),
        n_positive = sum(.data$correlation > 0, na.rm = TRUE),
        n_negative = sum(.data$correlation < 0, na.rm = TRUE),
        mean_correlation = mean(.data$correlation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(.data$n_genes))
    focal_by_comparison[[g]] <- by_comp
    readr::write_csv(by_comp, file.path(path_tbl, paste0(slug, "_by_comparison.csv")))
  } else {
    focal_by_comparison[[g]] <- NULL
  }
}

# --- SYMBOL → Entrez ----------------------------------------------------------
cat("\n=== Mapping gene symbols to Entrez (top 10 taxa gene set) ===\n")
entrez_mapping_top10 <- AnnotationDbi::mapIds(
  org.Dr.eg.db::org.Dr.eg.db,
  keys = unique_genes_top10,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
) %>%
  tibble::enframe(name = "SYMBOL", value = "ENTREZID") %>%
  dplyr::filter(!is.na(.data$ENTREZID))

readr::write_csv(entrez_mapping_top10, file.path(path_tbl, "gene_mapping_top10_taxa.csv"))

entrez_ids_top10 <- entrez_mapping_top10$ENTREZID
rate_top10 <- if (length(unique_genes_top10) > 0L) {
  nrow(entrez_mapping_top10) / length(unique_genes_top10)
} else {
  NA_real_
}
cat("  Mapped ", length(entrez_ids_top10), " / ", length(unique_genes_top10), " genes (", round(rate_top10 * 100, 1), "%)\n", sep = "")

focal_entrez_map <- list()
focal_entrez_ids <- list()
focal_mapping_rate <- list()

for (g in focal_genera) {
  slug <- focal_slug(g)
  g_genes <- focal_genes[[g]]

  cat("\n=== Mapping gene symbols to Entrez (", g, ") ===\n", sep = "")
  entrez_mapping <- AnnotationDbi::mapIds(
    org.Dr.eg.db::org.Dr.eg.db,
    keys = g_genes,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  ) %>%
    tibble::enframe(name = "SYMBOL", value = "ENTREZID") %>%
    dplyr::filter(!is.na(.data$ENTREZID))

  readr::write_csv(entrez_mapping, file.path(path_tbl, paste0("gene_mapping_", slug, ".csv")))

  entrez_ids <- entrez_mapping$ENTREZID
  rate <- if (length(g_genes) > 0L) nrow(entrez_mapping) / length(g_genes) else NA_real_

  cat("  Mapped ", length(entrez_ids), " / ", length(g_genes), " genes (", round(rate * 100, 1), "%)\n", sep = "")

  focal_entrez_map[[g]] <- entrez_mapping
  focal_entrez_ids[[g]] <- entrez_ids
  focal_mapping_rate[[g]] <- rate
}

# --- Enrichment: top 10 taxa ---------------------------------------------------
cat("\n=== Enrichment: genes for top 10 genera (union) ===\n")
enrich_top10 <- run_go_kegg(entrez_ids_top10, "top10")

readr::write_csv(enrich_top10$go_bp_df, file.path(path_tbl, "go_bp_top10_taxa.csv"))
readr::write_csv(enrich_top10$go_mf_df, file.path(path_tbl, "go_mf_top10_taxa.csv"))
readr::write_csv(enrich_top10$go_cc_df, file.path(path_tbl, "go_cc_top10_taxa.csv"))
readr::write_csv(enrich_top10$kegg_df, file.path(path_tbl, "kegg_top10_taxa.csv"))

gt_enrich_tbl <- function(df, title, subtitle) {
  gt_tbl <- df %>%
    gt::gt() %>%
    gt::tab_header(title = title, subtitle = subtitle)
  gt_tbl <- style_gt_significance(gt_tbl, df, alpha = 0.1)
  gt_tbl
}

gt_go_bp_top10 <- gt_enrich_tbl(
  enrich_top10$go_bp_df,
  title = "GO BP enrichment (top 10 genera gene set)",
  subtitle = "clusterProfiler enrichGO; BH-adjusted p < 0.1 highlighted"
)
gt_go_mf_top10 <- gt_enrich_tbl(
  enrich_top10$go_mf_df,
  title = "GO MF enrichment (top 10 genera gene set)",
  subtitle = "clusterProfiler enrichGO; BH-adjusted p < 0.1 highlighted"
)
gt_go_cc_top10 <- gt_enrich_tbl(
  enrich_top10$go_cc_df,
  title = "GO CC enrichment (top 10 genera gene set)",
  subtitle = "clusterProfiler enrichGO; BH-adjusted p < 0.1 highlighted"
)
gt_kegg_top10 <- gt_enrich_tbl(
  enrich_top10$kegg_df,
  title = "KEGG enrichment (top 10 genera gene set)",
  subtitle = "clusterProfiler enrichKEGG; BH-adjusted p < 0.1 highlighted"
)

gt::gtsave(gt_go_bp_top10, file.path(path_tbl, "go_bp_top10_taxa.html"))
gt::gtsave(gt_go_mf_top10, file.path(path_tbl, "go_mf_top10_taxa.html"))
gt::gtsave(gt_go_cc_top10, file.path(path_tbl, "go_cc_top10_taxa.html"))
gt::gtsave(gt_kegg_top10, file.path(path_tbl, "kegg_top10_taxa.html"))

saveRDS(enrich_top10$go_bp, file.path(path_stats, "go_bp_top10_taxa.rds"))
saveRDS(enrich_top10$go_mf, file.path(path_stats, "go_mf_top10_taxa.rds"))
saveRDS(enrich_top10$go_cc, file.path(path_stats, "go_cc_top10_taxa.rds"))
saveRDS(enrich_top10$kegg, file.path(path_stats, "kegg_top10_taxa.rds"))

save_gg_list(
  enrich_top10,
  prefix = "top10_taxa",
  title_prefix = "Functional enrichment",
  subtitle = "Genes with significant partial correlations to top 10 genera (by association count)"
)

# --- Enrichment: focal genera ---------------------------------------------------
focal_enrich <- list()
focal_gt <- list()

for (g in focal_genera) {
  slug <- focal_slug(g)
  cat("\n=== Enrichment: ", g, "-associated genes ===\n", sep = "")

  enrich_g <- run_go_kegg(focal_entrez_ids[[g]], slug)
  focal_enrich[[g]] <- enrich_g

  readr::write_csv(enrich_g$go_bp_df, file.path(path_tbl, paste0("go_bp_", slug, ".csv")))
  readr::write_csv(enrich_g$go_mf_df, file.path(path_tbl, paste0("go_mf_", slug, ".csv")))
  readr::write_csv(enrich_g$go_cc_df, file.path(path_tbl, paste0("go_cc_", slug, ".csv")))
  readr::write_csv(enrich_g$kegg_df, file.path(path_tbl, paste0("kegg_", slug, ".csv")))

  gt_go_bp <- gt_enrich_tbl(
    enrich_g$go_bp_df,
    title = paste0("GO BP enrichment (", g, " gene set)"),
    subtitle = "clusterProfiler enrichGO; BH-adjusted p < 0.1 highlighted"
  )
  gt_go_mf <- gt_enrich_tbl(
    enrich_g$go_mf_df,
    title = paste0("GO MF enrichment (", g, " gene set)"),
    subtitle = "clusterProfiler enrichGO; BH-adjusted p < 0.1 highlighted"
  )
  gt_go_cc <- gt_enrich_tbl(
    enrich_g$go_cc_df,
    title = paste0("GO CC enrichment (", g, " gene set)"),
    subtitle = "clusterProfiler enrichGO; BH-adjusted p < 0.1 highlighted"
  )
  gt_kegg <- gt_enrich_tbl(
    enrich_g$kegg_df,
    title = paste0("KEGG enrichment (", g, " gene set)"),
    subtitle = "clusterProfiler enrichKEGG; BH-adjusted p < 0.1 highlighted"
  )

  gt::gtsave(gt_go_bp, file.path(path_tbl, paste0("go_bp_", slug, ".html")))
  gt::gtsave(gt_go_mf, file.path(path_tbl, paste0("go_mf_", slug, ".html")))
  gt::gtsave(gt_go_cc, file.path(path_tbl, paste0("go_cc_", slug, ".html")))
  gt::gtsave(gt_kegg, file.path(path_tbl, paste0("kegg_", slug, ".html")))

  saveRDS(enrich_g$go_bp, file.path(path_stats, paste0("go_bp_", slug, ".rds")))
  saveRDS(enrich_g$go_mf, file.path(path_stats, paste0("go_mf_", slug, ".rds")))
  saveRDS(enrich_g$go_cc, file.path(path_stats, paste0("go_cc_", slug, ".rds")))
  saveRDS(enrich_g$kegg, file.path(path_stats, paste0("kegg_", slug, ".rds")))

  save_gg_list(
    enrich_g,
    prefix = slug,
    title_prefix = "Functional enrichment",
    subtitle = paste0("Genes with significant partial correlations to ", g)
  )

  focal_gt[[g]] <- list(go_bp = gt_go_bp, go_mf = gt_go_mf, go_cc = gt_go_cc, kegg = gt_kegg)
}

focal_four_combined_ggplots <- run_focal_four_combined_outputs(
  focal_enrich_lst = focal_enrich,
  display_order = focal_genera_display
)

# Project-relative paths (portable; work with knitr root.dir = project root)
path_res_rel <- file.path("Results", "07__FunctionalAnno")
path_fig_rel <- file.path(path_res_rel, "Figures")
path_tbl_rel <- file.path(path_res_rel, "Tables")
path_stats_rel <- file.path(path_res_rel, "Stats")

# --- Bundle --------------------------------------------------------------------
bundle <- list(
  meta = list(
    run_date = as.character(Sys.Date()),
    script = "Code/01__Analysis/07__FunctionalAnno.R",
    skipped = FALSE,
    n_sig_edges = nrow(correlation_data),
    n_genes_top10_unique = length(unique_genes_top10),
    n_entrez_top10 = length(entrez_ids_top10),
    mapping_rate_top10 = rate_top10,
    focal_genera = focal_genera,
    focal_gene_counts = stats::setNames(vapply(focal_genes, length, integer(1L)), names(focal_genes)),
    focal_mapping_rate = focal_mapping_rate,
    n_genes_culicoidibacter_unique = length(focal_genes[["Culicoidibacter"]]),
    n_entrez_culicoidibacter = length(focal_entrez_ids[["Culicoidibacter"]]),
    mapping_rate_culicoidibacter = focal_mapping_rate[["Culicoidibacter"]],
    top_n_genera = top_n_taxa,
    input_csv = path_csv_06
  ),
  tables = list(
    top_10_taxa_summary = top_10_taxa_summary,
    gene_counts_by_taxa = gene_counts_by_taxa,
    focal_by_comparison = focal_by_comparison,
    culicoidibacter_by_comparison = focal_by_comparison[["Culicoidibacter"]],
    go_bp_top10_taxa = enrich_top10$go_bp_df,
    go_mf_top10_taxa = enrich_top10$go_mf_df,
    go_cc_top10_taxa = enrich_top10$go_cc_df,
    kegg_top10_taxa = enrich_top10$kegg_df,
    focal_enrich_tables = lapply(focal_enrich, function(e) {
      list(go_bp = e$go_bp_df, go_mf = e$go_mf_df, go_cc = e$go_cc_df, kegg = e$kegg_df)
    }),
    go_bp_culicoidibacter = focal_enrich[["Culicoidibacter"]]$go_bp_df,
    go_mf_culicoidibacter = focal_enrich[["Culicoidibacter"]]$go_mf_df,
    go_cc_culicoidibacter = focal_enrich[["Culicoidibacter"]]$go_cc_df,
    kegg_culicoidibacter = focal_enrich[["Culicoidibacter"]]$kegg_df
  ),
  paths = list(
    figures = list(
      go_bp_top10 = file.path(path_fig_rel, "go_bp_top10_taxa_dotplot.png"),
      go_mf_top10 = file.path(path_fig_rel, "go_mf_top10_taxa_dotplot.png"),
      go_cc_top10 = file.path(path_fig_rel, "go_cc_top10_taxa_dotplot.png"),
      kegg_top10 = file.path(path_fig_rel, "kegg_top10_taxa_dotplot.png"),
      go_bp_culico = file.path(path_fig_rel, "go_bp_culicoidibacter_dotplot.png"),
      go_mf_culico = file.path(path_fig_rel, "go_mf_culicoidibacter_dotplot.png"),
      go_cc_culico = file.path(path_fig_rel, "go_cc_culicoidibacter_dotplot.png"),
      kegg_culico = file.path(path_fig_rel, "kegg_culicoidibacter_dotplot.png"),
      go_bp_focal_four = file.path(path_fig_rel, "go_bp_focal_four_dotplot.png"),
      go_mf_focal_four = file.path(path_fig_rel, "go_mf_focal_four_dotplot.png"),
      go_cc_focal_four = file.path(path_fig_rel, "go_cc_focal_four_dotplot.png"),
      kegg_focal_four = file.path(path_fig_rel, "kegg_focal_four_dotplot.png")
    ),
    tables = list(
      combined_edges_top10 = file.path(path_tbl_rel, "genes_for_top10_taxa_edges.csv"),
      focal_edges = stats::setNames(
        lapply(focal_genera, function(g) {
          file.path(path_tbl_rel, paste0(focal_slug(g), "_partial_correlation_edges.csv"))
        }),
        focal_genera
      ),
      culicoidibacter_edges = file.path(path_tbl_rel, "culicoidibacter_partial_correlation_edges.csv"),
      go_bp_top10 = file.path(path_tbl_rel, "go_bp_top10_taxa.csv"),
      go_mf_top10 = file.path(path_tbl_rel, "go_mf_top10_taxa.csv"),
      go_cc_top10 = file.path(path_tbl_rel, "go_cc_top10_taxa.csv"),
      kegg_top10 = file.path(path_tbl_rel, "kegg_top10_taxa.csv"),
      go_bp_culico = file.path(path_tbl_rel, "go_bp_culicoidibacter.csv"),
      go_mf_culico = file.path(path_tbl_rel, "go_mf_culicoidibacter.csv"),
      go_cc_culico = file.path(path_tbl_rel, "go_cc_culicoidibacter.csv"),
      kegg_culico = file.path(path_tbl_rel, "kegg_culicoidibacter.csv"),
      go_bp_top10_html = file.path(path_tbl_rel, "go_bp_top10_taxa.html"),
      go_mf_top10_html = file.path(path_tbl_rel, "go_mf_top10_taxa.html"),
      go_cc_top10_html = file.path(path_tbl_rel, "go_cc_top10_taxa.html"),
      kegg_top10_html = file.path(path_tbl_rel, "kegg_top10_taxa.html"),
      go_bp_culico_html = file.path(path_tbl_rel, "go_bp_culicoidibacter.html"),
      go_mf_culico_html = file.path(path_tbl_rel, "go_mf_culicoidibacter.html"),
      go_cc_culico_html = file.path(path_tbl_rel, "go_cc_culicoidibacter.html"),
      kegg_culico_html = file.path(path_tbl_rel, "kegg_culicoidibacter.html"),
      go_bp_focal_four_combined_source = file.path(path_tbl_rel, "go_bp_focal_four_combined_source.csv"),
      go_mf_focal_four_combined_source = file.path(path_tbl_rel, "go_mf_focal_four_combined_source.csv"),
      go_cc_focal_four_combined_source = file.path(path_tbl_rel, "go_cc_focal_four_combined_source.csv"),
      kegg_focal_four_combined_source = file.path(path_tbl_rel, "kegg_focal_four_combined_source.csv"),
      go_bp_focal_four_term_by_taxon_wide = file.path(path_tbl_rel, "go_bp_focal_four_term_by_taxon_wide.csv"),
      go_bp_focal_four_term_by_taxon_wide_html = file.path(path_tbl_rel, "go_bp_focal_four_term_by_taxon_wide.html"),
      go_mf_focal_four_term_by_taxon_wide = file.path(path_tbl_rel, "go_mf_focal_four_term_by_taxon_wide.csv"),
      go_mf_focal_four_term_by_taxon_wide_html = file.path(path_tbl_rel, "go_mf_focal_four_term_by_taxon_wide.html"),
      go_cc_focal_four_term_by_taxon_wide = file.path(path_tbl_rel, "go_cc_focal_four_term_by_taxon_wide.csv"),
      go_cc_focal_four_term_by_taxon_wide_html = file.path(path_tbl_rel, "go_cc_focal_four_term_by_taxon_wide.html"),
      kegg_focal_four_term_by_taxon_wide = file.path(path_tbl_rel, "kegg_focal_four_term_by_taxon_wide.csv"),
      kegg_focal_four_term_by_taxon_wide_html = file.path(path_tbl_rel, "kegg_focal_four_term_by_taxon_wide.html")
    ),
    stats_dir = path_stats_rel
  ),
  table_go_bp_top10_taxa = gt_go_bp_top10,
  table_go_mf_top10_taxa = gt_go_mf_top10,
  table_go_cc_top10_taxa = gt_go_cc_top10,
  table_kegg_top10_taxa = gt_kegg_top10,
  focal_gt = focal_gt,
  table_go_bp_culicoidibacter = focal_gt[["Culicoidibacter"]]$go_bp,
  table_go_mf_culicoidibacter = focal_gt[["Culicoidibacter"]]$go_mf,
  table_go_cc_culicoidibacter = focal_gt[["Culicoidibacter"]]$go_cc,
  table_kegg_culicoidibacter = focal_gt[["Culicoidibacter"]]$kegg,
  enrich = list(
    top10 = enrich_top10,
    focal = focal_enrich,
    culicoidibacter = focal_enrich[["Culicoidibacter"]]
  ),
  ggplots_focal_four_combined = focal_four_combined_ggplots
)

saveRDS(bundle, bundle_rds)
message("Saved bundle: ", bundle_rds)

sieler2026_sync_main_figures_from_manifest(driver_script = "07__FunctionalAnno.R")

message("07__FunctionalAnno complete.")

} # end nrow(correlation_data) >= 1
