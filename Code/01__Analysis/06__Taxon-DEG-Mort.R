# 06__Taxon-DEG-Mort.R
# Created by: Michael Sieler
# Date last updated: 2026-04-24
#
# CLI: optional `--focal-exports-only` refreshes focal-genus edge CSVs + focal-four Venns from
#   existing `combined_sig_partial_correlations.csv` without rerunning partial-correlation inference.
#
# Description: Partial correlation network between genus-level relative abundance and host gene
#   expression (Spearman pre-screen, then nonparametric partial correlation controlling for
#   Total.Worm.Count via nptest::np.cor.test), across exposure-regime pairwise contrasts vs
#   A- T- P-. Saves combined results, tables, figures, and taxon_deg_mort__bundle.rds under
#   Results/06__Taxon-DEG-Mort/.
#
# Expected input:  Run from Sieler2026 root; Results/04__DiffGeneExp/Stats/*, Results/03__DiffAbund/,
#   Data/r_objects/ps-list__*.rds (ps.list[["All"]]), Code/00__Setup/04__TaxonGeneNetworkHelpers.R
#   (perform_pairwise_analysis and related DEG×DAT helpers).
# Expected output:  Results/06__Taxon-DEG-Mort/{Figures,Tables,Stats} and bundle RDS
#   (figures include partial-correlation displays, focal-four gene-set Venns, and related CSV/HTML).

init_rel <- file.path("Code", "00__Setup", "00__InitializeEnvironment.R")
if (!file.exists(init_rel)) {
  stop(
    "Set working directory to the Sieler2026 repository root, or run:\n",
    "  Rscript Code/01__Analysis/06__Taxon-DEG-Mort.R [--focal-exports-only]\n",
    "from the project root so `", init_rel, "` resolves."
  )
}
args_trailing <- commandArgs(trailingOnly = TRUE)
flag_focal_exports_only <- isTRUE("--focal-exports-only" %in% args_trailing)

source(init_rel)
# --- Focal-four Venn helpers (formerly Code/01__Analysis/06b__FocalFourGeneVenn.R) ---
.focal_four_taxon_outline_colors <- c(
  Culicoidibacter = "#E41A1C",
  Shewanella = "#377EB8",
  Flavobacterium = "#4DAF4A",
  Cetobacterium = "#984EA3"
)

build_focal_gene_sets_for_venn <- function(focal_edges_lst, display_order) {
  all_l <- stats::setNames(vector("list", length(display_order)), display_order)
  pos_l <- all_l
  neg_l <- all_l
  for (g in display_order) {
    df_g <- focal_edges_lst[[g]]
    if (is.null(df_g) || nrow(df_g) < 1L) {
      all_l[[g]] <- character(0)
      pos_l[[g]] <- character(0)
      neg_l[[g]] <- character(0)
    } else {
      all_l[[g]] <- unique(as.character(df_g$gene_id))
      pos_l[[g]] <- unique(as.character(df_g$gene_id[df_g$correlation > 0]))
      neg_l[[g]] <- unique(as.character(df_g$gene_id[df_g$correlation < 0]))
    }
  }
  list(all = all_l, positive = pos_l, negative = neg_l)
}

export_focal_gene_sets_long_csv <- function(sets_three, path_csv) {
  rows <- list()
  for (nm in names(sets_three)) {
    L <- sets_three[[nm]]
    for (g in names(L)) {
      genes <- L[[g]]
      if (length(genes) < 1L) {
        next
      }
      rows[[paste(nm, g, sep = "__")]] <- tibble::tibble(
        direction = nm,
        Taxon = g,
        gene_id = genes
      )
    }
  }
  if (length(rows) < 1L) {
    return(invisible(NULL))
  }
  readr::write_csv(dplyr::bind_rows(rows), path_csv)
  invisible(path_csv)
}

focal_four_format_venn_subtitle <- function(txt, wrap_width = 50L) {
  txt <- as.character(txt)
  if (length(txt) != 1L || !nzchar(txt)) {
    return(txt)
  }
  paste(stringr::str_wrap(stringr::str_squish(txt), width = wrap_width), collapse = "\n")
}

save_focal_four_venn <- function(
    gene_sets_named,
    title,
    subtitle,
    stem_no_ext,
    sq_in = 6.5
) {
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    stop("Install ggVennDiagram (see Code/00__Setup/01__Libraries.R).")
  }
  n_nonempty <- sum(vapply(gene_sets_named, length, integer(1L)) > 0L)
  if (n_nonempty < 2L) {
    message("Skipping Venn (< 2 non-empty gene sets): ", stem_no_ext)
    return(invisible(NULL))
  }

  nm <- names(gene_sets_named)
  set_color_vec <- unname(.focal_four_taxon_outline_colors[nm])
  set_color_vec[is.na(set_color_vec)] <- "grey45"
  names(set_color_vec) <- nm

  fill_alpha <- 0.19
  fill_by_id <- vapply(nm, function(g) {
    grDevices::adjustcolor(set_color_vec[[g]], alpha.f = fill_alpha)
  }, character(1L))
  names(fill_by_id) <- as.character(seq_along(nm))

  venn <- ggVennDiagram::Venn(gene_sets_named)
  pdata <- ggVennDiagram::process_data(venn)
  setedge <- ggVennDiagram::venn_setedge(pdata)
  setmeta <- ggVennDiagram::venn_set(pdata)
  setedge_path <- setedge %>%
    dplyr::left_join(
      setmeta %>% dplyr::select("id", "name"),
      by = "id"
    )
  setedge_poly <- setedge_path %>%
    dplyr::mutate(poly_fill = fill_by_id[.data$id])
  setlbl <- ggVennDiagram::venn_setlabel(pdata) %>%
    dplyr::mutate(
      hj = dplyr::case_when(
        .data$name == "Culicoidibacter" ~ 0,
        .data$name == "Cetobacterium" ~ 1,
        TRUE ~ 0.5
      ),
      nx = dplyr::case_when(
        .data$name == "Culicoidibacter" ~ 0.014,
        .data$name == "Cetobacterium" ~ -0.014,
        TRUE ~ 0
      )
    )
  reglbl <- ggVennDiagram::venn_regionlabel(pdata)

  subtitle_use <- focal_four_format_venn_subtitle(subtitle)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      ggplot2::aes(x = .data$X, y = .data$Y, group = .data$id, fill = .data$poly_fill),
      data = setedge_poly,
      colour = NA
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_path(
      ggplot2::aes(
        x = .data$X,
        y = .data$Y,
        group = .data$id,
        colour = .data$name
      ),
      data = setedge_path,
      linewidth = 0.55
    ) +
    ggplot2::geom_label(
      ggplot2::aes(x = .data$X, y = .data$Y, label = .data$count),
      data = reglbl,
      family = "sans",
      size = 3,
      linewidth = 0,
      fill = ggplot2::alpha("white", 0.58),
      colour = "grey15"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = .data$X + .data$nx,
        y = .data$Y,
        label = .data$name,
        colour = .data$name,
        hjust = .data$hj
      ),
      data = setlbl,
      family = "sans",
      size = 3.2,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = set_color_vec[nm],
      breaks = nm,
      limits = nm,
      guide = "none"
    ) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = title, subtitle = subtitle_use) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13, margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        colour = "grey25",
        lineheight = 1.22,
        margin = ggplot2::margin(b = 6)
      ),
      plot.margin = ggplot2::margin(t = 10, r = 48, b = 14, l = 48)
    )

  w_in <- sq_in + 0.9
  h_in <- sq_in + 0.55
  ggplot2::ggsave(paste0(stem_no_ext, ".png"), p, width = w_in, height = h_in, units = "in", dpi = 300L)
  ggplot2::ggsave(paste0(stem_no_ext, ".pdf"), p, width = w_in, height = h_in, units = "in", device = "pdf")
  invisible(p)
}

summarize_focal_four_gene_overlap_one <- function(gene_sets_named, direction) {
  nm <- names(gene_sets_named)
  if (length(nm) < 1L) {
    return(tibble::tibble(direction = direction))
  }
  sets <- lapply(nm, function(t) {
    unique(as.character(gene_sets_named[[t]]))
  })
  names(sets) <- nm

  gene_taxon <- dplyr::bind_rows(lapply(nm, function(t) {
    g <- sets[[t]]
    if (length(g) < 1L) {
      return(tibble::tibble(gene_id = character(0), Taxon = character(0)))
    }
    tibble::tibble(gene_id = g, Taxon = t)
  })) %>%
    dplyr::distinct(.data$gene_id, .data$Taxon)

  by_g <- gene_taxon %>%
    dplyr::group_by(.data$gene_id) %>%
    dplyr::summarise(n_taxa = dplyr::n(), .groups = "drop")

  out <- tibble::tibble(direction = direction)
  out$n_union_genes <- if (nrow(by_g) < 1L) 0L else nrow(by_g)
  out$n_genes_in_all_four_taxa <- if (nrow(by_g) < 1L) {
    0L
  } else {
    sum(by_g$n_taxa == length(nm), na.rm = TRUE)
  }
  out$n_genes_in_exactly_one_taxon <- if (nrow(by_g) < 1L) {
    0L
  } else {
    sum(by_g$n_taxa == 1L, na.rm = TRUE)
  }
  out$n_genes_in_exactly_two_taxa <- if (nrow(by_g) < 1L) {
    0L
  } else {
    sum(by_g$n_taxa == 2L, na.rm = TRUE)
  }
  out$n_genes_in_exactly_three_taxa <- if (nrow(by_g) < 1L) {
    0L
  } else {
    sum(by_g$n_taxa == 3L, na.rm = TRUE)
  }
  out$n_genes_in_at_least_two_taxa <- if (nrow(by_g) < 1L) {
    0L
  } else {
    sum(by_g$n_taxa >= 2L, na.rm = TRUE)
  }

  for (t in nm) {
    col <- paste0("n_genes_input_", gsub("[^A-Za-z0-9]+", "_", t, perl = TRUE))
    out[[col]] <- length(sets[[t]])
  }

  if (length(nm) >= 2L) {
    for (i in seq_along(nm)) {
      for (j in seq_len(i - 1L)) {
        a <- nm[[j]]
        b <- nm[[i]]
        pair_lab <- paste0("n_intersect_", a, "__", b)
        out[[pair_lab]] <- length(intersect(sets[[a]], sets[[b]]))
      }
    }
  }

  out
}

summarize_focal_four_gene_overlap_table <- function(vsets_three) {
  dplyr::bind_rows(
    summarize_focal_four_gene_overlap_one(vsets_three$all, "all_significant_edges"),
    summarize_focal_four_gene_overlap_one(vsets_three$positive, "any_positive_partial_cor"),
    summarize_focal_four_gene_overlap_one(vsets_three$negative, "any_negative_partial_cor")
  )
}

gt_focal_four_overlap_summary <- function(overlap_tbl) {
  num_cols <- setdiff(names(overlap_tbl), "direction")
  gt::gt(overlap_tbl) %>%
    gt::tab_header(
      title = "Focal four taxa: overlap among partial-correlation gene sets",
      subtitle = paste0(
        "Gene symbols per taxon (BH FDR < 0.1 on partial correlations; module 06); ",
        "union = unique genes in any set; pairwise columns = |A ∩ B| (counts)."
      )
    ) %>%
    gt::cols_label(direction = "Gene set rule") %>%
    gt::fmt_number(columns = dplyr::all_of(num_cols), decimals = 0, sep_mark = ",") %>%
    gt::tab_options(table.font.size = gt::px(12))
}

#' Build focal-four Venns + overlap tables from significant partial-correlation edges.
#'
#' @param sig_partial Rows from `combined_sig_partial_correlations` (module 06).
#' @param path_fig Absolute path to `06__Taxon-DEG-Mort/Figures`.
#' @param path_tbl Absolute path to `06__Taxon-DEG-Mort/Tables`.
#' @param display_order Character vector of four taxon names (plot / Venn order).
#' @return Invisibly, a list with `venn_gene_overlap_summary` and `gt_venn_gene_overlap_summary`, or NULL.
run_focal_four_gene_venn_outputs <- function(sig_partial, path_fig, path_tbl, display_order) {
  if (nrow(sig_partial) < 1L) {
    message("Skipping focal-four gene Venns: no significant partial-correlation edges.")
    return(invisible(NULL))
  }

  focal_edges_lst <- stats::setNames(vector("list", length(display_order)), display_order)
  for (g in display_order) {
    focal_edges_lst[[g]] <- sig_partial %>% dplyr::filter(.data$TaxaID == g)
  }

  subt_venn_all <- "FDR-significant partial correlations to each genus (module 06)."
  subt_venn_pos <- "At least one positive partial correlation per genus (FDR < 0.1 edges)."
  subt_venn_neg <- "At least one negative partial correlation per genus (FDR < 0.1 edges)."

  vsets <- build_focal_gene_sets_for_venn(focal_edges_lst, display_order)
  export_focal_gene_sets_long_csv(
    vsets,
    file.path(path_tbl, "focal_four_gene_sets_venn_long.csv")
  )
  save_focal_four_venn(
    vsets$all,
    title = "Gene overlap: focal genera (all significant edges)",
    subtitle = subt_venn_all,
    stem_no_ext = file.path(path_fig, "focal_four_genes_venn_all")
  )
  save_focal_four_venn(
    vsets$positive,
    title = "Gene overlap: focal genera (positive partial correlation)",
    subtitle = subt_venn_pos,
    stem_no_ext = file.path(path_fig, "focal_four_genes_venn_positive")
  )
  save_focal_four_venn(
    vsets$negative,
    title = "Gene overlap: focal genera (negative partial correlation)",
    subtitle = subt_venn_neg,
    stem_no_ext = file.path(path_fig, "focal_four_genes_venn_negative")
  )

  overlap_tbl <- summarize_focal_four_gene_overlap_table(vsets)
  readr::write_csv(overlap_tbl, file.path(path_tbl, "focal_four_venn_gene_overlap_summary.csv"))
  gt_overlap <- gt_focal_four_overlap_summary(overlap_tbl)
  gt::gtsave(gt_overlap, file.path(path_tbl, "focal_four_venn_gene_overlap_summary.html"))

  invisible(
    list(
      venn_gene_overlap_summary = overlap_tbl,
      gt_venn_gene_overlap_summary = gt_overlap
    )
  )
}
if (isTRUE(flag_focal_exports_only)) {
  path_tbl <- file.path(path.results, "06__Taxon-DEG-Mort", "Tables")
  path_fig <- file.path(path.results, "06__Taxon-DEG-Mort", "Figures")
  dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
  in_csv <- file.path(path_tbl, "combined_sig_partial_correlations.csv")
  if (!file.exists(in_csv)) {
    stop("Required file missing: ", in_csv, "\nRun `Code/01__Analysis/06__Taxon-DEG-Mort.R` first.")
  }
  focal_genera <- c("Shewanella", "Culicoidibacter", "Flavobacterium", "Cetobacterium")
  focal_genera_display <- c("Culicoidibacter", "Shewanella", "Flavobacterium", "Cetobacterium")
  focal_slug <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
    x <- gsub("^_+|_+$", "", x, perl = TRUE)
    tolower(x)
  }
  sig_partial <- readr::read_csv(in_csv, show_col_types = FALSE)
  if (!"TaxaID" %in% names(sig_partial)) {
    stop("Expected column `TaxaID` not found in: ", in_csv)
  }
  for (g in focal_genera) {
    out_path <- file.path(path_tbl, paste0("partial_correlation_edges__", focal_slug(g), ".csv"))
    df_g <- sig_partial %>% dplyr::filter(.data$TaxaID == g)
    readr::write_csv(df_g, out_path)
  }
  message("Wrote focal-genus edge CSVs under: ", path_tbl)
  if (nrow(sig_partial) > 0L) {
    run_focal_four_gene_venn_outputs(sig_partial, path_fig, path_tbl, focal_genera_display)
    message("Updated focal-four gene Venns under: ", path_fig)
  } else {
    message("No rows in combined_sig_partial_correlations.csv; skipping Venns.")
  }
  quit(save = "no")
}

source(file.path(path.setup, "04__TaxonGeneNetworkHelpers.R"))

if (!exists("ps.list", inherits = TRUE) || is.null(ps.list) || !"All" %in% names(ps.list)) {
  stop(
    "ps.list not found or missing element 'All'. Run Code/00__Setup/04__DataPreProcess.R and ",
    "ensure ps-list__*.rds exists under Data/r_objects/."
  )
}

path_res <- file.path(path.results, "06__Taxon-DEG-Mort")
path_fig <- file.path(path_res, "Figures")
path_tbl <- file.path(path_res, "Tables")
path_stats <- file.path(path_res, "Stats")
dir.create(path_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(path_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(path_stats, recursive = TRUE, showWarnings = FALSE)

# Archive previous outputs (copy → verify → clear) so each run is clean but reproducible.
sieler2026_archive_module_outputs(
  path_res_module = path_res,
  module_name = "06__Taxon-DEG-Mort",
  subdirs = c("Figures", "Tables", "Stats"),
  clear_original = TRUE
)

bundle_rds <- file.path(path_stats, "taxon_deg_mort__bundle.rds")

path_deg_stats <- file.path(path.results, "04__DiffGeneExp", "Stats")
path_ma <- file.path(path.results, "03__DiffAbund", "Stats", "maaslin_ExposureRegimes_noTank", "significant_results.tsv")

top_n_genes <- 100L
top_n_taxa_plot <- 10L

# Focal taxa for reviewer-ready per-genus exports (overlap across modules)
focal_genera <- c("Shewanella", "Culicoidibacter", "Flavobacterium", "Cetobacterium")
# Display order for focal-four Venns and downstream plots (matches mortality palette)
focal_genera_display <- c("Culicoidibacter", "Shewanella", "Flavobacterium", "Cetobacterium")
focal_slug <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
  x <- gsub("^_+|_+$", "", x, perl = TRUE)
  tolower(x)
}

# --- Load inputs ----------------------------------------------------------------
deg_path <- file.path(path_deg_stats, "all_treatment_results.rds")
meta_path <- file.path(path_deg_stats, "metadata_final.rds")
counts_path <- file.path(path_deg_stats, "counts_int_filtered.rds")
for (p in c(deg_path, meta_path, counts_path, path_ma)) {
  if (!file.exists(p)) {
    stop("Required file missing: ", p)
  }
}

deg_results_all <- readRDS(deg_path)
metadata_final <- readRDS(meta_path)
counts_int_filt <- readRDS(counts_path)

dat_results_all <- readr::read_tsv(path_ma, show_col_types = FALSE) %>%
  dplyr::rename(taxa = "feature")

prep <- sieler2026_prepare_metadata_expr_for_pairwise(
  metadata_final = metadata_final,
  counts_int_filt = counts_int_filt,
  deg_results_all = deg_results_all
)

metadata <- prep$metadata
expr_counts <- prep$expr_counts

taxa_counts_all <- ps.list[["All"]] %>%
  phyloseq::otu_table() %>%
  base::as.data.frame()

set.seed(42)
results.all <- perform_pairwise_analysis(
  question_name = "All Treatments Analysis",
  dat_results = dat_results_all,
  deg_results = prep$deg_results,
  taxa_counts = taxa_counts_all,
  expr_counts = expr_counts,
  metadata = metadata,
  base_treatment = "A- T- P-",
  comparison_treatments = c(
    "A- T- P+",
    "A+ T- P-", "A+ T- P+",
    "A- T+ P-", "A- T+ P+",
    "A+ T+ P-", "A+ T+ P+"
  ),
  top_n = top_n_genes,
  top_by = "correlation"
)

# --- Save primary stats ---------------------------------------------------------
saveRDS(results.all, file.path(path_stats, "results_all.rds"))

sig_partial <- results.all$combined_sig_partial_correlations
readr::write_csv(
  sig_partial,
  file.path(path_tbl, "combined_sig_partial_correlations.csv")
)
readr::write_csv(
  results.all$combined_correlations,
  file.path(path_tbl, "combined_sig_spearman_correlations.csv")
)

per_taxon_summary <- if (nrow(sig_partial) > 0L) {
  sig_partial %>%
    dplyr::group_by(.data$TaxaID) %>%
    dplyr::summarise(
      n_gene_associations = dplyr::n(),
      n_positive = sum(.data$correlation > 0, na.rm = TRUE),
      n_negative = sum(.data$correlation < 0, na.rm = TRUE),
      mean_partial_r = mean(.data$correlation, na.rm = TRUE),
      mean_abs_partial_r = mean(abs(.data$correlation), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$n_gene_associations), dplyr::desc(.data$mean_abs_partial_r))
} else {
  tibble::tibble(
    TaxaID = character(),
    n_gene_associations = integer(),
    n_positive = integer(),
    n_negative = integer(),
    mean_partial_r = numeric(),
    mean_abs_partial_r = numeric()
  )
}
readr::write_csv(per_taxon_summary, file.path(path_tbl, "per_taxon_partial_correlation_summary.csv"))

# Per-focal-genus edge exports (for downstream functional annotation + reviewer requests)
if (nrow(sig_partial) > 0L) {
  for (g in focal_genera) {
    out_path <- file.path(path_tbl, paste0("partial_correlation_edges__", focal_slug(g), ".csv"))
    df_g <- sig_partial %>% dplyr::filter(.data$TaxaID == g)
    readr::write_csv(df_g, out_path)
  }
}

# Also export GT HTML for manuscript display.
gt_per_taxon_summary <- per_taxon_summary %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Genera: FDR-significant partial correlations with host genes",
    subtitle = "Infection burden (Total.Worm.Count) adjusted; summed across pairwise contrasts vs A- T- P-"
  )
num_cols <- intersect(c("mean_partial_r", "mean_abs_partial_r"), names(per_taxon_summary))
if (length(num_cols) > 0L) {
  gt_per_taxon_summary <- gt_per_taxon_summary %>%
    gt::fmt_number(columns = dplyr::all_of(num_cols), decimals = 3)
}
gt::gtsave(gt_per_taxon_summary, file.path(path_tbl, "per_taxon_partial_correlation_summary.html"))

# --- Focal four: gene-set Venns (same significant edges as partial-correlation exports) ----------
focal_four_venn_bundle <- NULL
if (nrow(sig_partial) > 0L) {
  focal_four_venn_bundle <- run_focal_four_gene_venn_outputs(
    sig_partial,
    path_fig,
    path_tbl,
    focal_genera_display
  )
}

# --- Figures ----------------------------------------------------------------------
x_gene <- 1
x_taxa <- 3

p_bipartite <- NULL
p_hist <- NULL
p_top10 <- NULL
p_all_sig <- NULL

if (nrow(sig_partial) > 0L) {
  top_taxa_ids <- sig_partial %>%
    dplyr::group_by(.data$TaxaID) %>%
    dplyr::summarise(
      n_correlations = dplyr::n(),
      mean_abs_r = mean(abs(.data$correlation), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$n_correlations), dplyr::desc(.data$mean_abs_r)) %>%
    dplyr::slice_head(n = top_n_taxa_plot) %>%
    dplyr::pull(.data$TaxaID)

  filt_edges <- sig_partial %>% dplyr::filter(.data$TaxaID %in% top_taxa_ids)

  p_hist <- filt_edges %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$correlation, fill = .data$comparison_name)) +
    ggplot2::geom_histogram(bins = 30L, alpha = 0.7, position = "identity") +
    theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
    ggplot2::labs(
      title = "Distribution of partial correlations (infection burden controlled)",
      subtitle = paste0("Top ", top_n_taxa_plot, " taxa by number of significant associations"),
      x = "Partial correlation",
      y = "Count",
      fill = "Comparison"
    ) +
    ggplot2::facet_wrap(~comparison_name, scales = "free_y")

  gene_levels <- filt_edges %>% dplyr::pull(.data$gene_id) %>% unique() %>% sort()
  taxa_levels <- filt_edges %>% dplyr::pull(.data$TaxaID) %>% unique() %>% sort()
  n_genes <- length(gene_levels)
  n_taxa <- length(taxa_levels)

  gene_pos <- tibble::tibble(gene_id = gene_levels, y = seq_len(n_genes))
  taxa_pos <- tibble::tibble(
    TaxaID = taxa_levels,
    y = seq(1, n_genes, length.out = n_taxa)
  )

  edges <- filt_edges %>%
    dplyr::left_join(gene_pos, by = "gene_id") %>%
    dplyr::rename(y_gene = y) %>%
    dplyr::left_join(taxa_pos, by = "TaxaID") %>%
    dplyr::rename(y_taxa = y)

  top_genes_by_cor <- filt_edges %>%
    dplyr::group_by(.data$gene_id) %>%
    dplyr::summarise(max_abs_correlation = max(abs(.data$correlation), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$max_abs_correlation)) %>%
    dplyr::slice_head(n = top_n_genes) %>%
    dplyr::pull(.data$gene_id)

  top_edges <- edges %>%
    dplyr::filter(.data$gene_id %in% top_genes_by_cor) %>%
    dplyr::filter(!is.na(.data$y_gene), !is.na(.data$y_taxa)) %>%
    dplyr::arrange(dplyr::desc(abs(.data$correlation))) %>%
    dplyr::slice_head(n = top_n_genes)

  gene_pos_f <- gene_pos %>%
    dplyr::filter(.data$gene_id %in% top_genes_by_cor) %>%
    dplyr::mutate(y = seq_len(dplyr::n()))
  taxa_conn <- top_edges %>% dplyr::pull(.data$TaxaID) %>% unique()
  n_gf <- nrow(gene_pos_f)
  taxa_pos_f <- taxa_pos %>%
    dplyr::filter(.data$TaxaID %in% taxa_conn) %>%
    dplyr::mutate(y = seq(1, n_gf, length.out = dplyr::n()))

  top_edges_u <- top_edges %>%
    dplyr::select(-dplyr::any_of(c("y_gene", "y_taxa"))) %>%
    dplyr::left_join(gene_pos_f, by = "gene_id") %>%
    dplyr::rename(y_gene = y) %>%
    dplyr::left_join(taxa_pos_f, by = "TaxaID") %>%
    dplyr::rename(y_taxa = y)

  p_bipartite <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = top_edges_u,
      ggplot2::aes(
        x = x_gene, xend = x_taxa,
        y = .data$y_gene, yend = .data$y_taxa,
        color = .data$correlation
      ),
      linewidth = 0.8
    ) +
    ggplot2::geom_text(
      data = gene_pos_f,
      ggplot2::aes(x = x_gene - 0.05, y = .data$y, label = .data$gene_id),
      hjust = 1, size = 2.5
    ) +
    ggplot2::geom_text(
      data = taxa_pos_f,
      ggplot2::aes(x = x_taxa + 0.05, y = .data$y, label = .data$TaxaID),
      hjust = 0, size = 2.5
    ) +
    ggplot2::scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1, 1),
      name = "Partial r"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = "Bipartite network: gene-taxa partial correlations",
      subtitle = "Controlling for infection burden (top gene-taxa pairs by |partial r| within top taxa)",
      caption = "Blue = negative, red = positive"
    ) +
    ggplot2::xlim(0.5, 3.5)

  top10_scatter <- per_taxon_summary %>%
    dplyr::arrange(
      dplyr::desc(.data$n_gene_associations),
      dplyr::desc(.data$mean_abs_partial_r)
    ) %>%
    dplyr::slice_head(n = 10L) %>%
    dplyr::mutate(
      is_focal = .data$TaxaID %in% focal_genera,
      # Keep focal outline ring sign-based (neg = blue, pos = red) but make label text taxon-based.
      outline_color = dplyr::case_when(
        !.data$is_focal ~ "black",
        .data$mean_partial_r < 0 ~ "#2166AC",
        .data$mean_partial_r > 0 ~ "#B2182B",
        TRUE ~ "gray35"
      ),
      label_color = dplyr::case_when(
        !.data$is_focal ~ "black",
        TRUE ~ unname(.focal_four_taxon_outline_colors[as.character(.data$TaxaID)])
      ),
      label_face = dplyr::if_else(
        .data$is_focal,
        "bold",
        "plain"
      ),
      label_size = dplyr::if_else(.data$is_focal, 4.25, 3)
    )

  set.seed(42)
  p_top10 <- ggplot2::ggplot(
    top10_scatter,
    ggplot2::aes(
      x = .data$mean_abs_partial_r,
      y = .data$n_gene_associations
    )
  ) +
    ggplot2::geom_point(
      shape = 21L,
      stroke = 0.6,
      color = "black",
      size = 6,
      ggplot2::aes(fill = .data$mean_partial_r),
      alpha = 1
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(top10_scatter, .data$is_focal),
      ggplot2::aes(
        x = .data$mean_abs_partial_r,
        y = .data$n_gene_associations,
        color = .data$outline_color
      ),
      shape = 21L,
      fill = NA,
      stroke = 1.15,
      size = 7,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Mean Partial Correlation",
      breaks = c(-1, -0.5, 0, 0.5, 1)
    ) +
    ggrepel::geom_label_repel(
      ggplot2::aes(
        label = .data$TaxaID,
        color = .data$label_color,
        fontface = .data$label_face,
        size = .data$label_size
      ),
      fill = ggplot2::alpha("white", 0.94),
      label.size = 0.32,
      label.padding = grid::unit(0.28, "lines"),
      point.padding = grid::unit(2, "lines"),
      box.padding = grid::unit(0.65, "lines"),
      min.segment.length = 0,
      segment.size = 0.25,
      segment.color = ggplot2::alpha("grey45", 0.55),
      force = 2,
      max.overlaps = 20L,
      seed = 42,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_identity() +
    theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(0.65), lineheight = 1.25)
    ) +
    ggplot2::labs(
      title = "Top 10 taxa by number of significant partial correlations",
      x = "Mean absolute partial correlation",
      y = "Number of significant gene associations",
      caption = "Blue = negative, red = positive\n(point fill = mean signed partial correlation)"
    )

  # All genera with ≥1 FDR-significant partial-correlation edge (same summary table); focal four labeled
  all_sig_scatter <- per_taxon_summary %>%
    dplyr::arrange(
      dplyr::desc(.data$n_gene_associations),
      dplyr::desc(.data$mean_abs_partial_r)
    ) %>%
    dplyr::mutate(
      is_focal = .data$TaxaID %in% focal_genera,
      outline_color = dplyr::case_when(
        !.data$is_focal ~ "black",
        .data$mean_partial_r < 0 ~ "#2166AC",
        .data$mean_partial_r > 0 ~ "#B2182B",
        TRUE ~ "gray35"
      ),
      label_color = dplyr::case_when(
        !.data$is_focal ~ "black",
        TRUE ~ unname(.focal_four_taxon_outline_colors[as.character(.data$TaxaID)])
      ),
      label_face = dplyr::if_else(
        .data$is_focal,
        "bold",
        "plain"
      ),
      label_size = dplyr::if_else(.data$is_focal, 4.25, 2.5)
    )

  # Room below y = 0 for ggrepel (e.g. low-count focal genera) without implying negative counts
  y_max_pts <- max(all_sig_scatter$n_gene_associations, na.rm = TRUE)
  y_hi <- max(25, ceiling(y_max_pts * 1.06 / 25) * 25)
  y_breaks <- seq(0, y_hi, by = 25)

  set.seed(42)
  p_all_sig <- ggplot2::ggplot(
    all_sig_scatter,
    ggplot2::aes(
      x = .data$mean_abs_partial_r,
      y = .data$n_gene_associations
    )
  ) +
    ggplot2::geom_point(
      data = dplyr::filter(all_sig_scatter, !.data$is_focal),
      ggplot2::aes(fill = .data$mean_partial_r),
      shape = 21L,
      color = "grey40",
      stroke = 0.35,
      size = 3,
      alpha = 0.72
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(all_sig_scatter, .data$is_focal),
      ggplot2::aes(fill = .data$mean_partial_r),
      shape = 21L,
      color = "black",
      stroke = 0.55,
      size = 5.8,
      alpha = 1
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(all_sig_scatter, .data$is_focal),
      ggplot2::aes(
        x = .data$mean_abs_partial_r,
        y = .data$n_gene_associations,
        color = .data$outline_color
      ),
      shape = 21L,
      fill = NA,
      stroke = 1.12,
      size = 6.6,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Mean Partial Correlation",
      breaks = c(-1, -0.5, 0, 0.5, 1),
      guide = ggplot2::guide_colorbar(order = 1L)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-20, y_hi),
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    ggrepel::geom_label_repel(
      data = dplyr::filter(all_sig_scatter, .data$is_focal),
      ggplot2::aes(
        label = .data$TaxaID,
        color = .data$label_color,
        fontface = .data$label_face,
        size = .data$label_size
      ),
      fill = ggplot2::alpha("white", 0.94),
      label.size = 0.35,
      label.padding = grid::unit(0.3, "lines"),
      point.padding = grid::unit(1.4, "lines"),
      box.padding = grid::unit(0.75, "lines"),
      min.segment.length = 0,
      segment.size = 0.28,
      segment.color = ggplot2::alpha("grey45", 0.55),
      force = 3,
      max.overlaps = 30L,
      seed = 42,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_size_identity() +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      data = tibble::tibble(
        TaxaID = factor(focal_genera_display, levels = focal_genera_display),
        mean_abs_partial_r = stats::median(all_sig_scatter$mean_abs_partial_r, na.rm = TRUE),
        n_gene_associations = stats::median(all_sig_scatter$n_gene_associations, na.rm = TRUE)
      ),
      mapping = ggplot2::aes(
        x = .data$mean_abs_partial_r,
        y = .data$n_gene_associations,
        color = .data$TaxaID
      ),
      inherit.aes = FALSE,
      alpha = 0,
      size = 0.05,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      name = "Focal Taxa",
      values = .focal_four_taxon_outline_colors[as.character(focal_genera_display)],
      breaks = focal_genera_display,
      limits = focal_genera_display,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        order = 2L,
        nrow = 1L,
        byrow = TRUE,
        override.aes = list(alpha = 1, size = 3.2)
      )
    ) +
    theme_sieler2026_publication_with_grid(base_size = 14, legend_position = "bottom") +
    ggplot2::theme(
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.15, "cm"),
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = "All taxa with significant partial correlations (FDR < 0.1)",
      subtitle = paste0(
        nrow(all_sig_scatter),
        " significant genera; focal genera (boxed labels)"
      ),
      x = "Mean absolute partial correlation",
      y = "Number of significant gene associations",
      caption = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off")
} else {
  empty_lab <- "No significant partial correlations (FDR < 0.1)"
  p_bipartite <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = empty_lab) +
    ggplot2::theme_void()
  p_hist <- p_bipartite
  p_top10 <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = empty_lab) +
    ggplot2::theme_void()
  p_all_sig <- p_top10
}

save_gg <- function(plot, pdf_path, png_path, w_in, h_in) {
  ggplot2::ggsave(pdf_path, plot, width = w_in, height = h_in, units = "in", device = "pdf")
  ggplot2::ggsave(png_path, plot, width = w_in, height = h_in, units = "in", dpi = 300L)
}

save_gg(
  p_bipartite,
  file.path(path_fig, "bipartite_partial_correlations.pdf"),
  file.path(path_fig, "bipartite_partial_correlations.png"),
  6, 12
)
save_gg(
  p_hist,
  file.path(path_fig, "partial_correlation_histogram.pdf"),
  file.path(path_fig, "partial_correlation_histogram.png"),
  14, 10
)
save_gg(
  p_top10,
  file.path(path_fig, "top10_taxa_partial_cor_scatter.pdf"),
  file.path(path_fig, "top10_taxa_partial_cor_scatter.png"),
  8, 6.5
)
save_gg(
  p_all_sig,
  file.path(path_fig, "all_sig_taxa_partial_cor_scatter.pdf"),
  file.path(path_fig, "all_sig_taxa_partial_cor_scatter.png"),
  9, 7
)

# Project-relative paths (portable; work with knitr root.dir = project root)
path_res_rel <- file.path("Results", "06__Taxon-DEG-Mort")
path_fig_rel <- file.path(path_res_rel, "Figures")
path_tbl_rel <- file.path(path_res_rel, "Tables")
path_stats_rel <- file.path(path_res_rel, "Stats")

# --- Bundle -----------------------------------------------------------------------
bundle <- list(
  meta = list(
    run_date = as.character(Sys.Date()),
    script = "Code/01__Analysis/06__Taxon-DEG-Mort.R",
    n_sig_partial_edges = nrow(sig_partial),
    n_comparisons = length(results.all$pairwise_results),
    focal_genera_display = focal_genera_display
  ),
  paths = list(
    figures = list(
      bipartite_pdf = file.path(path_fig_rel, "bipartite_partial_correlations.pdf"),
      bipartite_png = file.path(path_fig_rel, "bipartite_partial_correlations.png"),
      histogram_pdf = file.path(path_fig_rel, "partial_correlation_histogram.pdf"),
      histogram_png = file.path(path_fig_rel, "partial_correlation_histogram.png"),
      top10_pdf = file.path(path_fig_rel, "top10_taxa_partial_cor_scatter.pdf"),
      top10_png = file.path(path_fig_rel, "top10_taxa_partial_cor_scatter.png"),
      all_sig_pdf = file.path(path_fig_rel, "all_sig_taxa_partial_cor_scatter.pdf"),
      all_sig_png = file.path(path_fig_rel, "all_sig_taxa_partial_cor_scatter.png"),
      focal_four_venn_all = file.path(path_fig_rel, "focal_four_genes_venn_all.png"),
      focal_four_venn_positive = file.path(path_fig_rel, "focal_four_genes_venn_positive.png"),
      focal_four_venn_negative = file.path(path_fig_rel, "focal_four_genes_venn_negative.png")
    ),
    tables = list(
      combined_sig_partial = file.path(path_tbl_rel, "combined_sig_partial_correlations.csv"),
      focal_partial_edges = stats::setNames(
        lapply(focal_genera, function(g) {
          file.path(path_tbl_rel, paste0("partial_correlation_edges__", focal_slug(g), ".csv"))
        }),
        focal_genera
      ),
      per_taxon_summary = file.path(path_tbl_rel, "per_taxon_partial_correlation_summary.csv"),
      per_taxon_summary_html = file.path(path_tbl_rel, "per_taxon_partial_correlation_summary.html"),
      focal_four_gene_sets_venn_long = file.path(path_tbl_rel, "focal_four_gene_sets_venn_long.csv"),
      focal_four_venn_gene_overlap_summary = file.path(path_tbl_rel, "focal_four_venn_gene_overlap_summary.csv"),
      focal_four_venn_gene_overlap_summary_html = file.path(path_tbl_rel, "focal_four_venn_gene_overlap_summary.html")
    ),
    stats_dir = path_stats_rel
  ),
  tables = list(
    per_taxon_summary = per_taxon_summary,
    focal_four_venn_gene_overlap_summary = if (!is.null(focal_four_venn_bundle)) {
      focal_four_venn_bundle$venn_gene_overlap_summary
    } else {
      NULL
    }
  ),
  table_per_taxon_summary = gt_per_taxon_summary,
  table_focal_four_venn_gene_overlap = if (!is.null(focal_four_venn_bundle)) {
    focal_four_venn_bundle$gt_venn_gene_overlap_summary
  } else {
    NULL
  },
  ggplots = list(all_sig_partial_cor_scatter = p_all_sig)
)
saveRDS(bundle, bundle_rds)
message("Saved bundle: ", bundle_rds)

sieler2026_sync_main_figures_from_manifest(driver_script = "06__Taxon-DEG-Mort.R")

message("06__Taxon-DEG-Mort complete.")
