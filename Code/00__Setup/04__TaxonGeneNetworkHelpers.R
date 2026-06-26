# 04__TaxonGeneNetworkHelpers.R
# Created by: Michael Sieler
# Date last updated: 2026-04-26
#
# Description: DEG×DAT integration helpers (gene × genus partial correlations) formerly split across
#   04b__DEGxDAT_Functions.R, plus Sieler2026 adapters for Sample_RNA–based DESeq2 outputs and gut_id_num IDs.
#
# Expected input: path.setup set (via 00__InitializeEnvironment.R); libraries attached.
# Expected output: perform_pairwise_analysis, perform_integrated_analysis, gen_get_top_correlations,
#   sieler2026_prepare_metadata_expr_for_pairwise in the search path.

# =============================================================================
# DEG×DAT core (formerly Code/00__Setup/04b__DEGxDAT_Functions.R)
# =============================================================================

# Helper function to get top correlations by correlation, taxa, or gene
gen_get_top_correlations <- function(correlation_results, top_n = 100, top_by = c("correlation", "taxa", "gene")) {
  top_by <- match.arg(top_by)
  
  # Check if we have any significant correlations
  sig_correlations <- correlation_results %>%
    dplyr::filter(fdr < 0.05)
  
  if(nrow(sig_correlations) == 0) {
    cat("WARNING: No significant correlations found (FDR < 0.05). Returning empty data frame.\n")
    return(data.frame())
  }
  
  if (top_by == "correlation") {
    top_correlations <- sig_correlations %>%
      dplyr::arrange(dplyr::desc(abs_correlation)) %>%
      head(top_n)
  } else if (top_by == "taxa") {
    top_taxa <- sig_correlations %>%
      dplyr::group_by(TaxaID) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      head(top_n) %>%
      dplyr::pull(TaxaID)
    
    if(length(top_taxa) == 0) {
      cat("WARNING: No taxa found with significant correlations. Returning empty data frame.\n")
      return(data.frame())
    }
    
    top_correlations <- sig_correlations %>%
      dplyr::filter(TaxaID %in% top_taxa)
  } else if (top_by == "gene") {
    top_genes <- sig_correlations %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      head(top_n) %>%
      dplyr::pull(gene_id)
    
    if(length(top_genes) == 0) {
      cat("WARNING: No genes found with significant correlations. Returning empty data frame.\n")
      return(data.frame())
    }
    
    top_correlations <- sig_correlations %>%
      dplyr::filter(gene_id %in% top_genes)
  }
  
  return(top_correlations)
}

# Function to perform integrated analysis for a single research question
perform_integrated_analysis <- function(question_name, 
                                       dat_results, 
                                       deg_results, 
                                       taxa_counts, 
                                       expr_counts, 
                                       metadata, 
                                       treatment_comparison,
                                       top_n = 100,
                                       top_by = c("correlation", "taxa", "gene")) {
  # This function performs the full DEGxDAT analysis for a single research question (single treatment comparison).
  # INPUTS:
  #   - question_name: string, name of the research question
  #   - dat_results: data frame of DAT results
  #   - deg_results: data frame of DEG results
  #   - taxa_counts: OTU/taxa count matrix
  #   - expr_counts: gene expression count matrix
  #   - metadata: sample metadata
  #   - treatment_comparison: vector of treatment names to compare
  #   - top_n: (deprecated - kept for backwards compatibility) no longer used
  #   - top_by: (deprecated - kept for backwards compatibility) no longer used
  #   Note: Partial correlation analysis now uses ALL significant correlations (FDR < 0.05),
  #         not just the top N, to ensure complete coverage of all significant gene-taxa pairs.
  # OUTPUT: list with all results for this comparison (correlations, significant genes/taxa, etc.)

  cat("\n=== ANALYSIS FOR:", question_name, "===\n")
  cat("Treatment comparison:", paste(treatment_comparison, collapse = " vs "), "\n")
  
  # Set seed for reproducibility
  set.seed(42)
  
  # Filter for significant DEGs and DATs
  # Check column names - DEG results may use 'Gene' instead of 'gene'
  if("Gene" %in% colnames(deg_results)) {
    deg_results <- deg_results %>% dplyr::rename(gene = Gene)  # Standardize to 'gene' column name
  }
  
  # Filter DEG results for this specific treatment comparison
  # If Treatment column exists, filter by comparison treatment
  if("Treatment" %in% colnames(deg_results)) {
    # Filter for both treatments in the comparison (all DEGs that are significant in either treatment)
    deg_sig <- deg_results %>%
      dplyr::filter(Treatment %in% treatment_comparison, padj < 0.05)
  } else {
    # If no Treatment column, use all significant results
    deg_sig <- deg_results %>%
      dplyr::filter(padj < 0.05)
  }
  
  # Check column names - DAT results may use 'qval', 'q.value', 'qval_joint', or 'qval_individual'
  # Note: We prioritize qval_individual over qval_joint to match Maaslin2 behavior (qval_individual is less conservative)
  # Debug: print available columns
  cat("DAT results columns:", paste(colnames(dat_results), collapse = ", "), "\n")
  
  if("qval" %in% colnames(dat_results)) {
    cat("Using 'qval' column for filtering\n")
    dat_sig <- dat_results %>%
      dplyr::filter(qval < 0.05)
  } else if("qval_individual" %in% colnames(dat_results)) {
    cat("Using 'qval_individual' column for filtering (matches Maaslin2 qval behavior)\n")
    dat_sig <- dat_results %>%
      dplyr::filter(qval_individual < 0.05) %>%
      dplyr::rename(qval = qval_individual)
  } else if("qval_joint" %in% colnames(dat_results)) {
    cat("Using 'qval_joint' column for filtering (more conservative than qval_individual)\n")
    dat_sig <- dat_results %>%
      dplyr::filter(qval_joint < 0.05) %>%
      dplyr::rename(qval = qval_joint)
  } else if("q.value" %in% colnames(dat_results)) {
    cat("Using 'q.value' column for filtering\n")
    dat_sig <- dat_results %>%
      dplyr::filter(q.value < 0.05) %>%
      dplyr::rename(qval = q.value)
  } else if("pval_joint" %in% colnames(dat_results)) {
    # Fallback to p-value if no q-value available
    cat("Using 'pval_joint' column for filtering (will calculate FDR)\n")
    dat_sig <- dat_results %>%
      dplyr::filter(pval_joint < 0.05) %>%
      dplyr::mutate(qval = p.adjust(pval_joint, method = "BH"))
  } else if("pval_individual" %in% colnames(dat_results)) {
    # Fallback to p-value if no q-value available
    cat("Using 'pval_individual' column for filtering (will calculate FDR)\n")
    dat_sig <- dat_results %>%
      dplyr::filter(pval_individual < 0.05) %>%
      dplyr::mutate(qval = p.adjust(pval_individual, method = "BH"))
  } else {
    stop("No q-value or p-value column found in DAT results. Available columns: ", paste(colnames(dat_results), collapse = ", "))
  }
  
  cat("Significant genes:", nrow(deg_sig), "\n")
  cat("Significant taxa:", nrow(dat_sig), "\n")
  
  # Check if we have any significant genes or taxa
  if(nrow(deg_sig) == 0 || nrow(dat_sig) == 0) {
    cat("WARNING: No significant genes or taxa found. Returning empty results.\n")
    return(list(
      question_name = question_name,
      treatment_comparison = treatment_comparison,
      correlation_results = data.frame(),
      sig_correlations = data.frame(),
      partial_cor_results = data.frame(),
      sig_partial_correlations = data.frame(),
      failed_pairs = list(),
      deg_sig = deg_sig,
      dat_sig = dat_sig,
      samples = character(0)
    ))
  }
  
  # Get relevant samples for this comparison
  rna_samples <- metadata %>%
    dplyr::filter(Treatment %in% treatment_comparison) %>%
    dplyr::filter(Sample %in% colnames(expr_counts)) %>%
    dplyr::pull(Sample) %>%
    as.character()
  
  cat("Samples for analysis:", length(rna_samples), "\n")
  
  # Check if we have enough samples
  if(length(rna_samples) < 3) {
    cat("WARNING: Insufficient samples for analysis. Returning empty results.\n")
    return(list(
      question_name = question_name,
      treatment_comparison = treatment_comparison,
      correlation_results = data.frame(),
      sig_correlations = data.frame(),
      partial_cor_results = data.frame(),
      sig_partial_correlations = data.frame(),
      failed_pairs = list(),
      deg_sig = deg_sig,
      dat_sig = dat_sig,
      samples = rna_samples
    ))
  }
  
  # Filter expression data for significant genes and relevant samples
  # Check if gene column exists, if genes are in row names, or if we need to use Gene
  if("gene" %in% colnames(expr_counts)) {
    gene_col <- "gene"
  } else if("Gene" %in% colnames(expr_counts)) {
    gene_col <- "Gene"
    expr_counts <- expr_counts %>% dplyr::rename(gene = Gene)
  } else if(!is.null(rownames(expr_counts)) && length(rownames(expr_counts)) > 0) {
    # If genes are in row names, convert to column
    expr_counts <- expr_counts %>%
      tibble::rownames_to_column(var = "gene")
    gene_col <- "gene"
    cat("Gene identifiers found in row names. Converted to 'gene' column.\n")
  } else {
    stop("No gene column or row names found in expression counts")
  }
  
  # Check for gene_id and gene_name columns
  if(!"gene_id" %in% colnames(expr_counts)) {
    # If gene_id doesn't exist, create it from gene column
    expr_counts <- expr_counts %>%
      dplyr::mutate(gene_id = gene)
  }
  if(!"gene_name" %in% colnames(expr_counts)) {
    # If gene_name doesn't exist, create it from gene column
    expr_counts <- expr_counts %>%
      dplyr::mutate(gene_name = gene)
  }
  
  filtered_expr <- expr_counts %>%
    dplyr::filter(gene %in% deg_sig$gene) %>%
    dplyr::select(gene, gene_id, gene_name, dplyr::any_of(rna_samples))
  
  # Filter taxa data for significant taxa and relevant samples
  # Check if taxa column is named 'taxa' or 'feature'
  if(!"taxa" %in% colnames(dat_sig)) {
    if("feature" %in% colnames(dat_sig)) {
      dat_sig <- dat_sig %>% dplyr::rename(taxa = feature)
    } else {
      stop("No 'taxa' or 'feature' column found in DAT results")
    }
  }
  
  # CRITICAL: Get sample mapping from phyloseq object
  # The taxa data columns are phyloseq sample names (like "f16", "f20")
  # We need to map these to expression sample IDs (like "16", "20", "46")
  cat("=== CREATING SAMPLE MAPPING ===\n")
  if(exists("ps.list") && "All" %in% names(ps.list)) {
    # Get sample data from phyloseq object
    ps_sample_data <- ps.list[["All"]] %>%
      microViz::samdat_tbl() %>%
      tibble::rownames_to_column(var = "phyloseq_sample")
    
    # Check if fecal.sample.number exists
    if("fecal.sample.number" %in% colnames(ps_sample_data)) {
      # Create mapping: phyloseq sample name -> expression sample ID
      sample_mapping <- ps_sample_data %>%
        dplyr::select(phyloseq_sample, fecal.sample.number) %>%
        dplyr::mutate(
          # Strip "f" prefix from fecal.sample.number to get expression sample ID
          expr_sample_id = gsub("^f", "", fecal.sample.number)
        ) %>%
        dplyr::filter(!is.na(fecal.sample.number) & fecal.sample.number != "")
      
      # Get phyloseq sample names that correspond to our RNA samples
      rna_phyloseq_samples_mapping <- sample_mapping %>%
        dplyr::filter(expr_sample_id %in% rna_samples)
      
      # The OTU table columns should match fecal.sample.number
      taxa_cols <- colnames(taxa_counts)
      rna_phyloseq_samples <- rna_phyloseq_samples_mapping %>%
        dplyr::filter(fecal.sample.number %in% taxa_cols) %>%
        dplyr::pull(fecal.sample.number)
      
      if(length(rna_phyloseq_samples) == 0) {
        # Fallback: try matching by phyloseq_sample
        rna_phyloseq_samples <- rna_phyloseq_samples_mapping %>%
          dplyr::filter(phyloseq_sample %in% taxa_cols) %>%
          dplyr::pull(phyloseq_sample)
      }
    } else {
      # Fallback: assume column names are like "f1", "f2" and match directly
      taxa_cols <- colnames(taxa_counts)
      rna_phyloseq_samples <- taxa_cols[taxa_cols %in% paste0("f", rna_samples)]
      
      # Create a simple mapping
      sample_mapping <- data.frame(
        phyloseq_sample = rna_phyloseq_samples,
        expr_sample_id = gsub("^f", "", rna_phyloseq_samples),
        fecal.sample.number = rna_phyloseq_samples
      )
    }
  } else {
    stop("Phyloseq object not found. Please ensure ps.list is loaded.")
  }
  
  # CRITICAL: Map ASV IDs to genus names
  # dat_sig$taxa contains genus names, but taxa_counts has ASV IDs as row names
  # We need to get the taxonomy table from phyloseq to map ASV IDs to genus names
  cat("=== MAPPING ASV IDs TO GENUS NAMES ===\n")
  if(exists("ps.list") && "All" %in% names(ps.list)) {
    # Get taxonomy table from phyloseq object
    tax_table <- ps.list[["All"]] %>%
      phyloseq::tax_table() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "ASV_ID")
    
    # Check which column contains genus names
    genus_col <- NULL
    if("Genus" %in% colnames(tax_table)) {
      genus_col <- "Genus"
    } else if("genus" %in% colnames(tax_table)) {
      tax_table <- tax_table %>% dplyr::rename(Genus = genus)
      genus_col <- "Genus"
    }
    
    if(!is.null(genus_col)) {
      # Create mapping: ASV_ID -> Genus
      asv_to_genus <- tax_table %>%
        dplyr::select(ASV_ID, Genus) %>%
        dplyr::filter(!is.na(Genus) & Genus != "")
      
      # Find which ASVs correspond to significant genera
      significant_genera <- unique(dat_sig$taxa)
      
      # Get ASV IDs that belong to significant genera
      significant_asvs <- asv_to_genus %>%
        dplyr::filter(Genus %in% significant_genera) %>%
        dplyr::pull(ASV_ID)
    } else {
      cat("⚠️  WARNING: Cannot map ASV IDs to genera. Using all ASVs.\n")
      significant_asvs <- rownames(taxa_counts)
      asv_to_genus <- data.frame(ASV_ID = character(), Genus = character())
    }
  } else {
    stop("Phyloseq object not found.")
  }
  
  # Filter taxa counts for significant taxa (ASVs) and relevant samples
  # taxa_counts has: rows = ASVs/taxa, columns = samples (phyloseq sample names)
  filtered_taxa <- taxa_counts %>%
    as.data.frame() %>%
    # Keep only significant ASVs (filter rows)
    dplyr::filter(rownames(.) %in% significant_asvs) %>%
    # Keep only relevant samples (select columns)
    dplyr::select(dplyr::any_of(rna_phyloseq_samples)) %>%
    # Transpose: now rows = samples, columns = taxa
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "fecal_sample") %>%
    # Map fecal sample names to expression sample IDs
    # The row names after transpose are the fecal.sample.number values (like "f16", "f20")
    dplyr::left_join(sample_mapping, by = c("fecal_sample" = "fecal.sample.number")) %>%
    dplyr::filter(!is.na(expr_sample_id) & expr_sample_id %in% rna_samples) %>%
    dplyr::select(-fecal_sample, -phyloseq_sample) %>%
    dplyr::rename(Sample = expr_sample_id) %>%
    dplyr::mutate(Sample = as.character(Sample))
  
  # Check if we have data after filtering
  cat("After filtering - Expression data:", nrow(filtered_expr), "genes,", ncol(filtered_expr), "columns\n")
  cat("After filtering - Taxa data:", nrow(filtered_taxa), "samples,", ncol(filtered_taxa), "columns\n")
  cat("Taxa sample names:", paste(head(filtered_taxa$Sample, 5), collapse = ", "), "\n")
  cat("RNA sample names:", paste(head(rna_samples, 5), collapse = ", "), "\n")
  
  if(nrow(filtered_expr) == 0 || ncol(filtered_taxa) <= 1) {
    cat("WARNING: No expression or taxa data available after filtering. Returning empty results.\n")
    cat("filtered_expr rows:", nrow(filtered_expr), "\n")
    cat("filtered_taxa columns:", ncol(filtered_taxa), "\n")
    return(list(
      question_name = question_name,
      treatment_comparison = treatment_comparison,
      correlation_results = data.frame(),
      sig_correlations = data.frame(),
      partial_cor_results = data.frame(),
      sig_partial_correlations = data.frame(),
      failed_pairs = list(),
      deg_sig = deg_sig,
      dat_sig = dat_sig,
      samples = rna_samples
    ))
  }
  
  # Normalize gene expression data (z-score per gene)
  # Handle cases where sd = 0 (constant values)
  normalized_expr <- filtered_expr %>%
    tidyr::pivot_longer(
      cols = -c(gene, gene_id, gene_name),
      names_to = "Sample",
      values_to = "count"
    ) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(
      mean_count = mean(count, na.rm = TRUE),
      sd_count = sd(count, na.rm = TRUE),
      z_score = dplyr::case_when(
        is.na(sd_count) | sd_count == 0 ~ 0,  # Constant values get z-score of 0
        TRUE ~ (count - mean_count) / sd_count
      )
    ) %>%
    dplyr::select(-c(count, mean_count, sd_count)) %>%
    tidyr::pivot_wider(
      names_from = Sample,
      values_from = z_score
    )
  
  # Aggregate ASVs by Genus and normalize taxa data
  # First, aggregate ASVs by Genus before normalization
  # This ensures correlations are calculated at the Genus level, not ASV level
  cat("=== AGGREGATING ASVs BY GENUS AND NORMALIZING TAXA DATA ===\n")
  
  # Initialize taxa_with_genus for use in partial correlation later
  taxa_with_genus <- data.frame()
  
  # Check if we have ASV to Genus mapping
  if(nrow(asv_to_genus) > 0 && ncol(filtered_taxa) > 1) {
    # Pivot filtered_taxa to long format and map ASVs to genera
    taxa_with_genus <- filtered_taxa %>%
      tidyr::pivot_longer(
        cols = -Sample,
        names_to = "ASV_ID",
        values_to = "abundance"
      ) %>%
      # Map ASV IDs to Genus names
      dplyr::left_join(asv_to_genus, by = "ASV_ID") %>%
      # Filter out ASVs without genus assignment
      dplyr::filter(!is.na(Genus) & Genus != "") %>%
      # Aggregate abundances by Genus and Sample (sum all ASVs in same genus)
      dplyr::group_by(Sample, Genus) %>%
      dplyr::summarise(
        abundance = sum(abundance, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::rename(TaxaID = Genus) %>%
      dplyr::mutate(Sample = as.character(Sample))
    
    # Now normalize by Genus (z-score per genus)
    # Handle cases where sd = 0 (constant values)
    prepared_taxa <- taxa_with_genus %>%
      dplyr::group_by(TaxaID) %>%
      dplyr::mutate(
        mean_abundance = mean(abundance, na.rm = TRUE),
        sd_abundance = sd(abundance, na.rm = TRUE),
        z_score = dplyr::case_when(
          is.na(sd_abundance) | sd_abundance == 0 ~ 0,  # Constant values get z-score of 0
          TRUE ~ (abundance - mean_abundance) / sd_abundance
        )
      ) %>%
      dplyr::select(Sample, TaxaID, z_score) %>%
      dplyr::rename(abundance = z_score) %>%
      dplyr::mutate(Sample = as.character(Sample))
  } else {
    # Fallback: if no ASV mapping, use original approach (but still handle zero variance)
    prepared_taxa <- filtered_taxa %>%
      tidyr::pivot_longer(
        cols = -Sample,
        names_to = "TaxaID",
        values_to = "abundance"
      ) %>%
      dplyr::group_by(TaxaID) %>%
      dplyr::mutate(
        mean_abundance = mean(abundance, na.rm = TRUE),
        sd_abundance = sd(abundance, na.rm = TRUE),
        z_score = dplyr::case_when(
          is.na(sd_abundance) | sd_abundance == 0 ~ 0,  # Constant values get z-score of 0
          TRUE ~ (abundance - mean_abundance) / sd_abundance
        )
      ) %>%
      dplyr::select(Sample, TaxaID, z_score) %>%
      dplyr::rename(abundance = z_score) %>%
      dplyr::mutate(Sample = as.character(Sample))
  }
  
  # Prepare expression data for correlation
  prepared_expr <- normalized_expr %>%
    tidyr::pivot_longer(
      cols = -c(gene, gene_id, gene_name),
      names_to = "Sample",
      values_to = "z_score"
    ) %>%
    dplyr::mutate(Sample = as.character(Sample))
  
  # Calculate gene-taxa correlations (Spearman)
  # Note: Correlations are calculated between genes and GENERA (not individual ASVs)
  # TaxaID column contains Genus names from aggregated ASV data
  correlation_results <- prepared_expr %>%
    dplyr::left_join(
      prepared_taxa,
      by = "Sample"
    ) %>%
    dplyr::group_by(gene_id, gene_name, TaxaID) %>%
    dplyr::summarise(
      correlation = cor(z_score, abundance, method = "spearman"),
      n_pairs = sum(!is.na(z_score) & !is.na(abundance)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      abs_correlation = abs(correlation)
    ) %>%
    dplyr::arrange(dplyr::desc(abs_correlation))
  
  # Calculate significance (p-value, FDR) - use correct formula for Spearman
  # Use actual n_pairs from the correlation calculation
  correlation_results <- correlation_results %>%
    dplyr::filter(!is.na(correlation) & n_pairs >= 3) %>%
    dplyr::mutate(
      # Calculate p-value for Spearman correlation using asymptotic approximation
      # For Spearman: t = r * sqrt((n-2)/(1-r^2)), p = 2*pt(-abs(t), df=n-2)
      t_stat = dplyr::case_when(
        abs_correlation >= 0.9999 ~ Inf,  # Perfect correlation, t = Inf
        TRUE ~ correlation * sqrt((n_pairs - 2) / pmax(1 - correlation^2, 1e-10))
      ),
      p_value = dplyr::case_when(
        abs_correlation >= 0.9999 ~ 0,  # Perfect correlation, p = 0
        TRUE ~ 2 * pt(-abs(t_stat), df = n_pairs - 2)
      )
    ) %>%
    dplyr::select(-t_stat) %>%
    dplyr::filter(!is.na(correlation) & !is.na(p_value))
  
  # Calculate FDR correction
  correlation_results <- correlation_results %>%
    dplyr::mutate(
      fdr = p.adjust(p_value, method = "BH")
    )
  
  # Get significant correlations (FDR < 0.05)
  sig_correlations <- correlation_results %>%
    dplyr::filter(fdr < 0.05) %>%
    dplyr::arrange(dplyr::desc(abs_correlation))
  
  cat("Significant correlations:", nrow(sig_correlations), "\n")
  
  # Check if we have any significant correlations before proceeding with partial correlation
  if(nrow(sig_correlations) == 0) {
    cat("WARNING: No significant correlations found. Skipping partial correlation analysis.\n")
    return(list(
      question_name = question_name,
      treatment_comparison = treatment_comparison,
      correlation_results = correlation_results,
      sig_correlations = sig_correlations,
      partial_cor_results = data.frame(),
      sig_partial_correlations = data.frame(),
      failed_pairs = list(),
      deg_sig = deg_sig,
      dat_sig = dat_sig,
      samples = rna_samples
    ))
  }
  
  # Use ALL significant correlations for partial correlation analysis
  # (Previously limited to top_n by top_by, but now using all for complete coverage)
  cat("Using ALL significant correlations for partial correlation analysis (", nrow(sig_correlations), " pairs)\n")
  cat("This ensures all significant gene-taxa pairs are included, not just the top 100.\n")
  
  # Get all gene-taxa pairs from significant correlations
  all_pairs <- sig_correlations %>%
    dplyr::select(gene_id, TaxaID) %>%
    dplyr::distinct()
  
  # Check if we have valid gene-taxa pairs
  if(nrow(all_pairs) == 0 || any(is.na(all_pairs$gene_id)) || any(is.na(all_pairs$TaxaID))) {
    cat("WARNING: No valid gene-taxa pairs found for partial correlation analysis.\n")
    return(list(
      question_name = question_name,
      treatment_comparison = treatment_comparison,
      correlation_results = correlation_results,
      sig_correlations = sig_correlations,
      partial_cor_results = data.frame(),
      sig_partial_correlations = data.frame(),
      failed_pairs = list(),
      deg_sig = deg_sig,
      dat_sig = dat_sig,
      samples = rna_samples
    ))
  }
  
  cat("Total unique gene-taxa pairs to analyze:", nrow(all_pairs), "\n")
  
  # Perform partial correlation analysis (controlling for worm counts)
  cat("Performing partial correlation analysis...\n")
  
  # Convert filtered expression data to wide format with samples as rows
  # Use normalized expression data (already has z-scores)
  expr_wide <- normalized_expr %>%
    dplyr::select(-c(gene, gene_name)) %>%
    tidyr::pivot_longer(
      cols = -gene_id,
      names_to = "Sample",
      values_to = "z_score"
    ) %>%
    tidyr::pivot_wider(
      names_from = gene_id,
      values_from = z_score
    ) %>%
    dplyr::mutate(Sample = as.character(Sample)) %>%
    dplyr::arrange(as.numeric(Sample))
  
  # Convert taxa data to wide format with samples as rows
  # Use z-score normalized taxa data for consistency with simple correlations
  # This ensures both analyses use the same normalization approach for fair comparison
  if(exists("taxa_with_genus") && nrow(taxa_with_genus) > 0) {
    # Normalize aggregated taxa data (z-score per genus) for consistency
    # This matches the normalization used in simple correlation analysis
    taxa_wide <- taxa_with_genus %>%
      dplyr::group_by(TaxaID) %>%
      dplyr::mutate(
        mean_abundance = mean(abundance, na.rm = TRUE),
        sd_abundance = sd(abundance, na.rm = TRUE),
        z_score = dplyr::case_when(
          is.na(sd_abundance) | sd_abundance == 0 ~ 0,  # Constant values get z-score of 0
          TRUE ~ (abundance - mean_abundance) / sd_abundance
        )
      ) %>%
      dplyr::select(Sample, TaxaID, z_score) %>%
      tidyr::pivot_wider(
        names_from = TaxaID,
        values_from = z_score
      ) %>%
      dplyr::mutate(Sample = as.character(Sample)) %>%
      dplyr::arrange(as.numeric(Sample))
  } else {
    # Use prepared_taxa (already z-score normalized)
    taxa_wide <- prepared_taxa %>%
      tidyr::pivot_wider(
        names_from = TaxaID,
        values_from = abundance
      ) %>%
      dplyr::mutate(Sample = as.character(Sample)) %>%
      dplyr::arrange(as.numeric(Sample))
  }
  
  # Create control variable matrix (worm counts/infection burden)
  # Note: Adjust column name based on your metadata structure
  control_matrix <- metadata %>%
    dplyr::filter(Sample %in% rna_samples) %>%
    dplyr::arrange(Sample)
  
  # Check what infection burden columns are available
  infection_cols <- c("Total.Worm.Count", "total_worm_count", "InfectionBurden", "infection_burden", 
                      "WormCount", "worm_count", "Total_Worm_Count")
  available_col <- infection_cols[infection_cols %in% colnames(control_matrix)][1]
  
  if(is.na(available_col)) {
    cat("WARNING: No infection burden column found. Using 0 as control (no partial correlation).\n")
    control_matrix <- matrix(0, nrow = nrow(control_matrix), ncol = 1)
  } else {
    cat("Using column", available_col, "for infection burden control.\n")
    control_matrix <- control_matrix %>%
      dplyr::select(dplyr::all_of(available_col)) %>%
      dplyr::mutate(!!rlang::sym(available_col) := ifelse(is.na(!!rlang::sym(available_col)), 0, !!rlang::sym(available_col))) %>%
      as.matrix()
  }
  
  # Verify sample alignment
  if(!all(expr_wide$Sample == taxa_wide$Sample) || 
     nrow(expr_wide) != nrow(control_matrix)) {
    stop("Sample order mismatch between expression, taxa, and control data")
  }
  
  # Initialize results storage for partial correlations
  partial_cor_results <- list()
  failed_pairs <- list()
  
  # Set up parallel processing
  cl <- parallel::makeCluster(6)  # Create a cluster with 6 cores
  on.exit(parallel::stopCluster(cl))  # Ensure cluster is stopped when done
  
  # Calculate partial correlations for all significant pairs
  cat("Calculating partial correlations for", nrow(all_pairs), "gene-taxa pairs...\n")
  cat("This may take some time depending on the number of pairs.\n")
  
  for(i in 1:nrow(all_pairs)) {
    gene <- all_pairs$gene_id[i]
    taxa <- all_pairs$TaxaID[i]
    
    # Progress indicator every 50 pairs
    if(i %% 50 == 0) {
      cat("Progress: ", i, "/", nrow(all_pairs), " pairs processed (", 
          round(100 * i / nrow(all_pairs), 1), "%)\n")
    }
    
    # Skip if gene or taxa is NA
    if(is.na(gene) || is.na(taxa)) {
      failed_pairs[[paste(gene, taxa, sep = "_")]] <- "NA gene or taxa ID"
      next
    }
    
    # Check if gene and taxa exist in the data
    if(!gene %in% colnames(expr_wide) || !taxa %in% colnames(taxa_wide)) {
      failed_pairs[[paste(gene, taxa, sep = "_")]] <- "Gene or taxa not found in data"
      next
    }
    
    # Get data for this gene-taxa pair
    x <- expr_wide[[gene]]
    y <- taxa_wide[[taxa]]
    
    # Skip if any data is missing
    if(any(is.na(c(x, y)))) {
      failed_pairs[[paste(gene, taxa, sep = "_")]] <- "Missing data"
      next
    }
    
    set.seed(42)
    
    # Calculate partial correlation
    res_pcor <- try(nptest::np.cor.test(
      x = x,
      y = y,
      z = control_matrix,
      partial = TRUE,
      parallel = TRUE,
      cl = cl,
      R = 1000,
      na.rm = TRUE
    ), silent = TRUE)
    
    # Only store results if calculation was successful
    if(!inherits(res_pcor, "try-error")) {
      partial_cor_results[[paste(gene, taxa, sep = "_")]] <- list(
        gene_id = gene,
        TaxaID = taxa,
        correlation = res_pcor$estimate,
        p_value = res_pcor$p.value
      )
    } else {
      failed_pairs[[paste(gene, taxa, sep = "_")]] <- as.character(res_pcor)
    }
  }
  
  # Convert partial correlation results to data frame
  if(length(partial_cor_results) > 0) {
    partial_cor_df <- do.call(rbind, lapply(partial_cor_results, function(x) {
      data.frame(
        gene_id = x$gene_id,
        TaxaID = x$TaxaID,
        correlation = x$correlation,
        p_value = x$p_value
      )
    })) %>%
      dplyr::mutate(
        fdr = p.adjust(p_value, method = "BH"),
        abs_correlation = abs(correlation)
      ) %>%
      dplyr::arrange(dplyr::desc(abs_correlation))
    
    # Get significant partial correlations (FDR < 0.1)
    sig_partial_cor <- partial_cor_df %>%
      dplyr::filter(fdr < 0.1) %>%
      dplyr::arrange(dplyr::desc(abs_correlation))
    
    cat("Significant partial correlations (FDR < 0.1):", nrow(sig_partial_cor), "\n")
  } else {
    partial_cor_df <- data.frame()
    sig_partial_cor <- data.frame()
    cat("No successful partial correlations were calculated.\n")
  }
  
  # Return all results for this comparison
  return(list(
    question_name = question_name,
    treatment_comparison = treatment_comparison,
    correlation_results = correlation_results,
    sig_correlations = sig_correlations,
    partial_cor_results = partial_cor_df,
    sig_partial_correlations = sig_partial_cor,
    failed_pairs = failed_pairs,
    deg_sig = deg_sig,
    dat_sig = dat_sig,
    samples = rna_samples
  ))
}

# Function to perform all pairwise comparisons for a research question
perform_pairwise_analysis <- function(question_name, 
                                     dat_results, 
                                     deg_results, 
                                     taxa_counts, 
                                     expr_counts, 
                                     metadata, 
                                     base_treatment,
                                     comparison_treatments,
                                     top_n = 100,
                                     top_by = c("correlation", "taxa", "gene")) {
  # This function runs DEGxDAT analysis for all pairwise comparisons between a base treatment and a set of comparison treatments.
  # INPUTS:
  #   - question_name: string, name of the research question
  #   - dat_results, deg_results, taxa_counts, expr_counts, metadata: as above
  #   - base_treatment: string, the reference treatment
  #   - comparison_treatments: vector of treatments to compare to base
  # OUTPUT: list with all pairwise results and combined summary tables

  cat("\n=== PAIRWISE ANALYSIS FOR:", question_name, "===\n")
  cat("Base treatment:", base_treatment, "\n")
  cat("Comparison treatments:", paste(comparison_treatments, collapse = ", "), "\n")
  
  # Store results for each pairwise comparison
  pairwise_results <- list()
  
  # Loop over each comparison treatment
  for(i in seq_along(comparison_treatments)) {
    comparison_treatment <- comparison_treatments[i]
    comparison_name <- paste(base_treatment, "vs", comparison_treatment)
    
    cat("\n--- Comparison", i, ":", comparison_name, "---\n")
    
    # Create treatment comparison vector
    treatment_comparison <- c(base_treatment, comparison_treatment)
    
    # Run integrated analysis for this pair
    result <- perform_integrated_analysis(
      question_name = comparison_name,
      dat_results = dat_results,
      deg_results = deg_results,
      taxa_counts = taxa_counts,
      expr_counts = expr_counts,
      metadata = metadata,
      treatment_comparison = treatment_comparison,
      top_n = top_n,
      top_by = top_by
    )
    
    # Add comparison info to result
    result$comparison_name <- comparison_name
    result$base_treatment <- base_treatment
    result$comparison_treatment <- comparison_treatment
    
    pairwise_results[[comparison_name]] <- result
  }
  
  # Combine significant correlations from all pairs (handle empty data frames)
  combined_correlations <- do.call(rbind, lapply(pairwise_results, function(x) {
    if(nrow(x$sig_correlations) > 0) {
      x$sig_correlations %>%
        dplyr::mutate(
          comparison_name = x$comparison_name,
          base_treatment = x$base_treatment,
          comparison_treatment = x$comparison_treatment
        )
    } else {
      data.frame()
    }
  }))
  
  # Combine all correlations from all pairs (handle empty data frames)
  combined_all_correlations <- do.call(rbind, lapply(pairwise_results, function(x) {
    if(nrow(x$correlation_results) > 0) {
      x$correlation_results %>%
        dplyr::mutate(
          comparison_name = x$comparison_name,
          base_treatment = x$base_treatment,
          comparison_treatment = x$comparison_treatment
        )
    } else {
      data.frame()
    }
  }))
  
  # Combine partial correlation results from all pairs
  combined_partial_correlations <- do.call(rbind, lapply(pairwise_results, function(x) {
    if(nrow(x$partial_cor_results) > 0) {
      x$partial_cor_results %>%
        dplyr::mutate(
          comparison_name = x$comparison_name,
          base_treatment = x$base_treatment,
          comparison_treatment = x$comparison_treatment
        )
    } else {
      data.frame()
    }
  }))
  
  # Combine significant partial correlations from all pairs
  combined_sig_partial_correlations <- do.call(rbind, lapply(pairwise_results, function(x) {
    if(nrow(x$sig_partial_correlations) > 0) {
      x$sig_partial_correlations %>%
        dplyr::mutate(
          comparison_name = x$comparison_name,
          base_treatment = x$base_treatment,
          comparison_treatment = x$comparison_treatment
        )
    } else {
      data.frame()
    }
  }))
  
  # Get all unique samples used in any pairwise comparison
  all_samples <- unique(unlist(lapply(pairwise_results, function(x) x$samples)))
  
  # Print summary of results
  cat("\n=== PAIRWISE ANALYSIS SUMMARY ===\n")
  cat("Total comparisons completed:", length(pairwise_results), "\n")
  cat("Comparisons with significant correlations:", sum(sapply(pairwise_results, function(x) nrow(x$sig_correlations) > 0)), "\n")
  cat("Comparisons with significant partial correlations:", sum(sapply(pairwise_results, function(x) nrow(x$sig_partial_correlations) > 0)), "\n")
  cat("Total significant correlations:", nrow(combined_correlations), "\n")
  cat("Total significant partial correlations:", nrow(combined_sig_partial_correlations), "\n")
  
  return(list(
    question_name = question_name,
    base_treatment = base_treatment,
    comparison_treatments = comparison_treatments,
    pairwise_results = pairwise_results,
    combined_correlations = combined_correlations,
    combined_all_correlations = combined_all_correlations,
    combined_partial_correlations = combined_partial_correlations,
    combined_sig_partial_correlations = combined_sig_partial_correlations,
    all_samples = all_samples
  ))
}
# =============================================================================
# Sieler2026 adapters
# =============================================================================

#' Build metadata and expression table for DEGxDAT helpers (expects `Sample` = gut_id_num
#'   matching phyloseq fecal.sample.number without the "f" prefix, and worm burden columns).
#'
#' @param metadata_final From Results/04__DiffGeneExp/Stats/metadata_final.rds (Day 60 gut rows).
#' @param counts_int_filt Integer gene-by-sample matrix; colnames TS047_RoL_RNA_*.
#' @param deg_results_all Tibble from all_treatment_results.rds (contains gene_id and gene_name).
#'
sieler2026_prepare_metadata_expr_for_pairwise <- function(
    metadata_final,
    counts_int_filt,
    deg_results_all) {
  gut_nums <- stringr::str_remove(colnames(counts_int_filt), "^TS047_RoL_RNA_")
  expr_counts <- base::as.data.frame(counts_int_filt, stringsAsFactors = FALSE)
  expr_counts <- tibble::rownames_to_column(expr_counts, var = "gene_id")
  stopifnot(ncol(expr_counts) == length(gut_nums) + 1L)
  colnames(expr_counts) <- c("gene_id", gut_nums)

  gene_names_map <- deg_results_all %>%
    dplyr::group_by(.data$gene_id) %>%
    dplyr::slice_head(n = 1L) %>%
    dplyr::ungroup() %>%
    dplyr::select("gene_id", "gene_name")
  expr_counts <- expr_counts %>%
    dplyr::left_join(gene_names_map, by = "gene_id") %>%
    dplyr::mutate(
      gene = .data$gene_id,
      gene_name = dplyr::coalesce(as.character(.data$gene_name), as.character(.data$gene_id))
    )

  if ("136" %in% colnames(expr_counts) && !"138" %in% colnames(expr_counts)) {
    colnames(expr_counts)[colnames(expr_counts) == "136"] <- "138"
  }

  metadata_out <- metadata_final %>%
    dplyr::mutate(Sample = as.character(.data$gut_id_num))

  if (any(metadata_out$Sample == "136", na.rm = TRUE) &&
      !any(metadata_out$Sample == "138", na.rm = TRUE)) {
    metadata_out <- metadata_out %>%
      dplyr::mutate(Sample = dplyr::if_else(.data$Sample == "136", "138", .data$Sample))
  }

  deg_out <- deg_results_all %>%
    dplyr::mutate(gene = .data$gene_id)

  list(metadata = metadata_out, expr_counts = expr_counts, deg_results = deg_out)
}
