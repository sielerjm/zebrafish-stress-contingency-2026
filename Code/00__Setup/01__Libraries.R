# 01__Libraries.R
# Created by: Michael Sieler
# Date last updated: 2026-04-26
#
# Description: Attaches project-wide R packages (general, plotting, microbiome, statistics).
#   Loads `tidyverse` last so its namespaces take precedence over overlapping exports.
#
# Expected input:  none (install missing packages before sourcing).
# Expected output:  Packages on the search path; object `cacheing` set (legacy spelling).

suppressPackageStartupMessages({
  # General -------------------------------------------------------------------
  library(forcats) # order/collapse factor levels (treatments, timepoints, taxonomic ranks)
  library(furrr) # parallel purrr over heavy jobs (per-taxon models, resampling)
  library(here) # stable paths to metadata, phyloseq RDS, and results
  library(knitr) # knit Rmd reports and methods documentation
  library(parallel) # extra cores for permutations, bootstraps, cross-validation
  library(purrr) # map over sample splits, list-columns, nested tidy workflows
  library(rcompanion) # transforms (e.g. Tukey) and companion stats for assays
  library(scales) # axis formatting (log, percent, commas) for abundance and diversity

  # Plotting ------------------------------------------------------------------
  library(ComplexHeatmap) # annotated heatmaps (taxa × samples) with metadata tracks
  library(ggbeeswarm) # beeswarm jitter for alpha diversity and read-count distributions
  library(ggExtra) # marginal histograms/densities on scatter plots (e.g. ordination)
  library(ggplot2) # core figures: diversity, composition, ordination, boxplots
  library(ggVennDiagram) # four-set gene overlap diagrams (functional annotation module)
  library(ggnewscale) # second color/fill scale when layering taxon + treatment
  library(ggpubr) # multi-panel layouts and comparison brackets for group tests
  library(ggrepel) # non-overlapping labels for ASV/taxon or sample names
  library(ggraph) # ggplot2-based network layouts for co-occurrence graphs
  library(gridExtra) # arrange multi-panel publication figures
  library(gt) # polished tables for taxon summaries and model output
  library(igraph) # graph objects underpinning co-occurrence / network analyses
  library(pheatmap) # fast clustered heatmaps for abundance or correlation matrices
  library(RColorBrewer) # palettes for treatments and abundance gradients
  library(tidygraph) # tbl_graph workflows paired with ggraph for networks

  # Microbiome (order retained from legacy script) -----------------------------
  library(maaslin3) # multivariable differential abundance vs sample metadata
  library(phyloseq) # ASV/OTU tables, taxonomy, phylogeny, and sample_data in one object
  library(phyloseqCompanion) # extra helpers around phyloseq objects
  library(vegan) # diversity indices, distances, NMDS/PCoA, adonis/perMANOVA
  library(microbiome) # transformations, diversity, compositional tools for phyloseq
  library(picante) # phylogenetic diversity and tree-based metrics (PD, MNTD, etc.)
  library(microViz) # tidy plotting and stats wrappers for phyloseq (barplots, ordination)

  # Statistics & genomics (order retained from legacy script) -------------------
  library(broom) # tidy data frames from models (useful for taxon or gene tables)
  library(car) # ANOVA variants for multi-factor microbiome experiments
  library(glmmTMB) # GLMMs for counts/proportions with random effects (tanks, subjects)
  library(lme4) # linear mixed models for repeated measures on hosts or blocks
  library(emmeans) # estimated marginal means and contrasts across factors
  library(multcomp) # multiple comparisons across many taxa or outcomes
  library(DESeq2) # count-based models (host RNA-seq; some 16S pipelines use similar GLMs)
  library(SummarizedExperiment) # matrix assays + colData for transcriptomics alongside 16S
  library(clusterProfiler) # GO/KEGG enrichment on host genes tied to infection studies
  library(org.Dr.eg.db) # Danio rerio ID mapping for zebrafish host transcriptomics
  library(GO.db) # Gene Ontology term lookups for enrichment results
  library(AnnotationDbi) # convert gene IDs (Entrez, symbol) for annotation joins
  library(KEGGREST) # pathway queries for interpreting host or metagenome function
  library(ppcor) # partial correlations between taxa and traits while controlling covariates
  library(nptest) # nonparametric tests when abundance residuals are ill-behaved
  library(minpack.lm) # neutral model (Sloan / Burns et al. nls fits)
  library(Hmisc) # Wilson CIs for neutral-model curves
  library(stats4) # MLE / negloglik in neutral-model comparisons

  # Load last: prefer tidyverse exports when names overlap with packages above
  library(tidyverse) # dplyr/tidyr/readr/tibble/purrr/ggplot2 for end-to-end tidy workflows
})

# Legacy project flag (downstream scripts may reference this spelling)
cacheing <- TRUE
