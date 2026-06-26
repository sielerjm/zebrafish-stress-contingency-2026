# install_dependencies.R
# Created by: Michael Sieler
# Date last updated: 2026-06-26
#
# Description: Install CRAN and Bioconductor packages required by Sieler2026 analysis
#   drivers. Package list mirrors Code/00__Setup/01__Libraries.R.
#
# Expected input:  Internet access; R 4.5+ recommended.
# Expected output:  Packages installed for local use (no lockfile; see SessionInfo CSV).

message("Sieler2026 — installing analysis dependencies...")
message("Reference versions: Manuscript/Supplementary/Software_SessionInfo__submission__2026-06-26.csv")

# CRAN packages ----------------------------------------------------------------
cran_pkgs <- c(
  "forcats", "furrr", "here", "knitr", "purrr", "rcompanion", "scales",
  "ComplexHeatmap", "ggbeeswarm", "ggExtra", "ggplot2", "ggVennDiagram",
  "ggnewscale", "ggpubr", "ggrepel", "ggraph", "gridExtra", "gt", "igraph",
  "pheatmap", "RColorBrewer", "tidygraph",
  "broom", "car", "glmmTMB", "lme4", "emmeans", "multcomp",
  "ppcor", "nptest", "minpack.lm", "Hmisc",
  "BiocManager", "devtools", "remotes"
)

install_if_missing <- function(pkgs, repos = "https://cloud.r-project.org") {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    message("Installing from CRAN: ", paste(missing, collapse = ", "))
    utils::install.packages(missing, repos = repos)
  }
}

install_if_missing(cran_pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  utils::install.packages("BiocManager")
}

# Bioconductor -----------------------------------------------------------------
bioc_pkgs <- c(
  "DESeq2", "SummarizedExperiment", "clusterProfiler", "org.Dr.eg.db",
  "GO.db", "AnnotationDbi", "KEGGREST"
)

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_bioc) > 0) {
  message("Installing from Bioconductor: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

# GitHub / special installs ----------------------------------------------------
special <- list(
  maaslin3 = function() {
    if (!requireNamespace("maaslin3", quietly = TRUE)) {
      message("Installing maaslin3 from GitHub (Biobakery)...")
      remotes::install_github("biobakery/Maaslin3")
    }
  },
  phyloseqCompanion = function() {
    if (!requireNamespace("phyloseqCompanion", quietly = TRUE)) {
      message("phyloseqCompanion not found — install manually if your workflow requires it.")
    }
  },
  microViz = function() {
    if (!requireNamespace("microViz", quietly = TRUE)) {
      message("Installing microViz from CRAN...")
      utils::install.packages("microViz")
    }
  },
  microbiome = function() {
    if (!requireNamespace("microbiome", quietly = TRUE)) {
      utils::install.packages("microbiome")
    }
  },
  picante = function() {
    if (!requireNamespace("picante", quietly = TRUE)) {
      utils::install.packages("picante")
    }
  },
  vegan = function() {
    if (!requireNamespace("vegan", quietly = TRUE)) {
      utils::install.packages("vegan")
    }
  },
  phyloseq = function() {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      BiocManager::install("phyloseq", update = FALSE, ask = FALSE)
    }
  },
  tidyverse = function() {
    if (!requireNamespace("tidyverse", quietly = TRUE)) {
      utils::install.packages("tidyverse")
    }
  }
)

invisible(lapply(special, function(fn) fn()))

message("Done. Source Code/00__Setup/01__Libraries.R to verify all packages load.")
