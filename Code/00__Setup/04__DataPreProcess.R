# 04__DataPreProcess.R
#
# Description:
#   One-shot preprocessing driver for Sieler2026 (Rules-of-Life zebrafish gut 16S). It takes the
#   merged DADA2 phyloseq + sample sheet, applies screening filters and metadata harmonization, then
#   saves objects the analysis notebooks load via 00__InitializeEnvironment.R.
#
# Experiment context (succinct): factorial antibiotic (A) × temperature (T) × parasite (P) on adult
#   zebrafish; fecal samples across days (Time in sample_data; final day often 60). "History" counts
#   A/T stressors before parasite exposure; parasite exposure is encoded in Parasite and the P leg of
#   Treatment. See Data/Context/ExperimentalDesignContext.md for the full design narrative.
#
# Pipeline order:
#   1) Load uncleaned phyloseq → parasite column fix → tax_fix → apply_pre_analysis_filters →
#      augment_filtered_phyloseq_metadata → save cleaned phyloseq.
#   2) Build ps.list subsets (All, Unexposed, Pre/Post, etc.) for community analyses.
#   3) Re-load uncleaned phyloseq for Mortality/Infection tables (not tax-filtered; same metadata chain).
#   4) Alpha / beta diversity on ps.list (populate_ps_list_alpha_diversity, build_beta_dist_matrices_for_ps_list in 03).
#   5) Save dated ps.list, data.list, and beta distance RDS under Data/r_objects/.
#
# Expected input (README: Data/DADA2/):
#   - Data/DADA2/pseq_uncleaned_*.rds  (merged DADA2 phyloseq + metadata; latest match is used)
#
# Expected output:
#   - Data/DADA2/pseq_cleaned_filtered_YYYY-MM-DD.rds (run date; 00__InitializeEnvironment.R picks latest)
#   - Data/r_objects/ps-list__DD_MM_YYYY.rds, data-list__DD_MM_YYYY.rds
#   - Data/r_objects/Rds/beta.dist.mat__YYYY-MM-DD.rds (beta distance matrices from 03__HelperFunctions.R)
#
# Created by: Michael Sieler
# Last updated: 2026-04-06

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(here))

# Project root: Rscript --file= when run from CLI; otherwise here::i_am so interactive runs in RStudio find Data/
script_path_from_cli <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[startsWith(ca, "--file=")])
  if (length(f) == 1L && nzchar(f)) {
    return(normalizePath(f, winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

.sp <- script_path_from_cli()
if (!is.na(.sp)) {
  proj.path <- normalizePath(file.path(dirname(.sp), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  here::i_am("Code/00__Setup/04__DataPreProcess.R")
  proj.path <- here::here()
}

# Core packages + 03__HelperFunctions (section 2b: DADA2 post-processing filters + metadata augmentation)
path.setup <- file.path(proj.path, "Code", "00__Setup")
source(file.path(path.setup, "01__Libraries.R"))
source(file.path(path.setup, "02__PlotSettings.R"))
source(file.path(path.setup, "03__HelperFunctions.R"))

path.data <- file.path(proj.path, "Data")
path.dada2 <- file.path(path.data, "DADA2")
path.objects <- file.path(path.data, "r_objects")
path.rds <- file.path(path.objects, "Rds")
dir.create(path.dada2, recursive = TRUE, showWarnings = FALSE)
dir.create(path.objects, recursive = TRUE, showWarnings = FALSE)
dir.create(path.rds, recursive = TRUE, showWarnings = FALSE)

# Diversity settings (legacy MicrobiomeProcessing defaults; used by populate_ps_list_* / build_beta_* in 03)
diversity.method <- list()
diversity.method[["alpha"]] <- c("shannon", "inverse_simpson", "observed")
diversity.method[["beta"]] <- c("bray", "canberra")

analysis.ID <- paste0("Sieler2026__DataPreProcess__", Sys.Date())

# Eight factorial cells: A = antibiotics, T = elevated temperature, P = parasite challenge (± each).
treatment_order <- c(
  "A- T- P-", "A- T- P+", "A+ T- P-", "A+ T- P+",
  "A- T+ P-", "A- T+ P+", "A+ T+ P-", "A+ T+ P+"
)

# Canonical prior-stressor mapping (A/T only; parasite status does not change prior-stressor count):
#   A- T- P- = 0, A- T- P+ = 0
#   A+ T- P- = 1, A+ T- P+ = 1
#   A- T+ P- = 1, A- T+ P+ = 1
#   A+ T+ P- = 2, A+ T+ P+ = 2
treatment_prior_stressor_num <- c(
  "A- T- P-" = 0L,
  "A- T- P+" = 0L,
  "A+ T- P-" = 1L,
  "A+ T- P+" = 1L,
  "A- T+ P-" = 1L,
  "A- T+ P+" = 1L,
  "A+ T+ P-" = 2L,
  "A+ T+ P+" = 2L
)

# Screening: ASVs must appear in at least min_prevalence_taxa samples; samples need at least
#   min_sample_reads total counts (after tax filtering). Aligns with prior DataProcessing notebooks.
min_prevalence_taxa <- 3L
min_sample_reads <- 5000L

# Decontam settings (optional).
# If PCR/kit blanks exist (e.g. sample_names like "PCRblank1"), we can remove contaminant taxa
# using decontam prevalence method, then drop the blank samples from the cleaned phyloseq.
run_decontam <- TRUE
decontam_threshold <- 0.5
blank_name_regex <- "(?i)^(pcrblank|kitblank)|blank"

# =============================================================================
# Load uncleaned phyloseq → screen → augment metadata → save cleaned object
# =============================================================================
raw_ps_path <- resolve_latest_rds(path.dada2, "^pseq_uncleaned_.*\\.rds$", error_if_empty = TRUE)
message("Using uncleaned phyloseq (newest mtime): ", basename(raw_ps_path))

message("Loading ", raw_ps_path)
# ps_work: the phyloseq object for the microbiome branch — OTU/ASV table, sample_data (design +
#   phenotypes), and tax_table when present. It is reassigned at each step until written as cleaned RDS.
ps_work <- readRDS(raw_ps_path)
ps_work <- ensure_parasite_column(ps_work)
ps_work <- rename_asvs_sequential(ps_work)

# Standardize ambiguous/unknown taxon labels before dropping taxa (microViz).
if (!is.null(phyloseq::tax_table(ps_work, errorIfNULL = FALSE))) {
  ps_work <- microViz::tax_fix(
    ps_work,
    suffix_rank = "current",
    anon_unique = TRUE,
    unknown = NA
  )
} else {
  warning("No tax_table on phyloseq; skipping tax_fix().")
}

# Optional: blank-based contaminant removal (decontam prevalence).
# We detect blanks by sample name by default; if you later add a dedicated metadata flag,
# swap `neg_controls` to use that column instead.
neg_controls <- stringr::str_detect(phyloseq::sample_names(ps_work), stringr::regex(blank_name_regex))
if (isTRUE(run_decontam) && any(neg_controls, na.rm = TRUE)) {
  message("Running decontam prevalence (n blanks = ", sum(neg_controls, na.rm = TRUE), ") ...")
  decon <- sieler2026_decontam_prevalence(
    ps_work,
    neg = neg_controls,
    threshold = decontam_threshold,
    remove_neg_samples = TRUE
  )
  ps_work <- decon$ps

  # Save decontam summary for auditing.
  decontam_out <- file.path(path.objects, paste0("decontam__", Sys.Date(), ".rds"))
  saveRDS(
    list(
      blank_name_regex = blank_name_regex,
      decontam_threshold = decontam_threshold,
      n_blank_samples = sum(neg_controls, na.rm = TRUE),
      contam_tbl = decon$contam_tbl
    ),
    decontam_out
  )
  message("Saved decontam summary: ", basename(decontam_out))
} else {
  message("Decontam skipped (no blanks detected by regex, or run_decontam=FALSE).")
}

# At this point ps_work still includes all samples and ASVs from DADA2 (after renames + tax_fix).
# apply_pre_analysis_filters() drops rare ASVs, suspect organelle / unresolved-phylum rows, and
# low-read samples; see 03__HelperFunctions.R. Output remains a phyloseq assigned back to ps_work.
ps_work <- apply_pre_analysis_filters(
  ps_work,
  min_prevalence = min_prevalence_taxa,
  min_reads = min_sample_reads
)

# Adds Treatment, History, treatment_group, time_point, Exp_Type, etc. (see 03__HelperFunctions.R).
ps_work <- augment_filtered_phyloseq_metadata(ps_work, treatment_order)

# Dated filename per run (ISO date, same convention as beta.dist.mat__YYYY-MM-DD.rds). Older runs remain
#   on disk unless removed; 00__InitializeEnvironment.R loads the newest pseq_cleaned_filtered_*.rds.
pseq_process_date <- format(Sys.Date(), "%Y-%m-%d")
clean_ps_path <- file.path(path.dada2, paste0("pseq_cleaned_filtered_", pseq_process_date, ".rds"))
saveRDS(ps_work, clean_ps_path)
message("Saved cleaned phyloseq: ", basename(clean_ps_path))

# =============================================================================
# ps.list: named phyloseq subsets for downstream Rmd (microbiome community analyses)
#   All = full cleaned object; filters use Parasite, Antibiotics, Temperature, Time.
# =============================================================================
ps.list <- list()
ps.list[["All"]] <- ps_work

# No parasite challenge over the time course (Parasite==0), or parasite-assigned fish at day 0 only.
ps.list[["Unexposed"]] <- ps.list[["All"]] %>%
  microViz::ps_filter((Parasite == 0) | (Parasite == 1 & Time == 0))

# Any post-baseline sample after at least one stressor factor was applied (non-zero Time).
ps.list[["Exposed"]] <- ps.list[["All"]] %>%
  microViz::ps_filter((Antibiotics == 1 | Temperature == 1 | Parasite == 1) & Time != 0)

# Day 0 (pre-parasite) across all tanks.
ps.list[["PreExposed"]] <- ps.list[["All"]] %>%
  microViz::ps_filter(Time == 0)

# Day 0 excluded (post-baseline time course).
ps.list[["PostExposed"]] <- ps.list[["All"]] %>%
  microViz::ps_filter(Time != 0)

# Final scheduled sampling day (often used for end-state community comparisons).
ps.list[["TimeFinal"]] <- ps.list[["All"]] %>%
  microViz::ps_filter(Time == 60)

# =============================================================================
# data.list: host-outcome tables (Mortality / Infection) — built from *uncleaned* phyloseq
#   so sample counts match the full experiment, not ASV-filtered rows. Mortality uses the same
#   metadata derivation as cleaned microbiome data (augment_filtered_phyloseq_metadata in
#   03__HelperFunctions.R), then flattens to a tibble for dplyr summaries.
#   Infection / Infection_Tank: Day 60 only, parasite-exposed (P+) fish, one row per fish (Sample);
#   prevalence = infected survivors with worms / surviving sampled fish at final time point.
# =============================================================================
ps_raw <- readRDS(raw_ps_path)
ps_raw <- ensure_parasite_column(ps_raw)

data.list <- list()

data.list[["Mortality"]] <- augment_filtered_phyloseq_metadata(ps_raw, treatment_order) %>%
  microViz::samdat_tbl()

# Day-60 survivor slice: one microbiome sample row per fish (fecal Sample ID).
mortality_infection_day60_p1 <- data.list[["Mortality"]] %>%
  dplyr::filter(.data$Time == 60L, .data$Parasite == 1L) %>%
  dplyr::distinct(.data$Sample, .keep_all = TRUE)

data.list[["Infection"]] <- mortality_infection_day60_p1 %>%
  dplyr::group_by(.data$Treatment) %>%
  dplyr::summarise(
    n_survivors_sampled = dplyr::n(),
    n_infected = sum(.data$Total.Worm.Count > 0, na.rm = TRUE),
    mean_worm_burden = mean(.data$Total.Worm.Count, na.rm = TRUE),
    percent_infected = round((.data$n_infected / .data$n_survivors_sampled) * 100, 1),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Treatment = factor(.data$Treatment, levels = treatment_order))

data.list[["Infection_Tank"]] <- mortality_infection_day60_p1 %>%
  dplyr::group_by(.data$Treatment, .data$Tank.ID) %>%
  dplyr::summarise(
    n_survivors_sampled = dplyr::n(),
    n_infected = sum(.data$Total.Worm.Count > 0, na.rm = TRUE),
    mean_worm_burden = mean(.data$Total.Worm.Count, na.rm = TRUE),
    percent_infected = round((.data$n_infected / .data$n_survivors_sampled) * 100, 1),
    HistoryLevel = dplyr::first(.data$HistoryLevel),
    HistoryLevelNum = dplyr::first(.data$HistoryLevelNum),
    A = dplyr::first(.data$A),
    T = dplyr::first(.data$T),
    P = dplyr::first(.data$P),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Treatment = factor(.data$Treatment, levels = treatment_order))

# =============================================================================
# Alpha / beta diversity: enrich ps.list; build beta.dist.mat for saved RDS (03__HelperFunctions.R)
# =============================================================================
message("Alpha diversity metrics on ps.list subsets ...")
ps.list <- populate_ps_list_alpha_diversity(ps.list, rank = "Genus")

message("Beta distance matrices per ps.list subset ...")
beta.dist.mat <- build_beta_dist_matrices_for_ps_list(
  ps.list,
  beta_methods = diversity.method[["beta"]]
)

# =============================================================================
# Save dated outputs (analysis notebooks load the latest ps-list / data-list by date or symlink)
# =============================================================================
date_tag <- format(Sys.Date(), "%d_%m_%Y")
ps_list_file <- file.path(path.objects, paste0("ps-list__", date_tag, ".rds"))
data_list_file <- file.path(path.objects, paste0("data-list__", date_tag, ".rds"))

saveRDS(ps.list, ps_list_file)
saveRDS(data.list, data_list_file)
message("Saved: ", basename(ps_list_file))
message("Saved: ", basename(data_list_file))

beta_out <- file.path(path.rds, paste0("beta.dist.mat__", Sys.Date(), ".rds"))
saveRDS(beta.dist.mat, beta_out)
message("Saved: ", basename(beta_out))

message("Done. analysis.ID = ", analysis.ID)
