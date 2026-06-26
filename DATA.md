# Data dictionary

**Created by:** Michael Sieler  
**Last updated:** 2026-06-26

## Overview

| Path | Description |
|------|-------------|
| `Data/Metadata/metadata.tsv` | Sample IDs, tank, treatment flags, time points |
| `Data/DADA2/` | Phyloseq objects from DADA2 pipeline |
| `Data/r_objects/` | List objects and beta distance matrices used at init |
| `Data/DEG/` | Salmon gene counts and transcriptomics metadata |
| `Data/Context/` | Experimental design narrative and schematic figures |
| `Data/SRA/ACCESSIONS.md` | NCBI BioProject accession (public-safe) |

## Phyloseq (`Data/DADA2/`)

| File | Role |
|------|------|
| `pseq_uncleaned_05052025.rds` | Uncleaned phyloseq object |
| `pseq_cleaned_filtered_2026-04-22.rds` | Canonical cleaned/filtered object (newest) |

Init (`00__InitializeEnvironment.R`) loads the **newest** file matching each pattern by modification time.

## List objects (`Data/r_objects/`)

| File | Role |
|------|------|
| `ps-list__22_04_2026.rds` | Named list of phyloseq subsets (All, Gut, etc.) |
| `data-list__22_04_2026.rds` | Companion metadata / model-ready tables |
| `Rds/beta.dist.mat__2026-04-22.rds` | Bray–Curtis and Canberra distance matrices |

## Transcriptomics (`Data/DEG/`)

| File | Role |
|------|------|
| `salmon.merged.gene_counts_length_scale__Corrected_f136-f138.tsv` | Length-scaled gene counts (canonical for DESeq2) |
| `tx2gene.tsv` | Transcript-to-gene mapping |
| `ROL_MajorExperiment__MetadataSheet__Corrected_05092025.xlsx` | RNA sample sheet |

Duplicate Salmon summaries (unscaled TPM, transcript-level RDS, etc.) are omitted from the public export to reduce size.

## Results bundles (`Results/*/Stats/`)

Serialized lists consumed by analysis drivers and `Code/02__Results/*.Rmd`. Key bundles:

| Bundle | Module |
|--------|--------|
| `diversity__gut__bundle.rds` | 01 Diversity |
| `composition__gut__bundle.rds` | 02 Composition |
| `diffabund__gut__bundle.rds` | 03 Diff abundance |
| `diffgeneexp__bundle.rds` | 04 Differential expression |
| MaAsLin `all_results.tsv` per model | 03 (summaries only in public repo) |

## Zenodo — large DESeq2 checkpoints

These files are **not** in GitHub (~215 MB total):

| File | Approx. size |
|------|--------------|
| `Results/04__DiffGeneExp/Stats/dds_treatment.rds` | ~54 MB |
| `Results/04__DiffGeneExp/Stats/dds_history.rds` | ~54 MB |
| `Results/04__DiffGeneExp/Stats/dds_parasite_history.rds` | ~54 MB |
| `Results/04__DiffGeneExp/Stats/dds_history_num.rds` | ~54 MB |

### Download instructions

1. Open the Zenodo record (DOI to be added after upload — see [ZENODO.md](ZENODO.md)).
2. Download `Sieler2026_dds_checkpoints.zip`.
3. Unzip into `Results/04__DiffGeneExp/Stats/`.
4. Verify checksums against `ZENODO_MANIFEST.txt` in the zip.

Alternatively, re-run `Code/01__Analysis/04__DiffGeneExp.R` to rebuild checkpoints (requires Salmon TSV in `Data/DEG/`).

## NCBI raw reads

- BioProject: [PRJNA1482558](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1482558)
- Public release: 2027-09-01
- Scope: 235 day-60 16S libraries; 72 day-60 intestinal RNA-seq libraries

## Sample identifiers

Treatment coding uses factorial antibiotic (±), temperature (±), and parasite (±) with time points T0, T14, T29, T60. See `Data/Context/ExperimentalDesignContext.md` for the full design narrative.
