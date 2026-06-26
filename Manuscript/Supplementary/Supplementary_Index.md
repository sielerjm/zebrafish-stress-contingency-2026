# Supplementary materials index

**Created by:** Michael Sieler  
**Last updated:** 2026-06-26  
**Build command:** `Rscript Code/01__Analysis/99__SupplementBuild.R all submission`

## Submission packaging (default)

| Root artifact | Contents |
|---------------|----------|
| `Sieler2026_SupplementaryTables__<mode>__<date>.zip` | 8 module `.xlsx` + INDEX `.xlsx` + INDEX manifest + **`Software_SessionInfo__<mode>__<date>.csv`** |
| `Sieler2026_SupplementaryFigures__<mode>__<date>.zip` | 8 module figure PDFs + manifests |
| `Supplementary_Tables__INDEX__<mode>__<date>.xlsx` | Master TOC (`workbook_file` → module workbook + sheet); META includes `software_session_info_csv` |
| `Software_SessionInfo__<mode>__<date>.csv` | Loose copy at supplementary root (also inside tables zip) |

Loose per-module files live under `Manuscript/Supplementary/_build/<mode>__<date>/`.

**Software provenance:** `Software_SessionInfo__*.csv` is written automatically by `sieler2026_supp_finalize_packages()` (or `Rscript Code/01__Analysis/99__SessionInfoExport.R <mode>`). Unified long format: `environment`, `package`, and `analysis_driver` rows (modules 01–08).

Table IDs (`Table S1.5.1`, `Table S2.8.1`, …) are unchanged; only file packaging differs.

## Per-module workbooks (inside tables zip or `_build/`)

| Module workbook | Tables (approx.) | Results subsection |
|-----------------|------------------|--------------------|
| `Supplementary_Tables__01__Diversity__<mode>__<date>.xlsx` | 5 | Alpha diversity |
| `Supplementary_Tables__02__Composition__<mode>__<date>.xlsx` | 9 | Composition / PERMANOVA |
| `Supplementary_Tables__03__DiffAbund__<mode>__<date>.xlsx` | 18 | MaAsLin differential abundance |
| `Supplementary_Tables__04__DiffGeneExp__<mode>__<date>.xlsx` | 3 | DEG contrasts |
| `Supplementary_Tables__05__Mort-Inf__<mode>__<date>.xlsx` | 18 | Mortality / infection |
| `Supplementary_Tables__06__Taxon-DEG-Mort__<mode>__<date>.xlsx` | 9 | Partial-correlation integration |
| `Supplementary_Tables__07__FunctionalAnno__<mode>__<date>.xlsx` | 44 | GO/KEGG enrichment |
| `Supplementary_Tables__08__NeutralModel__<mode>__<date>.xlsx` | 25 | Sloan neutral model |

**Total:** 131 tables across 8 workbooks.

## Legacy combined bundle (optional)

Add `--combined` to `99__SupplementBuild.R all` for archival single-file outputs under `_build/` plus `Sieler2026_SupplementaryTables__ALL__*` / `Sieler2026_SupplementaryFigures__ALL__*` zips.

## Module map (S prefix)

| Module | Table IDs | Figure IDs |
|--------|-----------|------------|
| 01 Diversity | S1.* | S1.* |
| 02 Composition | S2.* | S2.* |
| 03 DiffAbund | S3.* | S3.* |
| 04 DiffGeneExp | S4.* | S4.* |
| 05 Mort-Inf | S5.* | S5.* |
| 06 Taxon-DEG-Mort | S6.* | S6.* |
| 07 FunctionalAnno | S7.* | S7.* |
| 08 NeutralModel | S8.* | S8.* |

## bioRxiv / ISME attachment

Attach `Sieler2026_SupplementaryTables__submission__<date>.zip` and `Sieler2026_SupplementaryFigures__submission__<date>.zip` from the supplementary root, or host on Zenodo with the INDEX manifest. Cite `Table S<supp_id>` in prose as before.

## YAML maps

`supplement_map__NN__*.yml` — figure order and grouping per module (curated; stay at supplementary root).
