# Supplementary materials (index)

This folder holds submission checklists and tracking tables for **The ISME Journal** manuscript (historical stressors, zebrafish gut microbiome, host transcriptomics).

- **`Manuscript_Figures_Tables_Map.md`** — Crosswalk: draft PDF ↔ `Code/01__Analysis/*.R` ↔ `Results/` bundles.
- **`Supplementary_Index.md`** — Per-module workbook index and build notes.
- **`ISME_PreSubmission_Checklist.md`** — Pre-submission tasks (ISME + Nature artwork).

## Automated supplement assembly (tables + figures)

After a full `Results/<module>/` refresh, rebuild submission bundles from the project root.

## Figure-map YAML safety (prevent accidental rewrites)

The supplementary figure PDF is driven by curated per-module YAML mapping files under `Manuscript/Supplementary/`:

- `supplement_map__01__Diversity.yml`
- `supplement_map__02__Composition.yml`
- `supplement_map__03__DiffAbund.yml`
- `supplement_map__04__DiffGeneExp.yml`
- `supplement_map__05__Mort-Inf.yml`
- `supplement_map__06__Taxon-DEG-Mort.yml`
- `supplement_map__07__FunctionalAnno.yml`
- `supplement_map__08__NeutralModel.yml`

These files should not change during routine supplement builds. Only `99__BuildSupplementAll.R` can write `supplement_map__*.yml`, and only when you explicitly enable map writes.

### Safe default (recommended)

Build **per-module** table workbooks, per-module figure PDFs, and a master INDEX (no map writes):

```bash
Rscript Code/01__Analysis/99__SupplementBuild.R all submission
# or:
Rscript Code/01__Analysis/99__BuildSupplementAll.R submission
```

**Outputs (default):**

- **Submission root** (`Manuscript/Supplementary/`):
  - `Sieler2026_SupplementaryTables__<mode>__<date>.zip` (8 module xlsx + INDEX xlsx + INDEX manifest)
  - `Sieler2026_SupplementaryFigures__<mode>__<date>.zip` (8 module PDFs + manifests)
  - `Supplementary_Tables__INDEX__<mode>__<date>.xlsx` (master TOC; also inside tables zip)
- **Build artifacts** (`Manuscript/Supplementary/_build/<mode>__<date>/`):
  - Per-module `Supplementary_Tables__*.xlsx`, `Supplementary_Figures__*.pdf`, and `__manifest.csv` files
  - Older loose files from the supplementary root are moved to `_build/_archive/` on the next `all` build

**Legacy single-file bundle (optional):**

```bash
Rscript Code/01__Analysis/99__SupplementBuild.R all submission --combined
```

Produces legacy `__ALL__` Excel/PDF under `_build/<mode>__<date>/` plus optional `Sieler2026_SupplementaryTables__ALL__*` / `Sieler2026_SupplementaryFigures__ALL__*` zips at root.

### When you want the scripts to write maps

- Generate missing maps (only creates a YAML when it does not exist):

```bash
Rscript Code/01__Analysis/99__SupplementBuild.R all submission --generate-missing-maps --allow-map-writes
```

- Sync maps to match every figure stem currently on disk (overwrites existing YAMLs; use with care):

```bash
Rscript Code/01__Analysis/99__SupplementBuild.R all submission --sync-figure-maps --allow-map-writes --force-map-overwrite
```

If a map-writing flag is provided without the required opt-in flags, the build will warn and skip map writes, then continue building per-module outputs.

### Optional “lock” workflow (filesystem-level guardrail)

Set curated `supplement_map__*.yml` files to read-only when you do not intend to edit them. See `Code/01__Analysis/99__SupplementMapLock.R` (`lock` / `unlock` / `status`).

## Per-module table workbook

Scans `Results/<module>/Tables/*.csv`, assigns stable `supp_id` values (`<module_index>.<group_index>.<within_group_index>`), and writes `Table S<supp_id>` labels on the TOC sheet. Each data worksheet is named like the `supp_id` (e.g. `1.2.3`). Paths in META/TOC/manifest use `/<project_folder>/Results/...`. Optional overrides: `supplement_overrides__<module>.yml` in this folder.

```bash
Rscript Code/01__Analysis/99__BuildSupplementExcel.R 01__Diversity submission
```

## Per-module figure PDF

Reads `supplement_map__<module>.yml` (ordered `figures:` list). PDF-first under `Results/<module>/Figures/`, else PNG rasterized to one-page PDF. Manifest mirrors table `supp_id` scheme.

```bash
Rscript Code/01__Analysis/99__BuildSupplementFiguresPdf.R 01__Diversity submission
```

## Combined legacy scripts

- [`99__BuildSupplementCombined.R`](../../Code/01__Analysis/99__BuildSupplementCombined.R) — optional `__ALL__` Excel and/or figure PDF (`--tables-only`, `--figures-only`).
- [`99__SupplementBuild.R`](../../Code/01__Analysis/99__SupplementBuild.R) — documented entry point (`all`, `combined`, `excel`, `figures-pdf`).

**R packages:** `openxlsx`, `readr`, `yaml`, `pdftools`; for PNG without `magick`, `png` + `grid`.

**Citation in prose:** Use **`Table S<supp_id>`** / **`Figure S<supp_id>`** where `supp_id` matches the workbook TOC or figure manifest (`supp_id_override` in YAML when editorial numbering must stay fixed).

Last updated: 2026-06-25 (Michael Sieler)
