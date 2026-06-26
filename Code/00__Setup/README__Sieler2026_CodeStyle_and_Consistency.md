# Sieler2026 code style & consistency guide (for humans + LLM agents)

**Created by:** Michael Sieler  
**Last updated:** 2026-04-25  

This document is the **single source of truth** for how code, tables, figures, and Results reports should be written in this repository so that:

- New analysis modules look and behave like existing modules.
- Tables/figures are consistent across modules (aesthetics + file organization).
- A future LLM/agent can **audit** or **extend** the project without guessing conventions.

Scope: `Code/00__Setup/` + `Code/01__Analysis/` + `Code/02__Results/` + `Results/*` + `Manuscript/Supplementary/` (supplement assembly metadata only; not analysis code) + `Manuscript/MainFigures/` (main-text figure manifest and ISME notes; not analysis code).

---

## Constants (single source of truth)

To prevent “almost the same green” / “almost the same alpha” drift, every module should reuse a small set of constants (defined in helpers and referenced here).

- **Significance alpha defaults**
  - `ALPHA_DEFAULT`: 0.05 (most modules)
  - `ALPHA_ENRICHMENT_DEFAULT`: 0.1 (GO/KEGG style tables, when that is the intended cutoff)
- **GT significance styling**
  - `SIG_FILL_GREEN`: light green fill (currently used in helpers, e.g. `#e6f4ea`)
  - `SIG_TEXT_GREEN`: green text (enforced via `style_gt_significance()`)
- **P-value formatting**
  - `PVALUE_PRINT_THRESHOLD`: 1e-4 (print as `<0.0001`)

Policy: if a module needs different constants (e.g., FDR < 0.1 in module 06), it must be recorded in `bundle$meta` and in GT subtitles.

- **Figure colors (do not fork hex values in drivers)**
  - **Exposure regimes** (eight `Treatment` levels): `treatment_order`, `treatment_colors`, and `treatment_color_scale` in `02__PlotSettings.R`. The same hex values are exposed as named levels via `exposure_regime_colors()` in `03__HelperFunctions.R` (must stay in sync with `treatment_color_scale`).
  - **Prior stressor history** (0 / 1 / 2 prior stressors, i.e. `HistoryLevelNum` or `history_order` labels): `prior_stressor_history_colors_numeric` in `02__PlotSettings.R` (`0` → `#1B9E77`, `1` → `#D95F02`, `2` → `#7570B3`). `history_colors` / `history_color_scale` are derived from that vector. Trend plots (`alpha_diversity_stress_history_trend_plot`, `glmm_binomial_tank_history_numeric_trend_plot` in `03__HelperFunctions.R`) use the same numeric vector so composition panels (PCoA, betadisper) and trend panels match.

---

## Canonical initialization pattern (must be consistent)

### Required bootstrap

- **Analysis drivers** (`Code/01__Analysis/*.R`) must `source()`:
  - `Code/00__Setup/00__InitializeEnvironment.R`
- **Results Rmds** (`Code/02__Results/*.Rmd`) must bootstrap via:
  - `Code/00__Setup/source_project_init.R` then `sieler2026_source_initialize_project(root = root)`

Rationale: these patterns ensure `proj.path`, `path.*`, palettes/themes, and shared helpers are consistent regardless of working directory.

### What init provides (do not re-implement ad hoc)

From `00__InitializeEnvironment.R` + setup chain (`01` → `02` → `03`):

- **Paths**: `path.results`, `path.data`, `path.objects`, `path.dada2`, `path.deg`, etc.
- **Canonical phyloseq RDS pointers**: `path.pseq.cleaned`, `path.pseq.uncleaned` (newest by modification time).
- **Optional loaded objects** (when present under `Data/r_objects/`): `ps.list`, `data.list`
- **Plot defaults**: palettes, `ggplot2::theme_update()` session defaults (moderate sizes; see Figures section), and manuscript ggplot themes in `03__HelperFunctions.R` (`theme_sieler2026_*`).
- **Shared helpers**: GT helpers, manuscript ggplot theme, preprocessing utilities, and module-specific helper stacks.

---

## Repository-wide output contract (files, directories, naming)

### Where outputs go

Every analysis module \(NN\) writes to:

- `Results/NN__ModuleName/Figures/`
- `Results/NN__ModuleName/Tables/`
- `Results/NN__ModuleName/Stats/`

### File naming (prefer stability and grep-ability)

- **Figures**: informative snake-case names, include the core analysis dimension.
  - Good: `pcoa_bray_stress_history_parasite.png`, `alpha_factorial_ATP_shannon.pdf`
- **Tables**:
  - machine-readable CSV: `*.csv`
  - human-readable HTML GT exports: `*.html`
- **Stats**:
  - one module “bundle” RDS, e.g. `composition__gut__bundle.rds`

### Supplementary aggregation naming (tables + figures; hybrid defaults + overrides)

This repo generates many machine-readable tables and vector/raster figures under `Results/NN__.../`. For manuscript submission, we also assemble **supplementary Excel** (CSV sheets) and **supplementary figure PDFs** under `Manuscript/Supplementary/`.

Policy: **prefer deterministic, joinable filenames** so supplementary assembly is mostly automatic, and reserve hand-edited YAML only for true exceptions.

#### Canonical pairing key (`stem`)

For each “supplementary item” (one analysis block / one manuscript-facing result group), choose a stable snake_case **`stem`** that identifies the object across modalities:

- **Tables (CSV)**: `Tables/<stem>__<kind>.csv`
  - Recommended kinds: `table` (primary tidy export), `coef`, `anova`, `contrasts`, etc.
  - Example: `alpha_factorial_atp__table.csv`
- **Figures (PDF for supplement assembly)**: `Figures/<stem>__fig_<variant>.pdf`
  - Example: `alpha_factorial_atp__fig_shannon.pdf`, `alpha_factorial_atp__fig_simpson.pdf`
- **GT HTML exports** (human-readable companion): `Tables/<stem>__<kind>_gt.html` (optional; not required in the supplement Excel)

Why: a shared `stem` makes it trivial to glob “all figures for this table group” and to keep supplementary numbering stable across reruns.

#### Hybrid overrides (when filenames cannot be made joinable yet)

Keep small, explicit, git-tracked YAML files under `Manuscript/Supplementary/`:

- `supplement_overrides__NN__ModuleName.yml` (optional): per-table ordering/grouping/exclusions/captions for the supplement **Excel** TOC.
- `supplement_map__NN__ModuleName.yml` (as needed): ordered list of figure basenames for the supplement **figure PDF** compiler; optional `group:` / `order_within_group:` per figure so manifest `supp_id` aligns with the **table registry** from `sieler2026_supp_build_table_registry()` in [`Code/01__Analysis/98__SupplementShared.R`](../01__Analysis/98__SupplementShared.R).

Precedence (high → low):

1. `exclude: true` in overrides → omit from supplement outputs
2. explicit `supp_id` override (use sparingly)
3. default deterministic ordering rules (module → inferred group → basename)
4. `group` / `order_within_group` tweaks in overrides

#### Supplement numbering convention (manuscript-facing)

Default `supp_id` format is hierarchical and module-ordered (example):

- `1.2.3` = module 1 (e.g. `01__Diversity`), group 2 within that module, table 3 within that group

Figure PDF pages should reuse the **same hierarchical `supp_id` scheme** as tables (`module_index.group_index.within_group_index`), cited as **`Figure S<supp_id>`** in the figure manifest (tables use **`Table S<supp_id>`**). The third index counts items **within that modality** within the shared analysis `group` (so a first figure and first table in the same group can both be `…1`).

Supplementary **Excel** output: each data worksheet tab is named after that row’s `supp_id` (e.g. `1.2.3`, after Excel-safe sanitization). **Paths** in the META sheet, TOC sheet, and the workbook’s manifest CSV are `/<project_folder>/Results/...` (project basename only; no machine-specific prefix such as home directory). The **figure manifest CSV** uses the same path convention (`figure_path`) and the same `group` / `group_index` ordering as the table registry built by `99__BuildSupplementExcel.R` / `99__BuildSupplementFiguresPdf.R` via shared helpers in [`98__SupplementShared.R`](../01__Analysis/98__SupplementShared.R). **Default rebuild:** per-module artifacts under `_build/<mode>__<date>/`, master INDEX xlsx + zip bundles at supplementary root via `Rscript Code/01__Analysis/99__SupplementBuild.R all <mode>` (see [`99__BuildSupplementAll.R`](../01__Analysis/99__BuildSupplementAll.R)). **Legacy optional:** add `--combined` for `__ALL__` outputs under `_build/` plus `Sieler2026_SupplementaryTables__ALL__*` zips.

#### Practical guidance for authors/maintainers

- When adding a new analysis block in a driver, **allocate a `stem` first**, then name **all** exported artifacts around it.
- Prefer **PDF** for vector figures intended for supplementary PDF compilation; keep PNG for web/HTML if useful, but do not rely on PNG for the supplement PDF pipeline unless explicitly supported later.
- If you must keep legacy filenames temporarily, add a **minimal** `supplement_map__*.yml` row for that block only, then remove it after the driver is refactored to canonical names.

#### How to build supplementary submission artifacts (tables workbook + figures PDF)

The manuscript submission artifacts live under:

- `Manuscript/Supplementary/` (root: submission zips + INDEX xlsx + YAML maps + docs)
- `Manuscript/Supplementary/_build/<mode>__YYYY-MM-DD/` (loose per-module xlsx/pdf + manifests)

**Root (submission-facing):**

- `Sieler2026_SupplementaryTables__<mode>__YYYY-MM-DD.zip`
- `Sieler2026_SupplementaryFigures__<mode>__YYYY-MM-DD.zip`
- `Supplementary_Tables__INDEX__<mode>__YYYY-MM-DD.xlsx`

**`_build/` (regenerable loose files):**

- `Supplementary_Tables__<module>__<mode>__YYYY-MM-DD.xlsx` (+ per-module `__manifest.csv`)
- `Supplementary_Figures__<module>__<mode>__YYYY-MM-DD.pdf` (+ manifest) when a figure map exists
- `Supplementary_Tables__INDEX__<mode>__YYYY-MM-DD__manifest.csv`

Optional legacy combined outputs (`--combined` on `all`) under `_build/` plus `Sieler2026_SupplementaryTables__ALL__*` / `Sieler2026_SupplementaryFigures__ALL__*` zips at root.

To rebuild from the repo root:

- **Default (per-module + INDEX)**:
  - `Rscript Code/01__Analysis/99__SupplementBuild.R all submission`
- **Legacy combined only**:
  - `Rscript Code/01__Analysis/99__SupplementBuild.R combined submission`

**Tables (Excel) build rules:**

- Inputs are all `Results/<module>/Tables/*.csv` across module directories matching `^\\d{2}__`.
  - Special-case: module `03__DiffAbund` may also include `Results/03__DiffAbund/Stats/maaslin_*/significant_results.tsv` rows in the combined workbook registry (see the tables manifest).
- The workbook includes:
  - `META` sheet (build date, git SHA if available, package versions, included modules)
  - `TOC` sheet (Table S-citations, module/group, source CSV path, target sheet name)
  - One worksheet per table, with tab name derived from `supp_id` (Excel-safe sanitized).
- The tables manifest CSV is the full registry used to build the workbook, and is the canonical reference for:
  - `supp_id` numbering (`module_index.group_index.within_group_idx`)
  - `citation_table` format (`Table S<supp_id>`)
  - `csv_path` stored as a **portable project-relative** path (e.g., `/Sieler2026/Results/...`, not `/Users/...`)

**Figures (combined PDF) build rules:**

- Inputs are per-module YAML maps under `Manuscript/Supplementary/`:
  - `supplement_map__NN__ModuleName.yml`
  - Each YAML provides an ordered `figures:` list of basenames that should exist under `Results/<module>/Figures/`.
- **PDF-first resolution**: for each basename, the build resolves `*.pdf` first; if a PDF is missing but a `*.png` exists, it is converted to a 1-page PDF and included.
  - Policy: still **prefer saving PDF** from drivers for supplement-intended figures; PNG is treated as a fallback.
- Figure numbering (`supp_id`) is aligned to the module’s **table registry** group indices when possible:
  - Each figure is assigned a `group` via best-effort matching to known table `group` keys for that module, unless a YAML `group:` override is provided per-figure.
  - Figures are then numbered within each group in YAML order (or `order_within_group` if provided).
- The figures manifest CSV is the canonical reference for:
  - `supp_id` and `citation_figure` (`Figure S<supp_id>`)
  - `pdf_start_page` and `pdf_n_pages` within the combined PDF
  - `figure_path` stored as a **portable project-relative** path

**Editing figure maps (`supplement_map__*.yml`) safely:**

- Keep `figures:` entries as **basenames without extensions** (e.g., `pcoa_bray_factorial_ATP`, not `pcoa_bray_factorial_ATP.pdf`).
- If the automatic group matching is wrong or ambiguous, add per-figure fields:
  - `group: <table_group_key>`
  - `order_within_group: <integer>`
- The YAML header should include provenance (recommended):
  - `# Created by: Michael Sieler`
  - `# Date last updated: YYYY-MM-DD`

**Troubleshooting (common issues + fixes):**

- **Combined Excel fails with “No module with Tables/*.csv found”**
  - Confirm at least one module has `Results/NN__*/Tables/*.csv`.
  - If you intentionally moved tables elsewhere, either restore the standard location or update the supplement builder to ingest the new location.

- **Combined figure PDF is skipped (“No supplementary figure maps resolved”)**
  - Confirm `Manuscript/Supplementary/supplement_map__NN__*.yml` exists for at least one module and contains a non-empty `figures:` list.
  - Confirm `include_cover: yes` / `no` is valid YAML (use `yes|no`, not `true|false`, to match current conventions).

- **A figure basename can’t be resolved**
  - Ensure the YAML uses **basenames without extensions**.
  - Ensure the file exists under `Results/<module>/Figures/` as either:
    - `<basename>.pdf` (preferred), or
    - `<basename>.png` (fallback; will be converted to a 1-page PDF)
  - If the figure was saved with a nonstandard suffix, rename it to the canonical basename (preferred) rather than encoding the suffix into YAML.

- **Build stops with “Missing group_index for module … figure group(s)”**
  - This means the figure-to-table `group` matching produced a group string that is not present in that module’s table registry group map.
  - Fix options (preferred → fallback):
    - Rename the figure basename so it matches an existing table `group` key deterministically.
    - Add a per-figure YAML override `group: <existing_table_group_key>`.
    - Add/adjust the module’s supplement grouping map in the shared helpers (only if the group concept truly differs between modalities).

- **Unexpected `supp_id` numbering drift after reruns**
  - Tables: ensure all supplement-intended CSVs follow the stable `stem` + group naming conventions; avoid ad hoc new CSV filenames that create new groups.
  - Figures: preserve YAML order within each module map; if a new figure must be inserted without renumbering, use `order_within_group` to control placement.

- **Manifest paths look “weird” (start with `/Sieler2026/...`)**
  - This is intended: manifests store **portable project-relative** paths, not machine absolute paths. Do not “fix” these to `/Users/...`.

### Bundle contract (required for all modules)

Each driver under `Code/01__Analysis/NN__*.R` must save a single RDS bundle to `Results/.../Stats/` that contains:

- `meta`: model text / provenance / constants used in the module
- `tables_tidy`: plain tibbles/data.frames used to generate GT tables (optional but recommended)
- `table_*`: `gt_tbl` objects for display (preferred)
- `modules`: nested lists mapping named sub-analyses → figures/tables + **relative paths** to saved files

The paired Rmd should load and display from the bundle by default (no re-fitting).

### Bundle schema (recommended keys)

To make bundles predictable for agents and for cross-module tooling, use these conventions:

- **Required**
  - `meta`: list (run_date, script path, inputs, parameters, model formulas, cutoffs)
  - `paths`: list of *portable* (repo-relative preferred) paths for `figures`, `tables`, and `stats_dir`
- **Strongly recommended**
  - `tables`: list of plain data frames/tibbles (the “truth” used to create GT)
  - `table_*`: individual GT objects (named for manuscript tables)
  - `tables_gt`: optional list of GT objects when there are many (instead of many top-level `table_*`)
  - `modules`: nested list keyed by analysis block name (e.g., `exposure_regime_factorial`, `stress_history`, `interaction`, `stratified`)
- **Avoid**
  - absolute paths in `paths` unless you also store a repo-relative version; Results Rmds must not assume absolutes exist.

---

## Standard file headers (templates)

### Analysis driver template (`Code/01__Analysis/NN__*.R`)

Every driver should start with a consistent, skim-friendly header:

- filename
- created by: Michael Sieler
- date last updated (YYYY-MM-DD)
- description (high-level)
- expected input(s) and where they come from
- expected output(s) and exact `Results/` locations
- key cutoffs (alpha/FDR), seeds, and any slow steps

### Results Rmd template (`Code/02__Results/NN__*.Rmd`)

Every Results Rmd should declare:

- that it is display-first and reads a bundle from `Results/<module>/Stats/`
- `params$rerun_analysis` behavior (off by default)
- the path-resolution helpers (`resolve_bundle_path`, `include_png`)

---

## Dependency capture (publication readiness)

For peer-reviewed publication and long-term reproducibility, add one of:

- **Preferred**: `renv` lockfile at the repo root (once dependencies stabilize)
- **Minimum**: each module bundle `meta` should include:
  - `R.version` (from `R.version.string`)
  - a `sessionInfo()` capture (as text) or at least key package versions

Policy: bundles should be self-describing enough that a reviewer can understand the software environment used to generate figures/tables.

---

## Tables (gt) — required conventions

### Non-negotiables

- **All human-facing tables must be GT** (`gt`).
- Do **not** use `knitr::kable()` for results tables in Results Rmds.
  - Exception: if a table is strictly a tiny debug/provenance object, still prefer GT unless there is a strong reason not to.

### Significant p-value styling (canonical rule)

For any “p-like” column (p-value / q-value / FDR / padj / etc.):

- If value < 0.05, the cell must be:
  - **light green fill**
  - **green text**
  - **bold**

This rule applies consistently across all modules and all GT tables.

### Multiple significance thresholds (module-appropriate, but must be explicit)

Some modules intentionally use thresholds other than 0.05. When this happens:

- Keep the **styling rule** consistent (green fill + green text + bold), but
- Apply it to the **module’s documented threshold**, and
- State the threshold in the **GT subtitle** and in the bundle `meta`.

Current intentional examples in this repo:

- Module `06__Taxon-DEG-Mort`: BH FDR on partial-correlation p-values uses **FDR < 0.1** (historical choice aligned with archived helper).
- Module `07__FunctionalAnno`: enrichment uses `p.adjust` cutoffs like **0.1** (and may also use `qvalueCutoff`); plots/tables should label these explicitly.

If a table includes both \(p\) and FDR-like columns, prefer styling the **FDR-adjusted** column, unless the table is explicitly about raw p-values.

### P-like column detection (policy)

Treat the following as p-like by name matching (case-insensitive):

- `p`, `p_value`, `p.value`, `p-value`, `Pr(>...)`, `qvalue`, `q_value`, `q.value`, `fdr`, `padj`, `adj_p`, `adjusted_p`

If a module produces nonstandard names, **rename** them to a conventional key (`p_value`, `q_value`, `fdr`, `padj`) before converting to GT.

For enrichment outputs, treat these as p-like too (common in GO/KEGG exports):

- `p.adjust`, `p.adjusted`, `p_adj`, `p_adj_bh`, and derived transforms like `log_padj` (style should apply to the underlying adjusted p-values, not necessarily the log transform).

### Formatting rules

- **Numbers**: apply consistent rounding (typically 3 decimals), but allow domain exceptions:
  - p-values: display as `"<0.0001"` when below threshold (default threshold: 1e-4)
- **Column labels**: bold column labels unless a table is intentionally minimal.
- **Wide tables**: use a fixed table width to prevent HTML overflow when needed.

### GT checklist (for audits and new modules)

- [ ] Table is created with `gt::gt()` (not `kable()`).
- [ ] All p/q/FDR/padj columns are detected and styled consistently (green fill + green text + bold for <0.05).
- [ ] p-values are formatted consistently (use `"<0.0001"` style for very small values).
- [ ] Numeric columns have consistent rounding/formatting.
- [ ] The table has a clear `tab_header()` title and an informative subtitle (model formula or analysis description).
- [ ] The HTML export uses `gt::gtsave()` into `Results/<module>/Tables/`.
- [ ] A CSV mirror exists for the same result (when appropriate) for downstream reuse.

---

## Figures — aesthetics and reproducibility

### Canonical ggplot2 themes (`03__HelperFunctions.R`)

Use these in analysis drivers (after `source("Code/00__Setup/00__InitializeEnvironment.R")`, which loads `02` then `03`).

| Theme helper | Purpose |
|--------------|---------|
| `theme_sieler2026_publication(base_size = 14, legend_position = "bottom")` | Default for main-text style plots: `theme_classic()`-like axes, consistent base size, bottom legend, subtitle line height and plot margins tuned for export. |
| `theme_sieler2026_composition_layers()` | Adds light horizontal y-grid (`grey(0.87)`), no vertical grid, **bold** legend titles. Compose with `theme_sieler2026_publication()` or add on top of a microViz plot. |
| `theme_sieler2026_composition_figure(...)` | `theme_sieler2026_publication()` + `theme_sieler2026_composition_layers()` — use for composition / beta-diversity figures that share the gut-manuscript strip style. |
| `theme_sieler2026_publication_with_grid(...)` | Publication theme plus light **major** panel grid on **both** continuous axes (`grey(0.88)`); minor grids blank. Use for scatter, histogram, and dotplot-style figures with two quantitative axes so they stay aligned with classic-axis figures while keeping grid context. |
| `theme_sieler2026_trend_panel_grid_major_y()` | Returns a `ggplot2::theme()` fragment (not a full theme): horizontal **major y** grid at scale breaks (`grey(0.87)`, linewidth `0.35`); blanks major x and all minors. Append after `theme_sieler2026_publication()` for stressor-history **trend** plots and **`glmm_binomial_tank_history_numeric_trend_plot`** outputs (mortality / infection trends). |

### Panel grid — how to choose

- **Composition / beta-diversity strip** (stacked relative abundance, default PCoA, betadisper violins): use `theme_sieler2026_composition_figure()` or add `theme_sieler2026_composition_layers()` on top of a microViz plot. Default pattern: light **major y** guides only (no vertical majors, no minors). On **`coord_flip()`** horizontal barplots (e.g. Fig 2.1), this manifests as guides along the **mean relative abundance** axis (the continuous axis in panel coordinates).
- **PCoA main-text Bray-Curtis × parasite (Fig 2.3, `p_pcoa_bray_hp`):** `theme_sieler2026_composition_layers()` **plus** explicit overrides: **major x** matching major y, and **minor x + minor y** (lighter, thinner lines) for both MDS axes. (A Canberra mirror exists as `p_pcoa_canberra_hp` for supplementary comparisons.)
- **Numeric stressor trend + tank-level binomial trend** (Simpson QR main panel, mortality / infection trends): `theme_sieler2026_publication()` **+** `theme_sieler2026_trend_panel_grid_major_y()`.
- **Two-axis scatter / dotplots** (partial-correlation scatter, MaAsLin scatter, KEGG dotplot): `theme_sieler2026_publication_with_grid()`.
- **Neutral model Sloan curve (`plot_neutral_fitted` in `08__NeutralModel.R`):** `theme_sieler2026_publication()` + `neutral_plot_legend_theme()` + `panel.grid.major` (light grey on **both** axes), minors blank.

**Policy for new work:** any **new** manuscript-facing ggplot should follow the patterns above. **Supplementary-only** or legacy figures may still differ; do not bulk-refactor supplement-only outputs unless asked. When promoting a figure to main-text or copying its style, align to this section and the main-text reference table below.

**Line weight (ISME-oriented):** `SIELER2026_MIN_LINEWIDTH_MM` (currently `0.35`) is documented in helpers as a floor near ~1 pt at export; use `max(desired, SIELER2026_MIN_LINEWIDTH_MM)` for strokes that must stay visible when downsized for print (e.g. dot outlines, significance bars).

**Session defaults (`02__PlotSettings.R`):** `ggplot2::theme_update()` sets moderate axis and legend sizes (~11–12 pt effective) so exploratory plots are not dominated by oversized inherited text. **Manuscript-facing drivers should still append a full `theme_sieler2026_*`** so exports are explicit and stable; do not rely on `theme_update()` alone for final figures.

**Special plot types:** `theme_void()` remains appropriate for network-like or axis-free diagrams (e.g. bipartite partial-correlation layout in `06__Taxon-DEG-Mort.R`). Document any other exception in the driver header.

**Neutral model (`08__NeutralModel.R`):** Sloan neutral-model panels use `theme_sieler2026_publication()` plus `neutral_plot_legend_theme()` (bottom, horizontal legend box) where applicable; additional `panel.grid.major` (both axes, light grey) and blank minors; raster exports use **300 dpi** for PNG alongside PDF.

### Palettes and factor ordering

- Exposure regimes (8-level `Treatment`) must be consistent with:
  - `treatment_order` and `treatment_color_scale` in `02__PlotSettings.R`
- Stress history labels and colors must be consistent with:
  - `history_order`, `history_color_scale`, and **`prior_stressor_history_colors_numeric`** in `02__PlotSettings.R` (numeric `0`/`1`/`2` keys for trend plots).

Avoid defining “near duplicates” of these palettes inside module scripts. If a new palette is required, add it to `02__PlotSettings.R` (and document it here).

### Main-text figure manifest (audit trail)

- **Path:** `Manuscript/MainFigures/main_figures_manifest.csv` — one row per main manuscript figure/panel: driver script, relative paths to PNG/PDF, dimensions, dpi, theme notes, ISME-related flags (e.g. composite figure counting).
- **Supporting text:** `Manuscript/MainFigures/isme_gta_figure_requirements.txt` (ISME guide excerpts), `Manuscript/MainFigures/external_asset_color_reference.txt` (hex checklist for non-R Fig 1 / Table 1 assets).

Regenerate or re-audit this manifest when main-text figure basenames or drivers change.

### Main-text figure reference (code as of 2026-04-27)

One row per `panel_id` in `main_figures_manifest.csv`. Summarizes **implemented** ggplot theme and panel grids (not every supplementary variant).

| panel_id | basename | driver | Primary ggplot / helper | Theme and panel grid |
|----------|----------|--------|-------------------------|----------------------|
| 1.1 | ExperimentalDesignSchematic | external | Non-R asset | N/A — align hex with `Manuscript/MainFigures/external_asset_color_reference.txt` |
| 1.0 | ExposureSchematic | external | Non-R asset | N/A |
| 2.1 | genus_relative_abundance_by_treatment | `02__Composition.R` | `p_rel_abund` — `microViz::comp_barplot` + `coord_flip` + `facet_grid(rows = PriorStressorHistory)` | Custom `theme()` sizes + **`theme_sieler2026_composition_layers()`** (major y in panel coordinates → guides on mean relative abundance after flip); black bar outlines |
| 2.2 | simpson_diversity_trend_quasirandom | `01__Diversity.R` | `alpha_diversity_stress_history_trend_plot` → **`p_simpson_qr`** | **`theme_sieler2026_publication` + `theme_sieler2026_trend_panel_grid_major_y()`** |
| 2.3 | pcoa_bray_stress_history_parasite | `02__Composition.R` | **`p_pcoa_bray_hp`** | **`theme_sieler2026_composition_layers()`** + **major x** + **minor x and minor y** |
| 2.4 | betadisper_bray_stress_history_parasite | `02__Composition.R` | **`composition_betadisper_boxplot_history_parasite()`** | **`theme_sieler2026_composition_figure`** (publication + composition_layers → major y only) |
| 3.1 | significant_genes_by_treatment_bar | `04__DiffGeneExp.R` | **`p_treatment_deg`** | **`theme_sieler2026_publication`** only — **no** panel grid (faceted `geom_col`) |
| 4.1 | mortality_prior_stressor_trend_square | `05__Mort-Inf.R` | **`glmm_binomial_tank_history_numeric_trend_plot`** | **`theme_sieler2026_publication` + `theme_sieler2026_trend_panel_grid_major_y()`** |
| 4.2 | infection_prevalence_trend_predicted_square_quasirandom | `05__Mort-Inf.R` | same helper | same as 4.1 |
| 5.1 | maaslin_mortality_top_taxa_scatter_combined__Tank | `03__DiffAbund.R` | driver ggplot | **`theme_sieler2026_publication_with_grid`** |
| 5.2 | all_sig_taxa_partial_cor_scatter | `06__Taxon-DEG-Mort.R` | driver ggplot | **`theme_sieler2026_publication_with_grid`** |
| 5.3 | kegg_focal_four_dotplot | `07__FunctionalAnno.R` | driver dotplot | **`theme_sieler2026_publication_with_grid`**; fill = focal taxon; x = \(-\log_{10}\)(adj. *p*); KEGG x capped at 5 with outlier labels |
| 6.1 | neutral_model_Time60__Genus__focal_four_genera | `08__NeutralModel.R` | **`plot_neutral_fitted()`** → `p_neutral_focal_four_genus` | **`theme_sieler2026_publication` + `neutral_plot_legend_theme()`**; **`panel.grid.major`** both axes (`grey88`), minors blank |

Optional columns `bundle_rds_rel` and `refresh_mode` document how each panel can be refreshed without a full driver (`bundle_ggplot`, `driver_figures_only`, `external`, or `rerun_driver`).

### Manuscript `Exports/` mirror

- `sieler2026_path_main_figure_exports()` and `sieler2026_sync_main_figures_from_manifest()` (in `03__HelperFunctions.R`) copy manifest-listed PNG/PDF into `Manuscript/MainFigures/Exports/` using **`FIG_{panel}__basename`** names (e.g. `FIG_02-1__genus_relative_abundance_by_treatment.png`; see `sieler2026_main_figure_export_filename()`). To drop a main-text figure, edit `main_figures_manifest.csv` first; if it was bundle-refreshable, update `98__MainFiguresRefresh.R` as well.
- Main-numbered analysis drivers call sync once at the end of a successful run for their manifest rows.
- **Fast refresh from bundles:** `Rscript Code/01__Analysis/98__MainFiguresRefresh.R --panels=...` (see `Manuscript/MainFigures/README_exports.md`). Re-exports ggplots **as stored** in each bundle; update bundles with a full run when plot code changes.
- **Partial drivers:** `03__DiffAbund.R`, `04__DiffGeneExp.R`, and `05__Mort-Inf.R` implement `--figures-only` (skips `sieler2026_archive_module_outputs()` so existing `Results/` checkpoints are not wiped).

### Export requirements

- Save both **`pdf`** and **`png`** for manuscript figures unless there is a reason not to; **prefer PDF** for vector submission and supplement PDF assembly.
- Use **300 dpi** for PNG rasters (journal graphical-abstract style minimum is a useful benchmark).
- Use consistent figure sizing within a module; prefer named `width` / `height` (or `fig_in`, etc.) variables at the top of the figure block and mirror them in `Manuscript/MainFigures/main_figures_manifest.csv` when the figure is main-text.

### Figure checklist

- [ ] Uses `theme_sieler2026_publication()`, `theme_sieler2026_composition_figure()`, or `theme_sieler2026_publication_with_grid()` as appropriate (or an explicitly documented exception such as `theme_void()`).
- [ ] If the plot is a **main-text** stressor-history **trend** (`alpha_diversity_stress_history_trend_plot`) or **tank binomial trend** (`glmm_binomial_tank_history_numeric_trend_plot`) on a continuous or percent y-axis, append **`theme_sieler2026_trend_panel_grid_major_y()`** after `theme_sieler2026_publication()` (unless the panel is intentionally grid-free and documented).
- [ ] If the plot matches **Fig 2.3** PCoA layout (`p_pcoa_bray_hp`), apply **`theme_sieler2026_composition_layers()`** plus **major x** and **minor x/y** grid overrides as in `02__Composition.R`.
- [ ] Uses canonical palette(s) for `Treatment`, stress history (`prior_stressor_history_colors_numeric` / `history_color_scale`), and parasite exposure where applicable.
- [ ] Saves both `.pdf` and `.png` to `Results/<module>/Figures/` when the figure is manuscript-facing.
- [ ] File names reflect the analysis dimension (distance metric, stratification, factorial terms, etc.).
- [ ] All randomness is seeded with `set.seed(42)` immediately before stochastic steps.
- [ ] Main-text figures: confirm the row in `Manuscript/MainFigures/main_figures_manifest.csv` still matches the driver’s `ggsave()` width, height, and dpi.

---

## Analysis drivers (`Code/01__Analysis/*.R`) — structure and style

### Module skeleton (recommended)

1. Header block: description, expected input, expected output, date updated.
2. Root assertion + init sourcing (`00__InitializeEnvironment.R`).
3. Path setup: define module `path_res_*` and subfolders; create dirs.
4. Archive previous outputs (copy → verify → clear).
5. Constants and model-text strings (for GT subtitles + bundle meta).
6. Load required objects (`ps.list`, `data.list`, etc.) and validate inputs early.
7. Analysis steps (fit models / compute summaries) with explicit seeds.
8. Create and export figures + tables.
9. Assemble and save module bundle RDS.
10. When the module owns main-text manifest rows, call `sieler2026_sync_main_figures_from_manifest(driver_script = "NN__....R")` so `Manuscript/MainFigures/Exports/` stays current.

### Code style rules (repository conventions)

- Use tidyverse conventions; prefer pipes and tibbles.
- Prefer explicit namespaces for non-base functions where ambiguity exists:
  - `dplyr::`, `tidyr::`, `ggplot2::`, `gt::`, `readr::`, `stringr::`, `emmeans::`, etc.
- Avoid mutable global state beyond what init provides; pass objects into helpers.
- Prefer small helpers in `03__HelperFunctions.R` rather than duplicating code across modules.

### Driver checklist

- [ ] `source("Code/00__Setup/00__InitializeEnvironment.R")` is used (with a helpful error message if missing).
- [ ] Output directories exist and follow the standard `Figures/`, `Tables/`, `Stats/` structure.
- [ ] Previous outputs are archived via `sieler2026_archive_module_outputs()`.
- [ ] Seeds are set (`set.seed(42)`) before stochastic procedures.
- [ ] Tables are GT and follow the p-like styling rules.
- [ ] Figures use canonical themes (`theme_sieler2026_publication`, `theme_sieler2026_composition_figure`, `theme_sieler2026_publication_with_grid`, or `theme_sieler2026_trend_panel_grid_major_y` layered as in the Figures section) and palettes; exported as `.pdf` and `.png` (300 dpi for PNG).
- [ ] A bundle RDS is saved and contains both human-facing GT objects and tidy/raw mirrors as appropriate.

### Table output rule (audit finding)

Several modules currently export **CSV-only** summary tables and rely on the Results Rmd to render them (sometimes with `kable()`), rather than saving an HTML GT export from the driver.

**Policy going forward:** for every manuscript-relevant table:

- Save a `*.csv` **and**
- Save a `*.html` (GT via `gt::gtsave()`)

Exceptions should be documented in the driver header (e.g., “CSV-only because table is used only as an intermediate input to the next module”).

---

## Results Rmds (`Code/02__Results/*.Rmd`) — display-only contract + narrative flow

### Display-first contract

Results Rmds should:

- Load the module bundle from `Results/<module>/Stats/*.rds`.
- Display GT tables and figures from the bundle.
- Only rerun the driver when `params$rerun_analysis: true`.

### Table policy in Rmds (repository decision)

- Remove “Tidy data (CSV mirrors)” preview sections.
- Do not use `knitr::kable()` for results tables.

### Bundle path resolution helpers (required pattern)

Across Results Rmds, there are currently two slightly different patterns for locating saved outputs:

- Some Rmds define `resolve_bundle_path()` + `include_png()`, which tolerates absolute paths saved in bundles *and* paths that should be treated as repo-relative.
- Some Rmds directly build `file.path(proj.path, <bundle_relative_path>)` and then `knitr::include_graphics()`.

**Policy going forward:** every Results Rmd should define and use a common, robust path pattern:

- `resolve_bundle_path(p)` that tries:
  1) `file.exists(p)` as given
  2) `here::here(p)` (treat as repo-relative)
  3) `file.path(root, p)` (treat as repo-relative using the discovered root)
- `include_png(p)` that calls `knitr::include_graphics(normalizePath(...))`

Rationale: bundles may store either absolute paths (from the machine that ran the driver) or relative paths; Results Rmds should be portable across collaborators.

### Narrative flow (standard order across modules)

When a module contains multiple analysis “layers,” use this standard order:

1. **Exposure regime** analyses (factorial A × T × P; 8 regimes)
2. **Stressor history** analyses (HistoryLevel / HistoryLevelNum)
3. **Interaction** analyses (History × Parasite; or analogous interaction)
4. **Stratified** analyses (Parasite within each History stratum; or analogous within-stratum contrasts)
5. Appendices: provenance / session info only if needed

### Results Rmd checklist

- [ ] Bootstraps project root via `source_project_init.R`.
- [ ] Loads the module bundle; errors with a helpful message if missing.
- [ ] Displays GT tables directly (no `kable()`).
- [ ] Uses the standard narrative order when sections exist.
- [ ] Uses `knitr::include_graphics()` for saved figure files, via a path resolver that tolerates relative/absolute paths.

---

## Consistency audit checklist (what an agent should verify)

### Setup layer

- [ ] `00__InitializeEnvironment.R` sources `01__Libraries.R` → `02__PlotSettings.R` → `03__HelperFunctions.R` in that order.
- [ ] `source_project_init.R` is the only root-finding bootstrap used in Rmds.

### Tables

- [ ] Every displayed results table is GT.
- [ ] Significant cells in p-like columns are green fill + green text + bold for <0.05.
- [ ] Formatting of small p-values is consistent (`<0.0001` threshold).

### Figures

- [ ] Theme and palettes are consistent (`theme_sieler2026_*` + canonical colors; numeric stress history uses `prior_stressor_history_colors_numeric`).
- [ ] Exports include both `.png` and `.pdf` where applicable; PNG dpi is 300 for manuscript-quality rasters.

### Analyses / organization

- [ ] Drivers write to the standardized `Results/<module>/Figures|Tables|Stats` structure.
- [ ] Results Rmds are display-first and follow the standard narrative order.
- [ ] Bundles exist and contain all objects the paired Rmd expects.

---

## Current status (as of 2026-04-25)

The audit standardization pass (tables, portability, and Results Rmd display-only policy) has been implemented across the repository.

### Resolved items (kept here as a changelog)

- **GT-only Results tables**: Results Rmds no longer use `knitr::kable()` for results tables.
- **Removed CSV-mirror previews**: “Tidy data (CSV mirrors)” preview sections are removed from Results Rmds per repository policy.
- **Display-only Results Rmds**: Results Rmds no longer write outputs (`dir.create()`, `readr::write_csv()`, etc.). Output writing lives in drivers.
- **Portable figure/table paths**: Results Rmds use the standard `resolve_bundle_path()` + `include_png()` helpers.
- **Driver-owned GT exports**: manuscript-relevant tables are exported by drivers as **CSV + GT HTML** (via `gt::gtsave()`), and bundles store display-ready `table_*` GT objects where applicable.
- **Canonical significance styling**: `style_gt_significance()` enforces **green fill + green text + bold**, with module-appropriate alpha (e.g., enrichment \(p.adjust < 0.1\)).
- **Figure themes and colors**: `theme_sieler2026_publication_with_grid()`, `theme_sieler2026_trend_panel_grid_major_y()` on trend and binomial-trend panels, unified `prior_stressor_history_colors_numeric`, moderated `theme_update()` defaults, neutral-model publication theme and 300 dpi PNG; main-text manifest and reference table under `Manuscript/MainFigures/` (see Figures section).

### Remaining “allowed exceptions” (documented; not treated as problems)

- **Figure-first Results layout**: `Code/02__Results/05__Mort-Inf.Rmd` is intentionally figure-first (plate-style). This is an allowed exception to the narrative ordering rule.
- **Theme exceptions**: network-like or axis-free diagrams may use `theme_void()` when appropriate (document in the driver). Scatter, histogram, dotplot, and neutral-model main panels use `theme_sieler2026_publication` or `theme_sieler2026_publication_with_grid` unless a new exception is explicitly approved and recorded here.

---

## Punch list (fix phase) — deterministic order + concrete file targets

This section is a **work plan** for bringing the repository into full alignment with the policies above. It is written so an LLM/agent can execute changes in a low-backtracking order.

*(All items below are now completed; retained as a historical execution log.)*

### 0) Unify canonical GT styling helpers (do this first; everything else depends on it)

- [x] **Update canonical significance styling to match the rule (green fill + green text + bold)**  
  - **File**: `Code/00__Setup/03__HelperFunctions.R`  
  - **Change**: modify the canonical p-like styling helper so it applies:
    - fill: light green (existing is fine)
    - text: green
    - weight: bold  
  - **End state**: every module that uses the canonical helper automatically meets the table styling requirement.

- [x] **Consolidate redundant helper variants into one canonical path**  
  - **File**: `Code/00__Setup/03__HelperFunctions.R`  
  - **Change**:
    - choose one “p-like column detector” (single regex policy)
    - choose one “formatter” for p-like columns (`<0.0001` threshold rule)
    - choose one “styler” for p-like significance cells  
  - **End state**: `style_gt_significance()`, `highlight_gt_significance()`, and `composition_gt_*` variants are either removed or turned into thin wrappers that delegate to the canonical implementation.

### 1) Make Results Rmds display-only (no file writing; no re-analysis by default)

- [x] **Remove output-writing from `03__DiffAbund.Rmd`**  
  - **File**: `Code/02__Results/03__DiffAbund.Rmd`  
  - **Change**: eliminate `dir.create(...)` and `readr::write_csv(...)` blocks; instead display tables/figures from the driver bundle.  
  - **End state**: Results Rmd only reads the bundle and prints GT tables / figures.

- [x] **Standardize all Results Rmds on the portable path helpers**  
  - **Files**:  
    - `Code/02__Results/05__Mort-Inf.Rmd`  
    - `Code/02__Results/06__Taxon-DEG-Mort.Rmd`  
    - `Code/02__Results/07__FunctionalAnno.Rmd`  
    - `Code/02__Results/08__NeutralModel.Rmd`  
  - **Change**: add a shared pattern:
    - `resolve_bundle_path(p)` (try absolute, `here::here(p)`, then `file.path(root, p)`)
    - `include_png(p)` wrapper around `knitr::include_graphics(normalizePath(...))`  
  - **End state**: every Results Rmd is portable across machines and does not assume bundle paths are absolute.

### 2) Eliminate `knitr::kable()` from Results Rmds (GT-only policy)

- [x] **Remove “Tidy data (CSV mirrors)” sections**  
  - **Files**:  
    - `Code/02__Results/01__Diversity.Rmd`  
    - `Code/02__Results/02__Composition.Rmd`  
  - **Change**: delete those sections entirely (per repository decision), leaving only GT tables and figures from the bundle.

- [x] **Convert module 04 Results Rmd summary tables away from `kable()`**  
  - **File**: `Code/02__Results/04__DiffGeneExp.Rmd`  
  - **Dependency**: requires the module 04 driver to generate GT tables in the bundle (see section 3).  
  - **End state**: module 04 Results Rmd prints GT tables directly from the bundle.

### 3) Ensure every module exports manuscript-relevant tables as both CSV + GT HTML

Aligned modules (reference implementations):
- `Code/01__Analysis/01__Diversity.R` exports `*.csv` + `*.html` via `gt::gtsave()`
- `Code/01__Analysis/02__Composition.R` exports `*.csv` + `*.html` via `gt::gtsave()`

Fix targets:

- [x] **Module 04 (`04__DiffGeneExp`)**  
  - **File**: `Code/01__Analysis/04__DiffGeneExp.R`  
  - **Change**: for each summary table currently exported as CSV, also:
    - build a GT table (using canonical helper + correct alpha)
    - `gt::gtsave()` to `Results/04__DiffGeneExp/Tables/*.html`
    - store GT objects under `bundle$table_*` (or `bundle$tables_gt`) for Results Rmd display  
  - **End state**: module 04 has GT HTML exports and display-ready GT objects.

- [x] **Module 05 (`05__Mort-Inf`)**  
  - **File**: `Code/01__Analysis/05__Mort-Inf.R`  
  - **Change**: identify manuscript-relevant tables (at minimum the ones currently rendered in the Results Rmd) and export GT HTML versions from the driver.  
  - **End state**: `Results/05__Mort-Inf/Tables/*.html` exists alongside the CSVs; Results Rmd displays driver-produced GT objects.

- [x] **Module 03 (`03__DiffAbund`)**  
  - **File**: `Code/01__Analysis/03__DiffAbund.R`  
  - **Change**: move the summary tables currently computed in `03__DiffAbund.Rmd` into the driver; export as CSV + GT HTML; store GT objects in the bundle.  
  - **End state**: Rmd is display-only and no longer writes outputs.

- [x] **Module 06 (`06__Taxon-DEG-Mort`)**  
  - **File**: `Code/01__Analysis/06__Taxon-DEG-Mort.R`  
  - **Change**: export GT HTML for manuscript-relevant summaries (e.g., per-taxon summary; any tables that include p/FDR should use the canonical styling with the module’s stated threshold).  
  - **End state**: `Results/06__Taxon-DEG-Mort/Tables/*.html` exists for key tables.

- [x] **Module 07 (`07__FunctionalAnno`)**  
  - **File**: `Code/01__Analysis/07__FunctionalAnno.R`  
  - **Change**: export GT HTML for GO/KEGG tables (at least top-N previews for manuscript use), using canonical p-like styling (likely alpha = 0.1) on `p.adjust`.  
  - **End state**: `Results/07__FunctionalAnno/Tables/*.html` exists and Results Rmd can display those GT tables without reformatting.

### 4) Apply canonical p-like styling across all GT tables (including enrichment)

- [x] **MaAsLin3 q-values**  
  - **Files**: `Code/01__Analysis/03__DiffAbund.R` and/or `Code/02__Results/03__DiffAbund.Rmd` (after display-only refactor)  
  - **Change**: ensure q-value columns are recognized as p-like and styled at alpha = 0.05.

- [x] **GO/KEGG `p.adjust`**  
  - **Files**: `Code/01__Analysis/07__FunctionalAnno.R` and Results Rmd display  
  - **Change**: style `p.adjust` cells at the module’s explicit threshold (often 0.1), and document that threshold in table subtitles + bundle meta.

### 5) Figure theme and export alignment (completed 2026-04-24)

- [x] **Publication themes across modules**  
  - **Files**: `03__HelperFunctions.R` (`theme_sieler2026_publication_with_grid`), `03__DiffAbund.R`, `06__Taxon-DEG-Mort.R`, `07__FunctionalAnno.R`, `08__NeutralModel.R`, `02__PlotSettings.R` (session `theme_update` sizes; `prior_stressor_history_colors_numeric`).  
  - **End state**: main-text style ggplot figures use `theme_sieler2026_*`; prior stressor hex values are single-sourced; neutral-model PNG exports use 300 dpi; manifest lives under `Manuscript/MainFigures/`.

