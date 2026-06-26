# Main-text figure exports (`Manuscript/MainFigures/Exports/`)

**Created by:** Michael Sieler  
**Last updated:** 2026-04-27  

This folder holds **copies** of main-text figure assets for easy sharing, submission packaging, or manuscript linking. **Canonical outputs** remain under `Results/<module>/Figures/` (and paths in `main_figures_manifest.csv` point there).

---

## 1. How files get here

1. Analysis drivers in `Code/01__Analysis/` (`01`–`08`) write PNG/PDF under `Results/.../Figures/`.
2. At the end of each run, the driver calls `sieler2026_sync_main_figures_from_manifest(driver_script = "NN__....R")` in `Code/00__Setup/03__HelperFunctions.R`.
3. That function reads `Manuscript/MainFigures/main_figures_manifest.csv`, skips rows with `driver_script == external_design_asset`, and **copies** each listed `path_png` / `path_pdf` into `Exports/` (only if the source file exists; missing sources produce warnings).

**Manual sync (all non-external manifest rows):** from the Sieler2026 repo root, after sourcing init:

```r
source("Code/00__Setup/00__InitializeEnvironment.R")
sieler2026_sync_main_figures_from_manifest()
```

**Sync only one module’s rows:** pass `driver_script = "02__Composition.R"` (basename is matched). **Sync only selected panels:** `panel_ids = c("2.1", "2.4")`.

---

## 2. Export filename pattern (`FIG_*`)

Each export uses:

`FIG_{NN}-{sub}__<original_basename>.<ext>`

- `NN` and `sub` come from the manifest **`panel_id`**, split on the dot (e.g. `2.1` → `02` and `1`).
- `<original_basename>` is the filename under `Results/.../Figures/` (e.g. `genus_relative_abundance_by_treatment.png`).

**Examples**

| `panel_id` | Example export name |
|------------|---------------------|
| `2.1` | `FIG_02-1__genus_relative_abundance_by_treatment.png` |
| `2.2` | `FIG_02-2__simpson_diversity_trend_quasirandom.pdf` |
| `2.3` | `FIG_02-3__pcoa_bray_stress_history_parasite.pdf` |
| `5.1` | `FIG_05-1__maaslin_mortality_top_taxa_scatter_combined__Tank.pdf` |
| `5.2` | `FIG_05-2__all_sig_taxa_partial_cor_scatter.pdf` |

Implementation: `sieler2026_main_figure_export_filename()` in `Code/00__Setup/03__HelperFunctions.R`.

---

## 3. Fast refresh without re-running full drivers (`98__MainFiguresRefresh.R`)

For panels whose ggplot objects are stored in **bundle** RDS files (see manifest columns `bundle_rds_rel` and `refresh_mode`):

```bash
# From repository root
Rscript Code/01__Analysis/98__MainFiguresRefresh.R --panels=2.2,2.4
Rscript Code/01__Analysis/98__MainFiguresRefresh.R --panels=all
```

This re-`ggsave`s into `Results/.../Figures/` from the bundle, then syncs the touched `panel_id`s into `Exports/`.

**Limitation:** you are re-exporting the **ggplot last saved in the bundle**. If you changed theme or layer **code** only, run the owning driver once (or use `--figures-only` where available) so the bundle is rebuilt, then run `98` again if you only need new dimensions/DPI.

---

## 4. Targeted rebuilds (`--figures-only`)

These drivers **skip** `sieler2026_archive_module_outputs()` so existing `Results/` checkpoints are not wiped:

| Script | Flag | Use case |
|--------|------|----------|
| `03__DiffAbund.R` | `--figures-only` | Fig 5.1 from saved tank CSV (see script header for required files). |
| `04__DiffGeneExp.R` | `--figures-only` | Fig 3.1 from `Stats/fig31_significant_genes_by_treatment_bar.rds`. |
| `05__Mort-Inf.R` | `--figures-only` | Fig 4.1–4.2 (main-text trend panels) from `Stats/mortinf_main_text_figure_ggplots.rds`. |

Example:

```bash
Rscript Code/01__Analysis/05__Mort-Inf.R --figures-only
```

---

## 5. Adding or removing a main-text figure

### Add a new main-text panel

1. Add a row to **`Manuscript/MainFigures/main_figures_manifest.csv`** (`figure_id`, `panel_id`, `basename`, `path_png`, `path_pdf`, `driver_script`, dimensions, `dpi_png`, notes, and optionally `bundle_rds_rel`, `refresh_mode`).
2. Implement or extend the figure in the listed **`driver_script`**, and call `sieler2026_sync_main_figures_from_manifest()` at the end of that driver (pattern matches other `01`–`08` drivers).
3. If the panel should be refreshable from a **bundle** via `98__MainFiguresRefresh.R`, extend **`Code/01__Analysis/98__MainFiguresRefresh.R`**: add the `panel_id` to `bundle_refreshable` and add a matching branch in the dispatch loop.

### Remove a main-text panel

1. **`Manuscript/MainFigures/main_figures_manifest.csv`** — remove (or stop maintaining) the row; this is the **source of truth** for sync and naming.
2. **`Code/01__Analysis/98__MainFiguresRefresh.R`** — if that `panel_id` appears in `bundle_refreshable` or in the `for (pid in panel_ids)` dispatch, remove it so `--panels=all` does not expect missing bundle fields.
3. **`Manuscript/MainFigures/Exports/`** — delete the corresponding `FIG_*` files (or rely on the next full sync only copying rows that still exist in the manifest).
4. **Manuscript / slides / other YAML** — remove any hard-coded references to the old export or `Results/` path if you no longer ship that figure.

---

## 6. Related documentation

- **`main_figures_manifest.csv`** — one row per panel; driver and `Results/` paths.
- **`Code/00__Setup/README__Sieler2026_CodeStyle_and_Consistency.md`** — figure conventions, manifest columns (`bundle_rds_rel`, `refresh_mode`), and driver checklist including sync.
