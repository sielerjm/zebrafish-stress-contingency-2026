# Reproduce paper results

**Created by:** Michael Sieler  
**Last updated:** 2026-06-26

## Recommended order

Run from the **repository root**. Each driver sources `Code/00__Setup/00__InitializeEnvironment.R`, which loads the newest dated objects in `Data/r_objects/` and `Data/DADA2/`.

| Step | Script | Depends on | Approx. runtime* |
|------|--------|------------|------------------|
| 0 | `Code/00__Setup/04__DataPreProcess.R` | `Data/DADA2/pseq_*.rds` | 10–30 min |
| 1 | `Code/01__Analysis/00__Overview.R` | preprocessing | < 1 min |
| 2 | `Code/01__Analysis/01__Diversity.R` | preprocessing | 5–15 min |
| 3 | `Code/01__Analysis/02__Composition.R` | preprocessing | 10–30 min |
| 4 | `Code/01__Analysis/03__DiffAbund.R` | preprocessing | 30–90 min |
| 5 | `Code/01__Analysis/04__DiffGeneExp.R` | Salmon counts in `Data/DEG/`; optional `dds_*.rds` from Zenodo | 15–45 min |
| 6 | `Code/01__Analysis/05__Mort-Inf.R` | metadata | 5–10 min |
| 7 | `Code/01__Analysis/06__Taxon-DEG-Mort.R` | modules 03–05 | 10–30 min |
| 8 | `Code/01__Analysis/07__FunctionalAnno.R` | module 04 | 10–20 min |
| 9 | `Code/01__Analysis/08__NeutralModel.R` | preprocessing | 5–15 min |

*Wall-clock times vary with CPU cores and whether MaAsLin3 / DESeq2 re-fit models.

## One-command driver

```r
source("run_all.R")
```

`run_all.R` runs steps 1–9 (assumes preprocessing already produced `ps.list`). It does **not** re-run `04__DataPreProcess.R` by default.

## Single module

```bash
Rscript Code/01__Analysis/02__Composition.R
```

## Results layout

Each module writes to `Results/{NN}__{Name}/`:

- `Figures/` — ggplot/pdf/png exports
- `Tables/` — gt HTML tables (regenerable)
- `Stats/` — `*__bundle.rds` objects consumed by `Code/02__Results/*.Rmd`

Main-text figure panels are copied to `Manuscript/MainFigures/Exports/` when you run `Code/01__Analysis/98__MainFiguresRefresh.R` (not part of `run_all.R`).

## Module 06 dependencies

`06__Taxon-DEG-Mort.R` integrates MaAsLin outputs (03), DESeq2 results (04), and mortality/infection models (05). Run 03–05 first or use bundled RDS already in `Results/`.

## Environment record

Package versions used for the submission build:

`Manuscript/Supplementary/Software_SessionInfo__submission__2026-06-26.csv`

`renv.lock` is planned for the ISME Journal resubmission track; not required for bioRxiv reproduction.

## Zenodo artifacts

If re-fitting DESeq2 from scratch without Zenodo files, module 04 will rebuild `dds_*.rds`. To skip re-fitting, download the Zenodo bundle (see [DATA.md](DATA.md)) before running module 04 or knitting its Rmd.
