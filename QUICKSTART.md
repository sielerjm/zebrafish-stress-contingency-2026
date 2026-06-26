# Quick start

**Created by:** Michael Sieler  
**Last updated:** 2026-06-26

## Prerequisites

- R **4.5.1** or newer (see `Manuscript/Supplementary/Software_SessionInfo__submission__2026-06-26.csv`)
- Git
- ~2 GB free disk for clone + intermediate outputs
- Optional: Zenodo download for DESeq2 checkpoint RDS (module 04)

## 1. Clone and set working directory

```bash
git clone https://github.com/sielerjm/zebrafish-stress-contingency-2026.git
cd zebrafish-stress-contingency-2026
```

In RStudio: **File → Open Project** → `Sieler2026.Rproj`.

## 2. Install packages (once per machine)

From the repo root in R:

```r
source("install_dependencies.R")
```

This installs CRAN and Bioconductor packages listed in `Code/00__Setup/01__Libraries.R`. Installation may take 15–30 minutes on a fresh system.

## 3. Run the full pipeline (optional)

```r
source("run_all.R")
```

Or run one module:

```bash
Rscript Code/01__Analysis/01__Diversity.R
```

## 4. Large DESeq2 files (module 04)

`dds_*.rds` files are not in GitHub (~215 MB total). Download from Zenodo (see [DATA.md](DATA.md)) and place in:

```text
Results/04__DiffGeneExp/Stats/
```

Then re-run or knit `Code/02__Results/04__DiffGeneExp.Rmd`.

## 5. Knit results reports (optional)

After drivers finish:

```r
rmarkdown::render("Code/02__Results/01__Diversity.Rmd")
```

## Common errors

| Error | Fix |
|-------|-----|
| `Set working directory to the Sieler2026 repository root` | `setwd()` to clone root or open `.Rproj` |
| `ps.list$All not found` | Run `Code/00__Setup/04__DataPreProcess.R` first (needs phyloseq RDS in `Data/DADA2/`) |
| `dds_*.rds` missing | Download Zenodo bundle into `Results/04__DiffGeneExp/Stats/` |
| MaAsLin3 install fails | Install `devtools` then follow MaAsLin3 GitHub instructions; see SessionInfo CSV |

## Next steps

- Module dependencies and runtimes: [REPRODUCE.md](REPRODUCE.md)
- File dictionary: [DATA.md](DATA.md)
