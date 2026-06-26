# Zenodo — DESeq2 checkpoint bundle

**Created by:** Michael Sieler  
**Last updated:** 2026-06-27

## Published record (v1)

| Field | Value |
|-------|-------|
| DOI | [10.5281/zenodo.20941630](https://doi.org/10.5281/zenodo.20941630) |
| Record | [zenodo.org/records/20941630](https://zenodo.org/records/20941630) |
| File | `Sieler2026_dds_checkpoints.zip` (224.8 MB) |
| License | CC-BY 4.0 (+ OSU custom text on record) |

GitHub links here for download instructions: [DATA.md](DATA.md).

## Contents

| File | Install path after unzip |
|------|--------------------------|
| `dds_treatment.rds` | `Results/04__DiffGeneExp/Stats/` |
| `dds_history.rds` | `Results/04__DiffGeneExp/Stats/` |
| `dds_parasite_history.rds` | `Results/04__DiffGeneExp/Stats/` |
| `dds_history_num.rds` | `Results/04__DiffGeneExp/Stats/` |
| `ZENODO_MANIFEST.txt` | MD5 checksums |

## ISME / future versions

Create **version 2** on the same Zenodo record (e.g. add `renv.lock`, optional composition bundle) — do not replace v1. See project ISME Track B notes in the private repo.

## Regenerate local zip

From private Sieler2026 repo root:

```bash
bash scripts/export_public_repo.sh ~/Projects/zebrafish-stress-contingency-2026-public
# Writes publication_export/zenodo/Sieler2026_dds_checkpoints.zip
```
